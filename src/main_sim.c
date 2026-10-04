/*
 * main_sim.c -- end to end simulation of the infrared core signal chain
 * --------------------------------------------------------------------
 *      blackbody scene
 *        -> detector array response (gain/offset/curvature/1-per-pixel noise)
 *        -> two point non uniformity correction
 *        -> blind pixel detection + 3x3 median replacement
 *        -> dynamic range compression (plateau histogram equalisation)
 *        -> BT.656 style frame output + CRC protected serial frame
 *
 * plus the two actuators that sit next to the detector in a real core:
 *      -> TEC temperature control loop (PID with anti-windup)
 *      -> VCSEL constant current driver with temperature compensation
 *
 * Every number printed here is produced by the code in src/; the same numbers
 * are written to results/ as CSV/TXT so they can be pasted into a report
 * without transcription errors.
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "agc.h"
#include "badpix.h"
#include "frame.h"
#include "ir_lab.h"
#include "nuc.h"
#include "tec_ctrl.h"
#include "vcs_pwr.h"

#define RESULT_DIR "results/"

/* A tiny helper so that every line goes to stdout and results/metrics.txt at
 * the same time (one run, two copies, no transcription). */
static FILE *g_log = NULL;

static void say(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    (void)vprintf(fmt, ap);
    va_end(ap);
    if (g_log != NULL) {
        va_start(ap, fmt);
        (void)vfprintf(g_log, fmt, ap);
        va_end(ap);
    }
}

/* ------------------------------------------------------------------ */
/* 1. TEC temperature loop                                             */
/* ------------------------------------------------------------------ */

static void run_tec_section(void)
{
    tec_plant_cfg_t pl;
    tec_pid_cfg_t pc;
    tec_metrics_t m[3];
    const int aw[3] = {IR_TEC_AW_NONE, IR_TEC_AW_CLAMP, IR_TEC_AW_BACKCALC};
    const char *aw_name[3] = {"no anti-windup", "clamping", "back-calculation"};
    const double dt = 0.02;
    const size_t n = 7500;            /* 150 s */
    const double sp = -20.0, t0 = 25.0;
    double *t_c, *u_c, *t_n, *t_b;
    size_t i;
    FILE *f;
    int k;

    tec_plant_cfg_default(&pl);
    tec_pid_cfg_default(&pc);

    say("\n== 1. TEC temperature loop (cold platform) ==\n");
    say("plant: Rth=%.2f K/W  Cth=%.2f J/K  tau=%.1f s  dead time=%.2f s  "
        "ambient=%.1f C  drive=[%.0f, %.0f] W  sensor LSB=%.3f K\n",
        pl.rth, pl.cth, pl.rth * pl.cth, pl.deadtime_s, pl.ambient,
        pl.u_min, pl.u_max, pl.q_step);
    say("step: %.1f C -> %.1f C (%.0f K)   PID Kp=%.2f W/K Ki=%.2f W/(K*s) "
        "Kd=%.2f  dead zone=%.3f K\n",
        t0, sp, sp - t0, pc.kp, pc.ki, pc.kd, pc.deadzone);

    for (k = 0; k < 3; ++k) {
        pc.aw_mode = aw[k];
        if (tec_run_step(&pl, &pc, sp, t0, dt, n, NULL, NULL, NULL, &m[k]) != IR_OK) {
            say("TEC: simulation failed\n");
            return;
        }
        say("  %-16s overshoot=%6.2f %%  settle(+-0.1 K)=%6.2f s  "
            "steady-state err=%+.4f K  actuator saturated=%6.2f s\n",
            aw_name[k], m[k].overshoot_pct, m[k].settle_time_s,
            m[k].ss_error_K, m[k].sat_time_s);
    }
    say("  anti-windup effect: overshoot %.2f %% -> %.2f %% (-%.1f %% relative), "
        "settling %.2f s -> %.2f s (-%.1f %%), "
        "time spent saturated %.2f s -> %.2f s\n",
        m[0].overshoot_pct, m[1].overshoot_pct,
        (m[0].overshoot_pct > 0.0)
            ? 100.0 * (m[0].overshoot_pct - m[1].overshoot_pct) / m[0].overshoot_pct
            : 0.0,
        m[0].settle_time_s, m[1].settle_time_s,
        (m[0].settle_time_s > 0.0)
            ? 100.0 * (m[0].settle_time_s - m[1].settle_time_s) / m[0].settle_time_s
            : 0.0,
        m[0].sat_time_s, m[1].sat_time_s);

    /* Raw traces for plotting. */
    t_c = (double *)malloc(n * sizeof(double));
    u_c = (double *)malloc(n * sizeof(double));
    t_n = (double *)malloc(n * sizeof(double));
    t_b = (double *)malloc(n * sizeof(double));
    if (t_c != NULL && u_c != NULL && t_n != NULL && t_b != NULL) {
        pc.aw_mode = IR_TEC_AW_CLAMP;
        (void)tec_run_step(&pl, &pc, sp, t0, dt, n, NULL, t_c, u_c, NULL);
        pc.aw_mode = IR_TEC_AW_NONE;
        (void)tec_run_step(&pl, &pc, sp, t0, dt, n, NULL, t_n, NULL, NULL);
        pc.aw_mode = IR_TEC_AW_BACKCALC;
        (void)tec_run_step(&pl, &pc, sp, t0, dt, n, NULL, t_b, NULL, NULL);

        f = ir_fopen_w(RESULT_DIR "tec_step.csv");
        ir_fprintf(f, "# TEC step response %.1f -> %.1f C, sampled every 200 ms\n",
                  t0, sp);
        ir_fprintf(f, "# t_s,T_no_aw_C,T_clamping_C,T_backcalc_C,u_clamping_W\n");
        for (i = 0u; i < n; i += 10u) {
            ir_fprintf(f, "%.2f,%.5f,%.5f,%.5f,%.5f\n", (double)i * dt,
                       t_n[i], t_c[i], t_b[i], u_c[i]);
        }
        fclose(f);
    }
    free(t_c); free(u_c); free(t_n); free(t_b);

    /* Gain sweep: Kp x Ki, with and without anti-windup. */
    {
        const double kps[4] = {2.5, 5.0, 8.0, 12.0};
        const double kis[4] = {0.625, 1.25, 2.5, 5.0};
        FILE *g = ir_fopen_w(RESULT_DIR "tec_gain_sweep.csv");
        int a, b;
        ir_fprintf(g, "# Kp_W_per_K,Ki_W_per_K_s,aw_mode,overshoot_pct,"
                      "settle_s,ss_err_K,sat_s\n");
        for (a = 0; a < 4; ++a) {
            for (b = 0; b < 4; ++b) {
                int w;
                for (w = 0; w < 2; ++w) {
                    tec_metrics_t mm;
                    pc.kp = kps[a];
                    pc.ki = kis[b];
                    pc.aw_mode = (w == 0) ? IR_TEC_AW_NONE : IR_TEC_AW_CLAMP;
                    if (tec_run_step(&pl, &pc, sp, t0, dt, n, NULL, NULL, NULL,
                                     &mm) == IR_OK) {
                        ir_fprintf(g, "%.3f,%.4f,%d,%.3f,%.3f,%.5f,%.3f\n",
                                   kps[a], kis[b], pc.aw_mode, mm.overshoot_pct,
                                   mm.settle_time_s, mm.ss_error_K, mm.sat_time_s);
                    }
                }
            }
        }
        fclose(g);
    }

    /* Robustness: does the loop still work when the thermal model is wrong?
     * NOTE this scan also exposes a physical limit, not a control one: the
     * drive needed to hold -20 C is (T_amb - T_set)/Rth = 45/Rth watts, so any
     * Rth below 45/30 = 1.5 K/W is simply out of the actuator's range. */
    {
        FILE *g = ir_fopen_w(RESULT_DIR "tec_robustness.csv");
        const double rths[5] = {1.2, 1.6, 2.0, 2.4, 2.8};
        const double cths[3] = {3.2, 4.0, 4.8};
        int a, b;
        ir_fprintf(g, "# Rth_K_per_W,Cth_J_per_K,required_W,overshoot_pct,"
                      "settle_s,ss_err_K\n");
        pc.kp = 5.0;
        pc.ki = 2.5;
        pc.aw_mode = IR_TEC_AW_CLAMP;
        for (a = 0; a < 5; ++a) {
            for (b = 0; b < 3; ++b) {
                tec_metrics_t mm;
                tec_plant_cfg_t p2 = pl;
                p2.rth = rths[a];
                p2.cth = cths[b];
                if (tec_run_step(&p2, &pc, sp, t0, dt, n, NULL, NULL, NULL,
                                 &mm) == IR_OK) {
                    ir_fprintf(g, "%.3f,%.3f,%.3f,%.3f,%.3f,%.5f\n",
                               p2.rth, p2.cth, (pl.ambient - sp) / p2.rth,
                               mm.overshoot_pct, mm.settle_time_s, mm.ss_error_K);
                }
            }
        }
        fclose(g);
        say("  robustness scan written (results/tec_robustness.csv); "
            "Rth <= %.1f K/W cannot reach the setpoint at all because it needs "
            "more than the %.0f W of drive available\n",
            (pl.ambient - sp) / (pl.u_max), -pl.u_min);
    }
}

/* ------------------------------------------------------------------ */
/* 2. NUC + blind pixels + AGC chain                                   */
/* ------------------------------------------------------------------ */

static void run_imaging_section(void)
{
    nuc_array_cfg_t ac;
    nuc_truth_t truth;
    nuc_coeff_t k2, k1;
    ir_rand_t rng;
    size_t n, i;
    const double phi_low = 100.0, phi_high = 300.0;
    const double phi_eval = phi_low + 0.68 * (phi_high - phi_low);
    double *s_low, *s_high, *s_scene, *c2, *c1, *resp, *resp_faulty, *ideal;
    double *replaced, *agc_out;
    uint8_t *flags, *injected, *frame_stream;
    badpix_cfg_t bc;
    badpix_replace_cfg_t rc;
    badpix_eval_t ev_naive, ev_clip;
    badpix_residual_t rs;
    agc_scene_t scene;
    agc_scene_cfg_t scc;
    agc_cfg_t gc;
    agc_analysis_t an;
    agc_metrics_t am;
    frame_timing_cfg_t fc;
    frame_report_t fr;
    size_t n_replaced;

    nuc_array_cfg_default(&ac);
    n = ac.w * ac.h;

    s_low    = (double *)malloc(n * sizeof(double));
    s_high   = (double *)malloc(n * sizeof(double));
    s_scene  = (double *)malloc(n * sizeof(double));
    c2       = (double *)malloc(n * sizeof(double));
    c1       = (double *)malloc(n * sizeof(double));
    resp     = (double *)malloc(n * sizeof(double));
    resp_faulty = (double *)malloc(n * sizeof(double));
    ideal    = (double *)malloc(n * sizeof(double));
    replaced = (double *)malloc(n * sizeof(double));
    agc_out  = (double *)malloc(n * sizeof(double));
    flags    = (uint8_t *)malloc(n);
    injected = (uint8_t *)malloc(n);
    frame_stream = (uint8_t *)malloc(23040u);

    if (nuc_truth_init(&truth, &ac, 20240501u) != IR_OK ||
        nuc_coeff_init(&k2, ac.w, ac.h) != IR_OK ||
        nuc_coeff_init(&k1, ac.w, ac.h) != IR_OK ||
        s_low == NULL || s_high == NULL || s_scene == NULL || c2 == NULL ||
        c1 == NULL || resp == NULL || resp_faulty == NULL || ideal == NULL ||
        replaced == NULL || agc_out == NULL || flags == NULL ||
        injected == NULL || frame_stream == NULL) {
        say("imaging section: allocation failed\n");
        return;
    }

    say("\n== 2. Non-uniformity correction (two point) ==\n");
    say("array %ux%u, gain spread %.1f %% RMS, offset spread %.1f counts, "
        "curvature %.2f, temporal noise %.2f counts\n",
        (unsigned)ac.w, (unsigned)ac.h, 100.0 * ac.gain_sigma,
        ac.offset_sigma, ac.curv_sigma, ac.noise_sigma);
    say("NU [%%] = 100 * sigma / |mean| over a uniform scene\n");

    ir_rand_seed(&rng, 4242u);
    (void)nuc_render(&truth, &ac, phi_low, 1, &rng, s_low);
    (void)nuc_render(&truth, &ac, phi_high, 1, &rng, s_high);
    (void)nuc_calib_two_point(&k2, s_low, s_high, phi_low, phi_high);
    (void)nuc_calib_one_point(&k1, s_low, phi_low);

    {
        FILE *f = ir_fopen_w(RESULT_DIR "nuc_nu_vs_flux.csv");
        int q;
        ir_fprintf(f, "# phi,NU_raw_pct,NU_2point_pct,NU_1point_pct,"
                      "mean_corrected\n");
        for (q = 0; q <= 20; ++q) {
            double phi = phi_low + (phi_high - phi_low) * (double)q / 20.0;
            double nb, n2, n1, mc;
            ir_rand_seed(&rng, 909u);
            (void)nuc_render(&truth, &ac, phi, 1, &rng, s_scene);
            (void)nuc_apply(&k2, s_scene, c2);
            (void)nuc_apply(&k1, s_scene, c1);
            nb = nuc_nu(s_scene, n, NULL, NULL);
            n2 = nuc_nu(c2, n, NULL, NULL);
            n1 = nuc_nu(c1, n, NULL, NULL);
            mc = ir_mean(c2, n);
            ir_fprintf(f, "%.1f,%.4f,%.4f,%.4f,%.4f\n", phi, nb, n2, n1, mc);
            if (q == 0 || q == 10 || q == 20 || q == 13) {
                say("  phi=%6.1f  NU raw=%7.3f %%  after 2-point=%6.3f %%  "
                    "after 1-point=%7.3f %%\n", phi, nb, n2, n1);
            }
        }
        fclose(f);
    }

    {
        double nb, n2, n1;
        ir_rand_seed(&rng, 909u);
        (void)nuc_render(&truth, &ac, phi_eval, 1, &rng, s_scene);
        (void)nuc_apply(&k2, s_scene, c2);
        (void)nuc_apply(&k1, s_scene, c1);
        nb = nuc_nu(s_scene, n, NULL, NULL);
        n2 = nuc_nu(c2, n, NULL, NULL);
        n1 = nuc_nu(c1, n, NULL, NULL);
        say("  headline at phi=%.1f (68 %% of the calibration span): "
            "raw %.3f %% -> two point %.3f %% -> single point %.3f %%\n",
            phi_eval, nb, n2, n1);
        say("  two point correction improves NU by %.1fx; single point "
            "(offset only) by %.1fx and it degrades as the scene moves away "
            "from the calibration flux\n",
            (n2 > 0.0) ? nb / n2 : 0.0, (n1 > 0.0) ? nb / n1 : 0.0);
    }

    /* ---- blind pixels on the differential response rate ---- */
    say("\n== 3. Blind pixel detection and replacement ==\n");
    for (i = 0u; i < n; ++i) {
        resp[i] = s_high[i] - s_low[i];
        resp_faulty[i] = resp[i];
        ideal[i] = phi_high - phi_low;   /* what a nominal pixel should read */
    }
    {
        size_t n_inj = 0u, target = n / 50u;   /* 2 % */
        ir_rand_seed(&rng, 31337u);
        memset(injected, 0, n);
        while (n_inj < target) {
            size_t idx = (size_t)(ir_rand_uniform(&rng) * (double)n);
            double factor;
            if (idx >= n || injected[idx]) {
                continue;
            }
            if (n_inj % 5u < 2u) {
                factor = 0.02;       /* dead */
            } else if (n_inj % 5u < 4u) {
                factor = 3.0;        /* hot */
            } else {
                factor = 0.5;        /* weak: 50 % responsivity, the hard case */
            }
            injected[idx] = 1u;
            resp_faulty[idx] = factor * resp[idx];
            ++n_inj;
        }
    }

    badpix_cfg_default(&bc);
    badpix_replace_cfg_default(&rc);
    bc.iterations = 1;
    (void)badpix_detect(resp_faulty, n, &bc, flags, NULL);
    badpix_eval(flags, injected, n, &ev_naive);

    badpix_cfg_default(&bc);
    (void)badpix_detect(resp_faulty, n, &bc, flags, NULL);
    badpix_eval(flags, injected, n, &ev_clip);
    (void)badpix_replace(resp_faulty, ac.w, ac.h, flags, &rc, replaced,
                         &n_replaced);
    badpix_residual(resp_faulty, replaced, ideal, injected, n, &rs);

    say("injected %u blind pixels (%.2f %%): 40 %% dead, 40 %% hot, 20 %% at "
        "50 %% responsivity\n", (unsigned)ev_clip.n_injected,
        100.0 * (double)ev_clip.n_injected / (double)n);
    say("  3-sigma, one pass over all pixels : detected %6.2f %%  "
        "false alarms %5.3f %%\n", ev_naive.detection_rate_pct,
        ev_naive.false_alarm_pct);
    say("  3-sigma, iterated sigma clipping  : detected %6.2f %%  "
        "false alarms %5.3f %%\n", ev_clip.detection_rate_pct,
        ev_clip.false_alarm_pct);
    say("  3x3 median replacement of %u pixels: |error| vs nominal "
        "mean %.2f -> %.2f counts, RMSE %.2f -> %.2f, worst %.2f -> %.2f\n",
        (unsigned)n_replaced, rs.mae_before, rs.mae_after, rs.rmse_before,
        rs.rmse_after, rs.max_before, rs.max_after);

    {
        FILE *f = ir_fopen_w(RESULT_DIR "badpix_summary.csv");
        ir_fprintf(f, "# metric,naive_1pass,iterated_clip\n");
        ir_fprintf(f, "detection_rate_pct,%.4f,%.4f\n",
                   ev_naive.detection_rate_pct, ev_clip.detection_rate_pct);
        ir_fprintf(f, "false_alarm_pct,%.4f,%.4f\n",
                   ev_naive.false_alarm_pct, ev_clip.false_alarm_pct);
        ir_fprintf(f, "detected_count,%u,%u\n", (unsigned)ev_naive.detected,
                   (unsigned)ev_clip.detected);
        ir_fprintf(f, "false_alarm_count,%u,%u\n",
                   (unsigned)ev_naive.false_alarm, (unsigned)ev_clip.false_alarm);
        ir_fprintf(f, "mae_before_counts,%.4f,%.4f\n", rs.mae_before, rs.mae_before);
        ir_fprintf(f, "mae_after_counts,%.4f,%.4f\n", rs.mae_after, rs.mae_after);
        fclose(f);
    }

    /* ---- AGC on a synthetic scene ---- */
    say("\n== 4. Dynamic range compression (AGC) ==\n");
    agc_scene_cfg_default(&scc);
    if (agc_scene_synth(&scene, &scc, 2024u) != IR_OK) {
        say("AGC scene allocation failed\n");
        return;
    }
    agc_cfg_default(&gc);
    agc_metrics(scene.img, n, scene.w, scene.h, scene.target, scene.bg, &am);
    {
        /* Saved before `am` is reused for the output metrics below. */
        const double in_sd = am.sd;
        const double in_local = am.local_contrast;
        const double in_mich = am.tb_michelson_pct;
        const double in_weber = am.tb_weber_pct;
        const double in_cnr = am.cnr;

        say("scene %ux%u: background %.1f with gradient (%.0f, %.0f), target "
            "+%.0f on %d x %d pixels, noise sigma %.1f\n",
            (unsigned)scene.w, (unsigned)scene.h, scc.bg_base, scc.bg_grad_x,
            scc.bg_grad_y, scc.target_delta, scc.target_size, scc.target_size,
            scc.noise_sigma);
        say("  input        sd=%7.3f  local contrast=%.4f  Michelson T/B=%6.2f %%  "
            "Weber T/B=%8.2f %%  CNR=%7.3f\n", in_sd, in_local,
            in_mich, in_weber, in_cnr);

        {
            const int modes[3] = {AGC_MODE_LINEAR, AGC_MODE_HEQ,
                                  AGC_MODE_PLATEAU};
            const char *names[3] = {"linear stretch", "plain HEQ",
                                    "plateau HEQ"};
            FILE *f = ir_fopen_w(RESULT_DIR "agc_modes.csv");
            int k;
            ir_fprintf(f, "# mode,input_sd,output_sd,sd_gain_pct,"
                          "local_contrast,michelson_pct,weber_pct,cnr,"
                          "entropy_bits\n");
            for (k = 0; k < 3; ++k) {
                gc.mode = modes[k];
                (void)agc_analyze(&gc, scene.img, n, &an);
                (void)agc_apply(&gc, &an, scene.img, n, agc_out, NULL);
                agc_metrics(agc_out, n, scene.w, scene.h, scene.target,
                            scene.bg, &am);
                say("  %-14s sd=%7.3f (%+6.1f %%)  local=%.4f (%+6.1f %%)  "
                    "Michelson=%6.2f %%  Weber=%8.2f %%  CNR=%7.3f\n",
                    names[k], am.sd, 100.0 * (am.sd / in_sd - 1.0),
                    am.local_contrast,
                    100.0 * (am.local_contrast / in_local - 1.0),
                    am.tb_michelson_pct, am.tb_weber_pct, am.cnr);
                ir_fprintf(f, "%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                           names[k], in_sd, am.sd,
                           100.0 * (am.sd / in_sd - 1.0), am.local_contrast,
                           am.tb_michelson_pct, am.tb_weber_pct, am.cnr,
                           am.entropy_bits);
            }
            fclose(f);
        }

        /* Plateau sweep: one knob, a whole family of mappings. */
        {
            FILE *f = ir_fopen_w(RESULT_DIR "agc_plateau_sweep.csv");
            const double fr[6] = {0.25, 0.5, 1.0, 1.5, 2.0, 3.0};
            int k;
            ir_fprintf(f, "# plateau_frac,output_sd,local_contrast,"
                          "michelson_pct,weber_pct,cnr\n");
            for (k = 0; k < 6; ++k) {
                agc_cfg_default(&gc);
                gc.mode = AGC_MODE_PLATEAU;
                gc.plateau_frac = fr[k];
                (void)agc_analyze(&gc, scene.img, n, &an);
                (void)agc_apply(&gc, &an, scene.img, n, agc_out, NULL);
                agc_metrics(agc_out, n, scene.w, scene.h, scene.target,
                            scene.bg, &am);
                ir_fprintf(f, "%.2f,%.4f,%.4f,%.4f,%.4f,%.4f\n", fr[k], am.sd,
                           am.local_contrast, am.tb_michelson_pct,
                           am.tb_weber_pct, am.cnr);
            }
            fclose(f);
        }
    }

    /* ---- frame output ---- */
    say("\n== 5. Frame output and timing ==\n");
    {
        uint8_t *img8 = (uint8_t *)malloc(n);
        size_t j;
        if (img8 != NULL) {
            for (j = 0u; j < n; ++j) {
                double v = agc_out[j];
                img8[j] = (uint8_t)((v < 0.0) ? 0.0 : (v > 255.0 ? 255.0 : v + 0.5));
            }
            frame_preset_ir128(&fc);
            if (frame_emit(&fc, img8, fc.h_active, fc.v_active, frame_stream,
                           (size_t)frame_total_clocks(&fc), &fr) == IR_OK) {
                say("  128x128 IR raster: %d active + %d blanking pixels/line, "
                    "%d active + %d blanking lines, %d clock/pixel\n",
                    fc.h_active, fc.h_blank, fc.v_active, fc.v_blank,
                    fc.clocks_per_pixel);
                say("    clocks/frame = %ld (expected %ld)  active pixel clocks "
                    "= %ld (expected %d)  lines = %ld  all counts ok = %s\n",
                    fr.clocks_emitted, fr.total_clocks, fr.active_pixel_clocks,
                    fc.h_active * fc.v_active, fr.lines,
                    fr.counts_ok ? "yes" : "NO");
            }
            free(img8);
        }
    }
    {
        frame_timing_cfg_t b;
        frame_preset_bt656_625(&b);
        say("  BT.656 625/50 reference raster: %d active + %d blanking "
            "pixels/line, %d + %d lines, %d clocks/pixel\n",
            b.h_active, b.h_blank, b.v_active, b.v_blank, b.clocks_per_pixel);
        say("    clocks/line = %ld   clocks/frame = %ld   "
            "= %.3f ms at 27 MHz\n", frame_line_clocks(&b),
            frame_total_clocks(&b), 1000.0 * (double)frame_total_clocks(&b) / 27.0e6);
        {
            frame_report_t br;
            (void)frame_emit(&b, NULL, 0, 0, NULL, 0u, &br);
            say("    walked with no image data: clocks=%ld active=%ld "
                "hsync=%ld vsync=%ld counts ok = %s\n",
                br.clocks_emitted, br.active_pixel_clocks, br.hsync_asserted,
                br.vsync_asserted, br.counts_ok ? "yes" : "NO");
        }
    }

    /* ---- serial frame ---- */
    {
        uint8_t pkt[32];
        size_t plen = 0u, consumed = 0u;
        irproto_frame_t pf;
        uint16_t crc_ok, crc_bad;
        uint8_t payload[16];
        size_t j;

        say("\n== 6. Serial frame protocol (CRC-16/CCITT) ==\n");
        {
            const uint8_t tv[9] = {'1','2','3','4','5','6','7','8','9'};
            say("  CRC-16/CCITT-FALSE self test \"123456789\" -> 0x%04X "
                "(standard vector 0x29B1)\n", ir_crc16_ccitt(tv, 9u));
        }
        for (j = 0u; j < sizeof(payload); ++j) {
            payload[j] = (uint8_t)(j * 7u + 3u);
        }
        (void)irproto_pack(0x21u, payload, sizeof(payload), pkt, sizeof(pkt),
                           &plen);
        crc_ok = ir_crc16_ccitt(&pkt[2], 19u);
        pkt[8] ^= 0x08u;                    /* single bit corruption */
        crc_bad = ir_crc16_ccitt(&pkt[2], 19u);
        pkt[8] ^= 0x08u;                    /* restore */
        {
            int rc = irproto_parse(pkt, plen, &pf, &consumed);
            say("  packed %u payload bytes into %u bytes; parse rc=%d, "
                "type=0x%02X, len=%u, consumed=%u, payload match=%s\n",
                (unsigned)sizeof(payload), (unsigned)plen, rc, pf.type, pf.len,
                (unsigned)consumed,
                (memcmp(pf.payload, payload, sizeof(payload)) == 0) ? "yes" : "no");
        }
        pkt[8] ^= 0x08u;
        say("  one flipped bit in the payload: CRC 0x%04X -> 0x%04X, "
            "parse rc=%d (expect %d = CRC error)\n", crc_ok, crc_bad,
            irproto_parse(pkt, plen, &pf, &consumed), IRPROTO_ERR_CRC);
    }

    nuc_truth_free(&truth);
    nuc_coeff_free(&k2);
    nuc_coeff_free(&k1);
    agc_scene_free(&scene);
    free(s_low); free(s_high); free(s_scene); free(c2); free(c1);
    free(resp); free(resp_faulty); free(ideal); free(replaced); free(agc_out);
    free(flags); free(injected); free(frame_stream);
}

/* ------------------------------------------------------------------ */
/* 3. VCSEL driver                                                     */
/* ------------------------------------------------------------------ */

static void run_vcs_section(void)
{
    vcs_cfg_t c;
    vcs_sweep_cfg_t sc;
    vcs_sweep_metrics_t sm;
    vcs_pi_cfg_t pi;
    vcs_step_metrics_t stm;
    vcs_line_t li, le;
    const char *names[3] = {"fixed current", "temperature feedforward",
                            "closed loop APC"};
    double *pt, *pp, *pi_arr;
    int mode;
    size_t np;

    vcs_cfg_default(&c);
    vcs_sweep_cfg_default(&sc);
    vcs_pi_cfg_default(&pi);
    np = (size_t)sc.n_points;
    pt = (double *)malloc(np * sizeof(double));
    pp = (double *)malloc(np * sizeof(double));
    pi_arr = (double *)malloc(np * sizeof(double));
    if (pt == NULL || pp == NULL || pi_arr == NULL) {
        say("VCS section: allocation failed\n");
        free(pt); free(pp); free(pi_arr);
        return;
    }

    say("\n== 7. VCSEL constant current driver and temperature drift ==\n");
    li = vcs_fit_ith(&c, sc.t_lo, sc.t_hi);
    le = vcs_fit_eta(&c, sc.t_lo, sc.t_hi);
    say("device: Ith(25 C)=%.4f mA, Ith(70 C)=%.4f mA (T0=%.0f K), "
        "eta(25)=%.3f W/A, eta(70)=%.3f W/A, target power %.2f mW\n",
        1000.0 * vcs_ith(&c, 25.0), 1000.0 * vcs_ith(&c, 70.0), c.t0_char,
        vcs_eta(&c, 25.0), vcs_eta(&c, 70.0), 1000.0 * c.p_target);
    say("current needed: %.4f mA at 25 C, %.4f mA at 70 C "
        "(driver limit %.0f mA)\n", 1000.0 * vcs_required_current(&c, 25.0),
        1000.0 * vcs_required_current(&c, 70.0), 1000.0 * c.i_max);
    say("firmware model (least squares over 25..70 C): "
        "Ith = %.6e + %.6e*T,  eta = %.6e + %.6e*T\n",
        li.a, li.b, le.a, le.b);

    {
        FILE *f = ir_fopen_w(RESULT_DIR "vcs_power_sweep.csv");
        ir_fprintf(f, "# t_C,mode,i_set_mA,p_out_mW,deviation_pct\n");
        for (mode = 0; mode < 3; ++mode) {
            if (vcs_sweep(&c, &sc, mode, pt, pp, pi_arr, np, &sm) != IR_OK) {
                continue;
            }
            say("  %-26s max |dP|=%7.3f %%   RMS |dP|=%7.3f %%   "
                "P in [%.3f, %.3f] mW   I in [%.3f, %.3f] mA\n",
                names[mode], sm.max_dev_pct, sm.rms_dev_pct,
                1000.0 * sm.p_min, 1000.0 * sm.p_max,
                1000.0 * sm.i_min, 1000.0 * sm.i_max);
            {
                int k;
                for (k = 0; k < sc.n_points; ++k) {
                    ir_fprintf(f, "%.1f,%s,%.5f,%.6f,%.5f\n", pt[k], names[mode],
                               1000.0 * pi_arr[k], 1000.0 * pp[k],
                               100.0 * (pp[k] - c.p_target) / c.p_target);
                }
            }
        }
        fclose(f);
    }

    /* Inner current loop step response. */
    {
        const size_t n = 2000;
        double dt = pi.dt;
        double *io = (double *)malloc(n * sizeof(double));
        if (io != NULL) {
            (void)vcs_current_step(&c, &pi, 6.0e-3, n, NULL, io, NULL, &stm);
            say("  current loop step to 6.000 mA: overshoot=%.3f %%, "
                "settle(+-1 %%)=%.2f us, steady-state error=%+.4f uA, "
                "final=%.6f mA\n", stm.overshoot_pct, 1e6 * stm.settle_time_s,
                1e6 * stm.ss_error_A, 1000.0 * io[n - 1u]);
            free(io);
        }
        say("  (driver tau = %.0f us, PI Kp=%.1f Ki=%.0f A/(A*s), "
            "control period %.2f us)\n", 1e6 * c.tau_drv, pi.kp, pi.ki, 1e6 * dt);
    }

    /* APC with a noisy monitor photodiode. */
    {
        double rms = 0.0, mx = 0.0;
        if (vcs_apc_noise_study(&c, 70.0, 20000u, 10.0e-6, 555u, &rms, &mx) == IR_OK) {
            say("  APC at 70 C with MPD noise (%.1f uA RMS, MPD %.1f A/W, "
                "%.1f uA LSB): optical power deviation RMS %.3f %%, peak %.3f %%\n",
                1e6 * c.mpd_noise, c.mpd_resp, 1e6 * c.mpd_lsb, rms, mx);
        }
    }

    free(pt); free(pp); free(pi_arr);
}

/* ------------------------------------------------------------------ */

int main(void)
{
    g_log = ir_fopen_w(RESULT_DIR "metrics.txt");
    if (g_log == NULL) {
        printf("warning: could not open results/metrics.txt "
               "(is the results/ directory present?)\n");
    }

    say("ir-core-lab -- end to end infrared core signal chain simulation\n");
    say("All numbers below are produced by this run; see results/*.csv for the "
        "raw traces.\n");

    run_tec_section();
    run_imaging_section();
    run_vcs_section();

    say("\n== summary ==\n");
    say("This is a pure numerical simulation: no detector, no TEC, no VCSEL, "
        "no MCU.\n");
    say("See README.md \"Boundaries\" for what is and is not modelled.\n");

    if (g_log != NULL) {
        fclose(g_log);
    }
    return 0;
}
