/*
 * test_ir.c -- self contained test suite for ir-core-lab
 * -----------------------------------------------------
 * No test framework, no dependencies: every check is a CHECK() macro that
 * increments a counter and records a failure with file:line.  The process
 * exits non-zero if anything failed, and the last line printed on success is
 *
 *      <N> checks passed
 *
 * Covered here:
 *   * ir_lab      RNG determinism, statistics, NU definition, file helpers
 *   * tec_ctrl    PID dead zone, output clamp, all three anti-windup modes,
 *                 FOPDT step response metrics, error paths
 *   * nuc         two point algebra exactness, NU before/after, one point
 *                 comparison, dead pixel guard, error paths
 *   * badpix      3-sigma detection, naive vs iterated clipping, borders,
 *                 median/mean replacement, residual, guard rails
 *   * agc         histogram, percentile clipping, monotonicity, plateau -> HEQ
 *                 and plateau -> linear limits, contrast metrics
 *   * frame       BT.656 reference rasters, state machine counters, stream
 *                 content, CRC-16/CCITT standard vectors, framing round trip,
 *                 corruption detection
 *   * vcs_pwr     threshold/slope temperature model, current loop step
 *                 response, fixed vs feedforward vs APC drift compensation
 */

#include <math.h>
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

static int g_checks = 0;
static int g_fails = 0;

#define CHECK(cond)                                                          \
    do {                                                                     \
        ++g_checks;                                                          \
        if (!(cond)) {                                                       \
            ++g_fails;                                                       \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);           \
        }                                                                    \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                \
    do {                                                                     \
        double va_ = (a);                                                    \
        double vb_ = (b);                                                    \
        ++g_checks;                                                          \
        if (!(fabs(va_ - vb_) <= (tol))) {                                   \
            ++g_fails;                                                       \
            printf("FAIL %s:%d  |%s - %s| = |%.12g - %.12g| = %.3g > %.3g\n",\
                   __FILE__, __LINE__, #a, #b, va_, vb_,                     \
                   fabs(va_ - vb_), (double)(tol));                          \
        }                                                                    \
    } while (0)

/* ================================================================== */
/* 1. ir_lab                                                           */
/* ================================================================== */

static void test_ir_lab(void)
{
    ir_rand_t a, b, c;
    double x[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[5] = {2.0, 4.0, 6.0, 8.0, 10.0};
    double mean = 0.0, sd = 0.0, nu;
    int i;

    /* RNG determinism */
    ir_rand_seed(&a, 12345u);
    ir_rand_seed(&b, 12345u);
    ir_rand_seed(&c, 54321u);
    CHECK(ir_rand_u32(&a) == ir_rand_u32(&b));
    CHECK(ir_rand_u32(&a) == ir_rand_u32(&b));
    CHECK(ir_rand_u32(&a) != ir_rand_u32(&c));
    ir_rand_seed(&a, 0u);   /* zero state must not lock up */
    CHECK(ir_rand_u32(&a) != 0u);

    /* Uniform range */
    ir_rand_seed(&a, 7u);
    {
        int ok = 1;
        double mn = 2.0, mx = -1.0;
        for (i = 0; i < 20000; ++i) {
            double u = ir_rand_uniform(&a);
            if (u < 0.0 || u >= 1.0) {
                ok = 0;
            }
            if (u < mn) { mn = u; }
            if (u > mx) { mx = u; }
        }
        CHECK(ok);
        CHECK(mn < 0.01);
        CHECK(mx > 0.99);
    }

    /* Gaussian: mean ~ 0, sd ~ 1 */
    ir_rand_seed(&a, 99u);
    {
        double *g = (double *)malloc(50000u * sizeof(double));
        CHECK(g != NULL);
        if (g != NULL) {
            double m, s;
            for (i = 0; i < 50000; ++i) {
                g[i] = ir_rand_gauss(&a);
            }
            s = ir_stddev(g, 50000u, &m);
            CHECK(fabs(m) < 0.03);
            CHECK(fabs(s - 1.0) < 0.03);
            free(g);
        }
    }

    /* Statistics on hand-checked data */
    CHECK_NEAR(ir_mean(x, 5u), 3.0, 1e-15);
    CHECK_NEAR(ir_stddev(x, 5u, &mean), sqrt(2.0), 1e-12);
    CHECK_NEAR(mean, 3.0, 1e-15);
    CHECK_NEAR(ir_median_copy(x, 5u), 3.0, 1e-15);
    {
        double q[4] = {4.0, 1.0, 3.0, 2.0};
        CHECK_NEAR(ir_median_copy(q, 4u), 2.5, 1e-15);
    }
    CHECK_NEAR(ir_percentile_copy(x, 5u, 0.0), 1.0, 1e-15);
    CHECK_NEAR(ir_percentile_copy(x, 5u, 50.0), 3.0, 1e-15);
    CHECK_NEAR(ir_percentile_copy(x, 5u, 100.0), 5.0, 1e-15);
    CHECK_NEAR(ir_min(x, 5u), 1.0, 1e-15);
    CHECK_NEAR(ir_max(x, 5u), 5.0, 1e-15);
    CHECK_NEAR(ir_rmse(x, y, 5u), sqrt((1.0 + 4.0 + 9.0 + 16.0 + 25.0) / 5.0), 1e-12);
    CHECK_NEAR(ir_mean_abs_diff(x, y, 5u), 3.0, 1e-15);

    /* Degenerate inputs must not crash or divide by zero */
    CHECK_NEAR(ir_mean(NULL, 5u), 0.0, 1e-15);
    CHECK_NEAR(ir_mean(x, 0u), 0.0, 1e-15);
    CHECK_NEAR(ir_stddev(NULL, 0u, NULL), 0.0, 1e-15);
    CHECK_NEAR(ir_median_copy(NULL, 0u), 0.0, 1e-15);
    CHECK_NEAR(ir_min(NULL, 0u), 0.0, 1e-15);

    /* NU definition: NU = 100*sigma/mean */
    nu = ir_nu_percent(x, 5u, NULL, NULL);
    CHECK_NEAR(nu, 100.0 * sqrt(2.0) / 3.0, 1e-10);
    CHECK_NEAR(nu, 47.140452, 1e-5);
    /* A constant array has zero non-uniformity. */
    {
        double k[64];
        for (i = 0; i < 64; ++i) {
            k[i] = 7.5;
        }
        CHECK_NEAR(ir_nu_percent(k, 64u, &mean, &sd), 0.0, 1e-15);
        CHECK_NEAR(mean, 7.5, 1e-15);
        CHECK_NEAR(sd, 0.0, 1e-15);
    }

    /* File helper: write then read back bytes (binary mode).
     * The temporary file lives in the current directory on purpose: the test
     * suite must not depend on any directory existing besides the one it is
     * run from. */
    {
        char text[32];
        FILE *f;
        size_t got;
        CHECK(ir_write_text("_test_write.tmp", "hello\x1a world") == IR_OK);
        f = fopen("_test_write.tmp", "rb");
        CHECK(f != NULL);
        got = 0u;
        if (f != NULL) {
            got = fread(text, 1u, sizeof(text) - 1u, f);
            text[got] = '\0';
            fclose(f);
        }
        /* The 0x1A byte must survive: this only holds because the writer used
         * "wb".  In text mode on Windows the stream would stop at 0x1A. */
        CHECK(got == 12u);
        CHECK(text[5] == (char)0x1A);
        remove("_test_write.tmp");
        CHECK(ir_write_text(NULL, "x") == IR_ERR_PARAM);
        CHECK(ir_write_text("no_such_dir_zz/x.txt", "x") == IR_ERR_IO);
    }
}

/* ================================================================== */
/* 2. tec_ctrl                                                         */
/* ================================================================== */

static void test_tec(void)
{
    tec_plant_cfg_t pl;
    tec_pid_cfg_t pc;
    tec_pid_t st;
    tec_metrics_t m_none, m_clamp, m_bc;
    const double dt = 0.02;
    const size_t n = 7500u;
    double u;
    int i;

    tec_plant_cfg_default(&pl);
    tec_pid_cfg_default(&pc);
    CHECK_NEAR(pl.rth * pl.cth, 8.0, 1e-12);   /* tau */
    CHECK(pl.u_min < 0.0);
    CHECK(pl.u_max > 0.0);
    CHECK(pl.ambient > 0.0);
    CHECK(pc.deadzone > 0.0);
    CHECK(pc.aw_mode == IR_TEC_AW_CLAMP);

    /* Sensor quantisation */
    CHECK_NEAR(tec_sensor_read(&pl, 1.234), 1.23, 1e-12);
    CHECK_NEAR(tec_sensor_read(&pl, -1.236), -1.24, 1e-12);
    pl.q_step = 0.0;
    CHECK_NEAR(tec_sensor_read(&pl, 1.234567), 1.234567, 1e-15);
    tec_plant_cfg_default(&pl);

    /* Reset clears the state */
    st.integ = 5.0;
    st.dfilt = 5.0;
    st.primed = 1;
    st.saturated = 1;
    tec_pid_reset(&st);
    CHECK_NEAR(st.integ, 0.0, 1e-15);
    CHECK_NEAR(st.dfilt, 0.0, 1e-15);
    CHECK(st.primed == 0);
    CHECK(st.saturated == 0);

    /* Output clamp */
    tec_pid_reset(&st);
    pc.kp = 100.0;
    pc.ki = 0.0;
    pc.deadzone = 0.0;
    u = tec_pid_step(&st, &pc, -20.0, 25.0, dt, -30.0, 30.0);
    CHECK_NEAR(u, -30.0, 1e-12);
    CHECK(st.last_u_raw < -30.0);
    u = tec_pid_step(&st, &pc, 100.0, 25.0, dt, -30.0, 30.0);
    CHECK_NEAR(u, 30.0, 1e-12);

    /* Dead zone: a small error contributes nothing to P or I */
    tec_pid_reset(&st);
    pc.kp = 5.0;
    pc.ki = 2.5;
    pc.deadzone = 0.02;
    u = tec_pid_step(&st, &pc, -20.0, -19.995, dt, -30.0, 30.0);  /* e = -5 mK */
    CHECK_NEAR(u, 0.0, 1e-15);
    CHECK_NEAR(st.integ, 0.0, 1e-15);
    u = tec_pid_step(&st, &pc, -20.0, -19.5, dt, -30.0, 30.0);    /* e = -0.5 K */
    CHECK(u < -2.0);

    /* Derivative on measurement: a ramp produces a negative correction */
    tec_pid_reset(&st);
    pc.kp = 0.0;
    pc.ki = 0.0;
    pc.kd = 1.0;
    pc.deadzone = 0.0;
    (void)tec_pid_step(&st, &pc, 0.0, 0.0, dt, -1e9, 1e9);
    (void)tec_pid_step(&st, &pc, 0.0, 1.0, dt, -1e9, 1e9);   /* +50 K/s */
    CHECK(st.dfilt > 0.0);
    u = tec_pid_step(&st, &pc, 0.0, 2.0, dt, -1e9, 1e9);
    CHECK(u < 0.0);
    /* Once the measurement stops moving, the filtered derivative decays. */
    {
        double d_prev = st.dfilt;
        (void)tec_pid_step(&st, &pc, 0.0, 2.0, dt, -1e9, 1e9);
        CHECK(st.dfilt < d_prev);
    }

    /* ---- anti-windup: integrator behaviour while pinned ---- */
    tec_pid_cfg_default(&pc);
    pc.deadzone = 0.0;
    tec_pid_reset(&st);
    pc.aw_mode = IR_TEC_AW_NONE;
    for (i = 0; i < 200; ++i) {
        u = tec_pid_step(&st, &pc, -20.0, 25.0, dt, -30.0, 30.0);
    }
    CHECK_NEAR(u, -30.0, 1e-12);
    CHECK(st.integ < -100.0);          /* charged far past what the plant needs */
    {
        double integ_none = st.integ;

        pc.aw_mode = IR_TEC_AW_CLAMP;
        tec_pid_reset(&st);
        for (i = 0; i < 200; ++i) {
            u = tec_pid_step(&st, &pc, -20.0, 25.0, dt, -30.0, 30.0);
        }
        CHECK_NEAR(u, -30.0, 1e-12);
        CHECK_NEAR(st.integ, 0.0, 1e-12);   /* never charged while pinned */
        CHECK(st.integ > integ_none);
        CHECK(st.saturated == 1);

        pc.aw_mode = IR_TEC_AW_BACKCALC;
        tec_pid_reset(&st);
        for (i = 0; i < 200; ++i) {
            u = tec_pid_step(&st, &pc, -20.0, 25.0, dt, -30.0, 30.0);
        }
        /* Back calculation drives the integrator to the value that makes the
         * *unclipped* output equal the rail, i.e. integ -> u_min here. */
        CHECK(st.integ > -35.0);
        CHECK(st.integ < -15.0);
        CHECK(st.integ > integ_none);
    }

    /* ---- closed loop step response ---- */
    tec_plant_cfg_default(&pl);
    tec_pid_cfg_default(&pc);
    pc.aw_mode = IR_TEC_AW_NONE;
    CHECK(tec_run_step(&pl, &pc, -20.0, 25.0, dt, n, NULL, NULL, NULL,
                       &m_none) == IR_OK);
    pc.aw_mode = IR_TEC_AW_CLAMP;
    CHECK(tec_run_step(&pl, &pc, -20.0, 25.0, dt, n, NULL, NULL, NULL,
                       &m_clamp) == IR_OK);
    pc.aw_mode = IR_TEC_AW_BACKCALC;
    CHECK(tec_run_step(&pl, &pc, -20.0, 25.0, dt, n, NULL, NULL, NULL,
                       &m_bc) == IR_OK);

    /* The recorded headline behaviour of this model, asserted so that a
     * regression in the controller cannot pass silently. */
    CHECK(m_none.overshoot_pct > 20.0);
    CHECK(m_clamp.overshoot_pct < 1.0);
    CHECK(m_clamp.overshoot_pct < m_none.overshoot_pct);
    CHECK(m_none.settle_time_s > 30.0);
    CHECK(m_clamp.settle_time_s < 20.0);
    CHECK(m_clamp.settle_time_s < m_none.settle_time_s);
    CHECK(m_none.sat_time_s > m_clamp.sat_time_s);
    CHECK(m_clamp.sat_time_s > 0.0);
    CHECK(m_bc.overshoot_pct > m_clamp.overshoot_pct);
    CHECK(m_bc.overshoot_pct < m_none.overshoot_pct);
    CHECK(fabs(m_clamp.ss_error_K) <= 0.03);
    CHECK(fabs(m_none.ss_error_K) <= 0.03);
    CHECK(fabs(m_bc.ss_error_K) <= 0.03);
    CHECK(m_clamp.rms_error_K < 0.05);
    CHECK(m_clamp.n_samples == n);
    CHECK(m_clamp.peak_K >= -21.0);      /* did not dip far past the setpoint */

    /* A heating step must work symmetrically. */
    {
        tec_metrics_t mh;
        pc.aw_mode = IR_TEC_AW_CLAMP;
        CHECK(tec_run_step(&pl, &pc, 40.0, 25.0, dt, 3000u, NULL, NULL, NULL,
                           &mh) == IR_OK);
        CHECK(mh.overshoot_pct < 5.0);
        CHECK(mh.peak_K <= 40.0 + 2.0);
        CHECK(fabs(mh.ss_error_K) <= 0.03);
        CHECK(mh.settle_time_s > 0.0);
    }

    /* An unreachable setpoint must be reported, not silently "settled". */
    {
        tec_metrics_t mu;
        tec_plant_cfg_t p2 = pl;
        p2.rth = 1.0;                    /* needs 45 W, only 30 W available */
        pc.aw_mode = IR_TEC_AW_CLAMP;
        CHECK(tec_run_step(&p2, &pc, -20.0, 25.0, dt, 2000u, NULL, NULL, NULL,
                           &mu) == IR_OK);
        CHECK(mu.settle_time_s < 0.0);
        CHECK(mu.ss_error_K > 5.0);
        CHECK(mu.overshoot_pct == 0.0);
    }

    /* Zero dead time is legal and must settle too. */
    {
        tec_metrics_t m0;
        tec_plant_cfg_t p2 = pl;
        p2.deadtime_s = 0.0;
        pc.aw_mode = IR_TEC_AW_CLAMP;
        CHECK(tec_run_step(&p2, &pc, -20.0, 25.0, dt, n, NULL, NULL, NULL,
                           &m0) == IR_OK);
        CHECK(m0.settle_time_s > 0.0);
        CHECK(fabs(m0.ss_error_K) <= 0.03);
    }

    /* Error paths */
    CHECK(tec_run_step(NULL, &pc, 0.0, 0.0, dt, 10u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
    CHECK(tec_run_step(&pl, NULL, 0.0, 0.0, dt, 10u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
    CHECK(tec_run_step(&pl, &pc, 0.0, 0.0, 0.0, 10u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
    CHECK(tec_run_step(&pl, &pc, 0.0, 0.0, dt, 0u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
    {
        tec_plant_cfg_t bad = pl;
        bad.rth = 0.0;
        CHECK(tec_run_step(&bad, &pc, 0.0, 0.0, dt, 10u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
        bad = pl;
        bad.u_min = 10.0;
        bad.u_max = -10.0;
        CHECK(tec_run_step(&bad, &pc, 0.0, 0.0, dt, 10u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
        bad = pl;
        bad.deadtime_s = -1.0;
        CHECK(tec_run_step(&bad, &pc, 0.0, 0.0, dt, 10u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
    }
    CHECK_NEAR(tec_pid_step(NULL, &pc, 0.0, 0.0, dt, -1.0, 1.0), 0.0, 1e-15);
    CHECK_NEAR(tec_pid_step(&st, NULL, 0.0, 0.0, dt, -1.0, 1.0), 0.0, 1e-15);
    tec_pid_reset(NULL);
    tec_plant_cfg_default(NULL);
    tec_pid_cfg_default(NULL);
    CHECK(1);   /* reaching here without a crash is itself a result */
}

/* ================================================================== */
/* 3. nuc                                                              */
/* ================================================================== */

static void test_nuc(void)
{
    nuc_array_cfg_t ac;
    nuc_truth_t tr;
    nuc_coeff_t k2, k1;
    ir_rand_t rng;
    size_t n, i;
    double *sl, *sh, *sm, *c2, *c1;
    const double phi_l = 100.0, phi_h = 300.0, phi_m = 200.0;

    nuc_array_cfg_default(&ac);
    CHECK(ac.w == 128u);
    CHECK(ac.h == 128u);
    CHECK(ac.gain_sigma > 0.0 && ac.gain_sigma < 1.0);
    CHECK(ac.noise_sigma > 0.0);
    n = ac.w * ac.h;

    sl = (double *)malloc(n * sizeof(double));
    sh = (double *)malloc(n * sizeof(double));
    sm = (double *)malloc(n * sizeof(double));
    c2 = (double *)malloc(n * sizeof(double));
    c1 = (double *)malloc(n * sizeof(double));
    CHECK(sl && sh && sm && c2 && c1);
    if (!sl || !sh || !sm || !c2 || !c1) {
        free(sl); free(sh); free(sm); free(c2); free(c1);
        return;
    }

    CHECK(nuc_truth_init(&tr, &ac, 20240501u) != IR_OK ? 0 : 1);
    CHECK(tr.w == ac.w && tr.h == ac.h);
    CHECK(nuc_coeff_init(&k2, ac.w, ac.h) == IR_OK);
    CHECK(nuc_coeff_init(&k1, ac.w, ac.h) == IR_OK);

    /* Generated pixel parameters obey the model's own clamps. */
    CHECK_NEAR(ir_mean(tr.gain, n), 1.0, 0.01);
    CHECK(ir_min(tr.gain, n) >= 0.5);
    CHECK(ir_max(tr.gain, n) <= 1.5);
    CHECK(ir_min(tr.curv, n) >= -0.3);
    CHECK(ir_max(tr.curv, n) <= 0.3);
    /* The gain spread must actually be there, otherwise the test is vacuous. */
    CHECK(ir_stddev(tr.gain, n, NULL) > 0.10);

    /* Rendering must implement S = O + G*phi*(1 + C*x) exactly. */
    CHECK(nuc_render(&tr, &ac, phi_m, 0, NULL, sm) == IR_OK);
    {
        ir_rand_seed(&rng, 1u);
        CHECK(nuc_render(&tr, &ac, phi_l, 0, NULL, sl) == IR_OK);
        CHECK(nuc_render(&tr, &ac, phi_h, 0, NULL, sh) == IR_OK);
        /* phi_l is half a reference span below phi_ref, so x = -0.5 */
        CHECK_NEAR(sl[0], tr.offset[0] + tr.gain[0] * phi_l *
                          (1.0 + tr.curv[0] * (-0.5)), 1e-9);
        CHECK_NEAR(sh[123u], tr.offset[123u] + tr.gain[123u] * phi_h *
                            (1.0 + tr.curv[123u] * 0.5), 1e-9);
        /* At phi_ref the curvature term vanishes by construction. */
        CHECK_NEAR(sm[7u], tr.offset[7u] + tr.gain[7u] * phi_m, 1e-9);
    }

    /* Two point calibration is exact for the linear model: with zero
     * curvature the correction recovers the flux at ANY scene level. */
    {
        nuc_array_cfg_t lin = ac;
        nuc_truth_t t2;
        nuc_coeff_t kk;
        double *a, *b, *s, *o;
        lin.curv_sigma = 0.0;
        lin.noise_sigma = 0.0;
        a = (double *)malloc(n * sizeof(double));
        b = (double *)malloc(n * sizeof(double));
        s = (double *)malloc(n * sizeof(double));
        o = (double *)malloc(n * sizeof(double));
        CHECK(a && b && s && o);
        if (a && b && s && o) {
            CHECK(nuc_truth_init(&t2, &lin, 5u) == IR_OK);
            CHECK(nuc_coeff_init(&kk, lin.w, lin.h) == IR_OK);
            CHECK(nuc_render(&t2, &lin, phi_l, 0, NULL, a) == IR_OK);
            CHECK(nuc_render(&t2, &lin, phi_h, 0, NULL, b) == IR_OK);
            CHECK(nuc_calib_two_point(&kk, a, b, phi_l, phi_h) == IR_OK);
            CHECK(kk.mode == NUC_MODE_TWO_POINT);
            CHECK_NEAR(kk.flux_low, phi_l, 1e-15);
            CHECK_NEAR(kk.flux_high, phi_h, 1e-15);
            {
                double worst = 0.0;
                int q;
                for (q = 0; q <= 10; ++q) {
                    double phi = 50.0 + 30.0 * (double)q;   /* 50 .. 350 */
                    CHECK(nuc_render(&t2, &lin, phi, 0, NULL, s) == IR_OK);
                    CHECK(nuc_apply(&kk, s, o) == IR_OK);
                    for (i = 0u; i < n; ++i) {
                        double e = fabs(o[i] - phi);
                        if (e > worst) {
                            worst = e;
                        }
                    }
                }
                /* Extrapolating 150 counts beyond each calibration point must
                 * still be exact for a purely linear pixel. */
                CHECK(worst < 1e-9);
            }
            nuc_truth_free(&t2);
            nuc_coeff_free(&kk);
        }
        free(a); free(b); free(s); free(o);
    }

    /* NU before / after on the real (non-linear) array. */
    {
        double nu_raw, nu_2pt, nu_1pt;
        ir_rand_seed(&rng, 909u);
        CHECK(nuc_render(&tr, &ac, phi_m, 1, &rng, sm) == IR_OK);
        ir_rand_seed(&rng, 4242u);
        CHECK(nuc_render(&tr, &ac, phi_l, 1, &rng, sl) == IR_OK);
        CHECK(nuc_render(&tr, &ac, phi_h, 1, &rng, sh) == IR_OK);
        CHECK(nuc_calib_two_point(&k2, sl, sh, phi_l, phi_h) == IR_OK);
        CHECK(nuc_calib_one_point(&k1, sl, phi_l) == IR_OK);
        CHECK(k1.mode == NUC_MODE_ONE_POINT);
        CHECK(nuc_apply(&k2, sm, c2) == IR_OK);
        CHECK(nuc_apply(&k1, sm, c1) == IR_OK);

        nu_raw = nuc_nu(sm, n, NULL, NULL);
        nu_2pt = nuc_nu(c2, n, NULL, NULL);
        nu_1pt = nuc_nu(c1, n, NULL, NULL);

        CHECK(nu_raw > 10.0);            /* the FPN must dominate before NUC */
        CHECK(nu_2pt < 1.0);             /* and be almost gone after it       */
        CHECK(nu_2pt < nu_raw / 10.0);
        CHECK(nu_1pt < nu_raw);
        CHECK(nu_1pt > nu_2pt * 5.0);    /* offset-only leaves the gain spread */

        /* Corrected mean must sit on the scene flux, not on an arbitrary
         * pedestal: this is what makes the output radiometrically usable. */
        CHECK_NEAR(ir_mean(c2, n), phi_m, 0.5);
        CHECK_NEAR(ir_mean(c1, n), phi_m, 5.0);

        /* Calibration gains must centre on 1.0 (the correction is a relative
         * rescaling, not a renormalisation to something else). */
        CHECK_NEAR(ir_mean(k2.gain, n), 1.0, 0.1);
        CHECK(ir_min(k2.gain, n) > 0.0);
    }

    /* Off-point NU for single point correction grows with distance from the
     * calibration flux; two point does not. */
    {
        double nu1_near, nu1_far, nu2_near, nu2_far;
        ir_rand_seed(&rng, 11u);
        CHECK(nuc_render(&tr, &ac, phi_l + 10.0, 1, &rng, sm) == IR_OK);
        CHECK(nuc_apply(&k1, sm, c1) == IR_OK);
        CHECK(nuc_apply(&k2, sm, c2) == IR_OK);
        nu1_near = nuc_nu(c1, n, NULL, NULL);
        nu2_near = nuc_nu(c2, n, NULL, NULL);
        ir_rand_seed(&rng, 11u);
        CHECK(nuc_render(&tr, &ac, phi_h, 1, &rng, sm) == IR_OK);
        CHECK(nuc_apply(&k1, sm, c1) == IR_OK);
        CHECK(nuc_apply(&k2, sm, c2) == IR_OK);
        nu1_far = nuc_nu(c1, n, NULL, NULL);
        nu2_far = nuc_nu(c2, n, NULL, NULL);
        CHECK(nu1_far > nu1_near);       /* single point degrades on the move */
        CHECK(nu2_far < nu2_near);       /* two point improves (noise scales) */
        CHECK(nu2_far < 0.5);
    }

    /* A pixel that does not respond at all must not divide by zero. */
    {
        nuc_coeff_t kd;
        double s_a[4] = {10.0, 10.0, 0.0, 20.0};
        double s_b[4] = {10.0, 30.0, 0.0, 40.0};
        double src[4] = {10.0, 30.0, 0.0, 40.0};
        double dst[4] = {0.0, 0.0, 0.0, 0.0};
        CHECK(nuc_coeff_init(&kd, 2u, 2u) == IR_OK);
        CHECK(nuc_calib_two_point(&kd, s_a, s_b, 100.0, 300.0) == IR_OK);
        CHECK_NEAR(kd.gain[0], 1.0, 1e-15);      /* flat pixel -> pass through */
        CHECK_NEAR(kd.gain[2], 1.0, 1e-15);      /* 0 -> 0 pixel   -> pass through */
        CHECK_NEAR(kd.gain[1], 200.0 / 20.0, 1e-9);
        CHECK_NEAR(kd.gain[3], 200.0 / 20.0, 1e-9);
        CHECK(nuc_apply(&kd, src, dst) == IR_OK);
        CHECK(isfinite(dst[0]) && isfinite(dst[2]));
        nuc_coeff_free(&kd);
    }

    /* Error paths */
    CHECK(nuc_calib_two_point(&k2, NULL, sh, phi_l, phi_h) == IR_ERR_PARAM);
    CHECK(nuc_calib_two_point(&k2, sl, sh, 100.0, 100.0) == IR_ERR_RANGE);
    CHECK(nuc_calib_one_point(NULL, sl, phi_l) == IR_ERR_PARAM);
    CHECK(nuc_apply(&k2, NULL, c2) == IR_ERR_PARAM);
    CHECK(nuc_apply(NULL, sm, c2) == IR_ERR_PARAM);
    CHECK(nuc_truth_init(NULL, &ac, 1u) == IR_ERR_PARAM);
    CHECK(nuc_coeff_init(&k1, 0u, 4u) == IR_ERR_PARAM);
    CHECK(nuc_render(&tr, &ac, 100.0, 1, NULL, sm) == IR_ERR_PARAM);
    {
        nuc_array_cfg_t bad = ac;
        nuc_truth_t t3;
        bad.w = 0u;
        CHECK(nuc_truth_init(&t3, &bad, 1u) == IR_ERR_PARAM);
    }
    nuc_truth_free(&tr);
    nuc_coeff_free(&k2);
    nuc_coeff_free(&k1);
    nuc_truth_free(NULL);
    nuc_coeff_free(NULL);
    free(sl); free(sh); free(sm); free(c2); free(c1);
}

/* ================================================================== */
/* 4. badpix                                                           */
/* ================================================================== */

static void test_badpix(void)
{
    badpix_cfg_t bc;
    badpix_replace_cfg_t rc;
    badpix_stats_t st;
    badpix_eval_t ev;
    badpix_residual_t rs;
    uint8_t flags[64];
    uint8_t injected[64];
    double resp[64];
    double out[64];
    size_t n_rep;
    int i;

    badpix_cfg_default(&bc);
    badpix_replace_cfg_default(&rc);
    CHECK_NEAR(bc.k_sigma, 3.0, 1e-15);
    CHECK(bc.iterations >= 1);
    CHECK(bc.max_flagged_frac > 0.0 && bc.max_flagged_frac <= 1.0);
    CHECK(rc.method == BADPIX_REPLACE_MEDIAN3);
    CHECK(rc.min_valid >= 1u);

    /* A single gross outlier in a uniform array must be found by one pass. */
    for (i = 0; i < 64; ++i) {
        double z = (double)(i % 8) * 0.01;      /* tiny deterministic texture */
        resp[i] = 100.0 + z;
    }
    resp[27] = 1000.0;
    CHECK(badpix_detect(resp, 64u, &bc, flags, &st) == IR_OK);
    CHECK(flags[27] == 1u);
    CHECK(st.flagged >= 1u);
    CHECK(st.accepted == 64u - st.flagged);
    CHECK(st.hi > st.lo);
    {
        int others = 0;
        for (i = 0; i < 64; ++i) {
            if (i != 27 && flags[i]) {
                ++others;
            }
        }
        CHECK(others == 0);
    }

    /* A larger threshold must flag no more than a smaller one. */
    {
        uint8_t f2[64];
        badpix_cfg_t loose = bc;
        size_t c1 = 0u, c2 = 0u;
        loose.k_sigma = 8.0;
        CHECK(badpix_detect(resp, 64u, &bc, flags, NULL) == IR_OK);
        CHECK(badpix_detect(resp, 64u, &loose, f2, NULL) == IR_OK);
        for (i = 0; i < 64; ++i) {
            c1 += flags[i] ? 1u : 0u;
            c2 += f2[i] ? 1u : 0u;
        }
        CHECK(c2 < c1);
        CHECK(c2 == 0u);
    }

    /* A constant array has zero spread: nothing may be flagged. */
    {
        double flat[16];
        for (i = 0; i < 16; ++i) {
            flat[i] = 42.0;
        }
        CHECK(badpix_detect(flat, 16u, &bc, flags, &st) == IR_OK);
        CHECK(st.flagged == 0u);
        CHECK_NEAR(st.sd, 0.0, 1e-15);
    }

    /* The clipping guard must refuse to reject more than max_flagged_frac.
     * A two point distribution with a fraction p of outliers at distance d
     * satisfies sigma = d*sqrt(p(1-p)), so the outliers fall outside 3 sigma
     * exactly when p < 0.1 -- which is why 5 % outliers is the interesting
     * case and 50 % (used in the first version of this test) is caught by
     * nothing at all. */
    {
        double skew[64];
        badpix_cfg_t g = bc;
        for (i = 0; i < 64; ++i) {
            skew[i] = 100.0;
        }
        skew[5] = skew[17] = skew[42] = 300.0;   /* 3/64 = 4.7 % outliers */
        g.max_flagged_frac = 0.02;               /* guard below the true rate */
        CHECK(badpix_detect(skew, 64u, &g, flags, &st) == IR_OK);
        CHECK(st.aborted == 1);
        CHECK(st.flagged >= 3u);
        /* With a realistic guard the same array is clipped cleanly. */
        g.max_flagged_frac = 0.20;
        CHECK(badpix_detect(skew, 64u, &g, flags, &st) == IR_OK);
        CHECK(st.flagged == 3u);
        CHECK(st.aborted == 0);
        CHECK(flags[5] == 1u && flags[17] == 1u && flags[42] == 1u);
        CHECK(flags[0] == 0u);
    }

    /* Replacement: a flagged interior pixel takes the median of its
     * neighbours from the *input* image, not from the output being built. */
    for (i = 0; i < 64; ++i) {
        resp[i] = 100.0 + (double)i;
    }
    memset(flags, 0, sizeof(flags));
    flags[27] = 1u;                 /* 3x3 neighbourhood is 18..44 minus 27 */
    CHECK(badpix_replace(resp, 8u, 8u, flags, &rc, out, &n_rep) == IR_OK);
    CHECK(n_rep == 1u);
    CHECK(out[27] > resp[26]);
    CHECK(out[27] < resp[28]);
    CHECK_NEAR(out[27], 100.0 + 27.0, 0.51);   /* neighbours straddle the value */
    for (i = 0; i < 64; ++i) {
        if (i != 27) {
            CHECK_NEAR(out[i], resp[i], 1e-15);
        }
    }

    /* Mean replacement gives a different but still sane answer. */
    {
        badpix_replace_cfg_t rm = rc;
        rm.method = BADPIX_REPLACE_MEAN3;
        CHECK(badpix_replace(resp, 8u, 8u, flags, &rm, out, &n_rep) == IR_OK);
        CHECK(n_rep == 1u);
        CHECK_NEAR(out[27], 127.0, 1e-12);
    }

    /* A corner pixel has only three valid neighbours: it must still be
     * replaced when min_valid <= 3, and left alone when min_valid is higher. */
    {
        uint8_t f[64];
        badpix_replace_cfg_t r3 = rc, r9 = rc;
        memset(f, 0, sizeof(f));
        f[0] = 1u;
        r3.min_valid = 3u;
        r9.min_valid = 9u;
        CHECK(badpix_replace(resp, 8u, 8u, f, &r3, out, &n_rep) == IR_OK);
        CHECK(n_rep == 1u);
        CHECK(out[0] != resp[0]);
        CHECK(badpix_replace(resp, 8u, 8u, f, &r9, out, &n_rep) == IR_OK);
        CHECK(n_rep == 0u);
        CHECK_NEAR(out[0], resp[0], 1e-15);
    }

    /* A blind pixel must never be used as a neighbour. */
    {
        uint8_t f[64];
        memset(f, 0, sizeof(f));
        f[27] = 1u;
        f[19] = 1u;                 /* one neighbour is also blind */
        resp[19] = -1e6;
        CHECK(badpix_replace(resp, 8u, 8u, f, &rc, out, &n_rep) == IR_OK);
        CHECK(n_rep == 2u);
        CHECK(out[27] > 0.0);
        CHECK(out[19] > 0.0);
    }

    /* evaluate(): counts must be self consistent. */
    memset(injected, 0, sizeof(injected));
    injected[27] = 1u;
    injected[19] = 1u;
    memset(flags, 0, sizeof(flags));
    flags[27] = 1u;
    badpix_eval(flags, injected, 64u, &ev);
    CHECK(ev.n_injected == 2u);
    CHECK(ev.detected == 1u);
    CHECK(ev.missed == 1u);
    CHECK(ev.false_alarm == 0u);
    CHECK(ev.detected + ev.missed == ev.n_injected);
    CHECK_NEAR(ev.detection_rate_pct, 50.0, 1e-12);
    CHECK_NEAR(ev.false_alarm_pct, 0.0, 1e-12);

    /* residual(): before/after against an ideal response. */
    {
        double a[4] = {0.0, 0.0, 10.0, 12.0};
        double b[4] = {9.0, 9.0, 11.0, 12.0};
        double id[4] = {10.0, 10.0, 10.0, 10.0};
        uint8_t mask[4] = {1u, 1u, 0u, 0u};
        badpix_residual(a, b, id, mask, 4u, &rs);
        CHECK(rs.n == 2u);
        CHECK_NEAR(rs.mae_before, 10.0, 1e-12);
        CHECK_NEAR(rs.mae_after, 1.0, 1e-12);
        CHECK_NEAR(rs.rmse_before, 10.0, 1e-12);
        CHECK_NEAR(rs.rmse_after, 1.0, 1e-12);
        CHECK_NEAR(rs.max_before, 10.0, 1e-12);
        badpix_residual(NULL, b, id, mask, 4u, &rs);
        CHECK(rs.n == 0u);
        badpix_eval(NULL, injected, 64u, &ev);
        CHECK(ev.n_injected == 0u);
    }

    /* Error paths */
    CHECK(badpix_detect(NULL, 64u, &bc, flags, NULL) == IR_ERR_PARAM);
    CHECK(badpix_detect(resp, 64u, NULL, flags, NULL) == IR_ERR_PARAM);
    CHECK(badpix_detect(resp, 64u, &bc, NULL, NULL) == IR_ERR_PARAM);
    CHECK(badpix_detect(resp, 0u, &bc, flags, NULL) == IR_ERR_PARAM);
    {
        badpix_cfg_t bad = bc;
        bad.k_sigma = 0.0;
        CHECK(badpix_detect(resp, 64u, &bad, flags, NULL) == IR_ERR_PARAM);
        bad = bc;
        bad.iterations = 0;
        CHECK(badpix_detect(resp, 64u, &bad, flags, NULL) == IR_ERR_PARAM);
    }
    CHECK(badpix_replace(NULL, 8u, 8u, flags, &rc, out, &n_rep) == IR_ERR_PARAM);
    CHECK(badpix_replace(resp, 0u, 8u, flags, &rc, out, &n_rep) == IR_ERR_PARAM);
    CHECK(badpix_replace(resp, 8u, 8u, flags, &rc, NULL, &n_rep) == IR_ERR_PARAM);
    /* NULL config must fall back to the defaults, not to garbage.  `flags`
     * currently holds exactly one flagged pixel (index 27). */
    CHECK(badpix_replace(resp, 8u, 8u, flags, NULL, out, &n_rep) == IR_OK);
    CHECK(n_rep == 1u);
    {
        double out2[64];
        size_t n2 = 0u;
        CHECK(badpix_replace(resp, 8u, 8u, flags, &rc, out2, &n2) == IR_OK);
        CHECK(n2 == n_rep);
        CHECK_NEAR(out[27], out2[27], 1e-15);
    }
    badpix_cfg_default(NULL);
    badpix_replace_cfg_default(NULL);
    CHECK(1);
}

/* ================================================================== */
/* 5. agc                                                              */
/* ================================================================== */

static void test_agc(void)
{
    agc_cfg_t gc;
    agc_analysis_t an;
    agc_metrics_t m;
    agc_scene_cfg_t scc;
    agc_scene_t sc;
    double img[256];
    double out[256];
    uint8_t out8[256];
    uint8_t target[256];
    uint8_t bg[256];
    double map[AGC_MAX_BINS];
    size_t n = 256u, i;
    int k;

    agc_cfg_default(&gc);
    CHECK(gc.bins == AGC_MAX_BINS);
    CHECK(gc.mode == AGC_MODE_PLATEAU);
    CHECK(gc.plateau_frac > 0.0);
    CHECK(gc.out_max > gc.out_min);

    /* A synthetic two level frame: 4 bright pixels, the rest at 50/52. */
    for (i = 0u; i < n; ++i) {
        img[i] = ((i % 2u) == 0u) ? 50.0 : 52.0;
        target[i] = 0u;
        bg[i] = 1u;
    }
    img[0] = img[1] = img[16] = img[17] = 100.0;
    target[0] = target[1] = target[16] = target[17] = 1u;
    bg[0] = bg[1] = bg[16] = bg[17] = 0u;

    CHECK(agc_analyze(&gc, img, n, &an) == IR_OK);
    CHECK(an.bins == AGC_MAX_BINS);
    CHECK(an.total == n);
    CHECK(an.hi > an.lo);
    {
        size_t sum = 0u;
        for (k = 0; k < an.bins; ++k) {
            sum += an.hist[k];
        }
        CHECK(sum == n);
        sum = 0u;
        for (k = 0; k < an.bins; ++k) {
            sum += an.hist_clipped[k];
        }
        CHECK(sum == n);       /* plateau redistribution conserves the count */
    }
    CHECK(an.lo >= 50.0 - 1e-9);
    CHECK(an.hi <= 100.0 + 1e-9);

    /* Linear map hits both rails. */
    gc.mode = AGC_MODE_LINEAR;
    CHECK(agc_build_map(&gc, &an, map) == IR_OK);
    CHECK_NEAR(map[0], gc.out_min, 1e-12);
    CHECK_NEAR(map[AGC_MAX_BINS - 1], gc.out_max, 1e-12);
    {
        int mono = 1;
        for (k = 1; k < AGC_MAX_BINS; ++k) {
            if (map[k] < map[k - 1]) {
                mono = 0;
            }
        }
        CHECK(mono);
    }

    /* HEQ and plateau maps are monotonic non-decreasing too. */
    for (k = 0; k < 2; ++k) {
        int mono = 1;
        gc.mode = (k == 0) ? AGC_MODE_HEQ : AGC_MODE_PLATEAU;
        CHECK(agc_build_map(&gc, &an, map) == IR_OK);
        for (i = 1u; i < (size_t)AGC_MAX_BINS; ++i) {
            if (map[i] < map[i - 1]) {
                mono = 0;
            }
        }
        CHECK(mono);
    }

    /* A huge plateau means "no clipping at all", so plateau must reduce to
     * plain HEQ bit for bit (the algorithm's asymptotic limit). */
    {
        double map_p[AGC_MAX_BINS], map_h[AGC_MAX_BINS];
        agc_analysis_t a2;
        agc_cfg_t g2 = gc;
        g2.mode = AGC_MODE_PLATEAU;
        g2.plateau_frac = 1.0e9;
        CHECK(agc_analyze(&g2, img, n, &a2) == IR_OK);
        CHECK(agc_build_map(&g2, &a2, map_p) == IR_OK);
        g2.mode = AGC_MODE_HEQ;
        CHECK(agc_analyze(&g2, img, n, &a2) == IR_OK);
        CHECK(agc_build_map(&g2, &a2, map_h) == IR_OK);
        {
            int same = 1;
            for (k = 0; k < AGC_MAX_BINS; ++k) {
                if (fabs(map_p[k] - map_h[k]) > 1e-9) {
                    same = 0;
                }
            }
            CHECK(same);
        }
    }

    /* A tiny plateau flattens the histogram, so the map approaches the
     * linear stretch (it cannot be exactly equal: the uniform redistribution
     * leaves a remainder of `excess % bins` counts on the lowest bins). */
    {
        double map_p[AGC_MAX_BINS], map_l[AGC_MAX_BINS];
        agc_analysis_t a2;
        agc_cfg_t g2 = gc;
        g2.mode = AGC_MODE_PLATEAU;
        g2.plateau_frac = 1.0e-9;
        CHECK(agc_analyze(&g2, img, n, &a2) == IR_OK);
        CHECK(agc_build_map(&g2, &a2, map_p) == IR_OK);
        g2.mode = AGC_MODE_LINEAR;
        CHECK(agc_analyze(&g2, img, n, &a2) == IR_OK);
        CHECK(agc_build_map(&g2, &a2, map_l) == IR_OK);
        CHECK_NEAR(map_p[0], map_l[0], 8.0);
        CHECK_NEAR(map_p[AGC_MAX_BINS - 1], map_l[AGC_MAX_BINS - 1], 1e-9);
        CHECK_NEAR(map_p[AGC_MAX_BINS / 2], map_l[AGC_MAX_BINS / 2], 8.0);
        /* And it must be very different from full HEQ. */
        g2.mode = AGC_MODE_HEQ;
        CHECK(agc_analyze(&g2, img, n, &a2) == IR_OK);
        CHECK(agc_build_map(&g2, &a2, map_l) == IR_OK);
        CHECK(fabs(map_l[AGC_MAX_BINS / 2] - map_p[AGC_MAX_BINS / 2]) > 10.0);
    }

    /* apply(): the 8 bit output must stay in range and keep the two levels
     * apart. */
    gc.mode = AGC_MODE_LINEAR;
    CHECK(agc_analyze(&gc, img, n, &an) == IR_OK);
    CHECK(agc_apply(&gc, &an, img, n, out, out8) == IR_OK);
    {
        int in_range = 1;
        for (i = 0u; i < n; ++i) {
            if (out[i] < gc.out_min - 1e-9 || out[i] > gc.out_max + 1e-9) {
                in_range = 0;
            }
        }
        CHECK(in_range);
        CHECK(out[0] > out[2]);                /* target above background */
        CHECK(out8[0] > out8[2]);
        CHECK(out8[0] <= 255u);
    }

    /* Metrics on the synthetic frame: Michelson = |100-51|/(151)*100,
     * Weber = 49/51*100 and CNR = 49/1 (background sd is exactly 1). */
    agc_metrics(img, n, 16u, 16u, target, bg, &m);
    CHECK_NEAR(m.mean, (4.0 * 100.0 + 252.0 * 51.0) / 256.0, 1e-9);
    CHECK_NEAR(m.tb_michelson_pct, 100.0 * 49.0 / 151.0, 1e-9);
    CHECK_NEAR(m.tb_weber_pct, 100.0 * 49.0 / 51.0, 1e-9);
    CHECK_NEAR(m.bg_sd, 1.0, 1e-12);
    CHECK_NEAR(m.cnr, 49.0, 1e-9);
    CHECK(m.entropy_bits > 0.0);
    CHECK(m.entropy_bits < 2.0);

    /* A constant frame: zero contrast, zero entropy, no NaN. */
    for (i = 0u; i < n; ++i) {
        img[i] = 5.0;
    }
    agc_metrics(img, n, 16u, 16u, NULL, NULL, &m);
    CHECK_NEAR(m.sd, 0.0, 1e-15);
    CHECK_NEAR(m.rms_contrast_pct, 0.0, 1e-15);
    CHECK_NEAR(m.entropy_bits, 0.0, 1e-15);
    CHECK_NEAR(m.tb_michelson_pct, 0.0, 1e-15);

    /* A constant frame must not break agc_analyze (the lo/hi fallback). */
    gc.mode = AGC_MODE_PLATEAU;
    CHECK(agc_analyze(&gc, img, n, &an) == IR_OK);
    CHECK(an.hi > an.lo);

    /* Synthetic scene sanity. */
    agc_scene_cfg_default(&scc);
    CHECK(scc.w == 128u && scc.h == 128u);
    CHECK(scc.target_size > 0);
    CHECK(agc_scene_synth(&sc, &scc, 2024u) == IR_OK);
    {
        size_t nt = 0u, nb = 0u;
        double mt = 0.0, mb = 0.0;
        for (i = 0u; i < sc.w * sc.h; ++i) {
            if (sc.target[i]) { mt += sc.img[i]; ++nt; }
            if (sc.bg[i])     { mb += sc.img[i]; ++nb; }
        }
        CHECK(nt == (size_t)(scc.target_size * scc.target_size));
        CHECK(nb > 0u);
        CHECK(mt / (double)nt > mb / (double)nb);
        agc_metrics(sc.img, sc.w * sc.h, sc.w, sc.h, sc.target, sc.bg, &m);
        CHECK(m.cnr > 5.0);                     /* the target is visible */
        CHECK(m.local_contrast > 0.0);
    }

    /* The whole point of the module: contrast must go UP on a real scene. */
    {
        double in_sd;
        agc_metrics(sc.img, sc.w * sc.h, sc.w, sc.h, sc.target, sc.bg, &m);
        in_sd = m.sd;
        gc.mode = AGC_MODE_PLATEAU;
        gc.plateau_frac = 1.0;
        CHECK(agc_analyze(&gc, sc.img, sc.w * sc.h, &an) == IR_OK);
        {
            double *o = (double *)malloc(sc.w * sc.h * sizeof(double));
            CHECK(o != NULL);
            if (o != NULL) {
                CHECK(agc_apply(&gc, &an, sc.img, sc.w * sc.h, o, NULL) == IR_OK);
                agc_metrics(o, sc.w * sc.h, sc.w, sc.h, sc.target, sc.bg, &m);
                CHECK(m.sd > in_sd * 3.0);
                CHECK(m.tb_michelson_pct > 30.0);
                free(o);
            }
        }
    }

    /* Error paths */
    CHECK(agc_analyze(NULL, img, n, &an) == IR_ERR_PARAM);
    CHECK(agc_analyze(&gc, NULL, n, &an) == IR_ERR_PARAM);
    CHECK(agc_analyze(&gc, img, 0u, &an) == IR_ERR_PARAM);
    {
        agc_cfg_t bad = gc;
        bad.bins = 1;
        CHECK(agc_analyze(&bad, img, n, &an) == IR_ERR_RANGE);
        bad.bins = AGC_MAX_BINS + 1;
        CHECK(agc_analyze(&bad, img, n, &an) == IR_ERR_RANGE);
    }
    CHECK(agc_build_map(NULL, &an, map) == IR_ERR_PARAM);
    CHECK(agc_apply(&gc, &an, img, n, NULL, NULL) == IR_ERR_PARAM);
    CHECK(agc_scene_synth(NULL, &scc, 1u) == IR_ERR_PARAM);
    {
        agc_scene_cfg_t bad = scc;
        bad.w = 0u;
        CHECK(agc_scene_synth(&sc, &bad, 1u) == IR_ERR_PARAM);
    }
    agc_scene_free(&sc);
    CHECK(sc.img == NULL);
    agc_scene_free(NULL);
    agc_cfg_default(NULL);
    agc_scene_cfg_default(NULL);
    CHECK(1);
}

/* ================================================================== */
/* 6. frame / protocol                                                 */
/* ================================================================== */

static void test_frame(void)
{
    frame_timing_cfg_t c;
    frame_timing_t st;
    frame_signal_t sig;
    frame_report_t rep;
    int i;

    /* --- geometry --- */
    frame_preset_ir128(&c);
    CHECK(c.h_active == 128);
    CHECK(c.v_active == 128);
    CHECK(c.clocks_per_pixel == 1);
    CHECK(frame_line_clocks(&c) == 160);
    CHECK(frame_total_clocks(&c) == 23040);

    frame_preset_bt656_625(&c);
    CHECK(c.h_active == 720);
    CHECK(c.h_blank == 144);
    CHECK(c.v_active == 576);
    CHECK(c.v_blank == 49);
    CHECK(c.clocks_per_pixel == 2);
    CHECK(frame_line_clocks(&c) == 1728);
    CHECK(frame_total_clocks(&c) == 1080000);
    /* 1080000 clocks at 27 MHz = exactly one 40 ms frame (25 Hz). */
    CHECK_NEAR((double)frame_total_clocks(&c) / 27.0e6, 0.040, 1e-12);
    CHECK(frame_total_clocks(NULL) == 0);
    CHECK(frame_line_clocks(NULL) == 0);

    /* --- state machine --- */
    frame_preset_ir128(&c);
    frame_timing_init(&st, &c);
    CHECK(st.total_clocks == 23040);
    CHECK(st.frames_done == 0);

    /* First line: 128 data clocks then 32 blanking clocks. */
    for (i = 0; i < 128; ++i) {
        int done = frame_timing_step(&st, &sig);
        CHECK(done == 0);
        if (i == 0) {
            CHECK(sig.x == 0 && sig.y == 0);
            CHECK(sig.de == 1);
            CHECK(sig.hsync == 0);
            CHECK(sig.vsync == 0);
            CHECK(sig.active_pixel == 1);
        }
    }
    CHECK(st.active_pixels == 128);
    CHECK(st.lines_done == 1);
    (void)frame_timing_step(&st, &sig);          /* first blanking clock */
    CHECK(sig.de == 0);
    CHECK(sig.hsync == 1);
    CHECK(sig.x == 128);
    CHECK(sig.active_pixel == 0);

    /* Walk the rest of the frame. */
    {
        int completions = 0;
        long clk;
        for (clk = 129; clk < 23040; ++clk) {
            completions += frame_timing_step(&st, &sig);
        }
        CHECK(completions == 1);
    }
    CHECK(st.frames_done == 1);
    CHECK(st.last_active_pixels == 128L * 128L);
    CHECK(st.last_lines == 144);
    CHECK(st.last_active_lines == 128);
    CHECK(st.frame_clocks_done == 0);            /* counters reset for the next */

    /* The next frame must be identical in structure. */
    {
        int completions = 0;
        long clk;
        long active = 0;
        for (clk = 0; clk < 23040; ++clk) {
            completions += frame_timing_step(&st, &sig);
            if (sig.active_pixel) {
                ++active;
            }
        }
        CHECK(completions == 1);
        CHECK(active == 128L * 128L);
        CHECK(st.frames_done == 2);
    }

    /* --- full frame emission and structural verification --- */
    frame_preset_ir128(&c);
    {
        uint8_t *img = (uint8_t *)malloc(16384u);
        uint8_t *stream = (uint8_t *)malloc(23040u);
        CHECK(img != NULL && stream != NULL);
        if (img != NULL && stream != NULL) {
            for (i = 0; i < 16384; ++i) {
                img[i] = (uint8_t)((i * 7) & 0xFF);
            }
            memset(stream, 0xEE, 23040u);
            CHECK(frame_emit(&c, img, 128, 128, stream, 23040u, &rep) == IR_OK);
            CHECK(rep.counts_ok == 1);
            CHECK(rep.dims_match == 1);
            CHECK(rep.total_clocks == 23040);
            CHECK(rep.clocks_emitted == 23040);
            CHECK(rep.active_pixel_clocks == 16384);
            CHECK(rep.blanking_clocks == 23040 - 16384);
            CHECK(rep.lines == 144);
            CHECK(rep.active_lines == 128);
            CHECK(rep.pixels_emitted == 16384);
            CHECK(rep.hsync_asserted == 32L * 144L);
            CHECK(rep.vsync_asserted == 160L * 16L);

            /* Pixel k sits at clock (k/128)*160 + (k%128). */
            CHECK(stream[0] == img[0]);
            CHECK(stream[127] == img[127]);
            CHECK(stream[160] == img[128]);
            CHECK(stream[127 * 160 + 127] == img[16383]);
            CHECK(stream[128] == 0u);          /* blanking level */
            CHECK(stream[23039] == 0u);        /* last clock is blanking */
        }
        free(img);
        free(stream);
    }

    /* A mismatched image is refused instead of producing garbage. */
    {
        uint8_t small[4] = {1u, 2u, 3u, 4u};
        CHECK(frame_emit(&c, small, 2, 2, NULL, 0u, &rep) == FRAME_ERR_DIM);
        CHECK(rep.counts_ok == 0);
    }
    /* Walking without an image still verifies the raster. */
    CHECK(frame_emit(&c, NULL, 0, 0, NULL, 0u, &rep) == IR_OK);
    CHECK(rep.counts_ok == 1);
    CHECK(rep.pixels_emitted == 0);
    CHECK(frame_emit(NULL, NULL, 0, 0, NULL, 0u, &rep) == IR_ERR_PARAM);
    {
        frame_timing_cfg_t bad = c;
        bad.clocks_per_pixel = 0;
        CHECK(frame_emit(&bad, NULL, 0, 0, NULL, 0u, &rep) == IR_ERR_PARAM);
        bad = c;
        bad.h_active = 0;
        CHECK(frame_emit(&bad, NULL, 0, 0, NULL, 0u, &rep) == IR_ERR_PARAM);
    }

    /* BT.656 reference raster walks to exactly 1 080 000 clocks. */
    frame_preset_bt656_625(&c);
    CHECK(frame_emit(&c, NULL, 0, 0, NULL, 0u, &rep) == IR_OK);
    CHECK(rep.clocks_emitted == 1080000L);
    CHECK(rep.active_pixel_clocks == 720L * 576L * 2L);
    CHECK(rep.lines == 625);
    CHECK(rep.active_lines == 576);
    CHECK(rep.counts_ok == 1);

    /* --- CRC --- */
    {
        const uint8_t tv[9] = {'1','2','3','4','5','6','7','8','9'};
        uint8_t one[1];
        uint8_t zeros[4];
        one[0] = 'A';
        zeros[0] = zeros[1] = zeros[2] = zeros[3] = 0u;
        CHECK(ir_crc16_ccitt(tv, 9u) == 0x29B1u);   /* the standard vector */
        CHECK(ir_crc16_ccitt(one, 1u) == 0xB915u);
        CHECK(ir_crc16_ccitt(zeros, 4u) == 0x84C0u);
        CHECK(ir_crc16_ccitt(NULL, 0u) == 0xFFFFu);  /* init value */
        CHECK(ir_crc16_ccitt(tv, 8u) != ir_crc16_ccitt(tv, 9u));
        /* Every single bit flip in a 4 byte message must change the CRC. */
        {
            uint8_t msg[4] = {0x12u, 0x34u, 0x56u, 0x78u};
            uint16_t base = ir_crc16_ccitt(msg, 4u);
            int bit, byte, distinct = 1;
            for (byte = 0; byte < 4; ++byte) {
                for (bit = 0; bit < 8; ++bit) {
                    uint8_t keep = msg[byte];
                    msg[byte] = (uint8_t)(keep ^ (uint8_t)(1u << bit));
                    if (ir_crc16_ccitt(msg, 4u) == base) {
                        distinct = 0;
                    }
                    msg[byte] = keep;
                }
            }
            CHECK(distinct);
        }
    }

    /* --- framing --- */
    {
        uint8_t payload[16];
        uint8_t pkt[64];
        uint8_t copy[64];
        size_t len = 0u, consumed = 0u;
        irproto_frame_t fr;

        for (i = 0; i < 16; ++i) {
            payload[i] = (uint8_t)(i * 11u + 1u);
        }
        CHECK(irproto_pack(0x42u, payload, 16u, pkt, sizeof(pkt), &len) == IR_OK);
        CHECK(len == 16u + IRPROTO_OVERHEAD);
        CHECK(pkt[0] == IRPROTO_SOF0);
        CHECK(pkt[1] == IRPROTO_SOF1);
        CHECK(pkt[2] == 0x42u);
        CHECK(pkt[3] == 16u);
        CHECK(pkt[4] == 0u);
        CHECK(irproto_parse(pkt, len, &fr, &consumed) == IR_OK);
        CHECK(fr.type == 0x42u);
        CHECK(fr.len == 16u);
        CHECK(fr.payload[0] == payload[0]);
        CHECK(fr.payload[15] == payload[15]);
        CHECK(memcmp(fr.payload, payload, 16u) == 0);
        CHECK(consumed == len);
        {
            uint16_t crc = ir_crc16_ccitt(&pkt[2], 19u);
            /* 5 header bytes + 16 payload bytes -> CRC at [21] and [22]. */
            CHECK(pkt[21] == (uint8_t)(crc & 0xFFu));
            CHECK(pkt[22] == (uint8_t)((crc >> 8) & 0xFFu));
            CHECK(len == 23u);
        }

        /* Zero length payload is legal. */
        CHECK(irproto_pack(0u, NULL, 0u, pkt, sizeof(pkt), &len) == IR_OK);
        CHECK(len == IRPROTO_OVERHEAD);
        CHECK(irproto_parse(pkt, len, &fr, &consumed) == IR_OK);
        CHECK(fr.len == 0u);
        CHECK(consumed == IRPROTO_OVERHEAD);

        /* Corruption detection: every single bit flip in the protected
         * region must be caught. */
        for (i = 0; i < 16; ++i) {
            payload[i] = (uint8_t)(i * 11u + 1u);
        }
        CHECK(irproto_pack(0x42u, payload, 16u, pkt, sizeof(pkt), &len) == IR_OK);
        memcpy(copy, pkt, len);
        {
            int all_caught = 1;
            size_t b;
            for (b = 2u; b < len; ++b) {
                int bit;
                for (bit = 0; bit < 8; ++bit) {
                    memcpy(pkt, copy, len);
                    pkt[b] ^= (uint8_t)(1u << bit);
                    if (irproto_parse(pkt, len, &fr, &consumed) == IR_OK) {
                        all_caught = 0;
                    }
                }
            }
            /* CRC-16 catches 100 % of single bit errors in messages of this
             * length; anything else would be a bug in the CRC, not bad luck. */
            CHECK(all_caught);
            memcpy(pkt, copy, len);
        }

        /* The magic bytes themselves are not protected (as in most framing
         * schemes) so a corrupted magic must be reported as a SOF error. */
        pkt[0] ^= 0xFFu;
        CHECK(irproto_parse(pkt, len, &fr, &consumed) == IRPROTO_ERR_SOF);
        memcpy(pkt, copy, len);

        /* Truncated frames are "not yet complete", not "corrupt". */
        CHECK(irproto_parse(pkt, 1u, &fr, &consumed) == IRPROTO_ERR_SHORT);
        CHECK(irproto_parse(pkt, 4u, &fr, &consumed) == IRPROTO_ERR_SHORT);
        CHECK(irproto_parse(pkt, len - 1u, &fr, &consumed) == IRPROTO_ERR_SHORT);
        /* A length field that promises more than is present. */
        pkt[3] = 0xFFu;
        pkt[4] = 0x00u;
        CHECK(irproto_parse(pkt, len, &fr, &consumed) == IRPROTO_ERR_SHORT);
        memcpy(pkt, copy, len);

        /* Packing into a buffer that is too small. */
        CHECK(irproto_pack(1u, payload, 16u, pkt, 10u, &len) == IRPROTO_ERR_CAP);
        CHECK(irproto_pack(1u, NULL, 4u, pkt, sizeof(pkt), &len) == IR_ERR_PARAM);
        CHECK(irproto_pack(1u, payload, 16u, NULL, 64u, &len) == IR_ERR_PARAM);
        CHECK(irproto_parse(NULL, 10u, &fr, &consumed) == IR_ERR_PARAM);
        CHECK(irproto_parse(pkt, 10u, NULL, &consumed) == IR_ERR_PARAM);
    }

    frame_timing_init(NULL, &c);
    frame_timing_step(NULL, &sig);
    frame_timing_init(&st, NULL);
    CHECK(st.total_clocks == 0);
    CHECK(1);
}

/* ================================================================== */
/* 7. vcs_pwr                                                          */
/* ================================================================== */

static void test_vcs(void)
{
    vcs_cfg_t c;
    vcs_sweep_cfg_t sc;
    vcs_sweep_metrics_t sm;
    vcs_pi_cfg_t pi;
    vcs_step_metrics_t stm;
    vcs_line_t li, le;
    double *io;
    size_t n = 2000u;

    vcs_cfg_default(&c);
    vcs_sweep_cfg_default(&sc);
    vcs_pi_cfg_default(&pi);
    CHECK(sc.t_hi > sc.t_lo);
    CHECK(sc.n_points >= 2);
    CHECK(pi.kp > 0.0 && pi.ki > 0.0);
    CHECK(pi.dt > 0.0);
    CHECK(c.tau_drv > 0.0);

    /* --- device model --- */
    CHECK_NEAR(vcs_ith(&c, c.t_ref), c.ith_ref, 1e-15);
    CHECK_NEAR(vcs_eta(&c, c.t_ref), c.eta_ref, 1e-15);
    /* T0 = 65 K means the threshold roughly doubles over 45 K. */
    CHECK_NEAR(vcs_ith(&c, 70.0) / vcs_ith(&c, 25.0), 2.0, 0.02);
    CHECK(vcs_ith(&c, 70.0) > vcs_ith(&c, 50.0));
    CHECK(vcs_ith(&c, 50.0) > vcs_ith(&c, 25.0));
    CHECK(vcs_eta(&c, 70.0) < vcs_eta(&c, 25.0));
    CHECK(vcs_eta(&c, 70.0) > 0.0);

    /* Below threshold the output is the spontaneous emission floor. */
    CHECK_NEAR(vcs_power(&c, 25.0, 0.0), c.p_floor, 1e-15);
    CHECK_NEAR(vcs_power(&c, 25.0, vcs_ith(&c, 25.0)), c.p_floor, 1e-15);
    CHECK_NEAR(vcs_power(&c, 25.0, vcs_ith(&c, 25.0) + 1.0e-3),
               c.p_floor + 1.0e-3 * vcs_eta(&c, 25.0), 1e-15);
    CHECK(vcs_power(&c, 25.0, 10.0e-3) > vcs_power(&c, 25.0, 5.0e-3));

    /* Units: 1 mW/mA of slope efficiency with a 5 mW target needs about
     * 6 mA at 25 degC, not 6 A and not 6 uA.  This assertion exists because
     * the first version of the defaults had eta = 1e-3 W/A by mistake. */
    CHECK_NEAR(vcs_required_current(&c, 25.0), 5.98e-3, 1e-4);
    CHECK_NEAR(vcs_required_current(&c, 70.0), 8.778e-3, 5e-4);
    CHECK(vcs_required_current(&c, 70.0) > vcs_required_current(&c, 25.0));
    /* And driving that current must actually hit the target power. */
    CHECK_NEAR(vcs_power(&c, 25.0, vcs_required_current(&c, 25.0)),
               c.p_target, 1e-9);
    CHECK_NEAR(vcs_power(&c, 70.0, vcs_required_current(&c, 70.0)),
               c.p_target, 1e-9);

    /* --- linear fit used by the feedforward --- */
    li = vcs_fit_ith(&c, 25.0, 70.0);
    le = vcs_fit_eta(&c, 25.0, 70.0);
    CHECK(li.b > 0.0);
    CHECK(le.b < 0.0);
    CHECK_NEAR(vcs_line_eval(&li, 47.5), 1.4142e-3, 0.05e-3);
    /* An exponential cannot be captured by a line: the fit must be imperfect,
     * and that imperfection is exactly the feedforward residual. */
    CHECK(fabs(vcs_line_eval(&li, 47.5) - vcs_ith(&c, 47.5)) > 1e-5);
    CHECK(fabs(vcs_line_eval(&le, 47.5) - vcs_eta(&c, 47.5)) > 1e-4);

    /* --- the three compensation strategies over 25..70 degC --- */
    CHECK(vcs_sweep(&c, &sc, VCS_MODE_FIXED, NULL, NULL, NULL, 0u, &sm) == IR_OK);
    {
        double fixed_max = sm.max_dev_pct;
        double fixed_rms = sm.rms_dev_pct;
        CHECK(fixed_max > 30.0);          /* a fixed bias really does drift */
        CHECK(fixed_max < 60.0);
        CHECK(fixed_rms > 10.0);
        CHECK(sm.p_min < 0.7 * c.p_target);   /* power collapses at 70 C */
        CHECK_NEAR(sm.i_min, sm.i_max, 1e-15); /* the current never changes */
    }
    CHECK(vcs_sweep(&c, &sc, VCS_MODE_FF, NULL, NULL, NULL, 0u, &sm) == IR_OK);
    {
        double ff_max = sm.max_dev_pct;
        CHECK(ff_max < 5.0);
        CHECK(ff_max < 2.5);
        CHECK(ff_max > 0.1);              /* a linear fit is not perfect */
        CHECK(sm.i_max > sm.i_min);       /* the setpoint follows temperature */
        CHECK_NEAR(sm.i_min, 5.9e-3, 0.2e-3);
        CHECK_NEAR(sm.i_max, 8.7e-3, 0.2e-3);
    }
    CHECK(vcs_sweep(&c, &sc, VCS_MODE_APC, NULL, NULL, NULL, 0u, &sm) == IR_OK);
    {
        double apc_max = sm.max_dev_pct;
        CHECK(apc_max < 0.1);
        CHECK(apc_max < 0.02 * 5.0);
        CHECK(sm.max_dev_pct < 0.1);
    }

    /* The headline comparison, stated as an assertion. */
    {
        double d_fixed, d_ff, d_apc;
        CHECK(vcs_sweep(&c, &sc, VCS_MODE_FIXED, NULL, NULL, NULL, 0u, &sm) == IR_OK);
        d_fixed = sm.max_dev_pct;
        CHECK(vcs_sweep(&c, &sc, VCS_MODE_FF, NULL, NULL, NULL, 0u, &sm) == IR_OK);
        d_ff = sm.max_dev_pct;
        CHECK(vcs_sweep(&c, &sc, VCS_MODE_APC, NULL, NULL, NULL, 0u, &sm) == IR_OK);
        d_apc = sm.max_dev_pct;
        CHECK(d_ff < d_fixed / 10.0);
        CHECK(d_apc < d_ff / 10.0);
    }

    /* A driver that cannot reach the required current must be reported. */
    {
        vcs_cfg_t lim = c;
        vcs_sweep_metrics_t s2;
        lim.i_max = 7.0e-3;               /* cannot reach 8.78 mA at 70 C */
        CHECK(vcs_sweep(&lim, &sc, VCS_MODE_FIXED, NULL, NULL, NULL, 0u, &s2) == IR_OK);
        CHECK(s2.n_shortfall == 0);       /* 5.98 mA still fits */
        CHECK(vcs_sweep(&lim, &sc, VCS_MODE_FF, NULL, NULL, NULL, 0u, &s2) == IR_OK);
        CHECK(s2.n_shortfall > 0);
        CHECK_NEAR(s2.i_max, 7.0e-3, 1e-12);
    }

    /* --- inner constant current loop --- */
    io = (double *)malloc(n * sizeof(double));
    CHECK(io != NULL);
    if (io != NULL) {
        size_t k;
        double mx = 0.0;
        CHECK(vcs_current_step(&c, &pi, 6.0e-3, n, NULL, io, NULL, &stm) == IR_OK);
        CHECK_NEAR(io[n - 1u], 6.0e-3, 1e-5);
        CHECK(fabs(stm.ss_error_A) < 2.0e-6);
        CHECK(stm.overshoot_pct < 5.0);
        CHECK(stm.settle_time_s > 0.0);
        CHECK(stm.settle_time_s < 200.0e-6);
        for (k = 0u; k < n; ++k) {
            if (io[k] > mx) {
                mx = io[k];
            }
            CHECK(io[k] >= 0.0);
        }
        CHECK(mx <= pi.i_limit + 1e-12);
        /* The current must never exceed the driver limit even for a crazy
         * setpoint: the clamp is on the command, and the plant is first order
         * so it cannot overshoot its own command. */
        CHECK(vcs_current_step(&c, &pi, 100.0e-3, n, NULL, io, NULL, &stm) == IR_OK);
        CHECK(io[n - 1u] <= pi.i_limit + 1e-12);
        free(io);
    }

    /* --- APC with a noisy monitor photodiode --- */
    {
        double rms = 0.0, mx = 0.0;
        CHECK(vcs_apc_noise_study(&c, 70.0, 20000u, 10.0e-6, 555u,
                                  &rms, &mx) == IR_OK);
        CHECK(rms > 0.0);
        CHECK(rms < 1.0);        /* still an order of magnitude better than FF */
        CHECK(mx < 2.0);
        /* Deterministic for a given seed. */
        {
            double rms2 = 0.0, mx2 = 0.0;
            CHECK(vcs_apc_noise_study(&c, 70.0, 20000u, 10.0e-6, 555u,
                                      &rms2, &mx2) == IR_OK);
            CHECK_NEAR(rms, rms2, 1e-15);
            CHECK_NEAR(mx, mx2, 1e-15);
        }
    }

    /* Error paths */
    CHECK(vcs_sweep(NULL, &sc, VCS_MODE_FF, NULL, NULL, NULL, 0u, &sm) == IR_ERR_PARAM);
    CHECK(vcs_sweep(&c, NULL, VCS_MODE_FF, NULL, NULL, NULL, 0u, &sm) == IR_ERR_PARAM);
    {
        vcs_sweep_cfg_t bad = sc;
        bad.n_points = 1;
        CHECK(vcs_sweep(&c, &bad, 0, NULL, NULL, NULL, 0u, &sm) == IR_ERR_PARAM);
        bad = sc;
        bad.t_hi = bad.t_lo;
        CHECK(vcs_sweep(&c, &bad, 0, NULL, NULL, NULL, 0u, &sm) == IR_ERR_PARAM);
    }
    CHECK(vcs_current_step(NULL, &pi, 1e-3, 10u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
    CHECK(vcs_current_step(&c, NULL, 1e-3, 10u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
    {
        vcs_pi_cfg_t bad = pi;
        bad.dt = 0.0;
        CHECK(vcs_current_step(&c, &bad, 1e-3, 10u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
    }
    CHECK(vcs_current_step(&c, &pi, 1e-3, 0u, NULL, NULL, NULL, NULL) == IR_ERR_PARAM);
    CHECK(vcs_apc_noise_study(NULL, 50.0, 10u, 1e-5, 1u, NULL, NULL) == IR_ERR_PARAM);
    CHECK(vcs_apc_noise_study(&c, 50.0, 0u, 1e-5, 1u, NULL, NULL) == IR_ERR_PARAM);
    CHECK(vcs_apc_noise_study(&c, 50.0, 10u, 0.0, 1u, NULL, NULL) == IR_ERR_PARAM);
    CHECK_NEAR(vcs_ith(NULL, 25.0), 0.0, 1e-15);
    CHECK_NEAR(vcs_eta(NULL, 25.0), 0.0, 1e-15);
    CHECK_NEAR(vcs_power(NULL, 25.0, 1.0), 0.0, 1e-15);
    CHECK_NEAR(vcs_line_eval(NULL, 25.0), 0.0, 1e-15);
    vcs_cfg_default(NULL);
    vcs_pi_cfg_default(NULL);
    vcs_sweep_cfg_default(NULL);
    CHECK(1);
}

/* ================================================================== */

int main(void)
{
    printf("ir-core-lab self test\n");
    printf("---------------------\n");

    test_ir_lab();
    test_tec();
    test_nuc();
    test_badpix();
    test_agc();
    test_frame();
    test_vcs();

    printf("---------------------\n");
    if (g_fails != 0) {
        printf("%d of %d checks FAILED\n", g_fails, g_checks);
        return 1;
    }
    printf("%d checks passed\n", g_checks);
    return 0;
}
