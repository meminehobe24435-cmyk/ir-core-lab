/* tec_ctrl.c -- TEC temperature loop, see tec_ctrl.h for the model. */

#include "tec_ctrl.h"
#include "ir_lab.h"

#include <math.h>
#include <stdlib.h>

#define TEC_SETTLE_BAND_K 0.1   /* +-0.1 K settling band used by the metrics */

void tec_plant_cfg_default(tec_plant_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->rth        = 2.0;    /* K/W  */
    c->cth        = 4.0;    /* J/K  -> tau = 8 s */
    c->ambient    = 25.0;   /* degC */
    c->deadtime_s = 0.4;    /* s    */
    c->q_step     = 0.01;   /* K    */
    c->u_min      = -30.0;  /* W    */
    c->u_max      = 30.0;   /* W    */
}

void tec_pid_cfg_default(tec_pid_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    /*
     * Kp comes from the IMC rule for a FOPDT process with lambda = theta, so
     * that the closed loop is as fast as the dead time allows:
     *      Kp = tau / (K * (lambda + theta)) = 8 / (2 * 0.8) = 5.0 W/K
     *      Ti = tau                          = 8 s
     * The IMC integral value (Ki = Kp/Ti = 0.625 W/(K*s)) turns out to be too
     * slow *for this step*: the loop leaves saturation while the error is
     * still ~6 K (30 W / 5 W/K) and then crawls to the setpoint with the open
     * loop time constant of 8 s.  Ki was therefore picked from a gain sweep
     * (see results/tec_gain_sweep.csv, produced by `make sim`): 2.5 W/(K*s),
     * i.e. Ti = 2 s, removes the tail (settling 39.0 s -> 14.0 s) while
     * keeping the overshoot under 0.1 % *because* of the anti-windup.
     *
     * Derivative is left at zero on purpose: the thermistor is quantised to
     * 10 mK, and differentiating a quantised signal injects a limit cycle
     * larger than the +-0.1 K band we want to demonstrate (the sweep shows
     * this: Ki = 6.0 W/(K*s) already develops a 57 mK peak-to-peak limit
     * cycle).  The derivative path is still implemented and is exercised by
     * the unit tests with an analytic (noise free) input.
     */
    c->kp       = 5.0;
    c->ki       = 2.5;
    c->kd       = 0.0;
    c->nd       = 10.0;
    c->deadzone = 0.02;     /* K */
    c->aw_mode  = IR_TEC_AW_CLAMP;
    c->tt       = 2.0;      /* s, back calculation tracking constant (= Ti) */
}

void tec_pid_reset(tec_pid_t *st)
{
    if (st == NULL) {
        return;
    }
    st->integ      = 0.0;
    st->dfilt      = 0.0;
    st->prev_meas  = 0.0;
    st->last_u     = 0.0;
    st->last_u_raw = 0.0;
    st->primed     = 0;
    st->saturated  = 0;
}

double tec_sensor_read(const tec_plant_cfg_t *c, double temp)
{
    if (c == NULL || c->q_step <= 0.0) {
        return temp;
    }
    return floor(temp / c->q_step + 0.5) * c->q_step;
}

double tec_pid_step(tec_pid_t *st, const tec_pid_cfg_t *cfg,
                    double setpoint, double meas, double dt,
                    double u_min, double u_max)
{
    double e, e_ctl, dmeas, tf, alpha, u_raw, u;

    if (st == NULL || cfg == NULL || dt <= 0.0) {
        return 0.0;
    }

    /* --- error and dead zone ------------------------------------- */
    e = setpoint - meas;
    e_ctl = (fabs(e) < cfg->deadzone) ? 0.0 : e;

    /* --- filtered derivative on the measurement ------------------- */
    if (!st->primed) {
        st->prev_meas = meas;
        st->dfilt     = 0.0;
        st->primed    = 1;
    } else {
        dmeas = (meas - st->prev_meas) / dt;
        st->prev_meas = meas;
        if (cfg->kd != 0.0 && cfg->kp > 0.0 && cfg->nd > 0.0) {
            tf = (cfg->kd / cfg->kp) / cfg->nd;   /* filter time constant */
            alpha = (tf > 0.0) ? (tf / (tf + dt)) : 0.0;
        } else {
            alpha = 0.0;
        }
        st->dfilt = alpha * st->dfilt + (1.0 - alpha) * dmeas;
    }

    /* --- positional PID ------------------------------------------- */
    u_raw = cfg->kp * e_ctl + st->integ - cfg->kd * st->dfilt;

    if (u_raw > u_max) {
        u = u_max;
        st->saturated = 1;
    } else if (u_raw < u_min) {
        u = u_min;
        st->saturated = 1;
    } else {
        u = u_raw;
        st->saturated = 0;
    }

    /* --- anti-windup ---------------------------------------------- */
    switch (cfg->aw_mode) {
    case IR_TEC_AW_NONE:
        /* Naive: the integrator keeps charging while the actuator is
         * pinned, so it has to be discharged through the opposite error
         * before the loop reacts again. */
        st->integ += cfg->ki * e_ctl * dt;
        break;

    case IR_TEC_AW_BACKCALC:
        if (cfg->tt > 0.0) {
            st->integ += cfg->ki * e_ctl * dt + (u - u_raw) * dt / cfg->tt;
        } else {
            st->integ += cfg->ki * e_ctl * dt;
        }
        break;

    case IR_TEC_AW_CLAMP:
    default:
        /* Conditional integration: only accumulate when doing so does not
         * push the (already clipped) output further into the rail. */
        if (!((u_raw > u_max && e_ctl > 0.0) || (u_raw < u_min && e_ctl < 0.0))) {
            st->integ += cfg->ki * e_ctl * dt;
        }
        break;
    }

    st->last_u     = u;
    st->last_u_raw = u_raw;
    return u;
}

int tec_run_step(const tec_plant_cfg_t *plant, const tec_pid_cfg_t *pid,
                 double setpoint, double t0_temp, double dt, size_t n,
                 double *t_out, double *temp_out, double *u_out,
                 tec_metrics_t *m)
{
    double tau, tss, decay, u, u_del, sat_time;
    double *ring = NULL;
    double *tbuf = NULL;
    size_t i, nd, ridx;
    tec_pid_t st;

    if (plant == NULL || pid == NULL || dt <= 0.0 || n == 0u) {
        return IR_ERR_PARAM;
    }
    if (plant->rth <= 0.0 || plant->cth <= 0.0) {
        return IR_ERR_PARAM;
    }
    if (plant->deadtime_s < 0.0 || plant->u_min >= plant->u_max) {
        return IR_ERR_PARAM;
    }

    tau   = plant->rth * plant->cth;
    nd    = (size_t)(plant->deadtime_s / dt + 0.5);   /* delay in samples */
    decay = exp(-dt / tau);

    ring = (double *)calloc(nd + 1u, sizeof(double));
    tbuf = (double *)malloc(n * sizeof(double));
    if (ring == NULL || tbuf == NULL) {
        free(ring);
        free(tbuf);
        return IR_ERR_IO;
    }

    tec_pid_reset(&st);
    ridx = 0u;
    sat_time = 0.0;

    for (i = 0u; i < n; ++i) {
        double temp = (i == 0u) ? t0_temp : tbuf[i - 1u];
        double meas = tec_sensor_read(plant, temp);

        u = tec_pid_step(&st, pid, setpoint, meas, dt,
                         plant->u_min, plant->u_max);
        if (st.saturated) {
            sat_time += dt;
        }

        /* push the new drive, read the delayed one (nd samples ago) */
        ring[ridx] = u;
        ridx = (ridx + 1u) % (nd + 1u);
        u_del = ring[ridx];

        /* exact zero order hold update of the first order plant */
        tss = plant->ambient + plant->rth * u_del;
        tbuf[i] = tss + (temp - tss) * decay;

        if (t_out != NULL) {
            t_out[i] = (double)i * dt;
        }
        if (u_out != NULL) {
            u_out[i] = u;
        }
    }

    if (temp_out != NULL) {
        for (i = 0u; i < n; ++i) {
            temp_out[i] = tbuf[i];
        }
    }

    if (m != NULL) {
        double step = setpoint - t0_temp;
        double ext_min = tbuf[0], ext_max = tbuf[0];
        double t_last_out = 0.0;
        double ss_sum = 0.0, rms_sum = 0.0, mn, mx;
        size_t ss_n = 0u, rms_from, last_out_idx = 0u;
        int has_out = 0;

        for (i = 0u; i < n; ++i) {
            double T = tbuf[i];
            if (T < ext_min) {
                ext_min = T;
            }
            if (T > ext_max) {
                ext_max = T;
            }
            if (fabs(T - setpoint) > TEC_SETTLE_BAND_K) {
                t_last_out = (double)(i + 1u) * dt;
                last_out_idx = i;
                has_out = 1;
            }
        }

        ss_n = n / 10u;
        if (ss_n == 0u) {
            ss_n = 1u;
        }
        for (i = n - ss_n; i < n; ++i) {
            ss_sum += tbuf[i] - setpoint;
        }

        rms_from = n - (n / 5u);
        if (rms_from >= n) {
            rms_from = 0u;
        }
        mn = tbuf[rms_from];
        mx = tbuf[rms_from];
        for (i = rms_from; i < n; ++i) {
            double d = tbuf[i] - setpoint;
            rms_sum += d * d;
            if (tbuf[i] < mn) {
                mn = tbuf[i];
            }
            if (tbuf[i] > mx) {
                mx = tbuf[i];
            }
        }

        m->n_samples    = n;
        m->ss_error_K   = ss_sum / (double)ss_n;
        m->rms_error_K  = sqrt(rms_sum / (double)(n - rms_from));
        m->bandwidth_K  = mx - mn;

        /*
         * Overshoot has to be measured in the direction the setpoint moved:
         * this plant is asked to *cool* from 25 degC to -20 degC, so the
         * interesting extreme is the minimum, not the maximum.  Percentages
         * are referenced to the size of the step.
         */
        if (step < 0.0) {
            double beyond = setpoint - ext_min;        /* undershoot, K */
            m->peak_K = ext_min;
            m->overshoot_pct = (beyond > 0.0) ? 100.0 * beyond / fabs(step) : 0.0;
        } else {
            double beyond = ext_max - setpoint;        /* overshoot, K */
            m->peak_K = ext_max;
            m->overshoot_pct = (beyond > 0.0) ? 100.0 * beyond / fabs(step) : 0.0;
        }

        /* settling time: last instant the response was outside the band */
        if (!has_out) {
            m->settle_time_s = 0.0;
        } else if (last_out_idx == n - 1u) {
            m->settle_time_s = -1.0;   /* never settled inside the window */
        } else {
            m->settle_time_s = t_last_out;
        }

        /* saturation time was accumulated while stepping the controller */
        m->sat_time_s = sat_time;
    }

    free(ring);
    free(tbuf);
    return IR_OK;
}
