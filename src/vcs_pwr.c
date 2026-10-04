/* vcs_pwr.c -- see vcs_pwr.h for the device and control model. */

#include "vcs_pwr.h"

#include <math.h>

#define VCS_SETTLE_BAND_FRAC 0.01   /* +-1 % of the setpoint */

void vcs_cfg_default(vcs_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->t_ref   = 25.0;      /* degC */
    c->ith_ref = 1.0e-3;    /* A    */
    c->t0_char = 65.0;      /* K    -> Ith doubles in 45 K */
    /*
     * Slope efficiency: 1 mW of optical power per mA of drive current is
     * 1 mW/mA = 1 W/A.  (It is easy to write 1e-3 here by reflex and get a
     * setpoint three orders of magnitude too large -- this typo was caught by
     * the "required current at 25 degC" assertion in test/test_ir.c.)
     */
    c->eta_ref = 1.0;       /* W/A  -> 1 mW per mA */
    c->eta_a   = -5.0e-3;   /* 1/K  -> -0.5 %/K */
    c->eta_b   = -2.0e-5;   /* 1/K^2 */
    c->p_target = 5.0e-3;   /* W    */
    c->i_max    = 15.0e-3;  /* A    */
    c->tau_drv  = 20.0e-6;  /* s    20 us driver settling */
    c->i_quant  = 3.66e-6;  /* A    15 mA / 4096 DAC LSB */
    c->mpd_resp = 0.5;      /* A/W  */
    c->mpd_lsb  = 1.2e-6;   /* A    12 bit TIA/ADC LSB */
    c->mpd_noise = 7.5e-6;  /* A    0.3 % of the ~2.5 mA MPD current */
    c->p_floor  = 20.0e-6;  /* W    spontaneous emission below threshold */
}

void vcs_pi_cfg_default(vcs_pi_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    /*
     * IMC PI for a first order plant 1/(tau*s+1) with closed loop time
     * constant lambda = tau/2:
     *      Kp = tau / (lambda * 1) = 2
     *      Ki = Kp / tau          = 1e5 A/(A*s)
     */
    c->kp = 2.0;
    c->ki = 1.0e5;
    c->dt = 0.5e-6;
    c->i_limit = 15.0e-3;
}

void vcs_sweep_cfg_default(vcs_sweep_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->t_lo = 25.0;
    c->t_hi = 70.0;
    c->n_points = 46;   /* 1 degC steps */
}

double vcs_ith(const vcs_cfg_t *c, double t)
{
    if (c == NULL) {
        return 0.0;
    }
    return c->ith_ref * exp((t - c->t_ref) / c->t0_char);
}

double vcs_eta(const vcs_cfg_t *c, double t)
{
    double u;
    if (c == NULL) {
        return 0.0;
    }
    u = t - c->t_ref;
    return c->eta_ref * (1.0 + c->eta_a * u + c->eta_b * u * u);
}

double vcs_power(const vcs_cfg_t *c, double t, double i)
{
    double ith, eta;
    if (c == NULL) {
        return 0.0;
    }
    ith = vcs_ith(c, t);
    if (i <= ith) {
        return c->p_floor;
    }
    eta = vcs_eta(c, t);
    if (eta < 0.0) {
        return c->p_floor;
    }
    return c->p_floor + eta * (i - ith);
}

double vcs_required_current(const vcs_cfg_t *c, double t)
{
    double eta;
    if (c == NULL) {
        return 0.0;
    }
    eta = vcs_eta(c, t);
    if (eta <= 0.0) {
        return c->i_max + 1.0;   /* unreachable */
    }
    /* The spontaneous emission floor p_floor is subtracted because the model
     * reports P = p_floor + eta*(I - Ith) above threshold. */
    return vcs_ith(c, t) + (c->p_target - c->p_floor) / eta;
}

double vcs_line_eval(const vcs_line_t *l, double t)
{
    if (l == NULL) {
        return 0.0;
    }
    return l->a + l->b * t;
}

/* Least squares fit of f over n_points evenly spaced samples in [t_lo, t_hi]. */
static vcs_line_t vcs_fit(double (*f)(const vcs_cfg_t *, double),
                          const vcs_cfg_t *c, double t_lo, double t_hi)
{
    vcs_line_t l;
    const int N = 64;
    int i;
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    double denom;

    l.a = 0.0;
    l.b = 0.0;
    if (c == NULL || t_hi <= t_lo) {
        return l;
    }
    for (i = 0; i < N; ++i) {
        double t = t_lo + (t_hi - t_lo) * (double)i / (double)(N - 1);
        double y = f(c, t);
        sx += t;
        sy += y;
        sxx += t * t;
        sxy += t * y;
    }
    denom = (double)N * sxx - sx * sx;
    if (fabs(denom) < 1e-12) {
        l.a = sy / (double)N;
        l.b = 0.0;
        return l;
    }
    l.b = ((double)N * sxy - sx * sy) / denom;
    l.a = (sy - l.b * sx) / (double)N;
    return l;
}

static double ith_wrapper(const vcs_cfg_t *c, double t)
{
    return vcs_ith(c, t);
}

static double eta_wrapper(const vcs_cfg_t *c, double t)
{
    return vcs_eta(c, t);
}

vcs_line_t vcs_fit_ith(const vcs_cfg_t *c, double t_lo, double t_hi)
{
    return vcs_fit(ith_wrapper, c, t_lo, t_hi);
}

vcs_line_t vcs_fit_eta(const vcs_cfg_t *c, double t_lo, double t_hi)
{
    return vcs_fit(eta_wrapper, c, t_lo, t_hi);
}

double vcs_ff_current(const vcs_cfg_t *c, const vcs_line_t *ith,
                      const vcs_line_t *eta, double t)
{
    double e;
    if (c == NULL || ith == NULL || eta == NULL) {
        return 0.0;
    }
    e = vcs_line_eval(eta, t);
    if (e <= 1e-9) {
        return c->i_max;
    }
    return vcs_line_eval(ith, t) + (c->p_target - c->p_floor) / e;
}

int vcs_current_step(const vcs_cfg_t *c, const vcs_pi_cfg_t *pi,
                     double i_set, size_t n,
                     double *t_out, double *i_out, double *d_out,
                     vcs_step_metrics_t *m)
{
    double i_act = 0.0, integ = 0.0, dt, sat_time = 0.0;
    size_t k, i;

    if (c == NULL || pi == NULL || n == 0u || pi->dt <= 0.0) {
        return IR_ERR_PARAM;
    }
    dt = pi->dt;

    for (k = 0u; k < n; ++k) {
        double meas, err, d_raw, d;
        /* Quantised current sense. */
        if (c->i_quant > 0.0) {
            meas = floor(i_act / c->i_quant + 0.5) * c->i_quant;
        } else {
            meas = i_act;
        }
        err = i_set - meas;
        d_raw = pi->kp * err + integ;
        if (d_raw > pi->i_limit) {
            d = pi->i_limit;
        } else if (d_raw < 0.0) {
            d = 0.0;
        } else {
            d = d_raw;
        }
        /* Conditional integration: do not charge the integrator while the
         * driver is already pinned at a rail and the error pushes further. */
        if (!((d_raw > pi->i_limit && err > 0.0) || (d_raw < 0.0 && err < 0.0))) {
            integ += pi->ki * err * dt;
        }
        if (fabs(d - d_raw) > 1e-15) {
            sat_time += dt;
        }

        /* First order driver: exact zero order hold solution. */
        i_act = d + (i_act - d) * exp(-dt / c->tau_drv);

        if (t_out != NULL) {
            t_out[k] = (double)k * dt;
        }
        if (i_out != NULL) {
            i_out[k] = i_act;
        }
        if (d_out != NULL) {
            d_out[k] = d;
        }
    }

    if (m != NULL) {
        double band = VCS_SETTLE_BAND_FRAC * fabs(i_set);
        double peak = 0.0, ss_sum = 0.0, t_last_out = 0.0;
        size_t ss_n = n / 5u, last_out_idx = 0u;
        int has_out = 0;

        if (ss_n == 0u) {
            ss_n = 1u;
        }
        if (i_out != NULL) {
            for (i = 0u; i < n; ++i) {
                double d = i_out[i];
                if (d > peak) {
                    peak = d;
                }
                if (fabs(d - i_set) > band) {
                    t_last_out = (double)(i + 1u) * dt;
                    last_out_idx = i;
                    has_out = 1;
                }
            }
            for (i = n - ss_n; i < n; ++i) {
                ss_sum += i_out[i] - i_set;
            }
        } else {
            /* Without a trace we can only fall back to zeroed metrics. */
            for (i = 0u; i < n; ++i) {
                (void)i;
            }
        }

        m->overshoot_pct = (i_set > 0.0 && peak > i_set)
                               ? 100.0 * (peak - i_set) / i_set
                               : 0.0;
        m->ss_error_A = ss_sum / (double)ss_n;
        if (!has_out) {
            m->settle_time_s = 0.0;
        } else if (last_out_idx == n - 1u) {
            m->settle_time_s = -1.0;
        } else {
            m->settle_time_s = t_last_out;
        }
        m->sat_time_s = sat_time;
    }
    return IR_OK;
}

double vcs_apc_current(const vcs_cfg_t *c, double t, double i_start)
{
    const int N = 4000;
    const double dt_out = 10.0e-6;   /* outer loop period, 10 us */
    /*
     * Integrator gain in A per (W*s).  The closed outer loop is first order
     * with tau = 1/(eta*gain): with eta = 1 W/A and gain = 1000 A/(W*s) that
     * is 1 ms, i.e. 100x the 10 us update period.  An earlier value of 1e6
     * gave tau = 1 us, *shorter* than the update period, and the loop
     * diverged into a 180 % limit cycle -- the classic "sampling too slowly
     * for the gain you asked for" mistake.
     */
    const double gain = 1000.0;
    double i_cmd, integ;

    if (c == NULL) {
        return 0.0;
    }
    i_cmd = i_start;
    integ = i_cmd;

    /*
     * Outer average power control loop.  It is an integrator on the power
     * error: dI/dt = gain * (P_target - P_measured), giving a first order
     * closed loop with tau = 1/(eta*gain) ~ 1 ms for eta = 1 mW/mA.  The
     * inner constant current loop (tau_drv = 20 us) is treated as settled
     * because it is 50x faster; keeping the two loops separated by more than
     * a decade is exactly the design rule a real driver would follow.
     */
    {
        int k;
        for (k = 0; k < N; ++k) {
            double p_meas = vcs_power(c, t, i_cmd);
            double i_raw;

            /* Monitor photodiode chain: quantiser only (its noise is handled
             * as an RMS figure separately, this helper returns the converged
             * setpoint). */
            if (c->mpd_resp > 0.0 && c->mpd_lsb > 0.0) {
                double a = c->mpd_resp * p_meas;
                a = floor(a / c->mpd_lsb + 0.5) * c->mpd_lsb;
                p_meas = a / c->mpd_resp;
            }

            integ += gain * (c->p_target - p_meas) * dt_out;
            i_raw = integ;
            if (i_raw > c->i_max) {
                i_raw = c->i_max;
            }
            if (i_raw < 0.0) {
                i_raw = 0.0;
            }
            i_cmd = i_raw;
        }
    }
    return i_cmd;
}

int vcs_apc_noise_study(const vcs_cfg_t *c, double t, size_t n_iter, double dt,
                        uint32_t seed, double *p_dev_rms_pct,
                        double *p_dev_max_pct)
{
    const double gain = 1000.0;
    ir_rand_t rng;
    double i_cmd, integ, sum_sq = 0.0, mx = 0.0;
    size_t k, n_eval = 0u;

    if (c == NULL || n_iter == 0u || dt <= 0.0) {
        return IR_ERR_PARAM;
    }
    i_cmd = vcs_required_current(c, t);
    if (i_cmd > c->i_max) {
        i_cmd = c->i_max;
    }
    integ = i_cmd;
    ir_rand_seed(&rng, seed);

    for (k = 0u; k < n_iter; ++k) {
        double p_true = vcs_power(c, t, i_cmd);
        double p_meas = p_true;
        double i_raw, dev;

        /* Monitor chain: quantiser plus RMS noise, referred back to watts. */
        if (c->mpd_resp > 0.0) {
            double a = c->mpd_resp * p_true;
            if (c->mpd_noise > 0.0) {
                a += c->mpd_noise * ir_rand_gauss(&rng);
            }
            if (c->mpd_lsb > 0.0) {
                a = floor(a / c->mpd_lsb + 0.5) * c->mpd_lsb;
            }
            p_meas = a / c->mpd_resp;
        }

        integ += gain * (c->p_target - p_meas) * dt;
        i_raw = integ;
        if (i_raw > c->i_max) {
            i_raw = c->i_max;
        }
        if (i_raw < 0.0) {
            i_raw = 0.0;
        }
        i_cmd = i_raw;

        if (k >= n_iter / 2u) {
            dev = 100.0 * (p_true - c->p_target) / c->p_target;
            sum_sq += dev * dev;
            if (fabs(dev) > mx) {
                mx = fabs(dev);
            }
            ++n_eval;
        }
    }

    if (p_dev_rms_pct != NULL) {
        *p_dev_rms_pct = (n_eval > 0u) ? sqrt(sum_sq / (double)n_eval) : 0.0;
    }
    if (p_dev_max_pct != NULL) {
        *p_dev_max_pct = mx;
    }
    return IR_OK;
}

int vcs_sweep(const vcs_cfg_t *c, const vcs_sweep_cfg_t *sc, int mode,
              double *t_out, double *p_out, double *i_out, size_t cap,
              vcs_sweep_metrics_t *m)
{
    vcs_line_t ith_fit, eta_fit;
    int k;

    if (c == NULL || sc == NULL || sc->n_points < 2 || sc->t_hi <= sc->t_lo) {
        return IR_ERR_PARAM;
    }

    ith_fit = vcs_fit_ith(c, sc->t_lo, sc->t_hi);
    eta_fit = vcs_fit_eta(c, sc->t_lo, sc->t_hi);

    if (m != NULL) {
        m->max_dev_pct = 0.0;
        m->rms_dev_pct = 0.0;
        m->p_min = 1e30;
        m->p_max = -1e30;
        m->i_min = 1e30;
        m->i_max = -1e30;
        m->n_shortfall = 0;
    }

    {
        double sum_sq = 0.0;
        double i_prev = vcs_required_current(c, sc->t_lo);

        for (k = 0; k < sc->n_points; ++k) {
            double t = sc->t_lo +
                       (sc->t_hi - sc->t_lo) * (double)k / (double)(sc->n_points - 1);
            double iset, p, dev;

            switch (mode) {
            case VCS_MODE_FIXED:
                /* Calibrated once at t_ref and never updated. */
                iset = vcs_required_current(c, c->t_ref);
                break;
            case VCS_MODE_FF:
                iset = vcs_ff_current(c, &ith_fit, &eta_fit, t);
                break;
            case VCS_MODE_APC:
            default:
                /* Warm start from the previous temperature point, as a real
                 * APC loop would do. */
                iset = vcs_apc_current(c, t, i_prev);
                i_prev = iset;
                break;
            }

            if (iset > c->i_max) {
                iset = c->i_max;
                if (m != NULL) {
                    m->n_shortfall++;
                }
            }
            if (iset < 0.0) {
                iset = 0.0;
            }

            p = vcs_power(c, t, iset);
            dev = 100.0 * fabs(p - c->p_target) / c->p_target;
            sum_sq += dev * dev;

            if (m != NULL) {
                if (dev > m->max_dev_pct) {
                    m->max_dev_pct = dev;
                }
                if (p < m->p_min) {
                    m->p_min = p;
                }
                if (p > m->p_max) {
                    m->p_max = p;
                }
                if (iset < m->i_min) {
                    m->i_min = iset;
                }
                if (iset > m->i_max) {
                    m->i_max = iset;
                }
            }
            if (t_out != NULL && (size_t)k < cap) {
                t_out[k] = t;
            }
            if (p_out != NULL && (size_t)k < cap) {
                p_out[k] = p;
            }
            if (i_out != NULL && (size_t)k < cap) {
                i_out[k] = iset;
            }
        }

        if (m != NULL) {
            m->rms_dev_pct = sqrt(sum_sq / (double)sc->n_points);
        }
    }
    return IR_OK;
}
