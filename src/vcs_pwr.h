/*
 * vcs_pwr.h -- VCSEL constant current driver with temperature drift compensation
 * -----------------------------------------------------------------------------
 * DEVICE MODEL
 * ------------
 * A VCSEL is a diode: below the threshold current Ith it emits only
 * spontaneous light, above it the optical power rises almost linearly with a
 * slope efficiency eta.  Both parameters move with junction temperature, and
 * Ith moves fast (it is the reason a fixed bias current cannot hold a
 * constant optical power):
 *
 *      Ith(T) = Ith_ref * exp((T - T_ref) / T0)          T0 = characteristic
 *                                                        temperature, ~60-70 K
 *      eta(T) = eta_ref * (1 + a*(T-T_ref) + b*(T-T_ref)^2)
 *      P(T,I) = eta(T) * (I - Ith(T))     for I > Ith, ~0 otherwise
 *
 * T0 = 65 K means the threshold doubles over 45 K, which is the right order of
 * magnitude for a 850-1000 nm oxide confined VCSEL.  It is a textbook model,
 * not a measurement of any part.
 *
 * CONTROL PROBLEM
 * ---------------
 * Three strategies are compared over a 25 -> 70 degC junction sweep:
 *
 *   VCS_MODE_FIXED  I_set is a constant, calibrated once at 25 degC.
 *                   The optical power collapses as the threshold climbs.
 *   VCS_MODE_FF     Feedforward: firmware stores a *linear* fit of Ith(T) and
 *                   eta(T) (the usual "two point calibration + linear
 *                   interpolation" that fits in a small MCU) and computes
 *                       I_set(T) = Ith_fit(T) + P_target / eta_fit(T)
 *                   The residual error is exactly the curvature that the
 *                   linear fit cannot represent.
 *   VCS_MODE_APC    Closed loop average power control: a PI controller on the
 *                   monitor photodiode (MPD) reading trims I_set until the
 *                   measured power equals the target.  Residual error is set
 *                   by the MPD quantisation and noise, not by temperature.
 *
 * The inner loop is a PI constant current regulator around a first order
 * driver model with time constant tau_drv, a current limit and a current sense
 * quantiser.  Its step response (overshoot / settling / static error) is
 * measured by vcs_current_step().
 */

#ifndef VCS_PWR_H
#define VCS_PWR_H

#include <stddef.h>

#include "ir_lab.h"

#ifdef __cplusplus
extern "C" {
#endif

#define VCS_MODE_FIXED 0
#define VCS_MODE_FF    1
#define VCS_MODE_APC   2

typedef struct {
    double t_ref;       /* degC  reference temperature                 */
    double ith_ref;     /* A     threshold current at t_ref            */
    double t0_char;     /* K     characteristic temperature            */
    double eta_ref;     /* W/A   slope efficiency at t_ref             */
    double eta_a;       /* 1/K   linear slope efficiency tempco        */
    double eta_b;       /* 1/K^2 quadratic slope efficiency tempco     */
    double p_target;    /* W     target optical power                  */
    double i_max;       /* A     driver compliance limit               */
    double tau_drv;     /* s     driver current loop time constant     */
    double i_quant;     /* A     current sense / DAC resolution        */
    double mpd_resp;    /* A/W   monitor photodiode responsivity       */
    double mpd_lsb;     /* A     MPD TIA + ADC resolution              */
    double mpd_noise;   /* A     MPD reading noise, RMS                */
    double p_floor;     /* W     spontaneous emission floor below Ith  */
} vcs_cfg_t;

typedef struct {
    double kp;          /* A/A       proportional gain of the current loop */
    double ki;          /* A/(A*s)   integral gain                        */
    double dt;          /* s         control period                       */
    double i_limit;     /* A         output clamp                         */
} vcs_pi_cfg_t;

typedef struct {
    double overshoot_pct;
    double settle_time_s;   /* last time outside a +-1 % band, -1 if never */
    double ss_error_A;      /* mean error over the last 20 % of the run    */
    double sat_time_s;
} vcs_step_metrics_t;

typedef struct {
    double t_lo, t_hi;
    int    n_points;
} vcs_sweep_cfg_t;

typedef struct {
    double max_dev_pct;     /* max |P - P_target| / P_target * 100        */
    double rms_dev_pct;
    double p_min, p_max;
    double i_min, i_max;
    int    n_shortfall;     /* points where I_set hit the driver limit     */
} vcs_sweep_metrics_t;

void vcs_cfg_default(vcs_cfg_t *c);
void vcs_pi_cfg_default(vcs_pi_cfg_t *c);
void vcs_sweep_cfg_default(vcs_sweep_cfg_t *c);

/* Device model. */
double vcs_ith(const vcs_cfg_t *c, double t);
double vcs_eta(const vcs_cfg_t *c, double t);
double vcs_power(const vcs_cfg_t *c, double t, double i);
/* Exact current needed for p_target at temperature t (0 if unreachable). */
double vcs_required_current(const vcs_cfg_t *c, double t);

/* Least squares linear model  y = a + b*T  over [t_lo, t_hi]. */
typedef struct { double a, b; } vcs_line_t;
vcs_line_t vcs_fit_ith(const vcs_cfg_t *c, double t_lo, double t_hi);
vcs_line_t vcs_fit_eta(const vcs_cfg_t *c, double t_lo, double t_hi);
double     vcs_line_eval(const vcs_line_t *l, double t);

/* Feedforward setpoint from the linear fits. */
double vcs_ff_current(const vcs_cfg_t *c, const vcs_line_t *ith,
                      const vcs_line_t *eta, double t);

/*
 * PI constant current loop, plant  I' = (d - I)/tau_drv  with
 * d = clip(kp*e + integ, 0, i_limit), e = i_set - quantise(I).
 */
int vcs_current_step(const vcs_cfg_t *c, const vcs_pi_cfg_t *pi,
                     double i_set, size_t n,
                     double *t_out, double *i_out, double *d_out,
                     vcs_step_metrics_t *m);

/*
 * Temperature sweep over the three strategies.  Arrays are optional (cap
 * gives their capacity); mode selects VCS_MODE_*.
 */
int vcs_sweep(const vcs_cfg_t *c, const vcs_sweep_cfg_t *sc, int mode,
              double *t_out, double *p_out, double *i_out, size_t cap,
              vcs_sweep_metrics_t *m);

/* Converged current of the average power control loop at temperature t. */
double vcs_apc_current(const vcs_cfg_t *c, double t, double i_start);

/*
 * APC loop *with* monitor photodiode noise, which is what limits it in
 * practice (the noiseless loop only stops at the MPD quantiser, so its
 * residual is an unrealistically clean 0.008 %).  Runs n_iter outer loop
 * iterations of dt seconds at temperature t and reports the RMS and peak
 * deviation of the *actual* optical power over the second half of the run.
 */
int vcs_apc_noise_study(const vcs_cfg_t *c, double t, size_t n_iter, double dt,
                        uint32_t seed, double *p_dev_rms_pct,
                        double *p_dev_max_pct);

#ifdef __cplusplus
}
#endif

#endif /* VCS_PWR_H */
