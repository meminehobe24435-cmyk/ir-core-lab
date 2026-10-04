/*
 * tec_ctrl.h -- thermoelectric cooler (TEC) temperature control loop
 * ---------------------------------------------------------------
 * A digital PID controller driving a first order plus dead time (FOPDT)
 * model of the detector cold platform.
 *
 * PLANT MODEL
 * -----------
 * The cold platform is modelled as a single thermal capacitance Cth [J/K]
 * losing heat to the ambient through a thermal resistance Rth [K/W]:
 *
 *      Cth * dT/dt = (Tamb - T)/Rth + u(t - theta)
 *
 * with u the TEC drive power in watts (negative = cooling), Tamb the ambient
 * temperature and theta a pure transport delay that lumps together the TEC
 * thermal mass, the thermistor lag and the ADC/control pipeline.
 *
 * The open loop is therefore first order:
 *
 *      tau = Rth * Cth      [s]        DC gain K = Rth  [K/W]
 *
 * The dead time is implemented as a ring buffer of controller outputs, and
 * the continuous part is integrated with the *exact* zero order hold solution
 *
 *      T[k+1] = Tss + (T[k] - Tss) * exp(-dt/tau),   Tss = Tamb + Rth*u_delayed
 *
 * which avoids Euler discretisation error entirely.  With the numbers used
 * here (dt = 20 ms, tau = 8 s) the per-step Euler error would be ~1e-4 K per
 * step, i.e. the exact form keeps the "steady state error" numbers honest.
 *
 * CONTROLLER
 * ----------
 * Positional (non-incremental) PID, derivative acting on the measurement so
 * that a setpoint change does not produce a derivative kick:
 *
 *      e[k]     = SP - y[k]
 *      df[k]    = filtered rate of change of y      (first order, Tf = Td/N)
 *      u_raw    = Kp*e + I - Kd*df
 *      u        = clip(u_raw, u_min, u_max)
 *
 * Anti-windup, selected by aw_mode:
 *   IR_TEC_AW_NONE      I += Ki*e*dt                       (naive, windup prone)
 *   IR_TEC_AW_CLAMP     integration is skipped while the output is saturated
 *                       *and* the error would push it further into saturation
 *                       (conditional integration / clamping)
 *   IR_TEC_AW_BACKCALC  I += Ki*e*dt + (u - u_raw)*dt/Tt   (back calculation)
 *
 * A dead zone (neutral zone) of +-deadzone kelvin zeroes the P and I terms so
 * that measurement quantisation noise around the setpoint does not make the
 * integrator wander; the price is a bounded residual steady state error of at
 * most one dead zone, which is reported by the metrics below.
 */

#ifndef TEC_CTRL_H
#define TEC_CTRL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define IR_TEC_AW_NONE      0
#define IR_TEC_AW_CLAMP     1
#define IR_TEC_AW_BACKCALC  2

typedef struct {
    double rth;         /* K/W   thermal resistance platform -> ambient   */
    double cth;         /* J/K   thermal capacitance of the platform      */
    double ambient;     /* degC  ambient / heat sink temperature          */
    double deadtime_s;  /* s     pure transport delay theta (>= 0)        */
    double q_step;      /* K     thermistor + ADC resolution, <=0 disables*/
    double u_min;       /* W     most negative drive (cooling)            */
    double u_max;       /* W     most positive drive (heating)            */
} tec_plant_cfg_t;

typedef struct {
    double kp;          /* W/K       proportional gain                    */
    double ki;          /* W/(K*s)   integral gain                        */
    double kd;          /* W*s/K     derivative gain                      */
    double nd;          /* -         derivative filter divisor (Tf = Td/nd)*/
    double deadzone;    /* K         neutral zone half width              */
    int    aw_mode;     /* IR_TEC_AW_*                                    */
    double tt;          /* s         back calculation tracking time const */
} tec_pid_cfg_t;

typedef struct {
    double integ;       /* W     integral accumulator                     */
    double dfilt;       /* K/s   filtered derivative of the measurement   */
    double prev_meas;   /* degC  previous measurement                     */
    double last_u;      /* W     last controller output                   */
    double last_u_raw;  /* W     last controller output before clipping   */
    int    primed;      /* 0 until the first sample has been stored       */
    int    saturated;   /* 1 if the last output was clipped               */
} tec_pid_t;

typedef struct {
    double overshoot_pct;  /* peak excursion beyond SP / |step| * 100      */
    double settle_time_s;  /* last time outside +-0.1 K, -1 if never       */
    double ss_error_K;     /* mean error over the last 10 % of the run     */
    double peak_K;         /* extreme reached in the direction of the step */
    double rms_error_K;    /* RMS error over the last 20 % of the run      */
    double sat_time_s;     /* time the actuator spent clipped              */
    double bandwidth_K;   /* observed pk-pk of the last 20 % (limit cycle)*/
    size_t n_samples;
} tec_metrics_t;

/* Library defaults: single stage TEC on a small cold platform.
 *   tau = Rth*Cth = 2.0 K/W * 4.0 J/K = 8 s, theta = 0.4 s
 *   +-30 W of drive is a typical small single stage TEC budget, which is
 *   enough to hold the platform ~45 K below ambient in steady state. */
void tec_plant_cfg_default(tec_plant_cfg_t *c);
void tec_pid_cfg_default(tec_pid_cfg_t *c);
void tec_pid_reset(tec_pid_t *st);

/* Quantised sensor reading (identity when q_step <= 0). */
double tec_sensor_read(const tec_plant_cfg_t *c, double temp);

/* One controller step.  Returns the clipped drive power in watts. */
double tec_pid_step(tec_pid_t *st, const tec_pid_cfg_t *cfg,
                    double setpoint, double meas, double dt,
                    double u_min, double u_max);

/*
 * Closed loop step response.
 *   t_out/temp_out/u_out are optional (may be NULL) arrays of n samples.
 *   Returns IR_OK, or IR_ERR_PARAM on a bad configuration.
 * Metrics are computed against a +-0.1 K settling band and a dead zone
 * induced residual error.
 */
int tec_run_step(const tec_plant_cfg_t *plant, const tec_pid_cfg_t *pid,
                 double setpoint, double t0_temp, double dt, size_t n,
                 double *t_out, double *temp_out, double *u_out,
                 tec_metrics_t *m);

#ifdef __cplusplus
}
#endif

#endif /* TEC_CTRL_H */
