/* nuc.c -- implementation of the two point / one point NUC, see nuc.h. */

#include "nuc.h"

#include <math.h>
#include <stdlib.h>

void nuc_array_cfg_default(nuc_array_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->w = 128u;
    c->h = 128u;
    /*
     * Numbers chosen to sit in the range quoted for uncooled microbolometer
     * and cooled MCT arrays in the open literature (responsivity spread of
     * order 5-20 % RMS, temporal NETD equivalent noise well under 1 % of the
     * signal).  They are order-of-magnitude placeholders, NOT measurements of
     * any vendor part -- see the boundary section of README.md.
     */
    c->gain_sigma   = 0.12;   /* 12 % RMS responsivity spread  */
    c->offset_sigma = 6.0;    /* ~2 % of the 300 count swing  */
    c->curv_sigma   = 0.02;   /* 2 % quadratic residue        */
    c->noise_sigma  = 0.6;    /* 0.3 % of the mid flux        */
    c->flux_ref     = 200.0;
}

int nuc_truth_init(nuc_truth_t *t, const nuc_array_cfg_t *c, uint32_t seed)
{
    ir_rand_t rng;
    size_t i, n;

    if (t == NULL || c == NULL || c->w == 0u || c->h == 0u) {
        return IR_ERR_PARAM;
    }
    n = c->w * c->h;
    t->w = c->w;
    t->h = c->h;
    t->flux_ref = c->flux_ref;
    t->gain   = (double *)malloc(n * sizeof(double));
    t->offset = (double *)malloc(n * sizeof(double));
    t->curv   = (double *)malloc(n * sizeof(double));
    if (t->gain == NULL || t->offset == NULL || t->curv == NULL) {
        nuc_truth_free(t);
        return IR_ERR_IO;
    }

    ir_rand_seed(&rng, seed);
    for (i = 0u; i < n; ++i) {
        double g = 1.0 + c->gain_sigma * ir_rand_gauss(&rng);
        double o = c->offset_sigma * ir_rand_gauss(&rng);
        double k = c->curv_sigma * ir_rand_gauss(&rng);
        /* Physical sanity clamps: a responsivity cannot be negative and a
         * 30 % curvature is already far outside any real detector. */
        if (g < 0.5) {
            g = 0.5;
        }
        if (g > 1.5) {
            g = 1.5;
        }
        if (k < -0.3) {
            k = -0.3;
        }
        if (k > 0.3) {
            k = 0.3;
        }
        t->gain[i]   = g;
        t->offset[i] = o;
        t->curv[i]   = k;
    }
    return IR_OK;
}

void nuc_truth_free(nuc_truth_t *t)
{
    if (t == NULL) {
        return;
    }
    free(t->gain);
    free(t->offset);
    free(t->curv);
    t->gain = NULL;
    t->offset = NULL;
    t->curv = NULL;
    t->w = 0u;
    t->h = 0u;
}

int nuc_coeff_init(nuc_coeff_t *k, size_t w, size_t h)
{
    size_t n;
    if (k == NULL || w == 0u || h == 0u) {
        return IR_ERR_PARAM;
    }
    n = w * h;
    k->w = w;
    k->h = h;
    k->flux_low = 0.0;
    k->flux_high = 0.0;
    k->mode = NUC_MODE_TWO_POINT;
    k->gain   = (double *)malloc(n * sizeof(double));
    k->offset = (double *)malloc(n * sizeof(double));
    if (k->gain == NULL || k->offset == NULL) {
        nuc_coeff_free(k);
        return IR_ERR_IO;
    }
    return IR_OK;
}

void nuc_coeff_free(nuc_coeff_t *k)
{
    if (k == NULL) {
        return;
    }
    free(k->gain);
    free(k->offset);
    k->gain = NULL;
    k->offset = NULL;
    k->w = 0u;
    k->h = 0u;
}

int nuc_render(const nuc_truth_t *t, const nuc_array_cfg_t *c, double flux,
               int with_noise, ir_rand_t *rng, double *out)
{
    size_t i, n;
    double x;

    if (t == NULL || c == NULL || out == NULL || t->w != c->w || t->h != c->h) {
        return IR_ERR_PARAM;
    }
    if (with_noise && rng == NULL) {
        return IR_ERR_PARAM;
    }
    n = c->w * c->h;
    x = (flux - t->flux_ref) / t->flux_ref;

    for (i = 0u; i < n; ++i) {
        double s = t->offset[i] + t->gain[i] * flux * (1.0 + t->curv[i] * x);
        if (with_noise) {
            s += c->noise_sigma * ir_rand_gauss(rng);
        }
        out[i] = s;
    }
    return IR_OK;
}

int nuc_calib_two_point(nuc_coeff_t *k,
                        const double *s_low, const double *s_high,
                        double phi_low, double phi_high)
{
    size_t i, n;

    if (k == NULL || s_low == NULL || s_high == NULL) {
        return IR_ERR_PARAM;
    }
    if (k->gain == NULL || k->offset == NULL) {
        return IR_ERR_PARAM;
    }
    if (phi_high == phi_low) {
        return IR_ERR_RANGE;   /* degenerate: no flux difference */
    }

    n = k->w * k->h;
    for (i = 0u; i < n; ++i) {
        double d = s_high[i] - s_low[i];
        double a;
        /*
         * A pixel whose two responses are (numerically) identical carries no
         * information about the flux difference: it is dead.  The calibration
         * cannot invert it, so it is flagged by storing a = 1 (i.e. "pass
         * through") and leaving the offset to the blind pixel stage.  This
         * mirrors what production firmware does instead of dividing by zero.
         */
        if (fabs(d) < 1e-9) {
            a = 1.0;
        } else {
            a = (phi_high - phi_low) / d;
        }
        k->gain[i]   = a;
        k->offset[i] = phi_low - a * s_low[i];
    }
    k->flux_low = phi_low;
    k->flux_high = phi_high;
    k->mode = NUC_MODE_TWO_POINT;
    return IR_OK;
}

int nuc_calib_one_point(nuc_coeff_t *k, const double *s_low, double phi_low)
{
    size_t i, n;

    if (k == NULL || s_low == NULL || k->gain == NULL || k->offset == NULL) {
        return IR_ERR_PARAM;
    }
    n = k->w * k->h;
    for (i = 0u; i < n; ++i) {
        k->gain[i]   = 1.0;
        k->offset[i] = phi_low - s_low[i];
    }
    k->flux_low = phi_low;
    k->flux_high = phi_low;
    k->mode = NUC_MODE_ONE_POINT;
    return IR_OK;
}

int nuc_apply(const nuc_coeff_t *k, const double *in, double *out)
{
    size_t i, n;

    if (k == NULL || in == NULL || out == NULL) {
        return IR_ERR_PARAM;
    }
    if (k->gain == NULL || k->offset == NULL) {
        return IR_ERR_PARAM;
    }
    n = k->w * k->h;
    for (i = 0u; i < n; ++i) {
        out[i] = k->gain[i] * in[i] + k->offset[i];
    }
    return IR_OK;
}

double nuc_nu(const double *frame, size_t n, double *mean_out, double *sd_out)
{
    return ir_nu_percent(frame, n, mean_out, sd_out);
}
