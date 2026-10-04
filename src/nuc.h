/*
 * nuc.h -- non-uniformity correction (NUC) for an infrared focal plane array
 * -------------------------------------------------------------------------
 * WHY
 * ---
 * Two pixels of a real IR FPA never see the same irradiance the same way.
 * A widely used first order (radiometric) model of a single detector is
 *
 *      S_ij(phi) = O_ij + G_ij * phi * (1 + C_ij * x) + n_ij(t)
 *
 *      phi   incident flux (arbitrary units, proportional to scene radiance)
 *      x     = (phi - phi_ref) / phi_ref, normalised flux
 *      O_ij  per pixel additive offset   (dark current, ROIC offset, 1/f...)
 *      G_ij  per pixel multiplicative gain (responsivity spread)
 *      C_ij  per pixel residual curvature: a cheap stand-in for the
 *            non-linearity that a two point calibration cannot remove
 *      n_ij  zero mean temporal noise
 *
 * The spread of O_ij and G_ij over the array is the *fixed pattern noise*
 * (FPN); it is what makes an uncorrected IR image look like it was shot
 * through frosted glass.
 *
 * FIGURE OF MERIT
 * ---------------
 * Non-uniformity is defined here exactly as the IR community does it, on a
 * *uniform* scene so that only spatial structure is measured:
 *
 *      NU [%] = 100 * sigma(S) / |mean(S)|          (population std dev)
 *
 * TWO POINT CORRECTION
 * --------------------
 * Expose the array to two uniform blackbody scenes, phi_low and phi_high,
 * and store the two frames S_L, S_H.  For every pixel solve the linear
 * system S = a*phi + b:
 *
 *      a_ij = (phi_high - phi_low) / (S_H,ij - S_L,ij)
 *      b_ij = phi_low - a_ij * S_L,ij
 *
 * Correction is then a per pixel affine map, corrected = a_ij*S + b_ij, which
 * is exact for the linear part of the model at *every* flux and leaves only
 * the temporal noise plus whatever the curvature term contributes.
 *
 * ONE POINT CORRECTION (for comparison)
 * -------------------------------------
 *      corrected = S - S_L,ij + phi_low      (a = 1 for every pixel)
 * It removes the offset spread only; the gain spread survives and grows with
 * the distance from phi_low, which is why the metric below evaluates at a
 * flux that is *not* the calibration point.
 */

#ifndef NUC_H
#define NUC_H

#include <stddef.h>
#include <stdint.h>

#include "ir_lab.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t w, h;
    double gain_sigma;    /* std dev of G_ij around 1.0            */
    double offset_sigma;  /* std dev of O_ij, in flux units        */
    double curv_sigma;    /* std dev of C_ij                       */
    double noise_sigma;   /* std dev of the temporal noise n_ij    */
    double flux_ref;      /* phi_ref used by the curvature term    */
} nuc_array_cfg_t;

/* Ground truth of the simulated array (known only to the generator). */
typedef struct {
    size_t w, h;
    double flux_ref;
    double *gain;    /* G_ij */
    double *offset;  /* O_ij */
    double *curv;    /* C_ij */
} nuc_truth_t;

/* Per pixel correction coefficients: corrected = gain*S + offset */
typedef struct {
    size_t w, h;
    double flux_low, flux_high;   /* the two calibration fluxes          */
    double *gain;                 /* a_ij                                */
    double *offset;               /* b_ij                                */
    int     mode;                 /* NUC_MODE_TWO_POINT / _ONE_POINT     */
} nuc_coeff_t;

#define NUC_MODE_TWO_POINT 2
#define NUC_MODE_ONE_POINT 1

void nuc_array_cfg_default(nuc_array_cfg_t *c);
/* 128x128 array, deterministic pseudo random pixel parameters. */
int  nuc_truth_init(nuc_truth_t *t, const nuc_array_cfg_t *c, uint32_t seed);
void nuc_truth_free(nuc_truth_t *t);

int  nuc_coeff_init(nuc_coeff_t *k, size_t w, size_t h);
void nuc_coeff_free(nuc_coeff_t *k);

/* Render one frame at a given uniform flux.
 * with_noise = 0 gives the noise free response (useful as an oracle). */
int nuc_render(const nuc_truth_t *t, const nuc_array_cfg_t *c, double flux,
               int with_noise, ir_rand_t *rng, double *out);

/* Calibration from two uniform blackbody frames. */
int nuc_calib_two_point(nuc_coeff_t *k,
                        const double *s_low, const double *s_high,
                        double phi_low, double phi_high);
/* Single point (offset only) calibration. */
int nuc_calib_one_point(nuc_coeff_t *k,
                        const double *s_low, double phi_low);

int nuc_apply(const nuc_coeff_t *k, const double *in, double *out);

/* Convenience metric wrapper: NU [%] of an array. */
double nuc_nu(const double *frame, size_t n, double *mean_out, double *sd_out);

#ifdef __cplusplus
}
#endif

#endif /* NUC_H */
