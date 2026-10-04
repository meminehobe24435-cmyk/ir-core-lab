/*
 * badpix.h -- blind (defective) pixel detection and replacement
 * ------------------------------------------------------------
 * DETECTION CRITERION
 * -------------------
 * The pixel response rate is the differential response to two blackbody
 * scenes, r_ij = S_H,ij - S_L,ij.  A pixel is declared blind when its
 * response rate falls outside the array's acceptance window:
 *
 *      |r_ij - mu| > k * sigma          (k = 3 by default)
 *
 * The catch is that mu and sigma must be estimated from an array that
 * *contains* the defective pixels.  A single pass over all pixels lets the
 * outliers inflate sigma, which widens the window and hides the marginal
 * defects (a pixel at half responsivity is only ~3.3 sigma away from the
 * good population but the inflated sigma pushes it inside the window).
 * The standard fix is iterated sigma clipping:
 *
 *      repeat
 *          mu, sigma <- mean and std dev of the *accepted* pixels only
 *          flag every pixel outside [mu - k*sigma, mu + k*sigma]
 *      until no new pixel is flagged (or max_iter reached)
 *
 * Both variants are implemented (`iterations = 1` gives the naive single
 * pass) so that the benefit can be measured instead of asserted.
 *
 * A hard guard stops the clipping when more than `max_flagged_frac` of the
 * array would be rejected: an over-aggressive clip on a nearly uniform array
 * can otherwise walk the window down until most of the array is "defective".
 *
 * REPLACEMENT
 * -----------
 * A flagged pixel is replaced by the median of its valid 3x3 neighbours
 * (median rather than mean so that a cluster of dead pixels does not bias the
 * estimate).  Neighbours that are themselves flagged are excluded; if fewer
 * than `min_valid` valid neighbours exist the original value is kept.
 * The neighbourhood is always read from the *input* image, never from the
 * partially rewritten output, so the result does not depend on scan order.
 */

#ifndef BADPIX_H
#define BADPIX_H

#include <stddef.h>
#include <stdint.h>

#include "ir_lab.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BADPIX_REPLACE_MEDIAN3 0
#define BADPIX_REPLACE_MEAN3   1

typedef struct {
    double k_sigma;          /* threshold in sigma, 3.0 typical          */
    int    iterations;       /* >=1; 1 = naive single pass               */
    double max_flagged_frac; /* abort clipping above this rejected frac  */
} badpix_cfg_t;

typedef struct {
    double mean;             /* accepted population mean                 */
    double sd;               /* accepted population std dev              */
    double lo, hi;           /* final acceptance window                  */
    size_t flagged;
    size_t accepted;
    int    iterations_used;
    int    aborted;          /* 1 if the max_flagged_frac guard tripped  */
} badpix_stats_t;

typedef struct {
    int    method;           /* BADPIX_REPLACE_*                         */
    size_t min_valid;        /* minimum valid neighbours, default 3      */
} badpix_replace_cfg_t;

typedef struct {
    size_t n_injected;
    size_t detected;         /* injected and flagged                     */
    size_t missed;
    size_t false_alarm;      /* flagged but not injected                 */
    size_t total;
    double detection_rate_pct;
    double false_alarm_pct;  /* false alarms / number of good pixels     */
} badpix_eval_t;

typedef struct {
    size_t n;                /* how many pixels were evaluated           */
    double mae_before, mae_after;
    double rmse_before, rmse_after;
    double max_before, max_after;
} badpix_residual_t;

void badpix_cfg_default(badpix_cfg_t *c);
void badpix_replace_cfg_default(badpix_replace_cfg_t *c);

/* flags must point at n bytes; they are zeroed by this function. */
int  badpix_detect(const double *resp, size_t n, const badpix_cfg_t *cfg,
                   uint8_t *flags, badpix_stats_t *st);

/*
 * Replaces flagged pixels of `in` into `out` (in and out may be the same).
 * *n_replaced counts only the pixels whose value was actually changed; a
 * flagged pixel with fewer than `min_valid` valid neighbours keeps its
 * original value and is not counted.
 */
int  badpix_replace(const double *in, size_t w, size_t h,
                    const uint8_t *flags, const badpix_replace_cfg_t *cfg,
                    double *out, size_t *n_replaced);

void badpix_eval(const uint8_t *flags, const uint8_t *injected, size_t n,
                 badpix_eval_t *ev);

/* Residuals against an ideal (nominal) response, over the pixels selected by
 * `mask` (pass NULL to use every pixel). */
void badpix_residual(const double *resp, const double *replaced,
                     const double *ideal, const uint8_t *mask, size_t n,
                     badpix_residual_t *r);

#ifdef __cplusplus
}
#endif

#endif /* BADPIX_H */
