/*
 * agc.h -- infrared image dynamic range compression / enhancement
 * --------------------------------------------------------------
 * A 14 bit IR frame has to become an 8 bit display frame without losing the
 * low contrast detail a human observer is looking for.  Three mappings are
 * implemented so that they can be compared on the same scene:
 *
 *  1. LINEAR STRETCH (percentile clipping)
 *         lo, hi = the pl_low / pl_high percentiles of the frame
 *         out    = out_min + (out_max-out_min) * (v-lo)/(hi-lo)
 *     Simple, monotonic, preserves radiometry, but a single hot object
 *     squeezes the rest of the scene into a few grey levels.
 *
 *  2. HISTOGRAM EQUALISATION (HEQ)
 *         cdf[i] = sum_{j<=i} h[j]
 *         out[i] = out_min + (out_max-out_min) * cdf[i] / cdf[last]
 *     Maximises the output entropy/contrast, but on an IR scene where 90 %
 *     of the pixels are background it over-amplifies background noise and can
 *     actually *reduce* the contrast of a small target.
 *
 *  3. PLATEAU HISTOGRAM EQUALISATION  (the usual choice in thermal imagers)
 *         T = plateau_frac * N / bins            (the "plateau")
 *         h'[i] = min(h[i], T);  excess = sum(h[i]-T) for h[i]>T
 *         then redistribute `excess` uniformly over all bins and equalise h'
 *     Clipping the dominant background mode limits its gain so the small
 *     targets keep a larger share of the output range.
 *
 * CONTRAST METRICS
 * ----------------
 *   global sd              sigma of the output frame (contrast energy)
 *   RMS contrast           sigma / mean
 *   local contrast         mean over the frame of sigma(3x3) / mean(3x3)
 *   Michelson T/B contrast |mt - mb| / (mt + mb) * 100, over a target mask and
 *                          a local background mask
 *   Weber T/B contrast     |mt - mb| / mb * 100
 *   CNR                    |mt - mb| / sigma_b -- the contrast to *noise*
 *                          ratio.  This is the figure that actually predicts
 *                          whether an observer can see the target: a mapping
 *                          that raises the target contrast while amplifying
 *                          background noise equally gains nothing.
 *   entropy                -sum p*log2(p) over a 256 bin output histogram
 */

#ifndef AGC_H
#define AGC_H

#include <stddef.h>
#include <stdint.h>

#include "ir_lab.h"

#ifdef __cplusplus
extern "C" {
#endif

#define AGC_MAX_BINS 256

#define AGC_MODE_LINEAR  0
#define AGC_MODE_HEQ     1
#define AGC_MODE_PLATEAU 2

typedef struct {
    int    bins;          /* histogram bins, <= AGC_MAX_BINS (256 typical) */
    int    mode;          /* AGC_MODE_*                                   */
    double pl_low;        /* lower percentile used to set the input range  */
    double pl_high;       /* upper percentile                              */
    double plateau_frac;  /* plateau level as a multiple of N/bins         */
    double out_min;       /* output black level                            */
    double out_max;       /* output white level                            */
} agc_cfg_t;

typedef struct {
    double lo, hi;             /* input range actually used              */
    int    bins;
    size_t hist[AGC_MAX_BINS]; /* input histogram over [lo, hi]          */
    size_t hist_clipped[AGC_MAX_BINS]; /* histogram after plateau clipping */
    size_t total;
} agc_analysis_t;

typedef struct {
    double mean, sd;
    double rms_contrast_pct;   /* 100*sd/mean                            */
    double local_contrast;     /* mean of local sd / local mean          */
    double entropy_bits;
    double tb_michelson_pct;
    double tb_weber_pct;
    double cnr;                /* |mt - mb| / sigma_b                    */
    double bg_sd;              /* sigma of the background mask           */
    double out_min, out_max;
} agc_metrics_t;

typedef struct {
    size_t w, h;
    double bg_base;        /* background pedestal                       */
    double bg_grad_x;      /* horizontal background gradient            */
    double bg_grad_y;      /* vertical background gradient              */
    double target_delta;   /* hot target contrast above the background  */
    int    target_size;    /* side of the square target, pixels         */
    int    target_x, target_y;
    double cold_delta;     /* a second, negative contrast target        */
    int    cold_size, cold_x, cold_y;
    double noise_sigma;
} agc_scene_cfg_t;

typedef struct {
    size_t w, h;
    double *img;
    uint8_t *target;   /* 1 for the hot target pixels          */
    uint8_t *bg;       /* 1 for the local background annulus   */
} agc_scene_t;

void agc_cfg_default(agc_cfg_t *c);
void agc_scene_cfg_default(agc_scene_cfg_t *c);

int  agc_analyze(const agc_cfg_t *cfg, const double *img, size_t n,
                 agc_analysis_t *a);

/* Builds a value -> display level mapping for bin index i: map[i]. */
int  agc_build_map(const agc_cfg_t *cfg, const agc_analysis_t *a,
                   double *map /* [bins] */);

/*
 * Maps the frame.  `outd` receives the un-quantised display value in
 * [out_min, out_max] (metrics are computed on it so that the 8 bit
 * quantisation does not pollute the numbers); out8 may be NULL.
 */
int  agc_apply(const agc_cfg_t *cfg, const agc_analysis_t *a,
               const double *img, size_t n, double *outd, uint8_t *out8);

void agc_metrics(const double *out, size_t n, size_t w, size_t h,
                 const uint8_t *target, const uint8_t *bg, agc_metrics_t *m);

int  agc_scene_synth(agc_scene_t *s, const agc_scene_cfg_t *c, uint32_t seed);
void agc_scene_free(agc_scene_t *s);

#ifdef __cplusplus
}
#endif

#endif /* AGC_H */
