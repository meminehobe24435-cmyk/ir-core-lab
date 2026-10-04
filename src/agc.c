/* agc.c -- see agc.h for the algorithms. */

#include "agc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void agc_cfg_default(agc_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->bins = AGC_MAX_BINS;
    c->mode = AGC_MODE_PLATEAU;
    c->pl_low = 1.0;
    c->pl_high = 99.0;
    /*
     * Plateau level = plateau_frac * (N / bins), i.e. a multiple of the mean
     * bin population.  1.0 clips the dominant background mode at the average
     * bin height and redistributes the excess uniformly.
     * The sweep in results/agc_plateau_sweep.csv shows what the knob does:
     * small values converge to the linear stretch (best target/background
     * contrast, lowest background noise gain) and large values converge to
     * plain HEQ (highest global sd and local contrast, worst CNR).  1.0 is the
     * compromise used by the simulation.
     */
    c->plateau_frac = 1.0;
    c->out_min = 0.0;
    c->out_max = 255.0;
}

void agc_scene_cfg_default(agc_scene_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->w = 128u;
    c->h = 128u;
    /* A ground/sky background with a gentle thermal gradient, a warm target
     * and a cold patch.  Values are in the same arbitrary flux units the NUC
     * module uses, so the whole chain is dimensionally consistent. */
    c->bg_base = 200.0;
    c->bg_grad_x = 15.0;
    c->bg_grad_y = -8.0;
    c->target_delta = 30.0;
    c->target_size = 20;
    c->target_x = 40;
    c->target_y = 40;
    c->cold_delta = -15.0;
    c->cold_size = 12;
    c->cold_x = 90;
    c->cold_y = 90;
    c->noise_sigma = 2.0;
}

int agc_analyze(const agc_cfg_t *cfg, const double *img, size_t n,
                agc_analysis_t *a)
{
    double lo, hi, span;
    size_t i;
    int bins;

    if (cfg == NULL || img == NULL || a == NULL || n == 0u) {
        return IR_ERR_PARAM;
    }
    bins = cfg->bins;
    if (bins < 2 || bins > AGC_MAX_BINS) {
        return IR_ERR_RANGE;
    }

    memset(a, 0, sizeof(*a));
    a->bins = bins;
    a->total = n;

    lo = ir_percentile_copy(img, n, cfg->pl_low);
    hi = ir_percentile_copy(img, n, cfg->pl_high);
    if (!(hi > lo)) {
        /* Degenerate (constant) frame: fall back to a unit wide window so the
         * mapping stays well defined instead of dividing by zero. */
        lo = lo - 0.5;
        hi = lo + 1.0;
    }
    a->lo = lo;
    a->hi = hi;
    span = hi - lo;

    for (i = 0u; i < n; ++i) {
        double v = img[i];
        int idx;
        if (v <= lo) {
            idx = 0;
        } else if (v >= hi) {
            idx = bins - 1;
        } else {
            idx = (int)((v - lo) / span * (double)bins);
            if (idx < 0) {
                idx = 0;
            }
            if (idx >= bins) {
                idx = bins - 1;
            }
        }
        a->hist[idx]++;
    }

    /* Plateau clipping + uniform redistribution of the excess. */
    {
        double t_d = cfg->plateau_frac * (double)n / (double)bins;
        size_t t, excess = 0u, share, rem, j;
        if (cfg->mode == AGC_MODE_HEQ) {
            t = (size_t)-1;      /* no clipping: plain histogram equalisation */
        } else {
            /*
             * The plateau is clamped to at least one count.  A plateau of zero
             * would discard the histogram entirely and produce a flat,
             * uniform distribution, i.e. a linear ramp -- so "T -> 0" is the
             * limit in which plateau equalisation degenerates into the linear
             * stretch.  (Getting this backwards, by treating a sub-count
             * plateau as "no clipping", turns the small-plateau limit into
             * full HEQ instead: the test suite caught exactly that.)
             */
            t = (t_d > 0.0) ? (size_t)(t_d + 0.5) : 0u;
            if (t < 1u) {
                t = 1u;
            }
        }
        for (j = 0u; j < (size_t)bins; ++j) {
            if (a->hist[j] > t) {
                excess += a->hist[j] - t;
                a->hist_clipped[j] = t;
            } else {
                a->hist_clipped[j] = a->hist[j];
            }
        }
        share = excess / (size_t)bins;
        rem = excess % (size_t)bins;
        for (j = 0u; j < (size_t)bins; ++j) {
            a->hist_clipped[j] += share;
            if (j < rem) {
                a->hist_clipped[j] += 1u;
            }
        }
    }

    return IR_OK;
}

int agc_build_map(const agc_cfg_t *cfg, const agc_analysis_t *a,
                  double *map)
{
    int i, bins;
    double span;

    if (cfg == NULL || a == NULL || map == NULL) {
        return IR_ERR_PARAM;
    }
    bins = a->bins;
    if (bins < 2) {
        return IR_ERR_RANGE;
    }
    span = cfg->out_max - cfg->out_min;

    if (cfg->mode == AGC_MODE_LINEAR) {
        for (i = 0; i < bins; ++i) {
            map[i] = cfg->out_min + span * (double)i / (double)(bins - 1);
        }
        return IR_OK;
    }

    {
        const size_t *h = (cfg->mode == AGC_MODE_HEQ) ? a->hist
                                                      : a->hist_clipped;
        size_t cdf = 0u, cdf_last = 0u;
        for (i = 0; i < bins; ++i) {
            cdf_last += h[i];
        }
        if (cdf_last == 0u) {
            for (i = 0; i < bins; ++i) {
                map[i] = cfg->out_min + span * (double)i / (double)(bins - 1);
            }
            return IR_OK;
        }
        for (i = 0; i < bins; ++i) {
            cdf += h[i];
            map[i] = cfg->out_min + span * (double)cdf / (double)cdf_last;
        }
        /* The last entry is exactly out_max by construction. */
        map[bins - 1] = cfg->out_max;
    }
    return IR_OK;
}

static int agc_value_to_bin(const agc_analysis_t *a, double v)
{
    int idx;
    double span = a->hi - a->lo;
    if (span <= 0.0) {
        return 0;
    }
    if (v <= a->lo) {
        return 0;
    }
    if (v >= a->hi) {
        return a->bins - 1;
    }
    idx = (int)((v - a->lo) / span * (double)a->bins);
    if (idx < 0) {
        idx = 0;
    }
    if (idx >= a->bins) {
        idx = a->bins - 1;
    }
    return idx;
}

int agc_apply(const agc_cfg_t *cfg, const agc_analysis_t *a,
              const double *img, size_t n, double *outd, uint8_t *out8)
{
    double *map;
    size_t i;

    if (cfg == NULL || a == NULL || img == NULL || outd == NULL) {
        return IR_ERR_PARAM;
    }
    map = (double *)malloc((size_t)a->bins * sizeof(double));
    if (map == NULL) {
        return IR_ERR_IO;
    }
    if (agc_build_map(cfg, a, map) != IR_OK) {
        free(map);
        return IR_ERR_PARAM;
    }

    for (i = 0u; i < n; ++i) {
        int idx = agc_value_to_bin(a, img[i]);
        double v = map[idx];
        outd[i] = v;
        if (out8 != NULL) {
            double q = v + 0.5;
            if (q < 0.0) {
                q = 0.0;
            }
            if (q > 255.0) {
                q = 255.0;
            }
            out8[i] = (uint8_t)q;
        }
    }

    free(map);
    return IR_OK;
}

void agc_metrics(const double *out, size_t n, size_t w, size_t h,
                 const uint8_t *target, const uint8_t *bg, agc_metrics_t *m)
{
    double sd, mean, loc_sum = 0.0;
    size_t loc_n = 0u, i, x, y;

    if (m == NULL) {
        return;
    }
    memset(m, 0, sizeof(*m));
    if (out == NULL || n == 0u) {
        return;
    }

    sd = ir_stddev(out, n, &mean);
    m->mean = mean;
    m->sd = sd;
    m->rms_contrast_pct = (fabs(mean) > 1e-12) ? 100.0 * sd / fabs(mean) : 0.0;
    m->out_min = ir_min(out, n);
    m->out_max = ir_max(out, n);

    /* Local contrast: 3x3 window, sd/mean, border pixels skipped. */
    if (w >= 3u && h >= 3u && w * h == n) {
        for (y = 1u; y + 1u < h; ++y) {
            for (x = 1u; x + 1u < w; ++x) {
                double buf[9];
                int k = 0, dx, dy;
                double lm, ls;
                for (dy = -1; dy <= 1; ++dy) {
                    for (dx = -1; dx <= 1; ++dx) {
                        size_t xx = (size_t)((long)x + dx);
                        size_t yy = (size_t)((long)y + dy);
                        buf[k++] = out[yy * w + xx];
                    }
                }
                ls = ir_stddev(buf, 9u, &lm);
                if (fabs(lm) > 1e-12) {
                    loc_sum += ls / fabs(lm);
                    ++loc_n;
                }
            }
        }
    }
    m->local_contrast = (loc_n > 0u) ? loc_sum / (double)loc_n : 0.0;

    /* Target / background contrast. */
    if (target != NULL) {
        double st = 0.0, sb = 0.0, sqb = 0.0;
        size_t nt = 0u, nb = 0u;
        for (i = 0u; i < n; ++i) {
            if (target[i]) {
                st += out[i];
                ++nt;
            }
            if (bg != NULL && bg[i]) {
                sb += out[i];
                ++nb;
            }
        }
        if (nt > 0u && nb > 0u) {
            double mt = st / (double)nt;
            double mb = sb / (double)nb;
            double d = fabs(mt - mb);
            double varb;
            if (mt + mb > 1e-12) {
                m->tb_michelson_pct = 100.0 * d / (mt + mb);
            }
            if (fabs(mb) > 1e-12) {
                m->tb_weber_pct = 100.0 * d / fabs(mb);
            }
            /* Background std dev is computed over the local annulus only. */
            for (i = 0u; i < n; ++i) {
                if (bg != NULL && bg[i]) {
                    double dd = out[i] - mb;
                    sqb += dd * dd;
                }
            }
            varb = sqb / (double)nb;
            m->bg_sd = sqrt(varb);
            if (m->bg_sd > 1e-12) {
                m->cnr = d / m->bg_sd;
            }
        }
    }

    /* Entropy over a 256 level output histogram. */
    {
        size_t hist[256];
        double lo = m->out_min, hi = m->out_max, span = hi - lo;
        double ent = 0.0;
        int b;
        memset(hist, 0, sizeof(hist));
        for (i = 0u; i < n; ++i) {
            int idx;
            if (span <= 0.0) {
                idx = 0;
            } else {
                idx = (int)((out[i] - lo) / span * 255.0 + 0.5);
                if (idx < 0) {
                    idx = 0;
                }
                if (idx > 255) {
                    idx = 255;
                }
            }
            hist[idx]++;
        }
        for (b = 0; b < 256; ++b) {
            if (hist[b] > 0u) {
                double p = (double)hist[b] / (double)n;
                ent -= p * (log(p) / log(2.0));
            }
        }
        m->entropy_bits = ent;
    }
}

int agc_scene_synth(agc_scene_t *s, const agc_scene_cfg_t *c, uint32_t seed)
{
    ir_rand_t rng;
    size_t n, i, x, y;

    if (s == NULL || c == NULL || c->w == 0u || c->h == 0u) {
        return IR_ERR_PARAM;
    }
    n = c->w * c->h;
    s->w = c->w;
    s->h = c->h;
    s->img    = (double *)malloc(n * sizeof(double));
    s->target = (uint8_t *)malloc(n);
    s->bg     = (uint8_t *)malloc(n);
    if (s->img == NULL || s->target == NULL || s->bg == NULL) {
        agc_scene_free(s);
        return IR_ERR_IO;
    }
    memset(s->target, 0, n);
    memset(s->bg, 0, n);

    ir_rand_seed(&rng, seed);
    for (y = 0u; y < c->h; ++y) {
        for (x = 0u; x < c->w; ++x) {
            size_t idx = y * c->w + x;
            double v = c->bg_base;
            if (c->w > 1u) {
                v += c->bg_grad_x * (double)x / (double)(c->w - 1u);
            }
            if (c->h > 1u) {
                v += c->bg_grad_y * (double)y / (double)(c->h - 1u);
            }
            s->img[idx] = v;
        }
    }

    /* Hot square target with a 3 pixel wide background annulus around it. */
    if (c->target_size > 0) {
        int tx0 = c->target_x, ty0 = c->target_y;
        int tx1 = tx0 + c->target_size, ty1 = ty0 + c->target_size;
        for (y = 0u; y < c->h; ++y) {
            for (x = 0u; x < c->w; ++x) {
                int ix = (int)x, iy = (int)y;
                size_t idx = y * c->w + x;
                if (ix >= tx0 && ix < tx1 && iy >= ty0 && iy < ty1) {
                    s->img[idx] += c->target_delta;
                    s->target[idx] = 1u;
                } else if (ix >= tx0 - 6 && ix < tx1 + 6 &&
                           iy >= ty0 - 6 && iy < ty1 + 6) {
                    s->bg[idx] = 1u;
                }
            }
        }
    }
    if (c->cold_size > 0) {
        int tx0 = c->cold_x, ty0 = c->cold_y;
        int tx1 = tx0 + c->cold_size, ty1 = ty0 + c->cold_size;
        for (y = 0u; y < c->h; ++y) {
            for (x = 0u; x < c->w; ++x) {
                int ix = (int)x, iy = (int)y;
                if (ix >= tx0 && ix < tx1 && iy >= ty0 && iy < ty1) {
                    s->img[y * c->w + x] += c->cold_delta;
                }
            }
        }
    }

    if (c->noise_sigma > 0.0) {
        for (i = 0u; i < n; ++i) {
            s->img[i] += c->noise_sigma * ir_rand_gauss(&rng);
        }
    }
    return IR_OK;
}

void agc_scene_free(agc_scene_t *s)
{
    if (s == NULL) {
        return;
    }
    free(s->img);
    free(s->target);
    free(s->bg);
    s->img = NULL;
    s->target = NULL;
    s->bg = NULL;
    s->w = 0u;
    s->h = 0u;
}
