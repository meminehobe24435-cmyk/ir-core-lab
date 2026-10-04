/* badpix.c -- blind pixel detection / replacement, see badpix.h. */

#include "badpix.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void badpix_cfg_default(badpix_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->k_sigma = 3.0;
    c->iterations = 8;          /* with early exit when nothing new is flagged */
    c->max_flagged_frac = 0.20; /* never reject more than 20 % of the array    */
}

void badpix_replace_cfg_default(badpix_replace_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->method = BADPIX_REPLACE_MEDIAN3;
    c->min_valid = 3u;
}

static int cmp_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) {
        return -1;
    }
    if (da > db) {
        return 1;
    }
    return 0;
}

int badpix_detect(const double *resp, size_t n, const badpix_cfg_t *cfg,
                  uint8_t *flags, badpix_stats_t *st)
{
    size_t i, flagged, accepted;
    double mu, sd, lo, hi;
    int it, used = 0, aborted = 0;

    if (resp == NULL || cfg == NULL || flags == NULL || n == 0u) {
        return IR_ERR_PARAM;
    }
    if (cfg->k_sigma <= 0.0 || cfg->iterations < 1) {
        return IR_ERR_PARAM;
    }

    memset(flags, 0, n);

    /* Seed the statistics from the whole array (this is what the naive
     * single pass would use as its final answer as well). */
    sd = ir_stddev(resp, n, &mu);
    (void)sd;

    for (it = 0; it < cfg->iterations; ++it) {
        size_t new_flags = 0u;
        double s = 0.0, s2 = 0.0;

        accepted = 0u;
        for (i = 0u; i < n; ++i) {
            if (!flags[i]) {
                s += resp[i];
                ++accepted;
            }
        }
        if (accepted < 2u) {
            /* Everything is already rejected: nothing left to estimate from. */
            aborted = 1;
            break;
        }
        mu = s / (double)accepted;
        for (i = 0u; i < n; ++i) {
            if (!flags[i]) {
                double d = resp[i] - mu;
                s2 += d * d;
            }
        }
        sd = sqrt(s2 / (double)accepted);
        lo = mu - cfg->k_sigma * sd;
        hi = mu + cfg->k_sigma * sd;

        for (i = 0u; i < n; ++i) {
            if (!flags[i] && (resp[i] < lo || resp[i] > hi)) {
                flags[i] = 1u;
                ++new_flags;
            }
        }

        /* Count the total flagged so far. */
        flagged = 0u;
        for (i = 0u; i < n; ++i) {
            if (flags[i]) {
                ++flagged;
            }
        }
        used = it + 1;

        if (cfg->max_flagged_frac > 0.0 &&
            (double)flagged > cfg->max_flagged_frac * (double)n) {
            aborted = 1;
            break;
        }
        if (new_flags == 0u) {
            break;   /* converged: the window no longer rejects anything */
        }
    }

    if (st != NULL) {
        size_t acc = 0u, fl = 0u;
        double s = 0.0, s2 = 0.0;
        for (i = 0u; i < n; ++i) {
            if (flags[i]) {
                ++fl;
            } else {
                s += resp[i];
                ++acc;
            }
        }
        st->flagged = fl;
        st->accepted = acc;
        st->iterations_used = used;
        st->aborted = aborted;
        if (acc > 0u) {
            mu = s / (double)acc;
            for (i = 0u; i < n; ++i) {
                if (!flags[i]) {
                    double d = resp[i] - mu;
                    s2 += d * d;
                }
            }
            sd = sqrt(s2 / (double)acc);
        } else {
            mu = 0.0;
            sd = 0.0;
        }
        st->mean = mu;
        st->sd = sd;
        st->lo = mu - cfg->k_sigma * sd;
        st->hi = mu + cfg->k_sigma * sd;
    }

    return IR_OK;
}

int badpix_replace(const double *in, size_t w, size_t h,
                   const uint8_t *flags, const badpix_replace_cfg_t *cfg,
                   double *out, size_t *n_replaced)
{
    badpix_replace_cfg_t local;
    size_t x, y, n, replaced = 0u;

    if (in == NULL || flags == NULL || out == NULL || w == 0u || h == 0u) {
        return IR_ERR_PARAM;
    }
    n = w * h;

    if (cfg != NULL) {
        local = *cfg;
    } else {
        badpix_replace_cfg_default(&local);
    }
    if (local.min_valid < 1u) {
        local.min_valid = 1u;
    }

    for (y = 0u; y < h; ++y) {
        for (x = 0u; x < w; ++x) {
            size_t idx = y * w + x;
            double buf[9];
            size_t nv = 0u, k;
            int dx, dy;

            if (!flags[idx]) {
                out[idx] = in[idx];
                continue;
            }

            for (dy = -1; dy <= 1; ++dy) {
                for (dx = -1; dx <= 1; ++dx) {
                    size_t nx, ny, nidx;
                    if (dx == 0 && dy == 0) {
                        continue;
                    }
                    /* Border pixels simply have fewer neighbours. */
                    if (dx < 0 && x == 0u) {
                        continue;
                    }
                    if (dy < 0 && y == 0u) {
                        continue;
                    }
                    nx = (size_t)((long)x + dx);
                    ny = (size_t)((long)y + dy);
                    if (nx >= w || ny >= h) {
                        continue;
                    }
                    nidx = ny * w + nx;
                    if (flags[nidx]) {
                        continue;   /* never borrow from another blind pixel */
                    }
                    buf[nv++] = in[nidx];
                }
            }

            if (nv < local.min_valid) {
                out[idx] = in[idx];   /* not enough information: keep it */
            } else {
                if (local.method == BADPIX_REPLACE_MEAN3) {
                    double s = 0.0;
                    for (k = 0u; k < nv; ++k) {
                        s += buf[k];
                    }
                    out[idx] = s / (double)nv;
                } else {
                    qsort(buf, nv, sizeof(double), cmp_double);
                    out[idx] = (nv % 2u == 1u)
                                   ? buf[nv / 2u]
                                   : 0.5 * (buf[nv / 2u - 1u] + buf[nv / 2u]);
                }
                /* Only count pixels whose value actually changed: a flagged
                 * pixel with too few valid neighbours is left alone, and
                 * reporting it as "replaced" hid that fact during development. */
                ++replaced;
            }
        }
    }

    if (n_replaced != NULL) {
        *n_replaced = replaced;
    }
    (void)n;
    return IR_OK;
}

void badpix_eval(const uint8_t *flags, const uint8_t *injected, size_t n,
                 badpix_eval_t *ev)
{
    size_t i, inj = 0u, det = 0u, fp = 0u, good = 0u;

    if (ev == NULL) {
        return;
    }
    memset(ev, 0, sizeof(*ev));
    if (flags == NULL || injected == NULL || n == 0u) {
        return;
    }

    for (i = 0u; i < n; ++i) {
        if (injected[i]) {
            ++inj;
            if (flags[i]) {
                ++det;
            }
        } else {
            ++good;
            if (flags[i]) {
                ++fp;
            }
        }
    }

    ev->total = n;
    ev->n_injected = inj;
    ev->detected = det;
    ev->missed = inj - det;
    ev->false_alarm = fp;
    ev->detection_rate_pct = (inj > 0u) ? 100.0 * (double)det / (double)inj : 0.0;
    ev->false_alarm_pct = (good > 0u) ? 100.0 * (double)fp / (double)good : 0.0;
}

void badpix_residual(const double *resp, const double *replaced,
                     const double *ideal, const uint8_t *mask, size_t n,
                     badpix_residual_t *r)
{
    size_t i, cnt = 0u;
    double sum_b = 0.0, sum_a = 0.0, sq_b = 0.0, sq_a = 0.0;
    double mx_b = 0.0, mx_a = 0.0;

    if (r == NULL) {
        return;
    }
    memset(r, 0, sizeof(*r));
    if (resp == NULL || replaced == NULL || ideal == NULL || n == 0u) {
        return;
    }

    for (i = 0u; i < n; ++i) {
        double db, da;
        if (mask != NULL && !mask[i]) {
            continue;
        }
        ++cnt;
        db = fabs(resp[i] - ideal[i]);
        da = fabs(replaced[i] - ideal[i]);
        sum_b += db;
        sum_a += da;
        sq_b += db * db;
        sq_a += da * da;
        if (db > mx_b) {
            mx_b = db;
        }
        if (da > mx_a) {
            mx_a = da;
        }
    }

    if (cnt == 0u) {
        return;
    }
    r->n = cnt;
    r->mae_before  = sum_b / (double)cnt;
    r->mae_after   = sum_a / (double)cnt;
    r->rmse_before = sqrt(sq_b / (double)cnt);
    r->rmse_after  = sqrt(sq_a / (double)cnt);
    r->max_before  = mx_b;
    r->max_after   = mx_a;
}
