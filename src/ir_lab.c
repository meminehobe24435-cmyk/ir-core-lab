/* ir_lab.c -- see ir_lab.h for the rationale of every routine. */

#include "ir_lab.h"

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* RNG                                                                 */
/* ------------------------------------------------------------------ */

void ir_rand_seed(ir_rand_t *r, uint32_t seed)
{
    if (r == NULL) {
        return;
    }
    /* xorshift32 locks up on a zero state; 0x9E3779B9 is a safe stand-in. */
    r->s = (seed == 0u) ? 0x9E3779B9u : seed;
    r->have_spare = 0;
    r->spare = 0.0;
}

uint32_t ir_rand_u32(ir_rand_t *r)
{
    uint32_t x;
    if (r == NULL) {
        return 0u;
    }
    x = r->s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    r->s = x;
    return x;
}

double ir_rand_uniform(ir_rand_t *r)
{
    /* 24 bits -> [0, 1) with a resolution of 2^-24. */
    return (double)(ir_rand_u32(r) >> 8) * (1.0 / 16777216.0);
}

double ir_rand_gauss(ir_rand_t *r)
{
    double u1, u2, mag;
    if (r == NULL) {
        return 0.0;
    }
    if (r->have_spare) {
        r->have_spare = 0;
        return r->spare;
    }
    /* Avoid log(0): u1 is drawn from (0, 1]. */
    do {
        u1 = ir_rand_uniform(r);
    } while (u1 <= 0.0);
    u2 = ir_rand_uniform(r);
    mag = sqrt(-2.0 * log(u1));
    r->spare = mag * sin(2.0 * IR_PI * u2);
    r->have_spare = 1;
    return mag * cos(2.0 * IR_PI * u2);
}

/* ------------------------------------------------------------------ */
/* Statistics                                                          */
/* ------------------------------------------------------------------ */

double ir_mean(const double *x, size_t n)
{
    double s = 0.0;
    size_t i;
    if (x == NULL || n == 0u) {
        return 0.0;
    }
    for (i = 0u; i < n; ++i) {
        s += x[i];
    }
    return s / (double)n;
}

double ir_stddev(const double *x, size_t n, double *mean_out)
{
    double m, s = 0.0;
    size_t i;
    if (x == NULL || n == 0u) {
        if (mean_out != NULL) {
            *mean_out = 0.0;
        }
        return 0.0;
    }
    /* Two pass: the mean is exact enough that the second pass does not
     * suffer catastrophic cancellation for the dynamic ranges used here. */
    m = ir_mean(x, n);
    for (i = 0u; i < n; ++i) {
        double d = x[i] - m;
        s += d * d;
    }
    if (mean_out != NULL) {
        *mean_out = m;
    }
    return sqrt(s / (double)n);
}

static int ir_cmp_double(const void *a, const void *b)
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

double ir_median_copy(const double *x, size_t n)
{
    double *tmp;
    double v;
    if (x == NULL || n == 0u) {
        return 0.0;
    }
    tmp = (double *)malloc(n * sizeof(double));
    if (tmp == NULL) {
        return 0.0;
    }
    memcpy(tmp, x, n * sizeof(double));
    qsort(tmp, n, sizeof(double), ir_cmp_double);
    v = (n % 2u == 1u) ? tmp[n / 2u]
                       : 0.5 * (tmp[n / 2u - 1u] + tmp[n / 2u]);
    free(tmp);
    return v;
}

double ir_percentile_copy(const double *x, size_t n, double p)
{
    double *tmp;
    double pos, frac, v;
    size_t lo, hi;
    if (x == NULL || n == 0u) {
        return 0.0;
    }
    if (p < 0.0) {
        p = 0.0;
    }
    if (p > 100.0) {
        p = 100.0;
    }
    tmp = (double *)malloc(n * sizeof(double));
    if (tmp == NULL) {
        return 0.0;
    }
    memcpy(tmp, x, n * sizeof(double));
    qsort(tmp, n, sizeof(double), ir_cmp_double);
    pos = (p / 100.0) * (double)(n - 1u);
    lo = (size_t)floor(pos);
    hi = (size_t)ceil(pos);
    if (hi >= n) {
        hi = n - 1u;
    }
    frac = pos - (double)lo;
    v = tmp[lo] * (1.0 - frac) + tmp[hi] * frac;
    free(tmp);
    return v;
}

double ir_min(const double *x, size_t n)
{
    double m;
    size_t i;
    if (x == NULL || n == 0u) {
        return 0.0;
    }
    m = x[0];
    for (i = 1u; i < n; ++i) {
        if (x[i] < m) {
            m = x[i];
        }
    }
    return m;
}

double ir_max(const double *x, size_t n)
{
    double m;
    size_t i;
    if (x == NULL || n == 0u) {
        return 0.0;
    }
    m = x[0];
    for (i = 1u; i < n; ++i) {
        if (x[i] > m) {
            m = x[i];
        }
    }
    return m;
}

double ir_rmse(const double *a, const double *b, size_t n)
{
    double s = 0.0;
    size_t i;
    if (a == NULL || b == NULL || n == 0u) {
        return 0.0;
    }
    for (i = 0u; i < n; ++i) {
        double d = a[i] - b[i];
        s += d * d;
    }
    return sqrt(s / (double)n);
}

double ir_mean_abs_diff(const double *a, const double *b, size_t n)
{
    double s = 0.0;
    size_t i;
    if (a == NULL || b == NULL || n == 0u) {
        return 0.0;
    }
    for (i = 0u; i < n; ++i) {
        s += fabs(a[i] - b[i]);
    }
    return s / (double)n;
}

double ir_nu_percent(const double *x, size_t n, double *mean_out, double *sd_out)
{
    double m, sd;
    if (x == NULL || n == 0u) {
        if (mean_out != NULL) {
            *mean_out = 0.0;
        }
        if (sd_out != NULL) {
            *sd_out = 0.0;
        }
        return 0.0;
    }
    sd = ir_stddev(x, n, &m);
    if (mean_out != NULL) {
        *mean_out = m;
    }
    if (sd_out != NULL) {
        *sd_out = sd;
    }
    if (fabs(m) < 1e-12) {
        return 0.0;
    }
    return 100.0 * sd / fabs(m);
}

/* ------------------------------------------------------------------ */
/* File helpers                                                        */
/* ------------------------------------------------------------------ */

FILE *ir_fopen_w(const char *path)
{
    if (path == NULL) {
        return NULL;
    }
    return fopen(path, "wb");
}

int ir_write_text(const char *path, const char *text)
{
    FILE *f;
    size_t len, wrote;
    if (path == NULL || text == NULL) {
        return IR_ERR_PARAM;
    }
    f = ir_fopen_w(path);
    if (f == NULL) {
        return IR_ERR_IO;
    }
    len = strlen(text);
    wrote = fwrite(text, 1u, len, f);
    if (fclose(f) != 0) {
        return IR_ERR_IO;
    }
    return (wrote == len) ? IR_OK : IR_ERR_IO;
}

void ir_fprintf(FILE *f, const char *fmt, ...)
{
    va_list ap;
    if (f == NULL || fmt == NULL) {
        return;
    }
    va_start(ap, fmt);
    (void)vfprintf(f, fmt, ap);
    va_end(ap);
}
