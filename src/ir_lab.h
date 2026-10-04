/*
 * ir_lab.h -- shared numerical helpers for ir-core-lab
 * ---------------------------------------------------
 * Pure C99, no third-party dependencies (libc + libm only).
 *
 * This header holds the small amount of infrastructure that every module
 * needs: a deterministic pseudo random generator, robust statistics and a
 * couple of file helpers.  Everything is written so that a run of this
 * project is bit-for-bit reproducible on every platform:
 *
 *   * the RNG is an integer xorshift, never rand()/srand();
 *   * statistics are computed with a plain two-pass algorithm (mean first,
 *     then the sum of squared deviations) which is accurate enough for the
 *     array sizes used here and keeps -ffast-math style surprises away;
 *   * no M_PI (strict C99 hides it), no non-standard functions.
 *
 * NOTE ON COMPILATION: all translation units are built with
 * -ffp-contract=off.  Without it clang on macOS may fuse a*b+c into a single
 * FMA instruction, which changes the last bit of many results and makes the
 * floating point assertions in test/test_ir.c fail on that platform only.
 */

#ifndef IR_LAB_H
#define IR_LAB_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/* Constants                                                           */
/* ------------------------------------------------------------------ */

#ifndef IR_PI
#define IR_PI 3.14159265358979323846
#endif

/* Return values used across the project. */
#define IR_OK           0
#define IR_ERR_PARAM   (-1)   /* NULL pointer / nonsensical argument   */
#define IR_ERR_RANGE   (-2)   /* value outside the supported range     */
#define IR_ERR_IO      (-3)   /* file could not be opened or written   */

/* ------------------------------------------------------------------ */
/* Deterministic pseudo random numbers                                 */
/* ------------------------------------------------------------------ */

/*
 * xorshift32: state must be non-zero.  Period 2^32-1.
 *   x ^= x << 13; x ^= x >> 17; x ^= x << 5;
 * It is *not* cryptographic, it is only used to build reproducible noise.
 */
typedef struct {
    uint32_t s;
    int      have_spare;   /* Box-Muller generates two normals at a time */
    double   spare;
} ir_rand_t;

void     ir_rand_seed(ir_rand_t *r, uint32_t seed);
uint32_t ir_rand_u32(ir_rand_t *r);
/* Uniform in [0, 1).  Uses the top 24 bits so the result is exactly
 * representable as a float-friendly double. */
double   ir_rand_uniform(ir_rand_t *r);
/* Standard normal, Box-Muller transform, spare value cached. */
double   ir_rand_gauss(ir_rand_t *r);

/* ------------------------------------------------------------------ */
/* Statistics                                                          */
/* ------------------------------------------------------------------ */

double ir_mean(const double *x, size_t n);
/* Population standard deviation (divide by n).  mean_out may be NULL. */
double ir_stddev(const double *x, size_t n, double *mean_out);
/* Sorts a private copy; x is left untouched. */
double ir_median_copy(const double *x, size_t n);
/* Sorts a private copy; p in [0,100].  Nearest-rank style interpolation. */
double ir_percentile_copy(const double *x, size_t n, double p);
double ir_min(const double *x, size_t n);
double ir_max(const double *x, size_t n);
/* RMS of (a[i]-b[i]). */
double ir_rmse(const double *a, const double *b, size_t n);
/* Mean absolute difference. */
double ir_mean_abs_diff(const double *a, const double *b, size_t n);

/*
 * Non-uniformity (NU) of an array, the figure of merit used by the
 * infrared community for fixed pattern noise:
 *
 *      NU [%] = 100 * sigma / |mean|          (sigma = population std dev)
 *
 * It is a pure spatial statistic: NU is evaluated on a *uniform* scene, so
 * anything left in the array is fixed pattern noise of the focal plane.
 * mean_out / sd_out may be NULL.
 */
double ir_nu_percent(const double *x, size_t n, double *mean_out, double *sd_out);

/* ------------------------------------------------------------------ */
/* Small file helpers                                                  */
/* ------------------------------------------------------------------ */

/*
 * All output is written in BINARY mode ("wb").  On Windows a stream opened in
 * text mode treats 0x1A (Ctrl-Z) as end of file and rewrites "\n" as "\r\n",
 * so a raw byte that happens to be 0x1A silently truncates the file.  Every
 * writer in this project therefore opens with "wb", including plain text.
 */
FILE *ir_fopen_w(const char *path);
/* Writes a NUL terminated string, returns IR_OK / IR_ERR_IO. */
int   ir_write_text(const char *path, const char *text);
/* Appends fmt... to an already opened binary stream, "\n" terminated. */
void  ir_fprintf(FILE *f, const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif /* IR_LAB_H */
