/*
 * smerfs_fenv.c
 *
 * Thin wrappers around the public SMERFS C functions that save and restore
 * the floating-point environment across each call.
 *
 * Background
 * ----------
 * The hypergeometric series in smerfs.c produces mathematically benign
 * intermediate NaN/Inf values that resolve before the function returns:
 *
 *   - Divide-by-zero in the AMS55 15.3.11 finite sum (denominator (n+1-m)
 *     hits zero at n=m-1, but the convergence check exits the loop before
 *     that term is ever computed for well-converged series).
 *   - Underflow of s^m (s=1-z near 0, m up to n_m-1 ~ 128).
 *   - Invalid (NaN) from Inf*0 mixing the above two.
 *
 * When called from a NAG Fortran program the NAG runtime installs a SIGFPE
 * handler that traps these hardware exceptions and prints warnings (or aborts
 * with -ieee=stop).  We use POSIX feholdexcept / fesetenv to mask all FP
 * exceptions for the duration of each C kernel call, then restore the
 * caller's FP environment on return.  This eliminates all three warnings:
 *
 *   Warning: Floating invalid operation occurred
 *   Warning: Floating divide by zero occurred
 *   Warning: Floating underflow occurred
 *
 * The wrappers are named with the suffix _w (for "wrapped") and are called
 * by the Fortran interface module via bind(C, name=...).
 */

#include <fenv.h>
#include <stdint.h>
#include <complex.h>

/* ------------------------------------------------------------------ */
/* Forward-declare the real implementations from the other C files     */
/* ------------------------------------------------------------------ */
int hyp_llp1  (const double, const double, const int, const int,
               const double *, double complex *);
double complex hyp_lmz (const double complex, const int, const double,
               const double complex, const double complex);
int update_cov(const int, const int, const int,
               const double, const double, const double, const double,
               const double complex *, const double complex *,
               const double *, const double *,
               double *, double *);
int inverse   (const int, const int, const double *, double *);
int cholesky  (const int, const int, const double *, double *);
int state_space(const int, const int,
               const double *, const double *,
               double *, double *);
int zigg      (const int, int, const uint32_t *, float *);
int cov_legendre(const int, const int, const int,
                 const double, const double,
                 const double *, double *, double *);
int state_space_1(const int,
                  const double *, const double *,
                  double *, double *);

/* ------------------------------------------------------------------ */
/* Convenience macro: save FP env, call, clear exceptions, restore    */
/* ------------------------------------------------------------------ */
#define FENV_WRAP_BEGIN  fenv_t _fenv; feholdexcept(&_fenv);
#define FENV_WRAP_END    feclearexcept(FE_ALL_EXCEPT); fesetenv(&_fenv);

/* ------------------------------------------------------------------ */
/* Wrapped entry points (bound to the original C names in the Fortran  */
/* interface via bind(C, name='...'))                                  */
/* ------------------------------------------------------------------ */

int hyp_llp1_w(const double r, const double im, const int m, const int nz,
               const double *zvals, double complex *out)
{
    FENV_WRAP_BEGIN
    int rc = hyp_llp1(r, im, m, nz, zvals, out);
    FENV_WRAP_END
    return rc;
}

double complex hyp_lmz_w(const double complex l, const int m, const double z,
                          const double complex gamma_ratio0,
                          const double complex psi0)
{
    FENV_WRAP_BEGIN
    double complex res = hyp_lmz(l, m, z, gamma_ratio0, psi0);
    FENV_WRAP_END
    return res;
}

int update_cov_w(const int m_max, const int N, const int M,
                 const double norm_re, const double norm_im,
                 const double llp1_re, const double llp1_im,
                 const double complex *F, const double complex *H,
                 const double *tau_power, const double *eta_ratio2,
                 double *cov, double *cross_cov)
{
    FENV_WRAP_BEGIN
    int rc = update_cov(m_max, N, M,
                        norm_re, norm_im, llp1_re, llp1_im,
                        F, H, tau_power, eta_ratio2,
                        cov, cross_cov);
    FENV_WRAP_END
    return rc;
}

int inverse_w(const int N, const int M,
              const double *matrices, double *out)
{
    FENV_WRAP_BEGIN
    int rc = inverse(N, M, matrices, out);
    FENV_WRAP_END
    return rc;
}

int cholesky_w(const int N, const int M,
               const double *matrices, double *out)
{
    FENV_WRAP_BEGIN
    int rc = cholesky(N, M, matrices, out);
    FENV_WRAP_END
    return rc;
}

int state_space_w(const int N, const int M,
                  const double *cross_cov, const double *cov,
                  double *innov, double *trans)
{
    FENV_WRAP_BEGIN
    int rc = state_space(N, M, cross_cov, cov, innov, trans);
    FENV_WRAP_END
    return rc;
}

int zigg_w(const int num_needed, int num_ints,
           const uint32_t *rand_ints, float *out)
{
    FENV_WRAP_BEGIN
    int rc = zigg(num_needed, num_ints, rand_ints, out);
    FENV_WRAP_END
    return rc;
}

int cov_legendre_w(const int m_max, const int N, const int lmax,
                   const double c0, const double c2,
                   const double *z_pts, double *cov, double *cross_cov)
{
    /* No FP exceptions expected (pure real arithmetic, no denormals for
     * well-behaved inputs), but wrap anyway for consistency.           */
    FENV_WRAP_BEGIN
    int rc = cov_legendre(m_max, N, lmax, c0, c2, z_pts, cov, cross_cov);
    FENV_WRAP_END
    return rc;
}

int state_space_1_w(const int N,
                    const double *cross_cov, const double *cov,
                    double *innov, double *trans)
{
    FENV_WRAP_BEGIN
    int rc = state_space_1(N, cross_cov, cov, innov, trans);
    FENV_WRAP_END
    return rc;
}
