#include <complex.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int update_cov(const int m_max, const int N, const int M, 
	       const double norm_re, const double norm_im, 
	       const double llp1_re, const double llp1_im,
	       const double complex *restrict F, const double complex *restrict H,
	       const double *restrict tau_power, const double *restrict eta_ratio2,
	       double *restrict cov, double *restrict cross_cov)
{
  /*
    Update the covariance (cov) and cross-covariance matrices with a single term 
    from the partial fraction decomposition (i.e. single norm and l(l+1)).

    m_max     - final m term (i.e. 0...m_max)
    N         - Number of z points
    M         - number of terms in partial fraction decomposition
    norm_re   - real(norm)
    norm_im   - imag(norm)
    llp1_re   - real(l(l+1))
    llp1_im   - imag(l(l+1))
    F         - (m_max+M+1, N) complex array of 2_F_1 (l,l+1;1+m;(1-z)/2)
    H         - (m_max+M+1, N) complex array of 2_F_1 (l,l+1;1+m;(1+z)/2)
    tau_power - (2M-1, N) double array of tau^p, with p 
                from 1-M...M-1 where tau_i = ((1-z_i)/(1+z_i))^(1/2)
    eta_ratio - (N-1,) double array of tau_{i+1} / tau_i
    updating:

    cov       - (m_max+1, N, M, M) array of covariances
    cross_cov - (m_max+1, N-1, M, M) array of cross-covariances

    returns 0 on success, -1 on out-of-memory
    
   */
  const int n_cross = N-1;
  double complex zeta[M], Rc[M], zetaA[M], zAH[M];
  double complex norm = norm_re + I*norm_im;
  const double complex llp1 = llp1_re + I*llp1_im;

  double *eta = (double *)malloc(n_cross * sizeof(double));

  if (!eta)
    return -1;

  for (int i=0;i<n_cross;++i)
    eta[i] = 1.0;

  for (int m=0,mN=0;m<=m_max;++m, mN+=N)
    {
	  
      // zeta symbols
      zeta[0] = 1.0;
      zeta[1] = m -llp1/(m+1); // First zeta

      for (int i=2,v=m+2;i<M;++i,++v)
	{
	  const double complex mult = (v*(v-1) - llp1)/v;
          zeta[i] = zeta[i-1]*mult;
	}



      for (int p=0;p<M;++p)
	zetaA[p] = (p&1) ? -zeta[p] : zeta[p];

      for (int p=0, t_idx = (M-1)*N;p<M;++p, t_idx-=N)
	zAH[p] = zetaA[p] * H[(m+p)*N] * tau_power[t_idx];

      
      // Pre-multiply LHS by the norm
      for (int p=0;p<M;++p)
	zeta[p] *= norm;

      // First covariance
      for (int p=0,F_idx=mN,t_idx=(M-1)*N;p<M;++p,F_idx+=N, t_idx+=N)
	{
	  const double complex Lp = tau_power[t_idx] * (zeta[p] * F[F_idx]); 
      
	  for (int q=0;q<M;++q)
	    {
	      const double complex Rq = zAH[q];
	      // Update real part
	      *cov++ += creal(Lp) * creal(Rq) - cimag(Lp)*cimag(Rq);
	    }
	}

      // Covariance and cross
      for (int n=0,tau0=(M-1)*N+1;n<n_cross;++n, tau0++)
	{

	  for (int p=0, t_idx = tau0;p<M;++p, t_idx-=N)
	    {
	      Rc[p] = eta[n] * zAH[p];
	      zAH[p] = zetaA[p] * H[(m+p)*N+n+1] * tau_power[t_idx];
	    }	    


	  for (int p=0,F_idx=mN+n+1,t_idx=tau0;p<M;++p,F_idx+=N, t_idx+=N)
	    {
	      const double complex Lp = tau_power[t_idx] * (zeta[p] * F[F_idx]); 
	      for (int q=0;q<M;++q)
		{
		  // Update with real parts
		  *cov++ += creal(Lp) * creal(zAH[q]) - cimag(Lp)*cimag(zAH[q]);
		  *cross_cov++ += creal(Lp)*creal(Rc[q]) - cimag(Lp)*cimag(Rc[q]);
		}
	    }
	  
	  eta[n] *= eta_ratio2[n]; //      eta *= tau[1:] / tau[:-1]
	}

      norm *= ((m+1)*m - llp1)/((m+1)*(m+1)); // norm for next m
      
    }

  free(eta);
  return 0;
}

int update_cov_range(const int m_lo, const int m_hi, const int N, const int M,
                     const double norm_re, const double norm_im,
                     const double llp1_re, const double llp1_im,
                     const double complex *restrict F, const double complex *restrict H,
                     const double *restrict tau_power, const double *restrict eta_ratio2,
                     double *restrict cov, double *restrict cross_cov)
{
  const int n_cross = N-1;
  double complex zeta[M], Rc[M], zetaA[M], zAH[M];
  double complex norm = norm_re + I*norm_im;
  const double complex llp1 = llp1_re + I*llp1_im;

  double *eta = (double *)malloc(n_cross * sizeof(double));
  if (!eta) return -1;

  for (int i=0; i<n_cross; ++i)
    eta[i] = 1.0;

  for (int m=0; m<=m_hi; ++m) {
    if (m >= m_lo) {
      zeta[0] = 1.0;
      zeta[1] = m - llp1/(m+1);

      for (int i=2, v=m+2; i<M; ++i, ++v) {
        const double complex mult = (v*(v-1) - llp1)/v;
        zeta[i] = zeta[i-1]*mult;
      }

      for (int p=0; p<M; ++p)
        zetaA[p] = (p&1) ? -zeta[p] : zeta[p];

      for (int p=0, t_idx = (M-1)*N; p<M; ++p, t_idx-=N)
        zAH[p] = zetaA[p] * H[(m+p-m_lo)*N] * tau_power[t_idx];

      for (int p=0; p<M; ++p)
        zeta[p] *= norm;

      for (int p=0, F_idx=(m-m_lo)*N, t_idx=(M-1)*N; p<M; ++p, F_idx+=N, t_idx+=N) {
        const double complex Lp = tau_power[t_idx] * (zeta[p] * F[F_idx]);
        for (int q=0; q<M; ++q) {
          const double complex Rq = zAH[q];
          *cov++ += creal(Lp) * creal(Rq) - cimag(Lp)*cimag(Rq);
        }
      }

      for (int n=0, tau0=(M-1)*N+1; n<n_cross; ++n, tau0++) {
        for (int p=0, t_idx = tau0; p<M; ++p, t_idx-=N) {
          Rc[p] = eta[n] * zAH[p];
          zAH[p] = zetaA[p] * H[(m+p-m_lo)*N+n+1] * tau_power[t_idx];
        }

        for (int p=0, F_idx=(m-m_lo)*N+n+1, t_idx=tau0; p<M; ++p, F_idx+=N, t_idx+=N) {
          const double complex Lp = tau_power[t_idx] * (zeta[p] * F[F_idx]);
          for (int q=0; q<M; ++q) {
            *cov++ += creal(Lp) * creal(zAH[q]) - cimag(Lp)*cimag(zAH[q]);
            *cross_cov++ += creal(Lp)*creal(Rc[q]) - cimag(Lp)*cimag(Rc[q]);
          }
        }
      }
    }

    for (int n=0; n<n_cross; ++n)
      eta[n] *= eta_ratio2[n];
    norm *= ((m+1)*m - llp1)/((m+1)*(m+1));
  }

  free(eta);
  return 0;
}

/*
  cov_legendre
  ============
  Compute Mord=1 scalar covariance and cross-covariance for all m-modes
  (0..m_max) at all z-points using a direct associated Legendre series.

  Background
  ----------
  The original SMERFS pipeline uses Mord=2: the state vector is (g_m, g_m')
  where g_m' encodes an angular-derivative companion.  The derivative
  component's covariance involves zeta[1] = m - l(l+1)/(m+1).  For the
  power spectrum C_l = 1/(c0 + c2*(l(l+1))^2) the partial-fraction poles are
  l(l+1) = ±i*sqrt(c0/c2), so zeta[1] ~ ±i*sqrt(c0/c2) for m=0.
  When c2 is small (short correlation lengths, e.g. xcorr ~ 1 deg where
  c2 ~ 7e-8 and |Im(l(l+1))| ~ 3635) the 2x2 covariance matrix has
  condition number ~ R^4 ~ 1.7e14, making its Cholesky decomposition
  numerically singular.  Simultaneously, hyp_llp1 needs O(R) ~ O(3635)
  series terms to converge, hitting MAX_ITERATIONS (10000) before finishing.

  The fix: use Mord=1 (scalar covariance) and compute it directly as

      K_m(z_i, z_j) = sum_{l=m}^{lmax} C_l * phi_lm(z_i) * phi_lm(z_j)

  where phi_lm are 4pi-normalised associated Legendre basis functions:
      phi_l0(z)  = sqrt((2l+1)/(4*pi))           * P_l^0(z)
      phi_lm(z)  = sqrt((2l+1)/(4*pi)*2*(l-m)!/(l+m)!) * P_l^m(z)  for m>0

  so that sum_{m=0}^{l} phi_lm(z)^2 = (2l+1)/(4*pi)  (addition theorem)
  and hence sum_m K_m(z,z) = sum_l C_l*(2l+1)/(4pi) = C(0).

  phi_lm is computed via the standard 3-term recurrence that is stable for
  all l, m in double precision.  Only real arithmetic is used.

  Mord=1 covariances are always positive (sum of squares) so the Cholesky
  in state_space always succeeds.

  Arguments
  ---------
  m_max    (in)  : max m-mode (0..m_max), equals n_m-1 = nphi/2
  N        (in)  : number of z-points (= uhalf = nz/2+1)
  lmax     (in)  : Legendre truncation (>= m_max, typically 3*nz)
  c0, c2   (in)  : power spectrum C_l = 1/(c0 + c2*(l(l+1))^2)
  z_pts    (in)  : (N) cos(theta) values, equator-first (z_pts[0]~0, z_pts[N-1]~1)
  cov      (out) : (m_max+1)*N  array; cov[m*N+i]      = K_m(z_i, z_i)
  cross_cov(out) : (m_max+1)*(N-1) array; cross_cov[m*(N-1)+i] = K_m(z_i, z_{i+1})

  Returns 0 on success, -1 on allocation failure.
*/
int cov_legendre(const int m_max, const int N, const int lmax,
                 const double c0, const double c2,
                 const double *restrict z_pts,
                 double *restrict cov,
                 double *restrict cross_cov)
{
  const int N1 = N - 1;
  const int nm = m_max + 1;

  /* Working arrays: two rows (prev, curr) per m-mode, each of length N.
   * phi_prev[m*N .. m*N+N-1] = phi_{l-1}^m(z_i)
   * phi_curr[m*N .. m*N+N-1] = phi_l^m(z_i)
   * These are updated in place as l advances.
   */
  double *phi_prev = (double *)calloc(nm * N, sizeof(double));
  double *phi_curr = (double *)calloc(nm * N, sizeof(double));
  if (!phi_prev || !phi_curr) { free(phi_prev); free(phi_curr); return -1; }

  /* sin(theta) at each z-point */
  double *sin_t = (double *)malloc(N * sizeof(double));
  if (!sin_t) { free(phi_prev); free(phi_curr); return -1; }
  for (int i = 0; i < N; ++i) {
    double s2 = 1.0 - z_pts[i] * z_pts[i];
    sin_t[i] = (s2 > 0.0) ? sqrt(s2) : 0.0;
  }

  /* Zero output */
  for (int k = 0; k < nm * N;  ++k) cov[k]       = 0.0;
  for (int k = 0; k < nm * N1; ++k) cross_cov[k] = 0.0;

  /* ------------------------------------------------------------------ *
   * Initialise phi_prev[:,m] = phi_m^m and phi_curr[:,m] = phi_{m+1}^m
   * for each m = 0..m_max.
   *
   * Fully-normalised phi_m^m(cos theta):
   *   phi_m^m = (-1)^m * sqrt((2m+1)/(4*pi)) * sqrt(prod_{k=1}^m (2k-1)/(2k)) * sin^m(theta)
   *
   * We work in log-space to avoid overflow for large m, then exponentiate.
   * The sign (-1)^m is dropped because it cancels in phi^2 and phi_i*phi_j.
   * ------------------------------------------------------------------ */
  for (int m = 0; m <= m_max; ++m) {
    /* log|phi_m^m| = 0.5*log((2m+1)/(4pi)) + sum_{k=1}^m 0.5*log((2k-1)/(2k))
     * + m * log(sin(theta))                                                    */
    double log_c = 0.5 * log((2*m + 1) / (4.0 * M_PI));
    for (int k = 1; k <= m; ++k)
      log_c += 0.5 * log((2.0*k - 1.0) / (2.0*k));

    double *Pmm   = &phi_prev[m * N];
    double *Pmp1m = &phi_curr[m * N];

    for (int i = 0; i < N; ++i) {
      double st = sin_t[i];
      Pmm[i] = (st > 0.0) ? exp(log_c + m * log(st)) : 0.0;
      /* phi_{m+1}^m = sqrt(2m+3) * z * phi_m^m  */
      Pmp1m[i] = sqrt(2.0*m + 3.0) * z_pts[i] * Pmm[i];
    }
  }

  /* ------------------------------------------------------------------ *
   * Walk l from 0 to lmax.  For each l, accumulate contributions from
   * all m in [0, min(l, m_max)].
   * ------------------------------------------------------------------ */
  for (int l = 0; l <= lmax; ++l) {
    const double ll1 = (double)l * (double)(l + 1);
    const double Cl  = 1.0 / (c0 + c2 * ll1 * ll1);
    const int    m_hi = (l < m_max) ? l : m_max;

    for (int m = 0; m <= m_hi; ++m) {
      /* For l > m+1: advance the recurrence
       *   phi_l^m = a_lm * z * phi_{l-1}^m - b_lm * phi_{l-2}^m
       * with
       *   a_lm = sqrt((4l^2-1) / (l^2-m^2))
       *   b_lm = sqrt(((l-1)^2-m^2) * (4l^2-1) / ((l^2-m^2) * (4(l-1)^2-1)))
       * phi_prev[m*N..] plays the role of phi_{l-2}^m after the update.
       */
      if (l > m + 1) {
        const double l2  = (double)l * (double)l;
        const double lm1 = (double)(l - 1);
        const double m2  = (double)m * (double)m;
        const double a   = sqrt((4.0*l2 - 1.0) / (l2 - m2));
        const double b   = sqrt((lm1*lm1 - m2) * (4.0*l2 - 1.0)
                                / ((l2 - m2) * (4.0*lm1*lm1 - 1.0)));
        double *pp = &phi_prev[m * N];
        double *pc = &phi_curr[m * N];
        for (int i = 0; i < N; ++i) {
          double new_val = a * z_pts[i] * pc[i] - b * pp[i];
          pp[i] = pc[i];   /* phi_{l-1}^m <- old phi_curr */
          pc[i] = new_val; /* phi_l^m     <- new value    */
        }
      }

      /* phi_l^m is now in phi_curr[m*N..] */
      const double *Plm = (l == m) ? &phi_prev[m * N] : &phi_curr[m * N];
      /* Note: when l==m, phi_m^m is in phi_prev (initialised there);
       *       when l==m+1, phi_{m+1}^m is in phi_curr (also initialised).
       *       when l>m+1, we just advanced phi_curr to phi_l^m above.       */
      /* Correct pointer: phi_prev holds phi_{l-2}^m (or phi_m^m for l==m),
       * phi_curr holds phi_l^m (or phi_{m+1}^m for l==m+1).                */

      const double weight = Cl;   /* shape only; field is renormalized anyway */

      double *cov_m       = &cov[m * N];
      double *cross_cov_m = &cross_cov[m * N1];
      for (int i = 0; i < N;  ++i) cov_m[i]       += weight * Plm[i] * Plm[i];
      for (int i = 0; i < N1; ++i) cross_cov_m[i]  += weight * Plm[i] * Plm[i+1];
    }
  }

  free(phi_prev); free(phi_curr); free(sin_t);
  return 0;
}
