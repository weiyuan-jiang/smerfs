! Fortran 90 interface module for the SMERFS C library
!
! SMERFS - Stochastic Markov Evaluation of Random Fields on the Sphere
! (Peter Creasey)
!
! Public C functions wrapped here:
!   cov.c      : update_cov
!   linalg.c   : inverse, cholesky, state_space
!   smerfs.c   : hyp_llp1, hyp_lmz
!   ziggurat.c : zigg

module smerfs_interface
  use, intrinsic :: iso_c_binding, only : &
      c_int, c_int32_t, c_double, c_float, c_double_complex
  implicit none

  private
  public :: update_cov_f, update_cov_range_f, inverse_f, cholesky_f, state_space_f, &
            hyp_llp1_f, hyp_lmz_f, zigg_f, &
            cov_legendre_f, state_space_1_f

  ! ------------------------------------------------------------------
  ! Raw C bindings (private; called only through the wrappers below)
  ! ------------------------------------------------------------------
  interface

    ! ------------------------------------------------------------------
    ! cov.c
    !
    ! int update_cov(
    !   const int m_max, const int N, const int M,
    !   const double norm_re, const double norm_im,
    !   const double llp1_re, const double llp1_im,
    !   const double complex *F, const double complex *H,
    !   const double *tau_power, const double *eta_ratio2,
    !   double *cov, double *cross_cov)
    ! ------------------------------------------------------------------
    function c_update_cov(m_max, N, M, &
                          norm_re, norm_im, &
                          llp1_re, llp1_im, &
                          F, H, tau_power, eta_ratio2, &
                          cov, cross_cov) &
        bind(C, name='update_cov_w') result(rc)
      import :: c_int, c_double, c_double_complex
      integer(c_int), value, intent(in) :: m_max, N, M
      real(c_double), value, intent(in) :: norm_re, norm_im
      real(c_double), value, intent(in) :: llp1_re, llp1_im
      complex(c_double_complex), intent(in)    :: F(*)
      complex(c_double_complex), intent(in)    :: H(*)
      real(c_double),            intent(in)    :: tau_power(*)
      real(c_double),            intent(in)    :: eta_ratio2(*)
      real(c_double),            intent(inout) :: cov(*)
      real(c_double),            intent(inout) :: cross_cov(*)
      integer(c_int) :: rc
    end function c_update_cov

    function c_update_cov_range(m_lo, m_hi, N, M, &
                                norm_re, norm_im, &
                                llp1_re, llp1_im, &
                                F, H, tau_power, eta_ratio2, &
                                cov, cross_cov) &
        bind(C, name='update_cov_range_w') result(rc)
      import :: c_int, c_double, c_double_complex
      integer(c_int), value, intent(in) :: m_lo, m_hi, N, M
      real(c_double), value, intent(in) :: norm_re, norm_im
      real(c_double), value, intent(in) :: llp1_re, llp1_im
      complex(c_double_complex), intent(in)    :: F(*)
      complex(c_double_complex), intent(in)    :: H(*)
      real(c_double),            intent(in)    :: tau_power(*)
      real(c_double),            intent(in)    :: eta_ratio2(*)
      real(c_double),            intent(inout) :: cov(*)
      real(c_double),            intent(inout) :: cross_cov(*)
      integer(c_int) :: rc
    end function c_update_cov_range

    ! ------------------------------------------------------------------
    ! linalg.c
    !
    ! int inverse(const int N, const int M,
    !             const double *matrices, double *out)
    ! ------------------------------------------------------------------
    function c_inverse(N, M, matrices, out) &
        bind(C, name='inverse_w') result(rc)
      import :: c_int, c_double
      integer(c_int), value, intent(in) :: N, M
      real(c_double),        intent(in)  :: matrices(*)
      real(c_double),        intent(out) :: out(*)
      integer(c_int) :: rc
    end function c_inverse

    ! ------------------------------------------------------------------
    ! int cholesky(const int N, const int M,
    !              const double *matrices, double *out)
    ! ------------------------------------------------------------------
    function c_cholesky(N, M, matrices, out) &
        bind(C, name='cholesky_w') result(rc)
      import :: c_int, c_double
      integer(c_int), value, intent(in) :: N, M
      real(c_double),        intent(in)  :: matrices(*)
      real(c_double),        intent(out) :: out(*)
      integer(c_int) :: rc
    end function c_cholesky

    ! ------------------------------------------------------------------
    ! int state_space(const int N, const int M,
    !                 const double *cross_cov, const double *cov,
    !                 double *innov, double *trans)
    ! ------------------------------------------------------------------
    function c_state_space(N, M, cross_cov, cov, innov, trans) &
        bind(C, name='state_space_w') result(rc)
      import :: c_int, c_double
      integer(c_int), value, intent(in) :: N, M
      real(c_double),        intent(in)  :: cross_cov(*)
      real(c_double),        intent(in)  :: cov(*)
      real(c_double),        intent(out) :: innov(*)
      real(c_double),        intent(out) :: trans(*)
      integer(c_int) :: rc
    end function c_state_space

    ! ------------------------------------------------------------------
    ! smerfs.c
    !
    ! int hyp_llp1(const double llp1_real, const double llp1_imag,
    !              const int m, const int nz,
    !              const double *zvals, double complex *out)
    ! ------------------------------------------------------------------
    function c_hyp_llp1(llp1_real, llp1_imag, m, nz, zvals, out) &
        bind(C, name='hyp_llp1_w') result(rc)
      import :: c_int, c_double, c_double_complex
      real(c_double), value, intent(in) :: llp1_real, llp1_imag
      integer(c_int), value, intent(in) :: m, nz
      real(c_double),            intent(in)  :: zvals(*)
      complex(c_double_complex), intent(out) :: out(*)
      integer(c_int) :: rc
    end function c_hyp_llp1

    ! ------------------------------------------------------------------
    ! double complex hyp_lmz(const double complex l, const int m,
    !                        const double z,
    !                        const double complex gamma_ratio0,
    !                        const double complex psi0)
    !
    ! C passes complex scalars by value; Fortran VALUE attribute matches.
    ! ------------------------------------------------------------------
    function c_hyp_lmz(l, m, z, gamma_ratio0, psi0) &
        bind(C, name='hyp_lmz_w') result(res)
      import :: c_int, c_double, c_double_complex
      complex(c_double_complex), value, intent(in) :: l
      integer(c_int),            value, intent(in) :: m
      real(c_double),            value, intent(in) :: z
      complex(c_double_complex), value, intent(in) :: gamma_ratio0
      complex(c_double_complex), value, intent(in) :: psi0
      complex(c_double_complex) :: res
    end function c_hyp_lmz

    ! ------------------------------------------------------------------
    ! ziggurat.c
    !
    ! int zigg(const int num_needed, int num_ints,
    !          const uint32_t *rand_ints, float *out)
    ! ------------------------------------------------------------------
    function c_zigg(num_needed, num_ints, rand_ints, out) &
        bind(C, name='zigg_w') result(rc)
      import :: c_int, c_int32_t, c_float
      integer(c_int),    value, intent(in) :: num_needed, num_ints
      integer(c_int32_t),       intent(in)  :: rand_ints(*)
      real(c_float),            intent(out) :: out(*)
      integer(c_int) :: rc
    end function c_zigg

    ! ------------------------------------------------------------------
    ! cov.c
    !
    ! int cov_legendre(const int m_max, const int N, const int lmax,
    !                  const double c0, const double c2,
    !                  const double *z_pts,
    !                  double *cov, double *cross_cov)
    ! ------------------------------------------------------------------
    function c_cov_legendre(m_max, N, lmax, c0, c2, z_pts, cov, cross_cov) &
        bind(C, name='cov_legendre') result(rc)
      import :: c_int, c_double
      integer(c_int), value, intent(in) :: m_max, N, lmax
      real(c_double), value, intent(in) :: c0, c2
      real(c_double),        intent(in)  :: z_pts(*)
      real(c_double),        intent(out) :: cov(*)
      real(c_double),        intent(out) :: cross_cov(*)
      integer(c_int) :: rc
    end function c_cov_legendre

    ! ------------------------------------------------------------------
    ! linalg.c
    !
    ! int state_space_1(const int N,
    !                   const double *cross_cov, const double *cov,
    !                   double *innov, double *trans)
    ! ------------------------------------------------------------------
    function c_state_space_1(N, cross_cov, cov, innov, trans) &
        bind(C, name='state_space_1_w') result(rc)
      import :: c_int, c_double
      integer(c_int), value, intent(in) :: N
      real(c_double),        intent(in)  :: cross_cov(*)
      real(c_double),        intent(in)  :: cov(*)
      real(c_double),        intent(out) :: innov(*)
      real(c_double),        intent(out) :: trans(*)
      integer(c_int) :: rc
    end function c_state_space_1

  end interface

contains

  ! ==================================================================
  ! update_cov_f
  ! ==================================================================
  ! Update covariance and cross-covariance matrices with a single term
  ! from the partial fraction decomposition.
  !
  !   m_max      : final m index (0 .. m_max)
  !   N          : number of z points
  !   M          : number of terms in partial fraction decomposition
  !   norm_re/im : real and imaginary parts of the norm
  !   llp1_re/im : real and imaginary parts of l(l+1)
  !   F          : (m_max+M+1, N)     complex  2F1(l,l+1;1+m;(1-z)/2)
  !   H          : (m_max+M+1, N)     complex  2F1(l,l+1;1+m;(1+z)/2)
  !   tau_power  : (2M-1, N)          real     tau^p, p in 1-M..M-1
  !   eta_ratio2 : (N-1)              real     tau_{i+1}/tau_i
  !   cov        : (m_max+1, N, M, M) real     covariances (in/out)
  !   cross_cov  : (m_max+1,N-1,M,M) real     cross-covariances (in/out)
  !   rc         : 0 = success, -1 = out-of-memory
  subroutine update_cov_f(m_max, N, M, &
                           norm_re, norm_im, &
                           llp1_re, llp1_im, &
                           F, H, tau_power, eta_ratio2, &
                           cov, cross_cov, rc)
    integer,                   intent(in)    :: m_max, N, M
    real(c_double),            intent(in)    :: norm_re, norm_im
    real(c_double),            intent(in)    :: llp1_re, llp1_im
    complex(c_double_complex), intent(in)    :: F((m_max+M+1)*N)
    complex(c_double_complex), intent(in)    :: H((m_max+M+1)*N)
    real(c_double),            intent(in)    :: tau_power((2*M-1)*N)
    real(c_double),            intent(in)    :: eta_ratio2(N-1)
    real(c_double),            intent(inout) :: cov((m_max+1)*N*M*M)
    real(c_double),            intent(inout) :: cross_cov((m_max+1)*(N-1)*M*M)
    integer,                   intent(out)   :: rc

    rc = int( c_update_cov(int(m_max, c_int), int(N, c_int), int(M, c_int), &
                           norm_re, norm_im, llp1_re, llp1_im, &
                           F, H, tau_power, eta_ratio2, cov, cross_cov) )
  end subroutine update_cov_f

  subroutine update_cov_range_f(m_lo, m_hi, N, M, &
                                norm_re, norm_im, &
                                llp1_re, llp1_im, &
                                F, H, tau_power, eta_ratio2, &
                                cov, cross_cov, rc)
    integer,                   intent(in)    :: m_lo, m_hi, N, M
    real(c_double),            intent(in)    :: norm_re, norm_im
    real(c_double),            intent(in)    :: llp1_re, llp1_im
    complex(c_double_complex), intent(in)    :: F((m_hi-m_lo+1+M)*N)
    complex(c_double_complex), intent(in)    :: H((m_hi-m_lo+1+M)*N)
    real(c_double),            intent(in)    :: tau_power((2*M-1)*N)
    real(c_double),            intent(in)    :: eta_ratio2(N-1)
    real(c_double),            intent(inout) :: cov((m_hi-m_lo+1)*N*M*M)
    real(c_double),            intent(inout) :: cross_cov((m_hi-m_lo+1)*(N-1)*M*M)
    integer,                   intent(out)   :: rc

    rc = int( c_update_cov_range(int(m_lo, c_int), int(m_hi, c_int), &
                                 int(N, c_int), int(M, c_int), &
                                 norm_re, norm_im, llp1_re, llp1_im, &
                                 F, H, tau_power, eta_ratio2, cov, cross_cov) )
  end subroutine update_cov_range_f

  ! ==================================================================
  ! inverse_f
  ! ==================================================================
  ! Inverse of N symmetric M x M matrices (M = 2 or 3).
  !
  !   N        : number of matrices
  !   M        : matrix dimension (2 or 3)
  !   matrices : (N*M*M) input  real array
  !   out      : (N*M*M) output real array (inverses)
  !   rc       : 0 = success; i+1 = matrix i singular; -1 = M unsupported
  subroutine inverse_f(N, M, matrices, out, rc)
    integer,        intent(in)  :: N, M
    real(c_double), intent(in)  :: matrices(N*M*M)
    real(c_double), intent(out) :: out(N*M*M)
    integer,        intent(out) :: rc

    rc = int( c_inverse(int(N, c_int), int(M, c_int), matrices, out) )
  end subroutine inverse_f

  ! ==================================================================
  ! cholesky_f
  ! ==================================================================
  ! Batched lower Cholesky factorisation of N symmetric M x M matrices.
  !
  !   N        : number of matrices
  !   M        : matrix dimension (2 supported; 3 not yet in C code)
  !   matrices : (N*M*M) input  real array
  !   out      : (N*M*M) output real array (lower triangular factors)
  !   rc       : 0 = success; i+1 = matrix i not pos-def; -1 = M unsupported
  subroutine cholesky_f(N, M, matrices, out, rc)
    integer,        intent(in)  :: N, M
    real(c_double), intent(in)  :: matrices(N*M*M)
    real(c_double), intent(out) :: out(N*M*M)
    integer,        intent(out) :: rc

    rc = int( c_cholesky(int(N, c_int), int(M, c_int), matrices, out) )
  end subroutine cholesky_f

  ! ==================================================================
  ! state_space_f
  ! ==================================================================
  ! Construct Kalman-filter innovation and transition matrices.
  !
  !   N         : number of state-space steps
  !   M         : state dimension (2 or 3)
  !   cross_cov : ((N-1)*M*M) input  real array of cross-covariances
  !   cov       : (N*M*M)     input  real array of covariances
  !   innov     : (N*M*M)     output real array of innovation (Cholesky) matrices
  !   trans     : ((N-1)*M*M) output real array of transition matrices
  !   rc        : 0 = success; i+1 = step i failed; -1 = M unsupported
  subroutine state_space_f(N, M, cross_cov, cov, innov, trans, rc)
    integer,        intent(in)  :: N, M
    real(c_double), intent(in)  :: cross_cov((N-1)*M*M)
    real(c_double), intent(in)  :: cov(N*M*M)
    real(c_double), intent(out) :: innov(N*M*M)
    real(c_double), intent(out) :: trans((N-1)*M*M)
    integer,        intent(out) :: rc

    rc = int( c_state_space(int(N, c_int), int(M, c_int), &
                            cross_cov, cov, innov, trans) )
  end subroutine state_space_f

  ! ==================================================================
  ! hyp_llp1_f
  ! ==================================================================
  ! Evaluate 2F1(-l, l+1; 1+m; z) over an array of z values for fixed
  ! complex l(l+1) and integer m.
  !
  !   llp1_re/im : real and imaginary parts of l(l+1)
  !   m          : integer order (>= 0)
  !   nz         : number of z values
  !   zvals      : (nz) real input  array
  !   out        : (nz) complex output array
  !   rc         : 0 = success; 1 = failure (too many iterations)
  subroutine hyp_llp1_f(llp1_re, llp1_im, m, nz, zvals, out, rc)
    real(c_double),            intent(in)  :: llp1_re, llp1_im
    integer,                   intent(in)  :: m, nz
    real(c_double),            intent(in)  :: zvals(nz)
    complex(c_double_complex), intent(out) :: out(nz)
    integer,                   intent(out) :: rc

    rc = int( c_hyp_llp1(llp1_re, llp1_im, &
                         int(m, c_int), int(nz, c_int), &
                         zvals, out) )
  end subroutine hyp_llp1_f

  ! ==================================================================
  ! hyp_lmz_f
  ! ==================================================================
  ! Evaluate 2F1(-l, l+1; 1+m; z) for large z (0.5 < z < 1) using the
  ! psi/logarithmic series expansion (AMS55 15.3.11).
  !
  !   l            : complex degree l
  !   m            : integer order (>= 0)
  !   z            : real argument (0.5 < z < 1)
  !   gamma_ratio0 : Gamma(m+1)*Gamma(m) / Gamma(m-l) / Gamma(m+l+1)
  !   psi0         : psi(1)+psi(1+m)-psi(m-l)-psi(m+l+1)
  !   returns      : complex value of 2F1
  function hyp_lmz_f(l, m, z, gamma_ratio0, psi0) result(res)
    complex(c_double_complex), intent(in) :: l
    integer,                   intent(in) :: m
    real(c_double),            intent(in) :: z
    complex(c_double_complex), intent(in) :: gamma_ratio0
    complex(c_double_complex), intent(in) :: psi0
    complex(c_double_complex)             :: res

    res = c_hyp_lmz(l, int(m, c_int), z, gamma_ratio0, psi0)
  end function hyp_lmz_f

  ! ==================================================================
  ! zigg_f
  ! ==================================================================
  ! Generate float-precision standard normal random variates using the
  ! Ziggurat method (Marsaglia & Tsang 2000).
  !
  !   num_needed : number of normals requested
  !   num_ints   : number of 32-bit random integers supplied
  !   rand_ints  : (num_ints) integer(c_int32_t) array of random values
  !   out        : (num_needed) real(c_float) output array
  !   rc         : 0 = all produced; >0 = number still outstanding
  !                (random integers exhausted before num_needed reached)
  subroutine zigg_f(num_needed, num_ints, rand_ints, out, rc)
    integer,            intent(in)  :: num_needed, num_ints
    integer(c_int32_t), intent(in)  :: rand_ints(num_ints)
    real(c_float),      intent(out) :: out(num_needed)
    integer,            intent(out) :: rc

    rc = int( c_zigg(int(num_needed, c_int), int(num_ints, c_int), &
                     rand_ints, out) )
  end subroutine zigg_f

  ! ==================================================================
  ! cov_legendre_f
  ! ==================================================================
  ! Compute Mord=1 scalar covariance and cross-covariance for all
  ! m-modes (0..m_max) at all z-points via direct associated Legendre
  ! series.  Bypasses hyp_llp1 entirely; always numerically stable.
  ! Used when |Im(l(l+1))| = sqrt(c0/c2) is large (short xcorr).
  !
  !   m_max     : max m-mode index (= n_m - 1 = nphi/2)
  !   N         : number of z-points (= uhalf = nz/2+1)
  !   lmax      : Legendre truncation degree (typically 3*nz)
  !   c0, c2    : power spectrum C_l = 1/(c0 + c2*(l(l+1))^2)
  !   z_pts     : (N) cos(theta) values, equator-first
  !   cov       : (m_max+1, N)   output scalar covariances
  !   cross_cov : (m_max+1, N-1) output scalar cross-covariances
  !   rc        : 0 = success, -1 = allocation failure
  subroutine cov_legendre_f(m_max, N, lmax, c0, c2, z_pts, cov, cross_cov, rc)
    integer,        intent(in)  :: m_max, N, lmax
    real(c_double), intent(in)  :: c0, c2
    real(c_double), intent(in)  :: z_pts(N)
    real(c_double), intent(out) :: cov((m_max+1)*N)
    real(c_double), intent(out) :: cross_cov((m_max+1)*(N-1))
    integer,        intent(out) :: rc

    rc = int( c_cov_legendre(int(m_max, c_int), int(N, c_int), int(lmax, c_int), &
                             c0, c2, z_pts, cov, cross_cov) )
  end subroutine cov_legendre_f

  ! ==================================================================
  ! state_space_1_f
  ! ==================================================================
  ! Construct Kalman-filter innovation (scalar sqrt) and transition
  ! (scalar) for the Mord=1 case.  Called after cov_legendre_f.
  !
  !   N         : number of z steps (= uhalf)
  !   cross_cov : (N-1) input scalar cross-covariances
  !   cov       : (N)   input scalar covariances
  !   innov     : (N)   output innovation = sqrt(cov - trans^2 * cov_prev)
  !   trans     : (N-1) output transition = cross_cov / cov
  !   rc        : 0 = success, i+1 = step i failed (non-PD)
  subroutine state_space_1_f(N, cross_cov, cov, innov, trans, rc)
    integer,        intent(in)  :: N
    real(c_double), intent(in)  :: cross_cov(N-1)
    real(c_double), intent(in)  :: cov(N)
    real(c_double), intent(out) :: innov(N)
    real(c_double), intent(out) :: trans(N-1)
    integer,        intent(out) :: rc

    rc = int( c_state_space_1(int(N, c_int), cross_cov, cov, innov, trans) )
  end subroutine state_space_1_f

end module smerfs_interface
