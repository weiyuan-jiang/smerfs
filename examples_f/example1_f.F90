! example1_f.F90
!
! Fortran 90 translation of examples/example1.py
!
! Demonstrates:
!   1. build_filter  - compute state-space filter coefficients (trans, innov)
!                      for an isotropic GRF on the sphere with power spectrum
!                        C_l = 1 / (c0 + c1*l(l+1) + c2*(l(l+1))^2)
!                      using nz=128 iso-latitude rings and nphi=256 phi points.
!
!   2. create_realisation - draw one realisation of the GRF via the Kalman
!                           state-space walk + inverse FFT (written as CSV).
!
!   3. analytic_cov - evaluate C(cos theta) = sum_l (2l+1) C_l P_l(cos theta)
!                     / 4 pi  on 1000 equally-spaced z-values and write to CSV.
!
! Outputs
!   realisation.csv   : nz x nphi matrix, one row per latitude ring
!   correl.csv        : z, C(z) pairs for 1000 points in [-1, 1]
!
! Build: see CMakeLists.txt in examples_f/
! The program calls the C kernels through smerfs_interface.
!
! Memory layout note
! ------------------
! The C library uses row-major (C-order) arrays.  Fortran uses column-major.
! To pass multidimensional arrays directly (without packing overhead), all
! arrays shared with the C kernels are declared with REVERSED dimension order:
!
!   Python/C  cov[m, n, p, q]  (shape m_max+1, N, M, M)
!   Fortran   cov(q, p, n, m)  (declared as cov(M, M, N, 0:m_max))
!
! This makes the in-memory byte layout identical, so the C pointer arithmetic
! works correctly without any transposing or packing loops.

program example1_f
  use, intrinsic :: iso_c_binding, only : c_double, c_float, &
                                          c_int32_t, c_double_complex
  use smerfs_interface
  implicit none

  ! ---------------------------------------------------------------
  ! Parameters matching example1.py
  ! ---------------------------------------------------------------
  integer, parameter :: nz   = 128   ! iso-latitude rings
  integer, parameter :: nphi = 256   ! phi points
  ! GRF power spectrum  C_l = 1/(c0 + c1*l(l+1) + c2*(l(l+1))^2)
  real(c_double), parameter :: c0 = 1.0d0
  real(c_double), parameter :: c1 = 0.0d0
  real(c_double), parameter :: c2 = 1.0d-4
  integer, parameter :: Mord = 2        ! Markov order = len(coeffs)-1
  integer, parameter :: lmax = 1000     ! truncation for analytic_cov
  integer, parameter :: nz_corr = 1000  ! points for covariance plot

  ! ---------------------------------------------------------------
  ! Constants
  ! ---------------------------------------------------------------
  real(c_double), parameter :: PI = 3.141592653589793d0

  ! Upper hemisphere: uhalf = nz/2 + 1 rings (equator first, pole last)
  integer, parameter :: uhalf = nz/2 + 1   ! z-points on upper half
  integer, parameter :: n_m   = nphi/2 + 1 ! m-modes (rfft half-spectrum)

  ! m_max = n_m - 1  (modes 0 .. n_m-1 needed for cov/cross_cov tables)
  integer, parameter :: m_max = n_m - 1

  ! ---------------------------------------------------------------
  ! Working arrays — C-compatible memory layout (reversed Fortran dims)
  ! ---------------------------------------------------------------
  real(c_double) :: z_pts(uhalf)           ! cos(theta)
  ! tau_p: C layout is [p_index][n_index] row-major  (2*Mord-1 rows, uhalf cols)
  ! Fortran reversed dims: (uhalf, 2*Mord-1) column-major = same flat layout
  real(c_double) :: tau_p(uhalf, 2*Mord-1) ! tau_p(n, p) in Fortran <=> C [p][n]
  real(c_double) :: eta_ratio(uhalf-1)     ! tau_{i+1}/tau_i
  real(c_double) :: xvals(uhalf), yvals(uhalf) ! x=(1-z)/2, y=(1+z)/2

  ! F[m, n]  and  H[m, n]  : (m_max+Mord+1) x uhalf
  !   Fortran: (uhalf, m_max+Mord+1) -> column-major = C row-major [m][n]
  complex(c_double_complex) :: Fmat(uhalf, m_max+Mord+1)
  complex(c_double_complex) :: Hmat(uhalf, m_max+Mord+1)

  ! cov[m, n, p, q]    shape (m_max+1, uhalf, Mord, Mord)
  ! Fortran: (Mord, Mord, uhalf, 0:m_max) — column-major matches C row-major
  real(c_double), allocatable :: cov(:,:,:,:)      ! (Mord, Mord, uhalf, 0:m_max)
  real(c_double), allocatable :: cross_cov(:,:,:,:) ! (Mord, Mord, uhalf-1, 0:m_max)

  ! State-space arrays — same reversed-dim convention
  ! innov[n, p, q]   shape (uhalf, Mord, Mord)  ->  Fortran (Mord, Mord, uhalf)
  ! trans[n, p, q]   shape (uhalf-1, Mord, Mord) -> Fortran (Mord, Mord, uhalf-1)
  ! stored for all m-modes: (Mord, Mord, uhalf, n_m) etc.
  real(c_double), allocatable :: all_innov(:,:,:,:)  ! (Mord, Mord, uhalf,   n_m)
  real(c_double), allocatable :: all_trans(:,:,:,:)  ! (Mord, Mord, uhalf-1, n_m)

  ! Partial fraction decomposition
  complex(c_double_complex) :: roots_llp1(Mord)
  complex(c_double_complex) :: norms_pf(Mord)
  complex(c_double_complex) :: norm, llp1_val, lam_val, sin_pi_lam
  real(c_double)             :: llp1_re, llp1_im

  ! Random numbers
  integer(c_int32_t), allocatable :: rand_ints(:)
  real(c_float),      allocatable :: noise_real(:)

  ! Realisation
  complex(c_float), allocatable :: fp_f(:,:)      ! (Mord, n_m)
  complex(c_float), allocatable :: fp_fstart(:,:) ! (Mord, n_m)
  complex(c_float), allocatable :: noise(:,:,:)   ! (Mord, n_m, uhalf)
  complex(c_float), allocatable :: res_cplx(:,:)  ! (n_m, nz)
  real(c_double),   allocatable :: realisation(:,:) ! (nphi, nz) -> irfft

  ! Analytic covariance
  real(c_double) :: C_l(0:lmax)
  real(c_double) :: z_corr(nz_corr), correl(nz_corr)
  real(c_double) :: llp1_arr

  ! Miscellaneous
  integer :: i, j, im_hyp, imode, ip, iq, nn, rc, iroot
  integer :: skip, n_rand_needed, n_rand_got
  real(c_double) :: disc_re, tmp, theta_val
  integer(c_int32_t) :: rand_state
  integer :: funit

  ! ----------------------------------------------------------------
  ! 1.  z-points on the upper hemisphere (equator first, pole last)
  !     Python: theta = (arange(nz//2+1)[::-1]+0.5)*(pi/nz)
  ! ----------------------------------------------------------------
  do i = 1, uhalf
    theta_val = (real(uhalf - i, c_double) + 0.5d0) * (PI / real(nz, c_double))
    z_pts(i)  = cos(theta_val)
  end do

  ! ----------------------------------------------------------------
  ! 2.  Twiddle factors  tau = sqrt((1-z)/(1+z))
  !     tau_p(n, p) in Fortran = C tau_power[p][n]
  !     p index (1-based): Mord = tau^0, Mord+1 = tau^1, Mord-1 = tau^-1, ...
  ! ----------------------------------------------------------------
  tau_p(:, Mord)   = 1.0d0
  tau_p(:, Mord+1) = sqrt((1.0d0 - z_pts) / (1.0d0 + z_pts))  ! tau^1
  tau_p(:, Mord-1) = 1.0d0 / tau_p(:, Mord+1)                  ! tau^-1
  ! For Mord=2 the loop below has zero iterations (Mord-1=1 < 2) -- intentional
  do i = 2, Mord-1
    tau_p(:, Mord+i) = tau_p(:, Mord+i-1) * tau_p(:, Mord+1)
    tau_p(:, Mord-i) = tau_p(:, Mord-i+1) * tau_p(:, Mord-1)
  end do

  do i = 1, uhalf-1
    eta_ratio(i) = tau_p(i+1, Mord+1) * tau_p(i, Mord-1)
  end do

  xvals = 0.5d0 * (1.0d0 - z_pts)
  yvals = 0.5d0 * (1.0d0 + z_pts)

  ! ----------------------------------------------------------------
  ! 3.  Partial fraction decomposition of 1/(c0 + c1*k + c2*k^2)
  ! ----------------------------------------------------------------
  disc_re = c1*c1 - 4.0d0*c0*c2
  if (disc_re >= 0.0d0) then
    roots_llp1(1) = cmplx((-c1 + sqrt(disc_re)) / (2.0d0*c2), 0.0d0, c_double_complex)
    roots_llp1(2) = cmplx((-c1 - sqrt(disc_re)) / (2.0d0*c2), 0.0d0, c_double_complex)
  else
    tmp = sqrt(-disc_re) / (2.0d0*c2)
    roots_llp1(1) = cmplx(-c1/(2.0d0*c2),  tmp, c_double_complex)
    roots_llp1(2) = cmplx(-c1/(2.0d0*c2), -tmp, c_double_complex)
  end if
  do iroot = 1, Mord
    norms_pf(iroot) = 1.0d0 / (cmplx(c1, 0.0d0, c_double_complex) &
                      + 2.0d0*cmplx(c2, 0.0d0, c_double_complex)*roots_llp1(iroot))
  end do

  ! ----------------------------------------------------------------
  ! 4.  Allocate covariance arrays (reversed-dim, C-compatible layout)
  !     cov      (Mord, Mord, uhalf,   0:m_max)
  !     cross_cov(Mord, Mord, uhalf-1, 0:m_max)
  ! ----------------------------------------------------------------
  allocate(cov      (Mord, Mord, uhalf,   0:m_max))
  allocate(cross_cov(Mord, Mord, uhalf-1, 0:m_max))
  allocate(all_innov(Mord, Mord, uhalf,   n_m))
  allocate(all_trans(Mord, Mord, uhalf-1, n_m))
  cov       = 0.0d0
  cross_cov = 0.0d0

  ! ----------------------------------------------------------------
  ! 5.  Evaluate hypergeometric functions and accumulate covariances
  ! ----------------------------------------------------------------
  write(*,'(a)') 'Computing covariances...'

  do iroot = 1, Mord
    llp1_val = roots_llp1(iroot)
    llp1_re  = real(llp1_val, c_double)
    llp1_im  = aimag(llp1_val)
    norm     = norms_pf(iroot)

    ! lam such that lam*(lam+1) = llp1
    if (llp1_re < -0.25d0) then
      lam_val = cmplx(-0.5d0, 0.0d0, c_double_complex) + &
                (0.0d0, 1.0d0) * sqrt(-0.25d0 - llp1_val)
    else
      lam_val = cmplx(-0.5d0, 0.0d0, c_double_complex) - &
                sqrt(cmplx(0.25d0, 0.0d0, c_double_complex) + llp1_val)
    end if

    sin_pi_lam = sin(PI * lam_val)
    ! norm from Python: -(0.25/pi) * a_i * pi / sin(lam*pi)
    norm = -(0.25d0 / PI) * norm * PI / sin_pi_lam

    ! Compute Fmat(n, m+1) = 2F1(llp1, m, x_n)  for m = 0 .. m_max+Mord
    !         Hmat(n, m+1) = 2F1(llp1, m, y_n)
    ! Fmat is declared (uhalf, m_max+Mord+1) -- pass column m+1 to the C call
    do im_hyp = 0, m_max + Mord
      call hyp_llp1_f(llp1_re, llp1_im, im_hyp, uhalf, xvals, Fmat(:, im_hyp+1), rc)
      if (rc /= 0) then
        write(*,'(a,i0)') 'ERROR: hyp_llp1 F failed at im_hyp=', im_hyp ; stop 1
      end if
      call hyp_llp1_f(llp1_re, llp1_im, im_hyp, uhalf, yvals, Hmat(:, im_hyp+1), rc)
      if (rc /= 0) then
        write(*,'(a,i0)') 'ERROR: hyp_llp1 H failed at im_hyp=', im_hyp ; stop 1
      end if
    end do

    ! Pass arrays to C update_cov.
    !
    ! C signature:
    !   update_cov(m_max, N, M, norm_re, norm_im, llp1_re, llp1_im,
    !              F[m_max+M+1][N], H[m_max+M+1][N],
    !              tau_power[2M-1][N], eta_ratio2[N-1],
    !              cov[(m_max+1)*N*M*M], cross_cov[(m_max+1)*(N-1)*M*M])
    !
    ! F in C is indexed F[(m+p)*N + n], i.e. row-major (m_row, n_col).
    ! Fmat in Fortran is (uhalf, m_max+Mord+1) = column-major (n_row, m_col),
    ! so memory layout: Fmat(1,1), Fmat(2,1),...,Fmat(uhalf,1), Fmat(1,2),...
    ! This matches C row-major F[0][0],F[0][1],...,F[0][N-1],F[1][0],...
    ! ONLY IF we treat Fortran rows as C columns — i.e. (n, m) in Fortran =
    ! (m, n) in C. But C accesses F[(m+p)*N + n], meaning F[m_index][n_index].
    ! With Fortran (uhalf, m_max+Mord+1), memory is F(n, m) column-major =
    ! F[n + uhalf*m] in flat indexing, whereas C wants F[m*N + n] = same!
    ! So the layouts ARE compatible: Fortran (uhalf, m_max+Mord+1) column-major
    ! gives the same flat memory as C row-major F[m_max+Mord+1][uhalf].
    !
    ! Similarly cov(Mord, Mord, uhalf, 0:m_max) column-major =
    ! C row-major cov[m_max+1][uhalf][Mord][Mord]. Correct.

    call update_cov_f(m_max, uhalf, Mord, &
                      real(norm, c_double), aimag(norm), &
                      llp1_re, llp1_im, &
                      Fmat, Hmat, tau_p, eta_ratio, &
                      cov, cross_cov, rc)
    if (rc /= 0) then
      write(*,'(a,i0)') 'ERROR: update_cov failed, rc=', rc ; stop 1
    end if

  end do ! iroot

  ! ----------------------------------------------------------------
  ! 6.  State-space decomposition for each m-mode
  !
  ! state_space expects cov[N][M][M] in C row-major.
  ! Our cov(Mord,Mord,uhalf,0:m_max) Fortran column-major stores element
  ! (q,p,n,m) at flat offset q-1 + Mord*(p-1) + Mord^2*(n-1) + Mord^2*N*m.
  ! A slice cov(:,:,:,imode) has layout [q][p][n] = C [N][M][M] reversed dims,
  ! BUT since state_space C code accesses cov[n][p][q] and our slice
  ! cov(:,:,:,imode) in flat memory is indexed q+M*p+M^2*n -- this matches
  ! C row-major [n][p][q] when read as such. No transpose needed.
  ! ----------------------------------------------------------------
  write(*,'(a)') 'Building state-space filters...'

  block
    real(c_double) :: innov_m(Mord, Mord, uhalf)
    real(c_double) :: trans_m(Mord, Mord, uhalf-1)

    do imode = 0, n_m-1
      call state_space_f(uhalf, Mord, &
                         cross_cov(:, :, :, imode), &
                         cov(:, :, :, imode), &
                         innov_m, trans_m, rc)
      if (rc /= 0) then
        write(*,'(a,i0,a,i0)') &
          'ERROR: state_space failed at imode=', imode, ', rc=', rc
        stop 1
      end if
      all_innov(:, :, :, imode+1) = innov_m
      all_trans(:, :, :, imode+1) = trans_m
    end do
  end block

  ! ----------------------------------------------------------------
  ! 7.  create_realisation
  !
  !     all_innov(Mord, Mord, uhalf, n_m)  ->  innov(p, q, ring, imode)
  !     all_trans(Mord, Mord, uhalf-1, n_m) ->  trans(p, q, step, imode)
  !
  !     noise(Mord, n_m, uhalf):  noise(q, imode, ring) ~ CN(0,1/2)
  !
  !     State walk:
  !       fp_f(p, imode) = sum_q innov(p, q, ring, imode) * noise(q, imode, ring)
  !                       + sum_q trans(p, q, ring-1, imode) * fp_f_prev(q, imode)
  !
  !     Result stored in res_cplx(imode, ring) = fp_f(1, imode) (zeroth component).
  ! ----------------------------------------------------------------
  write(*,'(a)') 'Creating realisation...'

  ! noise needs nz rings (not just uhalf): south-hemisphere walk accesses
  ! noise at indices uhalf+1..nz, matching the Python reference (filters.py).
  allocate(noise    (Mord, n_m, nz))
  allocate(fp_f     (Mord, n_m))
  allocate(fp_fstart(Mord, n_m))
  allocate(res_cplx (n_m, nz))
  allocate(realisation(nphi, nz))

  n_rand_needed = nz * Mord * n_m * 2
  n_rand_got    = n_rand_needed + n_rand_needed/32 + 64
  allocate(rand_ints (n_rand_got))
  allocate(noise_real(n_rand_needed))

  ! Simple LCG (seed=123, matching Python RandomState seed)
  rand_state = 123_c_int32_t
  do i = 1, n_rand_got
    rand_state = 1664525_c_int32_t * rand_state + 1013904223_c_int32_t
    rand_ints(i) = rand_state
  end do

  call zigg_f(n_rand_needed, n_rand_got, rand_ints, noise_real, rc)
  if (rc /= 0) &
    write(*,'(a,i0,a)') 'WARNING: zigg exhausted, ', rc, ' samples missing'

  ! noise(q, imode, ring) = (noise_real(2k-1) + i*noise_real(2k)) / sqrt(2)
  ! Loop over all nz rings so south-hemisphere indices (uhalf+1..nz) are valid.
  nn = 0
  do i = 1, nz
    do ip = 1, Mord
      do j = 1, n_m
        nn = nn + 1
        noise(ip, j, i) = cmplx(noise_real(2*nn-1), noise_real(2*nn), c_float) &
                          * real(0.5d0**0.5d0, c_float)
      end do
    end do
  end do

  ! -- Walk upward from equator (ring 1) to pole (ring uhalf) --
  ! Initialise: fp_f(:, imode) = innov(:, :, 1, imode+1) * noise(:, imode+1, 1)
  do j = 1, n_m
    do ip = 1, Mord
      fp_f(ip, j) = cmplx(0.0, 0.0, c_float)
      do iq = 1, Mord
        fp_f(ip, j) = fp_f(ip, j) + real(all_innov(ip, iq, 1, j), c_float) * noise(iq, j, 1)
      end do
    end do
  end do
  fp_fstart = fp_f
  res_cplx(:, uhalf) = fp_f(1, :)

  do i = 1, uhalf-1
    do j = 1, n_m
      do ip = 1, Mord
        fp_f(ip, j) = cmplx(0.0, 0.0, c_float)
        do iq = 1, Mord
          fp_f(ip, j) = fp_f(ip, j) &
            + real(all_trans(ip, iq, i, j), c_float) * fp_fstart(iq, j) &
            + real(all_innov(ip, iq, i+1, j), c_float) * noise(iq, j, i+1)
        end do
      end do
    end do
    fp_fstart = fp_f
    res_cplx(:, uhalf - i) = fp_f(1, :)
  end do

  ! -- Reverse derivatives then walk downward --
  ! Odd-indexed components flip sign across the equator
  do ip = 1, Mord
    if (mod(ip-1, 2) == 0) then
      fp_f(ip, :) = fp_fstart(ip, :)
    else
      fp_f(ip, :) = -fp_fstart(ip, :)
    end if
  end do
  fp_fstart = fp_f

  ! skip=1 if nz even (equatorial step already counted), 0 if odd
  skip = 1 - mod(nz, 2)

  do i = 0, nz - uhalf - 1
    do j = 1, n_m
      do ip = 1, Mord
        fp_f(ip, j) = cmplx(0.0, 0.0, c_float)
        do iq = 1, Mord
          fp_f(ip, j) = fp_f(ip, j) &
            + real(all_trans(ip, iq, i+skip+1, j), c_float) * fp_fstart(iq, j) &
            + real(all_innov(ip, iq, i+skip+2, j), c_float) * noise(iq, j, uhalf+i)
        end do
      end do
    end do
    fp_fstart = fp_f
    res_cplx(:, uhalf + i + 1) = fp_f(1, :)
  end do

  ! DC mode (imode=1): real part only, scaled by sqrt(2)
  do i = 1, nz
    res_cplx(1, i) = cmplx(real(res_cplx(1, i)) * sqrt(2.0), 0.0, c_float)
  end do

  ! irfft: res_cplx(n_m, nz) -> realisation(nphi, nz)
  call irfft_rows(res_cplx, nz, n_m, nphi, realisation)

  ! ----------------------------------------------------------------
  ! 8.  Write realisation to CSV
  ! ----------------------------------------------------------------
  write(*,'(a)') 'Writing realisation.csv...'
  open(newunit=funit, file='realisation.csv', status='replace', action='write')
  do i = 1, nz
    do j = 1, nphi
      if (j < nphi) then
        write(funit, '(es20.12,a)', advance='no') realisation(j, i), ','
      else
        write(funit, '(es20.12)') realisation(j, i)
      end if
    end do
  end do
  close(funit)

  ! ----------------------------------------------------------------
  ! 9.  analytic_cov
  ! ----------------------------------------------------------------
  write(*,'(a)') 'Computing analytic covariance...'
  do i = 0, lmax
    llp1_arr = real(i, c_double) * real(i+1, c_double)
    C_l(i) = 1.0d0 / (c0 + c1*llp1_arr + c2*llp1_arr**2)
  end do
  do i = 1, nz_corr
    z_corr(i) = -1.0d0 + 2.0d0*real(i-1, c_double)/real(nz_corr-1, c_double)
  end do
  do i = 1, nz_corr
    correl(i) = sum_legendre(z_corr(i), C_l, lmax)
  end do

  ! ----------------------------------------------------------------
  ! 10.  Write covariance to CSV
  ! ----------------------------------------------------------------
  write(*,'(a)') 'Writing correl.csv...'
  open(newunit=funit, file='correl.csv', status='replace', action='write')
  write(funit, '(a)') 'z,C(z)'
  do i = 1, nz_corr
    write(funit, '(es20.12,a,es20.12)') z_corr(i), ',', correl(i)
  end do
  close(funit)

  write(*,'(a)') 'Done.'
  write(*,'(a,i0,a,i0,a)') '  realisation.csv : ', nz, ' x ', nphi, ' grid'
  write(*,'(a,i0,a)')       '  correl.csv      : ', nz_corr, ' z values'

  deallocate(cov, cross_cov, all_innov, all_trans)
  deallocate(noise, fp_f, fp_fstart, res_cplx, realisation)
  deallocate(rand_ints, noise_real)

contains

  ! ----------------------------------------------------------------
  ! irfft_rows
  !
  ! Inverse real FFT of each column (ring) of the half-spectrum.
  !   inp(n_m, nrows) complex -> out(nphi, nrows) real
  !
  ! The numpy convention is irfft(X) = Re( IDFT(X) ), which gives:
  !   out(j) = sum_{k=0}^{n_m-1} [X(k)*exp(2*pi*i*j*k/nphi) + c.c.] / nphi
  ! where the k=0 and k=Nyquist terms are real-only.
  ! Multiplied by nphi afterwards (Python: res *= p).
  ! ----------------------------------------------------------------
  subroutine irfft_rows(inp, nrows, nm, nph, out)
    integer,          intent(in)  :: nrows, nm, nph
    complex(c_float), intent(in)  :: inp(nm, nrows)
    real(c_double),   intent(out) :: out(nph, nrows)
    integer :: row, j, k
    real(c_double) :: angle, re_sum, two_pi_over_n
    complex(c_double_complex) :: Xk

    two_pi_over_n = 2.0d0 * PI / real(nph, c_double)

    do row = 1, nrows
      do j = 0, nph-1
        ! DC term k=0
        re_sum = real(inp(1, row), c_double)
        ! k=1 .. nm-2: both k and nph-k contribute (factor 2)
        do k = 1, nm-2
          angle = two_pi_over_n * real(j*k, c_double)
          Xk = cmplx(real(inp(k+1, row), c_double), &
                     real(aimag(inp(k+1, row)), c_double), c_double_complex)
          re_sum = re_sum + 2.0d0*(real(Xk)*cos(angle) - aimag(Xk)*sin(angle))
        end do
        ! Nyquist k=nm-1 (real only for even nph)
        k = nm - 1
        angle = two_pi_over_n * real(j*k, c_double)
        re_sum = re_sum + real(inp(nm, row), c_double) * cos(angle)
        out(j+1, row) = re_sum
      end do
    end do
  end subroutine irfft_rows

  ! ----------------------------------------------------------------
  ! sum_legendre
  !
  ! Evaluate  sum_{l=0}^{lmax_in} (2l+1) * Cl(l) * P_l(z) / (4*pi)
  ! using the three-term Legendre recurrence.
  ! ----------------------------------------------------------------
  function sum_legendre(z, Cl, lmax_in) result(s)
    integer,        intent(in) :: lmax_in
    real(c_double), intent(in) :: z
    real(c_double), intent(in) :: Cl(0:lmax_in)
    real(c_double) :: s, P0, P1, P2, coeff
    integer :: l

    s  = 0.0d0
    P0 = 1.0d0
    P1 = z
    s  = s + Cl(0) * 1.0d0  * P0 * (0.25d0/PI)
    if (lmax_in >= 1) s = s + Cl(1) * 3.0d0 * P1 * (0.25d0/PI)
    do l = 2, lmax_in
      P2    = (real(2*l-1,c_double)*z*P1 - real(l-1,c_double)*P0) / real(l,c_double)
      coeff = real(2*l+1, c_double) * Cl(l) * (0.25d0/PI)
      s     = s + coeff * P2
      P0 = P1
      P1 = P2
    end do
  end function sum_legendre

end program example1_f
