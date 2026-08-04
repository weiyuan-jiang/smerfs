! Thin Fortran driver for the shared sphere-random-fields implementation.
program example1_f
  use, intrinsic :: iso_c_binding, only: c_double
  use sphere_random_fields_mod, only: sphere_random_fields
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_TileCoordType, only: grid_def_type
  implicit none

  integer, parameter :: nz = 128, nphi = 256, lmax = 1000, nz_corr = 1000
  real(c_double), parameter :: pi = 3.141592653589793d0
  type(sphere_random_fields) :: sf
  type(pert_param_type) :: pert_param
  type(grid_def_type) :: pert_grid_f
  real, allocatable :: field(:,:), field2(:,:)
  real(c_double) :: cl(0:lmax), z(nz_corr), cov(nz_corr)
  integer :: i, j, unit

  ! Equivalent Gaussian sigma for coeffs=(1,0,1e-4), computed from
  ! C(theta)/C(0)=exp(-0.5) using the same Legendre truncation as Fortran.
  pert_param%xcorr = 6.046390243223
  pert_grid_f%N_lon = nphi
  pert_grid_f%N_lat = nz
  pert_grid_f%dlon = 360.0 / real(nphi)
  pert_grid_f%dlat = 180.0 / real(nz)
  sf = sphere_random_fields(pert_param, pert_grid_f)
  allocate(field(nphi, nz), field2(nphi, nz))
  call sf%generate_2d_Random_field_seeded(123, field, field2, pert_param%xcorr, pert_param%xcorr, &
       pert_grid_f%dlon, pert_grid_f%dlat)

  open(newunit=unit, file='realisation.csv', status='replace', action='write')
  do i = 1, nz
    write(unit,'(*(es20.12,:,","))') (field(j,i), j=1,nphi)
  end do
  close(unit)

  do i = 0, lmax
    cl(i) = 1.0d0 / (sf%c0 + sf%c2 * real(i*(i+1),c_double)**2)
  end do
  do i = 1, nz_corr
    z(i) = -1.0d0 + 2.0d0*real(i-1,c_double)/real(nz_corr-1,c_double)
    cov(i) = legendre_covariance(z(i), cl, lmax)
  end do
  open(newunit=unit, file='correl.csv', status='replace', action='write')
  write(unit,'(a)') 'z,C(z)'
  do i = 1, nz_corr
    write(unit,'(es20.12,a,es20.12)') z(i), ',', cov(i)
  end do
  close(unit)
  call sf%finalize()

contains

  function legendre_covariance(x, coeffs, lm) result(value)
    real(c_double), intent(in) :: x, coeffs(0:lm)
    integer, intent(in) :: lm
    real(c_double) :: value, p0, p1, p2
    integer :: ell
    value = coeffs(0)/(4.0d0*pi)
    p0 = 1.0d0
    p1 = x
    if (lm >= 1) value = value + 3.0d0*coeffs(1)*p1/(4.0d0*pi)
    do ell = 2, lm
      p2 = (real(2*ell-1,c_double)*x*p1-real(ell-1,c_double)*p0) / real(ell,c_double)
      value = value + real(2*ell+1,c_double)*coeffs(ell)*p2/(4.0d0*pi)
      p0 = p1
      p1 = p2
    end do
  end function legendre_covariance
end program example1_f
