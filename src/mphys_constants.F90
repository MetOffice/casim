module mphys_constants
  use variable_precision, only: wp

  implicit none

  real, parameter :: rhow = 997.0

  real :: rho0 = 1.22 ! reference density

  real :: fixed_cloud_number = 150.0*1.0e6
  real :: fixed_rain_number = 5.0e3
  real :: fixed_rain_mu = 0.0

  real :: fixed_ice_number = 1000.0

  real :: fixed_snow_number = 500.0
  real :: fixed_snow_mu = 0.0

  real :: fixed_graupel_number = 500.0
  real :: fixed_graupel_mu = 0.0

  real :: fixed_aerosol_number  = 50.0*1.0e6
  real :: fixed_aerosol_rm      = 0.05*1.0e-6
  real :: fixed_aerosol_sigma   = 1.5
  real :: fixed_aerosol_density = 1777.0

  real :: cp = 1005.0

  !< Some of these should really be functions of temperature...
  !  AH - constants set to match UM.
  real(wp) :: visair = 1.44E-5  ! kinematic viscosity of air
  real(wp) :: Cwater = 4180.0  ! specific heat capacity of liquid water
  real(wp) :: Cice = 2100.0    ! specific heat capacity of ice

  real(wp) :: Mw = 0.18015e-1 ! Molecular weight of water  [kg mol-1].
  real(wp) :: zetasa = 0.8e-1 ! Surface tension at solution-air
  ! interface
  real(wp) :: Ru = 8.314472   ! Universal gas constant
  real(wp) :: Rd = 287.05     ! gas constant for dry air
  real(wp) :: Rv = 461.5100164      ! gas constant for water vapour (matched the UM)
  ! Do not use 'eps' below - this is a Fortran intrinsic!
  real(wp) :: mp_eps = 1.6077    ! (Rv/Rd)
  real(wp) :: repsilon = 0.62198 ! 1.0 / mp_eps (matched the UM)
  real(wp) :: Dv = 0.226e-4   ! diffusivity of water vapour in air
  real(wp), parameter :: Lv = 0.2501e7   ! Latent heat of vapourization
  real(wp), parameter :: Ls = 0.2835e7   ! Latent heat of sublimation
  real(wp), parameter :: Lf = Ls - Lv    ! Latent heat of fusion
  real(wp) :: ka = 0.243e-1   ! thermal conductivity of air
  real(wp) :: g = 9.8         ! gravitational acceleration ms-2
  real(wp) :: SIGLV=8.0e-02   ! Liquid water-air surface tension [N m-1].

  real(wp), parameter :: pi = 3.14159265358979323846 ! pi (as used in UM)

  real(wp) :: m3_to_cm3 = 1.0e-6 ! Conversion factor from [m-3] to [cm-3]

end module mphys_constants
