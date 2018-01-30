module shipway_parameters
  use variable_precision, only: wp

  Implicit None

  !-------------------------------------------------------------
  ! The following variables describe the different 
  ! modes of dry aerosol
  !------------------------------------------------------------
  integer, parameter :: max_nmodes=3
  integer :: nmodes ! number of modes in aerosol distribution.
                    ! This is the size for the following arrays
  integer :: imode  ! current mode being considered
  real(wp) :: Ndi(max_nmodes)      ! Number concentration (m-3)
  real(wp) :: rdi(max_nmodes)      ! Geometric mean radius (m)
  real(wp) :: sigmad(max_nmodes)   ! standard deviation
  real(wp) :: bi(max_nmodes)       ! Solubility parameters (see K&C)
  real(wp) :: betai(max_nmodes)    ! Solubility parameters (see K&C)
  logical :: use_mode(max_nmodes)  ! A flag that is set to true if there is significant aerosol number in a given mode
  logical :: l_aerosol_set     ! Has the aerosol been set?

  !-------------------------------------------------------------
  ! The following variables are the terms used in the K&C
  ! formulation of differential activation
  !-------------------------------------------------------------

  real(wp) :: sigmas(max_nmodes)   ! supersaturation dispersion
  real(wp) :: s0i(max_nmodes)      ! geometric mean supersaturation

  real(wp) :: ai(max_nmodes)        ! Coefficient for lookup method
  real(wp) :: logsigmas(max_nmodes) ! Log of sigmas

  !-------------------------------------------------------------
  ! Some additions which we are convenient to include here for
  ! Coupling to Dan's driver
  !-------------------------------------------------------------
  real(wp) :: vantHoff(max_nmodes)  ! van't Hoff factor
  real(wp) :: epsv(max_nmodes)      ! fraction of soluble mass
  real(wp) :: massMole(max_nmodes)  ! Molar mass
  real(wp) :: density(max_nmodes)  ! Molar mass

  ! Numerical issues
  real(wp), parameter :: Nd_min=1.e0 ! minimum number concentration use to set use_mode

end module shipway_parameters
