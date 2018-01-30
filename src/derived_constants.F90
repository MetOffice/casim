module derived_constants
  use special, only: pi
  use mphys_parameters, only: rain_params, a_r, b_r, f_r,               &
                              a2_r, b2_r, f2_r, ice_params,             &
                              nucleated_ice_radius, nucleated_ice_mass

  use mphys_switches, only: l_abelshipway

  implicit none

  character(len=*), parameter, private :: ModuleName='DERIVED_CONSTANTS'

  private

  public set_constants
contains

  ! set up some constants here
  ! NB generally we will need p1,p2,p3 to be consistent
  ! between species (for phase conversions) (IS THIS REALLY THE CASE?),
  ! but sp1,sp2,sp3 etc need not be the same
  subroutine set_constants()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Local variables
    character(len=*), parameter :: RoutineName='SET_CONSTANTS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (l_abelshipway) then ! override any other fallspeed settings

      rain_params%a_x=4854.1
      rain_params%b_x=1.0
      rain_params%f_x=195.0
      rain_params%a2_x=-446.009
      rain_params%b2_x=0.782127
      rain_params%f2_x=4085.35

      a_r=4854.1
      b_r=1.0
      f_r=195.0
      a2_r=-446.009
      b2_r=0.782127
      f2_r=4085.35

    end if

    nucleated_ice_radius = 10.0e-6
    nucleated_ice_mass   =  4.0/3.0*pi*ice_params%density*(nucleated_ice_radius)**3

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine set_constants
end module derived_constants
