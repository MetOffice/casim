module derived_constants
  use variable_precision, only: wp
  use special, only: pi
  use mphys_parameters, only: p1, p2, p3, sp1, sp2, sp3, rain_params, a_r, b_r, f_r, &
       a2_r, b2_r, f2_r, ice_params, nucleated_ice_radius, nucleated_ice_mass, cloud_params, snow_params, graupel_params
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

    implicit none

    character(len=*), parameter :: RoutineName='SET_CONSTANTS'

    cloud_params%p1=p1
    cloud_params%sp1=sp1
    cloud_params%p2=p2
    cloud_params%sp2=sp2
    cloud_params%p3=p3
    cloud_params%sp3=sp3

    rain_params%p1=p1
    rain_params%sp1=sp1
    rain_params%p2=p2
    rain_params%sp2=sp2
    rain_params%p3=p3
    rain_params%sp3=sp3

    ice_params%p1=p1
    ice_params%sp1=sp1
    ice_params%p2=p2
    ice_params%sp2=sp2
    ice_params%p3=p3
    ice_params%sp3=sp3

    snow_params%p1=p1
    snow_params%sp1=sp1
    snow_params%p2=p2
    snow_params%sp2=sp2
    snow_params%p3=p3
    snow_params%sp3=sp3

    graupel_params%p1=p1
    graupel_params%sp1=sp1
    graupel_params%p2=p2
    graupel_params%sp2=sp2
    graupel_params%p3=p3
    graupel_params%sp3=sp3

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
    nucleated_ice_radius=10.0e-6
    nucleated_ice_mass=4.0/3.0*pi*ice_params%density*(nucleated_ice_radius)**3
  end subroutine set_constants
end module derived_constants
