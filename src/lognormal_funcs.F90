module lognormal_funcs

use variable_precision, only: wp
use thresholds, only: ccn_tidy, aeromass_small, aeronumber_small
use mphys_constants, only: pi

implicit none

contains

  !
  ! Calculate mean radius of lognormal distribution
  ! given Mass and number
  !
  function MNtoRm(M, N, density, sigma)   

    implicit none

    real(wp), intent(in) :: M, N, density, sigma
    real(wp) :: MNtoRm

!    if (N==0 .or. M==0) then ! shouldn't really be here - DPG_bug_checks - M and N can be negative since gets called for increments from evaporation (removal of aerosol from activated mode).
    if (abs(N).lt.aeronumber_small .or. abs(M).lt.aeromass_small) then ! shouldn't really be here. Maybe should make thresholds smaller since will be divided by the timestep in some instances.
      MNtoRm=0.0
    else
      MNtoRm=( 3.0*M*exp(-4.5*log(sigma)**2)/(4.0*N*pi*density) )**(1.0/3.0)
    end if
  end function MNtoRm

end module lognormal_funcs
