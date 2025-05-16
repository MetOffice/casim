! Module for random parameters scheme
! The parameters are updated in the UM and
! passed through to CASIM via this module
module casim_stph

  use variable_precision, only: wp

  implicit none

  character(len=*), parameter, private :: ModuleName='CASIM_STPH'

  ! Stochastic physics settings
  logical :: l_rp2_casim = .false. ! Switch for RP scheme
  real(wp) :: snow_a_x_rp = 12.0 ! Snow fallspeed
  real(wp) :: ice_a_x_rp = 6000000.0 ! Ice fallspeed
  real(wp) :: mpof_casim_rp = 0.5 ! mixed-phase overlap factor
  real(wp) :: fixed_cloud_number_rp = 150.0*1.0e6 ! fixed cloud number

end module casim_stph
