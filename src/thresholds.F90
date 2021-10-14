module thresholds

  use variable_precision, only: wp

  implicit none

  ! This module contains threshold values which are used for various
  ! decisions in the code: Mainly for tidying or deciding whether
  ! processes should be active or not.

  ! Small values...
  ! (values above which may switch on processes)
  real(wp) :: th_small = 1.0e-9   ! small theta increment
  real(wp) :: qv_small = 1.0e-9   ! small vapour
  real(wp) :: ql_small = 1.0e-9   ! small cloud mass
  real(wp) :: qr_small = 1.0e-9   ! small rain mass
  real(wp) :: nl_small = 1e-6    ! small cloud number
  real(wp) :: nr_small = 10.0     ! small rain number
  real(wp) :: m3r_small = 1.0e-25 ! small rain moment 3! obsolete
  real(wp) :: qi_small = 1.0e-10   ! small ice mass
  real(wp) :: ni_small = 1.0e-6   ! small ice number
  real(wp) :: qs_small = 1.0e-10   ! small snow mass
  real(wp) :: ns_small = 1.0e-6   ! small snow number
  real(wp) :: m3s_small = 1.0e-25 ! small snow moment 3! obsolete
  real(wp) :: qg_small = 1.0e-10   ! small graupel mass
  real(wp) :: ng_small = 1.0e-6   ! small graupel number
  real(wp) :: m3g_small = 1.0e-25 ! small graupel moment 3

  real(wp) :: cfliq_small=1e-3 !small liquid cloudfraction

  ! Thresholds for tidying up small numbers...
  ! (values below which may want to remove and ignore)
  real(wp) :: th_tidy = 0.0     ! tidy theta increment
  real(wp) :: qv_tidy = 0.0     ! tidy vapour
  real(wp) :: ql_tidy = 1.0e-10  ! tidy cloud mass
  real(wp) :: qr_tidy = 1.0e-10  ! tidy rain mass
  real(wp) :: nl_tidy = 1.0e-6   ! tidy cloud number
  real(wp) :: nr_tidy = 1.0      ! tidy rain number
  real(wp) :: m3r_tidy = 1.0e-30 ! tidy rain moment 3! obsolete
  real(wp) :: qi_tidy = 1.0e-10   ! tidy ice mass
  real(wp) :: ni_tidy = 1.0e-6   ! tidy ice number
  real(wp) :: qs_tidy = 1.0e-10   ! tidy snow mass
  real(wp) :: ns_tidy = 1.0e-6   ! tidy snow number
  real(wp) :: m3s_tidy = 1.0e-30 ! tidy snow moment 3! obsolete
  real(wp) :: qg_tidy = 1.0e-10   ! tidy graupel mass
  real(wp) :: ng_tidy = 1.0e-6   ! tidy graupel number
  real(wp) :: m3g_tidy = 1.0e-30 ! tidy graupel moment 3 ! obsolete
  real(wp) :: ccn_tidy = 0.1e0  ! tidy ccn number
  !(NB if ccn_tidy is too small, this can result in
  ! very large values of rd (mean radius).)

  ! Significant values...
  ! (values above which may switch on processes)
  real(wp) :: qi_sig = 1.0e-5   ! sig ice mass
  real(wp) :: qg_sig = 1.0e-4   ! sig graupel mass
  real(wp) :: qs_sig = 1.0e-5   ! sig snow mass
  real(wp) :: ql_sig = 1.0e-4   ! sig cloud mass
  real(wp) :: qr_sig = 1.0e-4   ! sig rain mass).)

  ! Large values...
  ! (values above which may want to limit things)
  real(wp) :: qi_large = 1.0e-1   ! large ice mass
  real(wp) :: qg_large = 1.0e-1   ! large graupel mass
  real(wp) :: qs_large = 1.0e-1   ! large snow mass
  real(wp) :: ql_large = 1.0e-1   ! large cloud mass
  real(wp) :: qr_large = 1.0e-1   ! large rain mass
  real(wp) :: ni_large = 1.0e10   ! large ice mass
  real(wp) :: ng_large = 1.0e10   ! large graupel mass
  real(wp) :: ns_large = 1.0e10   ! large snow mass
  real(wp) :: nl_large = 1.0e10   ! large cloud mass
  real(wp) :: nr_large = 1.0e10   ! large rain mass

  ! other threshold quantities
  real(wp) :: ss_small = 1.0e-3   ! small sub/supersaturation (fraction)
  real(wp) :: w_small = 1.0e-3    ! small w for droplet activation
  real(wp) :: rn_min = 1.0e-6     ! minimum rain number

  real(wp) :: aeromass_small = 1.0e-25      ! small aerosol mass (kg/kg)
  real(wp) :: aeronumber_small = 1.0e-6     ! small aerosol number (/kg)

  !-----------------------------------------------------
  ! arrays containing all the threshold information
  ! for relevant microphysics options.
  ! Set in mphys_switches.
  !-----------------------------------------------------
  real(wp), allocatable :: thresh_tidy(:)     ! Tiny values which can be tidied away
  real(wp), allocatable :: thresh_small(:)    ! Small values which we require to bother with some processes
  real(wp), allocatable :: thresh_sig(:)      ! Significant values needed for some processes
  real(wp), allocatable :: thresh_large(:)    ! Large values - may be a potential problem

  real(wp), allocatable :: thresh_atidy(:)     ! Tiny values which can be tidied away for aerosol
end module thresholds
