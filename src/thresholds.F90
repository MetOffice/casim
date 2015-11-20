MODULE thresholds

USE variable_precision, ONLY: wp

IMPLICIT NONE

  ! This module contains threshold values which are used for various
  ! decisions in the code: Mainly for tidying or deciding whether
  ! processes should be active or not.

  ! Small values...
  ! (values above which may switch on processes)
REAL(wp) :: th_small = 1.0e-9   ! small theta increment
REAL(wp) :: qv_small = 1.0e-9   ! small vapour
REAL(wp) :: ql_small = 1.0e-9   ! small cloud mass
REAL(wp) :: qr_small = 1.0e-9   ! small rain mass
REAL(wp) :: nl_small = 100.0    ! small cloud number
REAL(wp) :: nr_small = 10.0     ! small rain number
REAL(wp) :: m3r_small = 1.0e-25 ! small rain moment 3! obsolete
REAL(wp) :: qi_small = 1.0e-10   ! small ice mass
REAL(wp) :: ni_small = 1.0e-6   ! small ice number
REAL(wp) :: qs_small = 1.0e-10   ! small snow mass
REAL(wp) :: ns_small = 1.0e-6   ! small snow number
REAL(wp) :: m3s_small = 1.0e-25 ! small snow moment 3! obsolete
REAL(wp) :: qg_small = 1.0e-10   ! small graupel mass
REAL(wp) :: ng_small = 1.0e-6   ! small graupel number
REAL(wp) :: m3g_small = 1.0e-25 ! small graupel moment 3

  ! Thresholds for tidying up small numbers...
  ! (values below which may want to remove and ignore)
REAL(wp) :: th_tidy = 0.0     ! tidy theta increment
REAL(wp) :: qv_tidy = 0.0     ! tidy vapour
REAL(wp) :: ql_tidy = 1.0e-10  ! tidy cloud mass
REAL(wp) :: qr_tidy = 1.0e-10  ! tidy rain mass
REAL(wp) :: nl_tidy = 1.0e-6   ! tidy cloud number
REAL(wp) :: nr_tidy = 1.0      ! tidy rain number
REAL(wp) :: m3r_tidy = 1.0e-30 ! tidy rain moment 3! obsolete
REAL(wp) :: qi_tidy = 1.0e-20   ! tidy ice mass
REAL(wp) :: ni_tidy = 1.0e-6   ! tidy ice number
REAL(wp) :: qs_tidy = 1.0e-10   ! tidy snow mass
REAL(wp) :: ns_tidy = 1.0e-6   ! tidy snow number
REAL(wp) :: m3s_tidy = 1.0e-30 ! tidy snow moment 3! obsolete
REAL(wp) :: qg_tidy = 1.0e-10   ! tidy graupel mass
REAL(wp) :: ng_tidy = 1.0e-6   ! tidy graupel number
REAL(wp) :: m3g_tidy = 1.0e-30 ! tidy graupel moment 3 ! obsolete
REAL(wp) :: ccn_tidy = 0.1e+6  ! tidy ccn number
            !(NB if ccn_tidy is too small, this can result in
            ! very large values of rd (mean radius).)

  ! Significant values...
  ! (values above which may switch on processes)
REAL(wp) :: qi_sig = 1.0e-5   ! sig ice mass
REAL(wp) :: qg_sig = 1.0e-4   ! sig graupel mass
REAL(wp) :: qs_sig = 1.0e-5   ! sig snow mass
REAL(wp) :: ql_sig = 1.0e-4   ! sig cloud mass
REAL(wp) :: qr_sig = 1.0e-4   ! sig rain mass).)

  ! Large values...
  ! (values above which may want to limit things)
REAL(wp) :: qi_large = 1.0e-1   ! large ice mass
REAL(wp) :: qg_large = 1.0e-1   ! large graupel mass
REAL(wp) :: qs_large = 1.0e-1   ! large snow mass
REAL(wp) :: ql_large = 1.0e-1   ! large cloud mass
REAL(wp) :: qr_large = 1.0e-1   ! large rain mass
REAL(wp) :: ni_large = 1.0e10   ! large ice mass
REAL(wp) :: ng_large = 1.0e10   ! large graupel mass
REAL(wp) :: ns_large = 1.0e10   ! large snow mass
REAL(wp) :: nl_large = 1.0e10   ! large cloud mass
REAL(wp) :: nr_large = 1.0e10   ! large rain mass

  ! other threshold quantities
REAL(wp) :: ss_small = 1.0e-3   ! small sub/supersaturation (fraction)
REAL(wp) :: w_small = 1.0e-3    ! small w for droplet activation
REAL(wp) :: rn_min = 1.0e-6     ! minimum rain number

REAL(wp) :: aeromass_small = 1.0e-25      ! small aerosol mass (kg/kg)
REAL(wp) :: aeronumber_small = 1.0e-6     ! small aerosol number (/kg)

  !-----------------------------------------------------
  ! arrays containing all the threshold information
  ! for relevant microphysics options.
  ! Set in mphys_switches.
  !-----------------------------------------------------
REAL(wp), ALLOCATABLE :: thresh_tidy(:)     ! Tiny values which can be tidied away
REAL(wp), ALLOCATABLE :: thresh_small(:)    ! Small values which we require to bother with some processes
REAL(wp), ALLOCATABLE :: thresh_sig(:)      ! Significant values needed for some processes
REAL(wp), ALLOCATABLE :: thresh_large(:)    ! Large values - may be a potential problem

REAL(wp), ALLOCATABLE :: thresh_atidy(:)     ! Tiny values which can be tidied away for aerosol


END MODULE thresholds
