! Dummy module for cloud fraction scheme to ensure that MONC and KiD build

MODULE cloud_frac_scheme

IMPLICIT NONE

CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='CLOUD_FRAC_SCHEME'

CONTAINS

SUBROUTINE cloud_frac_casim_mphys(k, pressure, T_in, T_l, rhcrit_lev, &
                                  qs, qv, cloud_mass, qfields, cloud_mass_new)

USE variable_precision, ONLY: wp

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER,  INTENT(IN)  :: k
REAL(wp), INTENT(IN)  :: pressure, T_in, T_l, rhcrit_lev, qs, qv, cloud_mass, qfields
REAL(wp), INTENT(OUT) :: cloud_mass_new

CHARACTER(LEN=*),   PARAMETER :: RoutineName='CLOUD_FRAC_CASIM_MPHYS'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

cloud_mass_new = cloud_mass

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE cloud_frac_casim_mphys

END MODULE cloud_frac_scheme
