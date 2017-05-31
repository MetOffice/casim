! Dummy module for cloud fraction scheme to ensure that MONC and KiD build

MODULE cloud_frac_scheme

IMPLICIT NONE

CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='CLOUD_FRAC_SCHEME'

CONTAINS

SUBROUTINE cloud_frac_casim_mphys(k, pressure, abs_liquid_t, rhcrit_lev, &
                                  qs, qv, cloud_mass, qfields, cloud_mass_new)

USE variable_precision, ONLY: wp

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='CLOUD_FRAC_CASIM_MPHYS'

INTEGER,  INTENT(IN)  :: k
REAL(wp), INTENT(IN)  :: pressure, abs_liquid_t, rhcrit_lev, qs, qv, cloud_mass, qfields
REAL(wp), INTENT(OUT) :: cloud_mass_new

cloud_mass_new = cloud_mass

END SUBROUTINE cloud_frac_casim_mphys

END MODULE cloud_frac_scheme
