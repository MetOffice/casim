MODULE passive_fields

!Fields which aren't updated during the microphysics calculation

USE mphys_die, ONLY: throw_mphys_error
USE variable_precision, ONLY: wp
USE qsat_funs, ONLY: qsaturation
USE mphys_switches, ONLY: i_th

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: i_here, j_here
#elif  DEF_MODEL==MODEL_UM
USE diaghelp_um, ONLY: i_here, j_here
#elif  DEF_MODEL==MODEL_LEM
USE diaghelp_lem, ONLY: i_here, j_here
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here
#endif

IMPLICIT NONE

REAL(wp), ALLOCATABLE :: rho(:), pressure(:), z(:), exner(:), rexner(:)
REAL(wp), ALLOCATABLE :: z_half(:), z_centre(:), dz(:)
REAL(wp), ALLOCATABLE :: qws(:), qws0(:), TdegC(:), TdegK(:), w(:), tke(:)

REAL(wp), ALLOCATABLE :: rdz_on_rho(:)
REAL(wp) :: dt

LOGICAL, SAVE :: first_passive = .TRUE.


CONTAINS

SUBROUTINE set_passive_fields(kl, ku, dt_in, rho_in, p_in, exner_in,   &
     z_half_in, z_centre_in, dz_in, w_in, tke_in, qfields)

INTEGER, INTENT(IN) :: kl, ku
REAL(wp), INTENT(IN) :: dt_in
REAL(wp), INTENT(IN) :: rho_in(kl:ku), p_in(kl:ku), exner_in(kl:ku)
REAL(wp), INTENT(IN) :: z_half_in(kl-1:ku),z_centre_in(kl:ku),dz_in(kl:ku)
REAL(wp), INTENT(IN) :: w_in(kl:ku), tke_in(kl:ku)
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
INTEGER :: nz, k

REAL(wp), POINTER :: theta(:)

theta => qfields(:,i_th)

nz = ku-kl+1
dt = dt_in
IF (first_passive) THEN
  ALLOCATE(rho(nz))
  ALLOCATE(pressure(nz))
  ALLOCATE(exner(nz))
  ALLOCATE(rexner(nz))
  ALLOCATE(dz(nz))
  ALLOCATE(z_half(0:nz))
  ALLOCATE(z_centre(nz))
  ALLOCATE(w(nz))
  ALLOCATE(tke(nz))
  ALLOCATE(rdz_on_rho(nz))
  ALLOCATE(qws(nz))
  ALLOCATE(qws0(nz))
  ALLOCATE(TdegK(nz))
  ALLOCATE(TdegC(nz))
  first_passive = .FALSE.
END IF


rho(:)        = rho_in(kl:ku)
pressure(:)   = p_in(kl:ku)
exner(:)      = exner_in(kl:ku)
rexner(:)     = 1.0/exner(:)
w(:)          = w_in(kl:ku)
tke(:)        = tke_in(kl:ku)
dz(:)         = dz_in(kl:ku)
z_half(0:)     = z_half_in(kl-1:ku)
z_centre(:)   = z_centre_in(kl:ku)
rdz_on_rho(:) = 1.0/(dz_in(kl:ku)*rho_in(kl:ku))
DO k=1,nz
  TdegK(k)    = theta(k)*exner(k)
  TdegC(k)    = TdegK(k) - 273.15
  qws(k)      = qsaturation(TdegK(k), pressure(k)/100.0)
  qws0(k)     = qsaturation(273.15_wp, pressure(k)/100.0)
END DO

qws(nz)=1.0e-8

IF (ANY(qws==0.0)) THEN
  DO k=1,nz
    PRINT *, 'DEBUG qws', i_here, j_here, k, theta(k), TdegK(k), qws(k),rdz_on_rho(k)
  END DO
  CALL throw_mphys_error(2, 'passive_fields', 'Error in saturation calculation')
END IF

NULLIFY(theta)

END SUBROUTINE set_passive_fields

END MODULE passive_fields
