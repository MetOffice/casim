MODULE which_mode_to_use

  ! Routine decides which modes to put back
  ! re-evaporated (or potentially any other) aerosol

USE variable_precision, ONLY: wp
#if DEF_MODEL==MODEL_UM
USE ukca_option_mod, ONLY: l_ukca
IMPLICIT NONE
#else
IMPLICIT NONE
LOGICAL :: l_ukca = .FALSE.
#endif


INTEGER :: imethod = 1 ! method to use
INTEGER, PARAMETER :: iukca_method = 2 ! method to use if UKCA is used


REAL(wp) :: max_accumulation_mean_radius = 0.25e-6
REAL(wp) :: min_coarse_mean_radius = 1.0e-6

CONTAINS

SUBROUTINE which_mode(dm, dn, r1_in, r2_in, density,               &
     dm1, dm2, dn1, dn2)

    ! Subroutine calculates how much of the total increments
    ! to mass and number should be sent to each of two modes.

REAL(wp), INTENT(IN) :: dm  ! total mass increment
REAL(wp), INTENT(IN) :: dn  ! total number increment

REAL(wp), INTENT(IN) :: r1_in  ! mean radius of mode 1
REAL(wp), INTENT(IN) :: r2_in  ! mean radius of mode 2

REAL(wp), INTENT(IN) :: density ! density of aerosol, assumed the same across all modes

REAL(wp), INTENT(OUT) :: dm1  ! mass increment to mode 1
REAL(wp), INTENT(OUT) :: dn1  ! number increment to mode 1
REAL(wp), INTENT(OUT) :: dm2  ! mass increment to mode 2
REAL(wp), INTENT(OUT) :: dn2  ! number increment to mode 2

    ! local variables
REAL(wp) :: rm    ! mean radius of increment
REAL(wp) :: gamma ! convenience variable

REAL(wp) :: r1, r2     ! r1 and r2 (possibly modified) 
REAL(wp) :: r1_3, r2_3 ! r1**3 and r2**3

REAL(wp) :: FTPI ! 4/3*pi
REAL(wp) :: RFTPI ! 1./(4/3*pi)

IF (l_ukca)imethod = iukca_method



IF (dm/=0.0_wp .and. dm*dn > 0.0) THEN
  
  
  r1=min(max_accumulation_mean_radius, r1)
  r2=max(min_coarse_mean_radius, r2)

  dm1=0.0
  dm2=0.0
  dn1=0.0
  dn2=0.0

  FTPI=3.14159*4./3.
  RFTPI=1./FTPI

  rm = (RFTPI*dm/dn/density)**(1.0/3.0)
  

  SELECT CASE(imethod)
    CASE (iukca_method)
          ! Call ukca subroutine, or add UKCA code here?
    CASE default
      IF (rm > r2) THEN
        dm1=0.0
        dn1=0.0
        dm2=dm
        dn2=dn
      ELSE IF (rm < r1) THEN
        dm1=dm
        dn1=dn
        dm2=0.0
        dn2=0.0
      ELSE
        r1_3 = r1*r1*r1
        r2_3 = r2*r2*r2
        gamma = (rm*rm*rm - r1_3)/(r2_3 - r1_3)
        dn2=dn*(gamma)
        dn1=dn - dn2
        dm2=RFTPI*density*r2_3*dn2
        dm1=dm - dm2
      END IF
  END SELECT

END IF


END SUBROUTINE which_mode

END MODULE which_mode_to_use

