MODULE m3_incs

USE variable_precision, ONLY: wp
USE lookup, ONLY: Gfunc
IMPLICIT NONE

CONTAINS

SUBROUTINE m3_inc_type4(dm3_y, c_x, c_y, p3, dm3_x)
    ! for changes of phase from category y to category x
    ! dM3_x = -(c_y/c_x)**(p3/3) * dM3_y
    ! c_x, c_y are densities of x and y respectively
    ! p3 is the third moment (must be the same for both x and y)
    ! This assumes both categories are spherical.

REAL(wp), INTENT(IN) :: dm3_y
REAL(wp), INTENT(IN) :: c_x,c_y
REAL(wp), INTENT(IN) :: p3
REAL(wp), INTENT(OUT) :: dm3_x

dm3_x = -(c_y/c_x)**(p3/3.0) * dm3_y

END SUBROUTINE m3_inc_type4

SUBROUTINE m3_inc_type3(p1, p2, p3,   &
     dm1, dm2, dm3, mu)
    ! tendency of m3 as a function of dm1, dm2
    ! assuming initial value for mu

REAL(wp), INTENT(IN) :: dm1, dm2
REAL(wp), INTENT(IN) :: p1, p2, p3
REAL(wp), INTENT(IN) :: mu
REAL(wp), INTENT(OUT) :: dm3
REAL(wp) :: k1, k2, k3

k3 = p2-p1
k1 = (p3-p2)/k3
k2 = (p1-p3)/k3

dm3 = (Gfunc(mu, p1, p2, p3)*dm1**(-k1)*dm2**(-k2))

END SUBROUTINE m3_inc_type3

SUBROUTINE m3_inc_type2(m1, m2, m3, p1, p2, p3,   &
     dm1, dm2, dm3, mu_init)

    ! tendency of m3 as a function of m1,m2,dm1,dm2
    ! assuming shape parameter does not vary
REAL(wp), INTENT(IN) :: m1, m2, m3, dm1, dm2
REAL(wp), INTENT(IN) :: p1, p2, p3
REAL(wp), INTENT(OUT) :: dm3
REAL(wp), INTENT(IN), optional :: mu_init ! initial mu for type3
REAL(wp) :: k1, k2, k3

REAL(wp) :: fac


IF (m1 > 0.0)then

  k1 = p3-p2
  k2 = p1-p3
  k3 = p2-p1

  IF (ABS(p3-6.0)< SPACING(p3)) THEN ! p3=6
    ! since we've hardwired p1=3 and p2=0, we can reduce the
    ! calculation to...
    fac = ((m1+dm1)/m1)*((m1+dm1)/m1)*(m2/(m2+dm2)) - 1.0
  ELSE IF (ABS(p3-4.5)< SPACING(p3)) THEN ! p3=6
    ! since we've hardwired p1=3 and p2=0, we can reduce the
    ! calculation to.
    fac = SQRT(((m1+dm1)/m1)*((m1+dm1)/m1)*((m1+dm1)/m1)*(m2/(m2+dm2))) - 1.0
  ELSE IF (ABS(p3-1.5)< SPACING(p3)) THEN ! p3=6
    ! since we've hardwired p1=3 and p2=0, we can reduce the
    ! calculation to...
    fac = SQRT(((m1+dm1)/m1)*((m2+dm2)/m2)) - 1.0
  ELSE
    !      fac=((m1/(m1+dm1))**(k1/k3)*(m2/(m2+dm2))**(k2/k3)-1.)
    fac=EXP((k1*LOG(m1/(m1+dm1)) + k2*LOG(m2/(m2+dm2)))/k3)-1.0
  END IF
  IF (fac < -1.0) THEN
    PRINT*, 'm3inc_2 ERROR', fac, m1,m2,m3,dm1,dm2,k1,k2,k3
    PRINT*, SQRT(fac)
  END IF

  dm3 = m3*fac

ELSE ! If there is no pre-existing mass, then use type 3

  call m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu_init)

END IF

END SUBROUTINE m3_inc_type2

SUBROUTINE m3_inc_type1(m1, m2, dm1, dm2)
    ! tendency of m2 as a function of m1,dm1
    ! assuming shape parameter and slope do not vary

REAL(wp), INTENT(IN) :: m1, m2, dm1
REAL(wp), INTENT(OUT) :: dm2

dm2 = dm1*m1/m2

END SUBROUTINE m3_inc_type1

END MODULE m3_incs
