MODULE lookup

USE mphys_die, ONLY: throw_mphys_error
USE variable_precision, ONLY: wp
USE mphys_switches, ONLY: diag_mu_option, l_passive3m, max_mu, l_sm_fix_n0
USE mphys_parameters, ONLY: hydro_params, c_r, d_r, p1, p2, p3
USE mphys_constants, ONLY: fixed_rain_number,   &
     fixed_rain_mu
USE special, ONLY: GammaFunc

#if DEF_MODEL==MODEL_KiD
  ! Kid modules
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, j_here, k_here
USE runtime, ONLY: time
#elif DEF_MODEL==MODEL_LEM
USE diaghelp_lem, ONLY: i_here, j_here, k_here
USE com_params, ONLY: time
#elif DEF_MODEL==MODEL_UM
USE diaghelp_um, ONLY: i_here, j_here, k_here
USE timestep_mod, ONLY: time => timestep_number
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here, time
#endif

IMPLICIT NONE

INTEGER, PARAMETER :: nmu=501
REAL(wp) :: min_mu=0.0!, max_mu=35.
REAL(wp), ALLOCATABLE :: mu_g(:), mu_i(:)
REAL(wp), ALLOCATABLE :: mu_g_sed(:), mu_i_sed(:)
!  real(wp) :: j1, j2
!  real(wp) :: k1, k2
!  real(wp) :: l2, l3

INTERFACE get_slope
  MODULE PROCEDURE get_slope_3M, get_slope_2M, get_slope_1M
END INTERFACE

INTERFACE get_lam_n0
  MODULE PROCEDURE get_lam_n0_3M, get_lam_n0_2M, get_lam_n0_1M
END INTERFACE
CONTAINS

FUNCTION Gfunc(mu, p1, p2, p3)

REAL(wp), INTENT(IN) :: mu, p1, p2, p3
REAL(wp) :: Gfunc, GfuncL
REAL(wp) :: k1, k2, k3

k3 = p2-p1
k1 = (p3-p2)/k3
k2 = (p1-p3)/k3

!    GfuncL=exp(k1*log(GammaFunc(1.+mu+p1)) &
!         +k2*log(GammaFunc(1.+mu+p2)) &
!         +log(GammaFunc(1.+mu+p3)) &
!         )

Gfunc = GammaFunc(1.0+mu+p1)**k1    &
       *GammaFunc(1.0+mu+p2)**k2    &
       *GammaFunc(1.0+mu+p3)


END FUNCTION Gfunc

FUNCTION Hfunc(m1,m2,m3, p1, p2, p3)

REAL(wp), INTENT(IN) :: m1,m2,m3
REAL(wp), INTENT(IN) :: p1,p2,p3
REAL(wp) :: Hfunc
REAL(wp) :: k1, k2, k3

k3 = p2-p1
k1 = (p3-p2)/k3
k2 = (p1-p3)/k3

Hfunc = EXP(k1*LOG(m1) + k2*LOG(m2) + LOG(m3))

END FUNCTION Hfunc

SUBROUTINE set_mu_lookup(p1, p2, p3, index, VALUE)

REAL(wp), INTENT(IN) :: p1, p2, p3
REAL(wp), INTENT(INOUT) :: INDEX(:), VALUE(:)
INTEGER :: i, lb, ub

lb=LBOUND(index,1)
ub=UBOUND(index,1)

DO i=lb,ub
  INDEX(i) = min_mu + (max_mu - min_mu)/(nmu-1.0)*(i-1)
  VALUE(i) = Gfunc(INDEX(i), p1, p2, p3)
END DO

END SUBROUTINE set_mu_lookup


SUBROUTINE get_slope_generic(k, params, n0, lam, mu, mass, number, m3)

INTEGER, INTENT(IN) :: k
TYPE(hydro_params), INTENT(IN) :: params
REAL(wp), INTENT(OUT) :: n0, lam, mu
REAL(wp), INTENT(IN) :: mass
REAL(wp), INTENT(IN), OPTIONAL :: number, m3

REAL(wp) :: m1, m2, p1, p2, p3

m1=mass/params%c_x
p1=params%p1
p2=params%p2
p3=params%p3
IF (number <= 0) THEN
!  PRINT*, 'ERROR in lookup', time, k, i_here, j_here,     &
!         params%id, params%i_2m, m1, number, m3
!  CALL throw_mphys_error(2, 'get_slope_generic', 'ERROR in lookup')
  return
END IF
IF (params%l_3m) THEN
  IF (l_passive3m) THEN
    m2 = number
    mu = params%fix_mu
    CALL get_lam_n0(m1, m2, p1, p2, mu, lam, n0)
  ELSE
    m2 = number
    CALL get_mu(m1, m2, m3, p1, p2, p3, mu)
    CALL get_lam_n0(m1, m2, m3, p1, p2, p3, mu, lam, n0)
  END IF
ELSE IF (params%l_2m) THEN
  m2 = number
  mu = params%fix_mu
  CALL get_lam_n0(m1, m2, p1, p2, mu, lam, n0)
ELSE
  mu = params%fix_mu
  n0 = params%fix_n0
  CALL get_lam_n0(m1, p1, mu, lam, n0)
END IF
IF (lam <= 0) THEN
  PRINT*, 'ERROR2 in lookup', params%id, params%i_2m, m1, number, lam
  CALL throw_mphys_error(2, 'get_slope_generic', 'ERROR2 in lookup')
END IF
END SUBROUTINE get_slope_generic

SUBROUTINE get_slope_3M(k, mass, number, moment3, p1, p2, p3,   &
     n0, lam, mu)
    ! 3 moment version

INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN) :: mass, number, moment3
REAL(wp), INTENT(IN) :: p1, p2, p3
REAL(wp), INTENT(OUT) :: n0, lam, mu
REAL(wp) :: m1,m2,m3

m1=mass/c_r
m2=number
m3=moment3

CALL get_mu(m1, m2, m3, p1, p2, p3, mu)
CALL get_lam_n0(m1, m2, m3, p1, p2, p3, mu, lam, n0)

END SUBROUTINE get_slope_3M


SUBROUTINE get_slope_2M(k, mass, number, p1, p2, n0, lam, mu)
    ! 2 moment version

INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN) :: mass, number
REAL(wp), INTENT(IN) :: p1, p2
REAL(wp), INTENT(OUT) :: n0, lam, mu
REAL(wp) :: m1,m2

REAL(wp) :: Dm, Deq

m1=mass/c_r
m2=number

SELECT CASE (diag_mu_option)
  CASE default
    mu = 0.0!params%fix_mu
  CASE (1)
    mu = 11.8*(1000.0*(m2/m1)**(1.0/(p2-p1))-.7)**2 + 2.0
  CASE (2)
    Dm = 1.0e3*(m1/m2/c_r)**(1.0/p1) ! in mm
    Deq = 1.1
    IF (Dm <= Deq) THEN
      mu = 6.0*TANH((4.0*(Dm - Deq))**2) + 1.0
    ELSE
      mu = 30.0*TANH((1.0*(Dm - Deq))**2) + 1.0
    END IF
END SELECT

CALL get_lam_n0(m1, m2, p1, p2, mu, lam, n0)

END SUBROUTINE get_slope_2M

SUBROUTINE get_slope_1M(k, mass, p1, n0, lam, mu)
    ! 1 moment version

INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN)  :: p1
REAL(wp), INTENT(IN) :: mass
REAL(wp), INTENT(OUT) :: n0, lam, mu
REAL(wp) :: m1

m1=mass/c_r

mu = fixed_rain_mu

CALL get_lam_n0(m1, p1, mu, lam, n0)

END SUBROUTINE get_slope_1M

SUBROUTINE get_mu(m1, m2, m3, p1, p2, p3, mu, mu_g_o, mu_i_o)

REAL(wp), INTENT(IN) :: m1, m2, m3
REAL(wp), INTENT(IN) :: p1, p2, p3
REAL(wp), INTENT(IN), OPTIONAL :: mu_g_o(:), mu_i_o(:)
REAL(wp), INTENT(OUT) :: mu

REAL(wp) :: G
REAL(wp) :: muG(nmu), muI(nmu)
INTEGER  :: i
REAL(wp) :: k1, k2, k3, pos

IF (PRESENT(mu_g_o) .AND. PRESENT(mu_i_o)) THEN
  muG=mu_g_o
  muI=mu_i_o
ELSE
  muG=mu_g
  muI=mu_i
END IF

k1 = p3-p2
k2 = p1-p3
k3 = p2-p1

pos=SIGN(1.0_wp,k1*k2) ! scaling to +-1 preserves precision for < condition

G = Hfunc(m1,m2,m3,p1,p2,p3)

IF (pos*G < pos*muG(1)) THEN
  mu = -1.0e-10  ! set to be small and negative so that it is picked up
                     ! in checks later on
ELSE
  DO i=1,nmu-1
    IF ((muG(i)-G)*(muG(i+1)-G) <= 0.0)EXIT
  END DO

  IF (i==nmu) THEN
    mu=max_mu
  ELSE
    mu = muI(i) +     &
         (G-muG(i))/(muG(i+1)-muG(i))*(muI(i+1)-muI(i))
  END IF
END IF

END SUBROUTINE get_mu

SUBROUTINE get_lam_n0_3M(m1, m2, m3, p1, p2, p3, mu, lam, n0)
    ! 3M version

REAL(wp), INTENT(IN) :: m1, m2, m3, mu
REAL(wp), INTENT(IN) :: p1, p2, p3
REAL(wp), INTENT(OUT) :: lam, n0
REAL(wp) :: p, m
REAL(wp) :: l2

!    l2 = 1./(p2-p3)
!
!    lam = ((GammaFunc(1.+mu+p2)/GammaFunc(1.+mu+p3)) &
!         *(m3/m2))**l2
    ! Make sure we use m1 and m2 to calculate lambda, in case we've gone
    ! out of bounds with mu and need to modify m3.
l2 = 1.0/(p2-p1)

lam = ((GammaFunc(1.0+mu+p2)/GammaFunc(1.0+mu+p1))     &
         *(m1/m2))**l2

m=m2
p=p2

n0 = lam**(p)*m*GammaFunc(1.0+mu)/GammaFunc(1.0+mu+p)

END SUBROUTINE get_lam_n0_3M

SUBROUTINE get_lam_n0_2M(m1, m2, p1, p2, mu, lam, n0)
    ! 2M version

REAL(wp), INTENT(IN) :: m1, m2, mu
REAL(wp), INTENT(IN) :: p1, p2
REAL(wp), INTENT(OUT) :: lam, n0
REAL(wp) :: p, m
REAL(wp) :: j1

j1=1.0/(p1-p2)

lam = ((GammaFunc(1.0+mu+p1)/GammaFunc(1.0+mu+p2))     &
       *(m2/m1))**(j1)

m=m2
p=p2

n0 = lam**(p)*m*GammaFunc(1.0+mu)/GammaFunc(1.0+mu+p)

END SUBROUTINE get_lam_n0_2M

SUBROUTINE get_lam_n0_1M(m1, p1, mu, lam, n0)
    ! 1M version

REAL(wp), INTENT(IN) :: m1, mu, n0
REAL(wp), INTENT(IN) :: p1
REAL(wp), INTENT(OUT) :: lam

! 
! Fixing Nx is equivalent to having n0=na*lam**(1+mu)
! (c.f. LEM formulation, na, nb)
! if we want to fix na and nb, we can do this and lambda is
! given by:
! lam=(na*gamma(1+mu+p1)/m1)**(1/(1+mu+p1-nb))
!

lam = (n0*GammaFunc(1.0+mu+p1)/GammaFunc(1.0+mu)     &
       *m1**(-1.0))**(1.0/p1)

END SUBROUTINE get_lam_n0_1M


SUBROUTINE get_n0(m, p, mu, lam, n0)
    ! Get n0 given a moment and lamda and mu

REAL(wp), INTENT(IN) :: m, mu, lam
REAL(wp), INTENT(IN) :: p
REAL(wp), INTENT(OUT) :: n0

n0 = m/(GammaFunc(1+mu+p)*lam**(-p)/GammaFunc(1+mu))

END SUBROUTINE get_n0

FUNCTION moment(n0,lam,mu,p)

REAL(wp) :: n0, lam, mu, p
REAL(wp) :: moment

moment=n0*GammaFunc(1+mu+p)*lam**(-p)/GammaFunc(1+mu)

END FUNCTION moment

END MODULE lookup
