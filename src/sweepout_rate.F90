MODULE sweepout_rate

USE variable_precision, ONLY: wp
USE special, ONLY: pi, Gammafunc
USE mphys_constants, ONLY: rho0
USE mphys_parameters, ONLY: hydro_params

IMPLICIT NONE

CONTAINS

FUNCTION sweepout(n0, lam, mu, params, rho, mass_weight)
    ! Calculate the sweepout rate given the distribution
    ! and fallspeed parameters

REAL(wp), INTENT(IN) :: mu, n0, lam
TYPE(hydro_params), INTENT(IN) :: params
REAL(wp), INTENT(IN) :: rho ! air density
    ! if present and true use mass-weighted sweepout
LOGICAL, INTENT(IN), OPTIONAL :: mass_weight

REAL(wp) :: sweepout

    ! local variables
REAL(wp) :: G3_b_mu, G1_mu !< gamma functions (these may be taken directly from
                               !< params once code has been made more efficent)
REAL(wp) :: arg3 !< argument for gamma function
REAL(wp) :: coef !< coefficient = c_x if mass-weighting

INTEGER :: k ! loop index

arg3 = 3.0 + params%b_x + mu
coef = 1.0
IF (PRESENT(mass_weight)) THEN
  IF (mass_weight) THEN
    arg3 = arg3 + params%d_x
    coef = params%c_x
  END IF
END IF

G3_b_mu = GammaFunc(arg3)
G1_mu   = GammaFunc(1.0+mu)
sweepout = coef*(pi*n0*params%a_x/4.0)*G3_b_mu/G1_mu       &
       * (1.0 + params%f_x/lam)**(-arg3)                   &
       * lam**(1 + mu - arg3)                              &
       * (rho/rho0)**(params%g_x)


END FUNCTION sweepout

FUNCTION binary_collection(n0_X, lam_X, mu_X, n0_Y, lam_Y, mu_Y,   &
     params_X, params_Y, rho, mass_weight)
    ! Calculate the number of species Y collected by species X through
    ! binary collisions as both species sediment

REAL(wp), INTENT(IN) :: mu_X, n0_X, lam_X
REAL(wp), INTENT(IN) :: mu_Y, n0_Y, lam_Y
TYPE(hydro_params), INTENT(IN) :: params_X
TYPE(hydro_params), INTENT(IN) :: params_Y
REAL(wp), INTENT(IN) :: rho ! air density
    ! if present and true use mass-weighted value (i.e. total mass accreted)
LOGICAL, INTENT(IN), OPTIONAL :: mass_weight

REAL(wp) :: binary_collection

    ! local variables
REAL(wp) :: G1_X, G2_X, G3_X  !< gamma functions (these may be taken directly from
                                           !< params once code has been made more efficent)
REAL(wp) :: G1_Y, G2_Y, G3_Y  !< gamma functions (these may be taken directly from
                                           !< params once code has been made more efficent)
REAL(wp) :: arg1_X !< argument for gamma function
REAL(wp) :: arg2_X !< argument for gamma function
REAL(wp) :: arg3_X !< argument for gamma function
REAL(wp) :: arg1_Y !< argument for gamma function
REAL(wp) :: arg2_Y !< argument for gamma function
REAL(wp) :: arg3_Y !< argument for gamma function
REAL(wp) :: l_X1 !< lam_X^-arg1_X
REAL(wp) :: l_X2 !< lam_X^-arg2_X
REAL(wp) :: l_X3 !< lam_X^-arg3_X
REAL(wp) :: l_Y1 !< lam_Y^-arg1_Y
REAL(wp) :: l_Y2 !< lam_Y^-arg2_Y
REAL(wp) :: l_Y3 !< lam_Y^-arg3_Y

REAL(wp) :: coef !< coefficient = c_x if mass-weighting
REAL(wp) :: V_X  !< mass-weighted fall velocity for X
REAL(wp) :: V_Y  !< mass-weighted fall velocity for Y
REAL(wp) :: delV !< bulk fall-speed differential

INTEGER :: k ! loop index

arg1_X = 1.0 + mu_X
arg2_X = 2.0 + mu_X
arg3_X = 3.0 + mu_X
arg1_Y = 1.0 + mu_Y
arg2_Y = 2.0 + mu_Y
arg3_Y = 3.0 + mu_Y
coef = 1.0
IF (PRESENT(mass_weight)) THEN
  IF (mass_weight) THEN
    arg1_Y = arg1_Y + params_Y%d_x
    arg2_Y = arg2_Y + params_Y%d_x
    arg3_Y = arg3_Y + params_Y%d_x
    coef = params_Y%c_x
  END IF
END IF

G1_X = GammaFunc(arg3_X)
G2_X = GammaFunc(arg3_X)
G3_X = GammaFunc(arg3_X)
G1_Y = GammaFunc(arg3_Y)
G2_Y = GammaFunc(arg3_Y)
G3_Y = GammaFunc(arg3_Y)

l_X1  = lam_X**(-arg1_X)
l_X2  = lam_X**(-arg2_X)
l_X3  = lam_X**(-arg3_X)
l_Y1  = lam_Y**(-arg1_Y)
l_Y2  = lam_Y**(-arg2_Y)
l_Y3  = lam_Y**(-arg3_Y)

V_X   = params_X%a_x * lam_X**(-params_X%b_x)*(rho/rho0)**(params_X%g_x)     &
       *GammaFunc(1.0 + mu_X + params_X%d_x + params_X%b_x)/GammaFunc(1.0 + mu_X + params_X%d_x)

V_Y   = params_Y%a_x * lam_Y**(-params_Y%b_x)*(rho/rho0)**(params_Y%g_x)     &
       *GammaFunc(1.0 + mu_Y + params_Y%d_x + params_Y%b_x)/GammaFunc(1.0 + mu_Y + params_Y%d_x)

delV = MAX(MAX(V_X,V_Y)/4.0, ABS(V_X-V_Y))

binary_collection = coef*0.25*pi*n0_X*n0_Y*delV*     &
       (l_X3*l_Y1*G1_Y*G3_X + 2.0*l_X2*l_Y2*G2_Y*G2_X + l_X1*l_Y3*G3_Y*G1_X)

END FUNCTION binary_collection

END MODULE sweepout_rate
