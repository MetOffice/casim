MODULE ventfac

USE variable_precision, ONLY: wp
USE passive_fields, ONLY: rho
USE special, ONLY: pi, GammaFunc
USE mphys_parameters, ONLY: hydro_params, vent_1, vent_2
USE mphys_constants, ONLY: visair, rho0, Dv

IMPLICIT NONE

CONTAINS

SUBROUTINE ventilation(k, V, n0, lam, mu, params)

INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN) :: n0, lam, mu
REAL(wp), INTENT(OUT) :: V  ! bulk ventilation factor

TYPE(hydro_params), INTENT(IN) :: params

REAL(wp) :: T1, T2
REAL(wp) :: Sc ! Schmidt number

REAL(wp) :: a_x, b_x, f_x

a_x = params%a_x
b_x = params%b_x
f_x = params%f_x

Sc = visair/Dv

T1 = vent_1*(mu+1.0)/lam
T2 = vent_2*Sc**(1.0/3.0)     &
       * (a_x*(rho(k)/rho0)**(0.5)*rho(k)/visair)**(0.5)

V = 2.0*pi*n0/rho(k)*(T1     &
       + T2*GammaFunc(.5*b_x + mu + 2.5)/GammaFunc(1.0+mu) &
       *(1.0 + .5*f_x/lam)**(-(.5*b_x + mu + 2.5))*lam**(-.5*b_x - 1.5))
    ! for computational efficiency/accuracy changed from...
![' ']    !*(lam + .5*f_x)**(-(.5*b_x + mu + 2.5))*lam**(1.+mu)) &

END SUBROUTINE ventilation

END MODULE ventfac
