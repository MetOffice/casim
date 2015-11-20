MODULE gauss_casim_micro

  !< Calculate the integral needed in the aggregation calculations
  !<
  !< OPTIMISATION POTENTIAL - LOOKUP

USE variable_precision, ONLY: wp
USE lookup, ONLY: max_mu
IMPLICIT NONE

INTEGER, PARAMETER :: maxq = 5
REAL(wp) :: gaussfunc_save(maxq) ! max 5 values - I.e. this is only used with 2m schemes
                                ! should be extended
INTEGER, PARAMETER, PRIVATE :: nbins_a = 50
REAL(wp) :: gaussfunc_save_2D(nbins_a,maxq)
LOGICAL :: l_save_2D(nbins_a,maxq) = .FALSE.

CONTAINS

SUBROUTINE gaussfunclookup(iq, VALUE, a, b, init)

INTEGER, INTENT(IN) :: iq !< parameter index relating to variable we're considering
REAL(wp), INTENT(OUT) :: VALUE !< returned value
REAL(wp), INTENT(IN), OPTIONAL :: a, b  !< Value of a and b to use. Only required if initializing
LOGICAL, INTENT(IN), OPTIONAL :: init !< if present and true then initialize values
                                          !< otherwise use precalculated values

LOGICAL :: vinit

vinit=.FALSE.
IF (PRESENT(init) .AND. PRESENT(a) .AND. PRESENT(b))vinit=init
IF (vinit) THEN
  VALUE = gauss_casim_func(a,b)
  gaussfunc_save(iq) = gauss_casim_func(a,b)
ELSE
  VALUE = gaussfunc_save(iq)
END IF

END SUBROUTINE gaussfunclookup

SUBROUTINE gaussfunclookup_2d(iq, VALUE, a, b)

INTEGER, INTENT(IN) :: iq !< parameter index relating to variable we're considering
REAL(wp), INTENT(OUT) :: VALUE !< returned value
REAL(wp), INTENT(IN) :: a, b  !< Value of a and b to use.  (a is mu, b is b_x)

INTEGER :: ibin ! mu(a) bin in which we sit

ibin = INT((a/max_mu)*(nbins_a-1)) + 1
IF (.NOT. l_save_2D(ibin,iq)) THEN
  VALUE = gauss_casim_func(a,b)
  gaussfunc_save_2d(ibin,iq) = gauss_casim_func(a,b)
  l_save_2D(ibin,iq)=.TRUE.
ELSE
  VALUE = gaussfunc_save_2d(ibin,iq)
END IF

END SUBROUTINE gaussfunclookup_2d

FUNCTION gauss_casim_func(a, b)

REAL(wp), INTENT(IN) :: a, b !< function arguments
REAL(wp), PARAMETER ::   tmax = 18.0    !< Limit of integration
REAL(wp), PARAMETER ::   dt   = 0.08   !< step size
REAL(wp) :: sum, t1, t2
REAL(wp) :: Gauss_casim_Func

sum = 0.0
t1  = 0.5*dt
Gaus_t1: DO WHILE (t1  <=  tmax)
  t2 = 0.5*dt
Gaus_t2: DO WHILE (t2  <=  tmax)
    sum = sum + (t1 + t2)**2*ABS( (t1**b) - (t2**b))          &
           *(t1**a)*(t2**a)*EXP(-(t1 + t2))
    t2  = t2 + dt
  END DO Gaus_t2
  t1 = t1 + dt
END DO Gaus_t1

Gauss_casim_Func = sum*dt*dt
END FUNCTION Gauss_casim_Func

END MODULE gauss_casim_micro

