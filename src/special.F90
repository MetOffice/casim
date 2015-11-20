MODULE special

USE variable_precision, ONLY: wp
!  Use solvers, only: brent

IMPLICIT NONE

REAL(wp), PARAMETER :: euler=0.57721566
REAL(wp), PARAMETER :: pi=3.141592654

INTERFACE erfinv
  MODULE PROCEDURE erfinv1
END INTERFACE

INTERFACE gammafunc
!     module procedure gammafunc1
  MODULE PROCEDURE gammalookup
!     module procedure intrinsic_gamma
END INTERFACE

REAL(wp) :: set_a ! value to use for inversion of f(x)=a

REAL(wp), ALLOCATABLE :: gammalookup_arg(:)
REAL(wp), ALLOCATABLE :: gammalookup_val(:)
REAL(wp) :: gammalookup_xmin, gammalookup_xmax, gammalookup_dx
LOGICAL :: l_gammalookup_set = .FALSE.

CONTAINS

SUBROUTINE set_gammalookup(xmin, xmax, dx)

    ! NB The following should provide sufficient range
    ! and density of points for linear interpolation
    ! to provide appropriate accuracy for any values required.
REAL(wp), INTENT(IN) :: xmin !< Minimum value of argument
REAL(wp), INTENT(IN) :: xmax !< Maximum value of argument
REAL(wp), INTENT(IN) :: dx   !< spacing of argument calculations

    ! Local variables
REAL(wp) :: arg
INTEGER :: nargs
INTEGER :: i
gammalookup_xmin = xmin
gammalookup_xmax = xmax
gammalookup_dx = dx

nargs = CEILING((xmax - xmin)/dx + 1)
ALLOCATE(gammalookup_arg(nargs))
ALLOCATE(gammalookup_val(nargs))

arg = xmin
DO i = 1,nargs
  gammalookup_arg(i) = arg
  gammalookup_val(i) = gammaFunc1(arg)
  arg = arg + dx
END DO

END SUBROUTINE set_gammalookup

FUNCTION gammalookup(x)

REAL(wp), INTENT(IN) :: x
REAL(wp) :: gammalookup
REAL(wp) :: xmin=1e-12, xmax=100.0, dx=.0001

INTEGER :: i_minus
IF (.NOT. l_gammalookup_set) THEN
  CALL set_gammalookup(xmin, xmax, dx)
  l_gammalookup_set=.TRUE.
END IF
    ! Locate x in table
i_minus = INT((x - gammalookup_xmin)/gammalookup_dx) + 1
gammalookup = 0.5*(gammalookup_val(i_minus) + gammalookup_val(i_minus+1))

END FUNCTION gammalookup

  !================!
  ! Gamma function !
  !================!
FUNCTION gammafunc1(x)
REAL(wp), INTENT(IN) :: x
REAL(wp) :: f,g,z
REAL(wp) :: gammafunc1

f=HUGE(x)
g=1
z=x
IF ((z+INT(ABS(z))) /= 0) THEN
  DO WHILE(z<3)
          ! Lets use a recursion relation for Gamma functions
          ! to get a large argument and use Stirlings formula
    g=g*z
    z=z+1
  END DO

       ! This is just stirlings formula...
  f=(1.0-2.0*(1-2.0/(3.0*z*z))/(7.0*z*z))/(30.0*z*z)
  f=(1.0-f)/(12.0*z)+z*(LOG(z)-1)
  f=(EXP(f)/g)*SQRT(2.0*pi/z)

END IF

gammafunc1=f

END FUNCTION gammafunc1

!   function intrinsic_gamma(x)
!     real(wp), intent(IN) :: x
!     real(wp) :: intrinsic_gamma

!     intrinsic_gamma=gamma(x)

!   end function intrinsic_gamma

FUNCTION erfg(x,c)

REAL(wp), INTENT(IN) :: x
INTEGER, INTENT(IN) :: c ! 0 gives erf(x)
                             ! 1 gives erfc(x)
REAL(wp) :: erfg, f, z
INTEGER :: j, cc

REAL(wp) :: t, a1, a2, a3, a4, a5, p
INTEGER :: sign

z=x
cc=c

IF (ABS(z) < 1e-10) THEN
  f=0.0
ELSE IF (ABS(z) < 1.5) THEN
  j=3 + INT(9*ABS(z))
  f=1
  DO WHILE(j /= 0)
    f=1.0+f*z**2*(.5-j)/(j*(.5+j))
    j=j-1
  END DO
  f=cc+f*z*(2.0-4.0*cc)/SQRT(pi)
ELSE
  cc=cc*ABS(z)/z
  j=3 + INT(32/ABS(z))
  f=0.0
  DO WHILE(j /= 0)
    f=1.0/(f*j + SQRT(2.0*z*z))
    j=j-1
  END DO
  f=f*(cc*cc+cc-1.0)*SQRT(2.0/pi)*EXP(-z*z)+(1.0-cc)
END IF

    ! quick fix, but should do this properly...
f=f*((1-c)*ABS(z)/z +c)
erfg=f

END FUNCTION erfg

! Now  use the fortran intrinsic erf

FUNCTION erf(x)

REAL(wp), INTENT(IN) :: x
INTEGER, PARAMETER :: c=0
REAL(wp) :: erf

erf=erfg(x,c)

END FUNCTION erf

FUNCTION erfc(x)

REAL(wp), INTENT(IN) :: x
INTEGER, PARAMETER :: c=1
REAL(wp) :: erfc

erfc=erfg(x,c)

END FUNCTION erfc

!   function erfc(x)

!     real(wp), intent(IN) :: x
!     real(wp) :: erfc

!     erfc=1-erf(x)

!   end function erfc

FUNCTION erfinv1(x)
    ! Inverse of error function
    !
    ! This needs more work to get good accuracy
    !

REAL(wp), INTENT(IN) :: x
REAL(wp) :: erfinv1

erfinv1 = .5*SQRT(pi)*(x + pi/12.0*x*x*x + 7.0/480.0*pi*pi*x**5     &
       + 127.0/40320*pi**3*x**7 + 4369.0/5806080*pi**4*x**9 &
       + 34807.0/182476800.0*pi**5*x**11)

END FUNCTION erfinv1

FUNCTION erfinv2(x)
    ! Inverse of error function
    !
    ! Alternative version solves equation
    !

REAL(wp), INTENT(IN) :: x
REAL(wp) :: erfinv2
REAL(wp) :: work, work_old, diff, erfx, erfx_old

diff=9999.0
work_old=.2
work=work+1.0
DO WHILE(ABS(diff) > 1e-3)
  erfx=erf(work)-x
  erfx_old=erf(work_old)-x
  diff= -erfx*(work_old-work)/(erfx_old-erfx)
  work_old=work
  work = work + diff
END DO

erfinv2 = work

END FUNCTION erfinv2

FUNCTION erfinv3(x, tol)
    ! Inverse of error function
    !
    ! Alternative version solves equation
    !

REAL(wp), INTENT(IN) :: x
REAL(wp), OPTIONAL :: tol
REAL(wp) :: erfinv3
REAL(wp) :: work, diff, erfx, derfx, tolval

tolval=1e-3
IF (PRESENT(tol)) tolval=tol

IF (ABS(x) > .95) THEN
       ! don't converge well for abs(x)->1
       ! should treat this properly
       ! (i.e. more sophisticated solver)
       ! but don't really care too much
       ! about these values for now
  work = erfinv1(x)
ELSE

  diff=9999.0
  work=erfinv1(x)
  DO WHILE(ABS(diff) > tolval)
    erfx=erf(work)-x
    derfx=2.0*EXP(-work*work)/SQRT(pi)
    diff= -erfx/derfx
    work = work + diff
  END DO
END IF

erfinv3 = work

END FUNCTION erfinv3


  ! function erfinv4(x, tol)
  !   ! Inverse of error function
  !   !
  !   ! Alternative version solves equation
  !   !

  !   real(wp), intent(IN) :: x
  !   real(wp), optional :: tol
  !   real(wp) :: erfinv4
  !   real(wp) :: work, z, diff, erfx, derfx

  !   integer :: iflag, ibrent
  !   real(wp) :: tolx, tolf, a, b
  !   logical :: verbose=.false.

  !   ! seems to be better at solving for positive x
  !   ! so exploit anti-symmetry
  !   z=x
  !   if (x<0)z=-x

  !   set_a=z ! Set the value of a for the external function

  !   tolf=1e-5
  !   if (present(tol)) tolf=tol
  !   tolx=tolf

  !   work = erfinv1(z) ! initial guess
  !   a=work-.5*abs(work)
  !   b=work+1.5+abs(work)
  !   call brent(erf_a, a, b, work, iflag, ibrent &
  !      ,tolx, tolf, verbose)

  !   if (x<0)work=-work

  !   erfinv4 = work

  !   write(123,*) x, work

  ! end function erfinv4

FUNCTION erf_a(x)

REAL(wp), INTENT(IN) :: x
REAL(wp) :: erf_a

erf_a = erf(x) - set_a

END FUNCTION erf_a


END MODULE special
