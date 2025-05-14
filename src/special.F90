module special
  use variable_precision, only: wp
  use mphys_constants, only: pi
  !  Use solvers, only: brent

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='SPECIAL'

  real(wp), parameter :: euler=0.57721566

  ! pi is set to the same value as that used in the UM 

  interface erfinv
     module procedure erfinv1
  end interface erfinv

  interface gammafunc
     !     module procedure gammafunc1
     module procedure gammalookup
     !     module procedure intrinsic_gamma
  end interface gammafunc

  real(wp), allocatable :: gammalookup_arg(:)
  real(wp), allocatable :: gammalookup_val(:)
  real(wp) :: gammalookup_xmin, gammalookup_xmax, gammalookup_dx
  logical :: l_gammalookup_set=.false.
  
  public pi, Gammafunc, casim_erfc, erfinv
contains
  ! NB The following should provide sufficient range
  ! and density of points for linear interpolation
  ! to provide appropriate accuracy for any values required.
  subroutine set_gammalookup(xmin, xmax, dx)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SET_GAMMALOOKUP'

    real(wp), intent(in) :: xmin !< Minimum value of argument
    real(wp), intent(in) :: xmax !< Maximum value of argument
    real(wp), intent(in) :: dx   !< spacing of argument calculations

    ! Local variables
    real(wp) :: arg
    integer :: nargs
    integer :: i

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    gammalookup_xmin=xmin
    gammalookup_xmax=xmax
    gammalookup_dx=dx

    nargs=ceiling((xmax - xmin)/dx + 1)
    allocate(gammalookup_arg(nargs))
    allocate(gammalookup_val(nargs))

    arg=xmin
    do i=1, nargs
      gammalookup_arg(i)=arg
      gammalookup_val(i)=gammaFunc1(arg)
      arg=arg+dx
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine set_gammalookup

  function gammalookup(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GAMMALOOKUP'

    real(wp), intent(in) :: x
    real(wp) :: gammalookup

    real(wp) :: xmin=1e-12, xmax=100.0, dx=.0001
    integer :: i_minus

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (.not. l_gammalookup_set) then
      call set_gammalookup(xmin, xmax, dx)
      l_gammalookup_set=.true.
    end if
    ! Locate x in table
    i_minus=int((x - gammalookup_xmin)/gammalookup_dx)+1
    gammalookup=0.5*(gammalookup_val(i_minus)+gammalookup_val(i_minus+1))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function gammalookup

  !================!
  ! Gamma function !
  !================!
  function gammafunc1(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GAMMAFUNC1'

    real(wp), intent(in) :: x
    real(wp) :: gammafunc1

    real(wp) :: f,g,z

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    f=huge(x)
    g=1
    z=x
    if ((z+int(abs(z))) /= 0) then
      do while(z < 3)
        ! Lets use a recursion relation for Gamma functions
        ! to get a large argument and use Stirlings formula
        g=g*z
        z=z+1
      end do

      ! This is just stirlings formula...
      f=(1.0-2.0*(1-2.0/(3.0*z*z))/(7.0*z*z))/(30.0*z*z)
      f=(1.0-f)/(12.0*z)+z*(log(z)-1)
      f=(exp(f)/g)*sqrt(2.0*pi/z)
    end if
    gammafunc1=f

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function gammafunc1

  function erfg(x,c)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ERFG'

    real(wp), intent(in) :: x
    integer, intent(in) :: c ! 0 gives erf(x)
    ! 1 gives erfc(x)
    real(wp) :: erfg, f, z
    integer :: j, cc

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    z=x
    cc=c

    if (abs(z) < 1e-10) then
      f=0.0
    else if (abs(z) < 1.5) then
      j=3+int(9*abs(z))
      f=1
      do while(j /= 0)
        f=1.0+f*z**2*(.5-j)/(j*(.5+j))
        j=j-1
      end do
      f=cc+f*z*(2.0-4.0*cc)/sqrt(pi)
    else
      cc=cc*int(abs(z)/z)
      j=3+int(32/abs(z))
      f=0.0
      do while(j /= 0)
        f=1.0/(f*j + sqrt(2.0*z*z))
        j=j-1
      end do
      f=f*(cc*cc+cc-1.0)*sqrt(2.0/pi)*exp(-z*z)+(1.0-cc)
    end if

    ! quick fix, but should do this properly...
    f=f*((1-c)*abs(z)/z +c)
    erfg=f

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function erfg

  function casim_erfc(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='CASIM_ERFC'

    real(wp), intent(in) :: x
    integer, parameter :: c=1
    real(wp) :: casim_erfc

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    casim_erfc=erfg(x,c)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function casim_erfc

  ! Inverse of error function
  !
  ! This needs more work to get good accuracy
  function erfinv1(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ERFINV1'

    real(wp), intent(in) :: x
    real(wp) :: erfinv1

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    erfinv1=.5*sqrt(pi)*(x+pi/12.0*x*x*x+7.0/480.0*pi*pi*x**5 &
         +127.0/40320*pi**3*x**7+4369.0/5806080*pi**4*x**9 &
         +34807.0/182476800.0*pi**5*x**11)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function erfinv1

  ! Inverse of error function
  !
  ! Alternative version solves equation
  function erfinv2(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ERFINV2'

    real(wp), intent(in) :: x
    real(wp) :: erfinv2

    real(wp) :: work, work_old, diff, erfx, erfx_old

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    diff=9999.0
    work_old=.2
    work=1.0
    do while(abs(diff) > 1e-3)
      erfx=erf(work)-x
      erfx_old=erf(work_old)-x
      diff=-erfx*(work_old-work)/(erfx_old-erfx)
      work_old=work
      work=work + diff
    end do

    erfinv2=work

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function erfinv2

  ! Inverse of error function
  !
  ! Alternative version solves equation
  function erfinv3(x, tol)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ERFINV3'

    real(wp), intent(in) :: x
    real(wp), optional, intent(in) :: tol
    real(wp) :: erfinv3

    real(wp) :: work, diff, erfx, derfx, tolval

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    tolval=1e-3
    if (present(tol)) tolval=tol

    if (abs(x) > .95) then
      ! don't converge well for abs(x)->1
      ! should treat this properly
      ! (i.e. more sophisticated solver)
      ! but don't really care too much
      ! about these values for now
      work=erfinv1(x)
    else

      diff=9999.0
      work=erfinv1(x)
      do while(abs(diff) > tolval)
        erfx=erf(work)-x
        derfx=2.0*exp(-work*work)/sqrt(pi)
        diff=-erfx/derfx
        work=work+diff
      end do
    end if
    erfinv3=work

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function erfinv3

end module special
