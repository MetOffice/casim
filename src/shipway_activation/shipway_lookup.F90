module shipway_lookup

  use variable_precision, only: wp

  implicit none

  private

  character(len=*), parameter, private :: ModuleName='SHIPWAY_LOOKUP'

  real(wp), parameter :: tol=1.e-9

  ! parameters for lookup tables
  real(wp), parameter :: xmin=1.e-4  ! May need to adjust if parameters go out of range
  real(wp), parameter :: xmax=20000.

  real(wp), parameter :: ymin=.1
  real(wp), parameter :: ymax=3.

  integer, parameter :: nx=500
  integer, parameter :: ny=10

  integer, parameter :: nterms=10000 ! We can afford to be generous

  logical :: l_generate_tables=.true.

  real(wp) :: x_param, y_param
  real(wp) :: xvalues(nx)
  real(wp) :: yvalues(ny)
  real(wp) :: J_table(nx,ny)
  
  public lookup_I, xmax, xmin, ymax, ymin, generate_tables

contains

  subroutine generate_tables(filename)
    ! NB The lookup table is stored as I(x,y)/x^2, then the x^2 
    ! term is restored after the interpolation is done in lookup_I  

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(*), optional, intent(in) :: filename
    real(wp) :: dlx,dly
    logical :: fexist
    logical :: l_calculate=.true., l_writeout=.false.
    integer :: i,j 

    character(len=*), parameter :: RoutineName='GENERATE_TABLES'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    if (present(filename))then
      INQUIRE(FILE=filename, EXIST=fexist)
      if (fexist)then
        open(200, file=trim(filename), STATUS='old')
        read(200, *) xvalues
        read(200, *) yvalues
        read(200, *) J_table
        l_calculate=.false.
      else
        l_writeout=.true.
      end if
    end if

    if (l_calculate)then

      dlx=(log(xmax)-log(xmin))/(nx-1)
      dly=(log(ymax)-log(ymin))/(ny-1)
      do i=1,nx
        xvalues(i) = xmin*exp((i-1)*dlx)
      end do

      do j=1,ny
        yvalues(j) = ymin*exp((j-1)*dly)
      end do

      do j=1,ny
        do i=1,nx
          J_table(i,j) = J1(xvalues(i),yvalues(j),nterms)/(xvalues(i)*xvalues(i))
        end do
      end do

      ! Quick hack to ensure J1 is monotonic increasing.
      do j=1,ny
        do i=2,nx
          if (J_table(i,j) <= J_table(i-1,j)) J_table(i,j) = J_table(i-1,j) + tiny(xmin)
        end do
      end do

      if (l_writeout)then
        open(200, file=trim(filename), STATUS='new')
        write(200, *) xvalues
        write(200, *) yvalues
        write(200, *) J_table
      end if
    end if

    l_generate_tables=.false.

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine generate_tables

  real(wp) function J1_integrand(t)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! See equation 21 of Shipway (2015)    
    real(wp), intent(in) :: t

    real(wp) :: x, y, expnt

    ! local variables

    character(len=*), parameter :: RoutineName='J1_INTEGRAND'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    x=x_param
    y=y_param

    ! Limit exponent to prevent underflow - if it's small enough
    ! we don't care about the loss of accuracy
    ! This allows for the factor multiplying the exponential term to be
    ! a minimum of e**(-20.) = 2.e-9
    ! Not that we still expect some floating underflow in the
    ! part of the parameter space where the integral < tiny(1.0)
    ! This doesn't really matter since the left had side will presumabely be > tiny(1.0)!
!    expnt=-0.5*log(t)*log(t)/(y*y)
    expnt=max(-0.5*log(t)*log(t)/(y*y), (minexponent(t))*log(2.0) + 20.0)
    J1_integrand = x*sqrt(x*x-t*t)/t*exp(expnt) &
       /sqrt(((x**3-t**3)/x**3)**0.6)/sqrt(0.5)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


  end function J1_integrand

  real(wp) function J1(x, y, nterms)
    
    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(in) :: x, y
    integer, intent(in) :: nterms
    character(len=*), parameter :: RoutineName='J1'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    x_param=x
    y_param=y
    
    J1 = simpson(J1_integrand, tol, x - tol, nterms)

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function J1

  subroutine lookup_I(x,y,J1)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! NB The lookup table is stored as I(x,y)/x^2, then the x^2 
    ! term is restored after the interpolation is done in lookup_I  

    real(wp), intent(in) :: x, y
    real(wp), intent(out) :: J1
    integer :: m(1)
    integer :: ix, iy
    real(wp) :: dx1, dy1, dx2, dy2

    character(len=*), parameter :: RoutineName='LOOKUP_I'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    if (l_generate_tables)then
      call generate_tables()
    end if

    m=minloc((xvalues(:) - x)*(xvalues(:) - x))
    ix=m(1)
    m=minloc((yvalues(:) - y)*(yvalues(:) - y))
    iy=m(1)

    if (x < xvalues(ix)) ix=ix-1
    if (y < yvalues(iy)) iy=iy-1
    
    dx1=x-xvalues(ix)
    dx2=xvalues(ix+1)-x
    dy1=y-yvalues(iy)
    dy2=yvalues(iy+1)-y

    J1=x*x*(1./((dx1+dx2)*(dy1+dy2))) &
       *(J_table(ix,iy)*dx2*dy2 &
       + J_table(ix+1,iy)*dx1*dy2 &
       + J_table(ix,iy+1)*dx2*dy1 &
       + J_table(ix+1,iy+1)*dx1*dy1)

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine lookup_I


  function simpson(f,a,b,N)
    
    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp) :: f
    real(wp), intent(in) :: a,b
    integer, intent(in) :: N
    
    real(wp) :: s, simpson, fk1,fk2,fk3
    real(wp) :: h
    integer  :: i
    real(wp) :: tol ! to mitigate rounding issues at single precision

    character(len=*), parameter :: RoutineName='SIMPSON'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    s=0
    
    tol=(b-a)*1e-7
    
    h = (b-a-tol)/(2*N)
    fk3=f(a)
    do i=0,N-1
       fk1=fk3
       fk2=f(a+(2*i+1)*h)
       fk3=f(a+2*(i+1)*h)
       s=s+fk1+4*fk2+fk3
    end do
    
    simpson=h*s/3.

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function simpson

end module shipway_lookup
