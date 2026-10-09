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
  
  ! Used by gamma_p / inverse_gamma_p
  integer, parameter :: max_terms = 1000   ! cap on series / continued fraction
  real(wp), parameter :: rel_tol = 1.0e-14_wp

  public pi, Gammafunc, casim_erfc, erfinv, gamma_p, inverse_gamma_p
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

  !-----------------------------------------------------------------------
  ! Regularised lower incomplete gamma function and its inverse, used to
  ! truncate the gamma size distributions in the Phillips et al. secondary
  ! ice production collision integrals (ice_multiplication).
  ! References: NIST DLMF (https://dlmf.nist.gov/) 3.10, 8.2, 8.7, 8.9;
  ! Abramowitz and Stegun (1964) 26.2.23 and 26.4.17.
  !-----------------------------------------------------------------------
  !> Regularised lower incomplete gamma function P(a,x), a > 0, x >= 0.
  !> For x < a+1 the power series DLMF 8.7.1 is summed directly; otherwise
  !> P = 1 - Q with Q(a,x) from the continued fraction DLMF 8.9.2,
  !> evaluated with the modified Lentz algorithm (DLMF 3.10).
  function gamma_p(a, x) result(p)

    implicit none

    real(wp), intent(in) :: a, x
    real(wp) :: p

    real(wp), parameter :: tiny_value = 1.0e-300_wp
    real(wp) :: log_prefactor, term, series, denom
    real(wp) :: b_n, a_n, c, d, delta, fraction
    integer :: n

    if (x <= 0.0_wp .or. a <= 0.0_wp) then
      p = 0.0_wp
      return
    end if

    log_prefactor = a*log(x) - x - log_gamma(a)

    if (x < a + 1.0_wp) then
      ! P(a,x) = x^a e^-x / Gamma(a+1) * sum_n x^n / ((a+1)...(a+n))
      term = 1.0_wp/a
      series = term
      denom = a
      do n = 1, max_terms
        denom = denom + 1.0_wp
        term = term*x/denom
        series = series + term
        if (abs(term) < rel_tol*abs(series)) exit
      end do
      p = min(1.0_wp, series*exp(log_prefactor))
    else
      ! Q(a,x) = x^a e^-x / Gamma(a) * 1/(x+1-a- 1(1-a)/(x+3-a- 2(2-a)/(...)))
      b_n = x + 1.0_wp - a
      c = 1.0_wp/tiny_value
      d = 1.0_wp/b_n
      fraction = d
      do n = 1, max_terms
        a_n = -real(n, wp)*(real(n, wp) - a)
        b_n = b_n + 2.0_wp
        d = b_n + a_n*d
        if (abs(d) < tiny_value) d = tiny_value
        c = b_n + a_n/c
        if (abs(c) < tiny_value) c = tiny_value
        d = 1.0_wp/d
        delta = c*d
        fraction = fraction*delta
        if (abs(delta - 1.0_wp) < rel_tol) exit
      end do
      p = max(0.0_wp, 1.0_wp - fraction*exp(log_prefactor))
    end if

  end function gamma_p

  !> Inverse of P(a,x) in x: returns x such that P(a,x) = prob,
  !> for a > 0 and 0 < prob < 1.  The starting value is the Wilson-Hilferty
  !> approximation (A&S 26.4.17) with the normal quantile from A&S 26.2.23;
  !> the root is then bracketed and refined by Newton's method, falling
  !> back to bisection whenever a Newton step would leave the bracket.
  function inverse_gamma_p(prob, a) result(x)

    implicit none

    real(wp), intent(in) :: prob, a
    real(wp) :: x

    integer, parameter :: max_iter = 200
    real(wp) :: q, t, z, x_low, x_high, residual, density, x_new
    integer :: iter

    if (prob <= 0.0_wp .or. a <= 0.0_wp) then
      x = 0.0_wp
      return
    end if

    ! Standard normal quantile z with Phi(z) = prob
    q = min(prob, 1.0_wp - prob)
    q = max(q, 1.0e-300_wp)
    t = sqrt(-2.0_wp*log(q))
    z = t - (2.515517_wp + t*(0.802853_wp + t*0.010328_wp)) /                  &
            (1.0_wp + t*(1.432788_wp + t*(0.189269_wp + t*0.001308_wp)))
    if (prob < 0.5_wp) z = -z

    ! Wilson-Hilferty starting value
    x = a*(1.0_wp - 1.0_wp/(9.0_wp*a) + z/(3.0_wp*sqrt(a)))**3
    if (x <= 0.0_wp) x = 0.5_wp*a

    ! Bracket the root
    x_low = 0.0_wp
    x_high = max(x, a, 1.0_wp)
    do iter = 1, max_iter
      if (gamma_p(a, x_high) >= prob) exit
      x_low = x_high
      x_high = 2.0_wp*x_high
    end do
    if (prob >= 1.0_wp) then
      x = x_high
      return
    end if
    x = min(max(x, x_low), x_high)

    ! Safeguarded Newton iteration
    do iter = 1, max_iter
      residual = gamma_p(a, x) - prob
      if (residual > 0.0_wp) then
        x_high = x
      else
        x_low = x
      end if

      density = 0.0_wp
      if (x > 0.0_wp) density = exp((a - 1.0_wp)*log(x) - x - log_gamma(a))

      if (density > 0.0_wp) then
        x_new = x - residual/density
      else
        x_new = 0.5_wp*(x_low + x_high)
      end if
      if (x_new <= x_low .or. x_new >= x_high) x_new = 0.5_wp*(x_low + x_high)

      if (abs(x_new - x) <= 1.0e-12_wp*max(x_new, 1.0e-30_wp)) then
        x = x_new
        exit
      end if
      x = x_new
    end do

  end function inverse_gamma_p

end module special
