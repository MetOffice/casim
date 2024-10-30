module gauss_casim_micro

  !< Calculate the integral needed in the aggregation calculations
  !<
  !< OPTIMISATION POTENTIAL - LOOKUP

  use variable_precision, only: wp
  use mphys_switches, only: max_mu
  implicit none
  private

  integer, parameter :: maxq = 5
  real(wp) :: gaussfunc_save(maxq) ! max 5 values - I.e. this is only used with 2m schemes
  ! should be extended
  integer, parameter, private :: nbins_a = 50
  real(wp) :: gaussfunc_save_2D(nbins_a,maxq)
  logical :: l_save_2D(nbins_a,maxq) = .false.

!$OMP THREADPRIVATE(gaussfunc_save, gaussfunc_save_2D, l_save_2D)

  character(len=*), parameter, private :: ModuleName='GAUSS_CASIM_MICRO'

  public gauss_casim_func, gaussfunclookup, gaussfunclookup_2d
contains

  subroutine gaussfunclookup(iq, gauss_value, a, b)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: iq !< parameter index relating to variable we're considering
    real(wp), intent(out) :: gauss_value !< returned value
    real(wp), intent(in), optional :: a, b  !< Value of a and b to use. Only required if initializing

    ! Local variables
    character(len=*), parameter :: RoutineName='GAUSSFUNCLOOKUP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if ( present(a) .and. present(b)) then
      gauss_value=gauss_casim_func(a,b)
      gaussfunc_save(iq)=gauss_casim_func(a,b)
    else
      gauss_value=gaussfunc_save(iq)
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine gaussfunclookup

  subroutine gaussfunclookup_2d(iq, gauss_value, a, b)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: iq !< parameter index relating to variable we're considering
    real(wp), intent(out) :: gauss_value !< returned value
    real(wp), intent(in) :: a, b  !< Value of a and b to use.  (a is mu, b is b_x)

    ! Local variables
    integer :: ibin ! mu(a) bin in which we sit

    character(len=*), parameter :: RoutineName='GAUSSFUNCLOOKUP_2D'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ibin=int((a/max_mu)*(nbins_a-1))+1
    if (.not. l_save_2D(ibin,iq)) then
      gauss_value=gauss_casim_func(a,b)
      gaussfunc_save_2d(ibin,iq)=gauss_casim_func(a,b)
      l_save_2D(ibin,iq)=.true.
    else
      gauss_value=gaussfunc_save_2d(ibin,iq)
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine gaussfunclookup_2d

  function gauss_casim_func(a, b)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    real(wp), intent(in) :: a, b !< function arguments

    ! Local variables
    real(wp), parameter ::   tmax = 18.0    !< Limit of integration
    real(wp), parameter ::   dt   = 0.08   !< step size
    real(wp) :: gauss_sum, t1, t2
    real(wp) :: Gauss_casim_Func

    character(len=*), parameter :: RoutineName='GAUSS_CASIM_FUNC'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    gauss_sum=0.0
    t1=0.5*dt
    Gaus_t1: do while (t1 <= tmax)
      t2=0.5*dt
      Gaus_t2: do while (t2 <= tmax)
        gauss_sum=gauss_sum+(t1+t2)**2*abs((t1**b)-(t2**b))*(t1**a)*(t2**a)*exp(-(t1+t2))
        t2=t2+dt
      end do Gaus_t2
      t1=t1+dt
    end do Gaus_t1
    Gauss_casim_Func=gauss_sum*dt*dt

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function Gauss_casim_Func
end module gauss_casim_micro
