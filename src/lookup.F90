module lookup
  use mphys_die, only: throw_mphys_error, bad_values, std_msg
  use variable_precision, only: wp
  use mphys_switches, only: diag_mu_option, l_passive3m, max_mu, l_sm_fix_n0
  use mphys_parameters, only: hydro_params, c_r, d_r, p1, p2, p3
  use mphys_constants, only: fixed_rain_number, fixed_rain_mu
  use special, only: GammaFunc

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='LOOKUP'

  integer, parameter :: nmu=501
  real(wp) :: min_mu = 0.0  ! ( max_mu = 35 )
  real(wp), allocatable :: mu_g(:), mu_i(:)
  real(wp), allocatable :: mu_g_sed(:), mu_i_sed(:)

  interface get_lam_n0
     module procedure get_lam_n0_3M, get_lam_n0_2M, get_lam_n0_1M
  end interface get_lam_n0

  public Gfunc, get_slope_generic, moment, get_n0, get_mu, get_lam_n0, set_mu_lookup, mu_i, mu_g, mu_i_sed, mu_g_sed, nmu
contains

  function Gfunc(mu, p1, p2, p3)

    implicit none

    character(len=*), parameter :: RoutineName='GFUNC'

    real(wp), intent(in) :: mu, p1, p2, p3

    real(wp) :: Gfunc, GfuncL
    real(wp) :: k1, k2, k3

    k3=p2-p1
    k1=(p3-p2)/k3
    k2=(p1-p3)/k3

    !    GfuncL=exp(k1*log(GammaFunc(1.+mu+p1)) &
    !         +k2*log(GammaFunc(1.+mu+p2)) &
    !         +log(GammaFunc(1.+mu+p3)) &
    !         )

    Gfunc=GammaFunc(1.0+mu+p1)**k1*GammaFunc(1.0+mu+p2)**k2*GammaFunc(1.0+mu+p3)
  end function Gfunc

  function Hfunc(m1,m2,m3, p1, p2, p3)

    implicit none

    character(len=*), parameter :: RoutineName='HFUNC'

    real(wp), intent(in) :: m1,m2,m3
    real(wp), intent(in) :: p1,p2,p3
    real(wp) :: Hfunc
    real(wp) :: k1, k2, k3

    k3=p2-p1
    k1=(p3-p2)/k3
    k2=(p1-p3)/k3

    Hfunc=exp(k1*log(m1)+k2*log(m2)+log(m3))
  end function Hfunc

  subroutine set_mu_lookup(p1, p2, p3, index, value)

    implicit none

    character(len=*), parameter :: RoutineName='SET_MU_LOOKUP'

    real(wp), intent(in) :: p1, p2, p3
    real(wp), intent(inout) :: index(:), value(:)

    integer :: i, lb, ub

    lb=lbound(index,1)
    ub=ubound(index,1)

    do i=lb, ub
      index(i)=min_mu+(max_mu-min_mu)/(nmu-1.0)*(i-1)
      value(i)=Gfunc(index(i), p1, p2, p3)
    end do
  end subroutine set_mu_lookup

  subroutine get_slope_generic(k, params, n0, lam, mu, mass, number, m3)

    implicit none

    character(len=*), parameter :: RoutineName='GET_SLOPE_GENERIC'

    integer, intent(in) :: k
    type(hydro_params), intent(in) :: params
    real(wp), intent(out) :: n0, lam, mu
    real(wp), intent(in) :: mass
    real(wp), intent(in), optional :: number, m3

    real(wp) :: m1, m2, p1, p2, p3

    m1=mass/params%c_x
    p1=params%p1
    p2=params%p2
    p3=params%p3
    
    if (params%l_3m) then
      if (l_passive3m) then
        m2=number
        mu=params%fix_mu
        call get_lam_n0(m1, m2, p1, p2, mu, lam, n0)
      else
        m2=number
        call get_mu(m1, m2, m3, p1, p2, p3, mu)
        call get_lam_n0(m1, m2, m3, p1, p2, p3, mu, lam, n0)
      end if
    else if (params%l_2m) then
      m2=number
      mu=params%fix_mu
      call get_lam_n0(m1, m2, p1, p2, mu, lam, n0)
    else
      mu=params%fix_mu
      n0=params%fix_n0
      call get_lam_n0(m1, p1, mu, lam, n0)
    end if
    if (lam <= 0) then
      write(std_msg, *) 'ERROR in lookup', params%id, params%i_2m, m1, number, lam
      call throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)
    end if
  end subroutine get_slope_generic

  subroutine get_mu(m1, m2, m3, p1, p2, p3, mu, mu_g_o, mu_i_o)

    implicit none

    character(len=*), parameter :: RoutineName='GET_MU'

    real(wp), intent(in) :: m1, m2, m3
    real(wp), intent(in) :: p1, p2, p3
    real(wp), intent(in), optional :: mu_g_o(:), mu_i_o(:)
    real(wp), intent(out) :: mu

    real(wp) :: G
    real(wp) :: muG(nmu), muI(nmu)
    integer  :: i
    real(wp) :: k1, k2, k3, pos

    if (present(mu_g_o) .and. present(mu_i_o)) then
      muG=mu_g_o
      muI=mu_i_o
    else
      muG=mu_g
      muI=mu_i
    end if

    k1=p3-p2
    k2=p1-p3
    k3=p2-p1

    pos=sign(1.0_wp,k1*k2) ! scaling to +-1 preserves precision for < condition

    G=Hfunc(m1,m2,m3,p1,p2,p3)

    if (pos*G < pos*muG(1)) then
      mu=-1.0e-10  ! set to be small and negative so that it is picked up
      ! in checks later on
    else
      do i=1, nmu-1
        if ((muG(i)-G)*(muG(i+1)-G) <= 0.0) exit
      end do

      if (i==nmu) then
        mu=max_mu
      else
        mu=muI(i)+(G-muG(i))/(muG(i+1)-muG(i))*(muI(i+1)-muI(i))
      end if
    end if
  end subroutine get_mu

  ! 3M version
  subroutine get_lam_n0_3M(m1, m2, m3, p1, p2, p3, mu, lam, n0)

    implicit none

    character(len=*), parameter :: RoutineName='GET_LAM_N0_3M'
    
    real(wp), intent(in) :: m1, m2, m3, mu
    real(wp), intent(in) :: p1, p2, p3
    real(wp), intent(out) :: lam, n0

    real(wp) :: p, m
    real(wp) :: l2

    !    l2 = 1./(p2-p3)
    !
    !    lam = ((GammaFunc(1.+mu+p2)/GammaFunc(1.+mu+p3)) &
    !         *(m3/m2))**l2
    ! Make sure we use m1 and m2 to calculate lambda, in case we've gone
    ! out of bounds with mu and need to modify m3.
    l2=1.0/(p2-p1)

    lam=((GammaFunc(1.0+mu+p2)/GammaFunc(1.0+mu+p1))*(m1/m2))**l2

    m=m2
    p=p2

    n0=lam**(p)*m*GammaFunc(1.0+mu)/GammaFunc(1.0+mu+p)
  end subroutine get_lam_n0_3M

  ! 2M version
  subroutine get_lam_n0_2M(m1, m2, p1, p2, mu, lam, n0)

    implicit none

    character(len=*), parameter :: RoutineName='GET_LAM_N0_2M'
    
    real(wp), intent(in) :: m1, m2, mu
    real(wp), intent(in) :: p1, p2
    real(wp), intent(out) :: lam, n0

    real(wp) :: p, m
    real(wp) :: j1

    j1=1.0/(p1-p2)

    lam=((GammaFunc(1.0+mu+p1)/GammaFunc(1.0+mu+p2))*(m2/m1))**(j1)

    m=m2
    p=p2

    n0=lam**(p)*m*GammaFunc(1.0+mu)/GammaFunc(1.0+mu+p)
  end subroutine get_lam_n0_2M

  ! 1M version
  subroutine get_lam_n0_1M(m1, p1, mu, lam, n0)

    implicit none

    character(len=*), parameter :: RoutineName='GET_LAM_N0_1M'
    
    real(wp), intent(in) :: m1, mu, n0
    real(wp), intent(in) :: p1
    real(wp), intent(out) :: lam

    ! 
    ! Fixing Nx is equivalent to having n0=na*lam**(1+mu)
    ! (c.f. LEM formulation, na, nb)
    ! if we want to fix na and nb, we can do this and lambda is
    ! given by:
    ! lam=(na*gamma(1+mu+p1)/m1)**(1/(1+mu+p1-nb))
    !

    lam=(n0*GammaFunc(1.0+mu+p1)/GammaFunc(1.0+mu)*m1**(-1.0))**(1.0/p1)
  end subroutine get_lam_n0_1M

  ! Get n0 given a moment and lamda and mu
  subroutine get_n0(m, p, mu, lam, n0)

    implicit none

    character(len=*), parameter :: RoutineName='GET_N0'
    
    real(wp), intent(in) :: m, mu, lam
    real(wp), intent(in) :: p
    real(wp), intent(out) :: n0

    n0=m/(GammaFunc(1+mu+p)*lam**(-p)/GammaFunc(1+mu))
  end subroutine get_n0

  function moment(n0,lam,mu,p)

    implicit none

    character(len=*), parameter :: RoutineName='MOMENT'
  
    real(wp), intent(in) :: n0, lam, mu, p
    real(wp) :: moment

    moment=n0*GammaFunc(1+mu+p)*lam**(-p)/GammaFunc(1+mu)
  end function moment
end module lookup
