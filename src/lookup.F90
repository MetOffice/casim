module lookup
  use mphys_die, only: throw_mphys_error, bad_values, std_msg
  use variable_precision, only: wp
  use mphys_switches, only: max_mu, l_kfsm
! use mphys_switches, only: l_passive3m
  use mphys_parameters, only: hydro_params, a_i, b_i, &
                              a_s, b_s
  use special, only: GammaFunc
  use passive_fields, only: rho
!#if DEF_MODEL==UM
  use casim_moments_mod, only: casim_moments
!#endif
  implicit none
  private

  character(len=*), parameter, private :: ModuleName='LOOKUP'

  integer, parameter :: nmu=501
  real(wp) :: min_mu = 0.0  ! ( max_mu = 35 )
  real(wp), allocatable :: mu_g(:), mu_i(:)
  real(wp), allocatable :: mu_g_sed(:), mu_i_sed(:)

  interface get_lam_n0
     module procedure get_lam_n0_3M, get_lam_n0_2M, get_lam_n0_1M, get_lam_n0_1M_KF
  end interface get_lam_n0

  public Gfunc, get_slope_generic, get_slope_generic_kf, moment, get_n0,        &
         get_mu, get_lam_n0, set_mu_lookup, mu_i, mu_g, mu_i_sed, mu_g_sed, nmu
contains

  function Gfunc(mu, p1, p2, p3)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GFUNC'

    real(wp), intent(in) :: mu, p1, p2, p3

    real(wp) :: Gfunc
!   real(wp) :: GfuncL
    real(wp) :: k1, k2, k3

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    k3=p2-p1
    k1=(p3-p2)/k3
    k2=(p1-p3)/k3

    !    GfuncL=exp(k1*log(GammaFunc(1.+mu+p1)) &
    !         +k2*log(GammaFunc(1.+mu+p2)) &
    !         +log(GammaFunc(1.+mu+p3)) &
    !         )

    Gfunc=GammaFunc(1.0+mu+p1)**k1*GammaFunc(1.0+mu+p2)**k2*GammaFunc(1.0+mu+p3)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function Gfunc

  function Hfunc(m1,m2,m3, p1, p2, p3)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='HFUNC'

    real(wp), intent(in) :: m1,m2,m3
    real(wp), intent(in) :: p1,p2,p3
    real(wp) :: Hfunc
    real(wp) :: k1, k2, k3

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    k3=p2-p1
    k1=(p3-p2)/k3
    k2=(p1-p3)/k3

    Hfunc=exp(k1*log(m1)+k2*log(m2)+log(m3))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function Hfunc

  subroutine set_mu_lookup(p1, p2, p3, ind, val)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SET_MU_LOOKUP'

    real(wp), intent(in) :: p1, p2, p3
    real(wp), intent(inout) :: ind(:), val(:)

    integer :: i, lb, ub

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    lb=lbound(ind,1)
    ub=ubound(ind,1)

    do i=lb, ub
      ind(i)=min_mu+(max_mu-min_mu)/(nmu-1.0)*(i-1)
      val(i)=Gfunc(ind(i), p1, p2, p3)
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine set_mu_lookup

  subroutine get_slope_generic(params, n0, lam, mu, mass, num, m3)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GET_SLOPE_GENERIC'

    type(hydro_params), intent(in) :: params
    real(wp), intent(out) :: n0, lam, mu
    real(wp), intent(in) :: mass
    real(wp), intent(in), optional :: num, m3
                                    !! m3 NEEDED for 3rd moment code

    real(wp) :: m1, m2, p1, p2, p3

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    m1=mass/params%c_x
    ! p1=params%p1
    ! p2=params%p2
    ! p3=params%p3

    ! if (params%l_3m) then
    !   if (l_passive3m) then
    !     m2=num
    !     mu=params%fix_mu
    !      call get_lam_n0(m1, m2, params, lam, n0, params%id)
    !   else
    !     m2=num
    !     call get_mu(m1, m2, m3, p1, p2, p3, mu)
    !     call get_lam_n0(m1, m2, m3, params, mu, lam, n0)
    !   end if
    ! elseif (params%l_2m) then
    if (params%l_2m) then
      m2=num
      mu=params%fix_mu
      !call get_lam_n0(m1, m2, p1, p2, mu, lam, n0)
      call get_lam_n0(m1, m2, params, lam, n0, params%id)
    else
      mu=params%fix_mu
      n0=params%fix_n0
      !call get_lam_n0(m1, p1, mu, lam, n0)
      call get_lam_n0(m1, params, lam, n0)
    end if
    if (lam <= 0) then
      write(std_msg, *) 'ERROR in lookup', params%id, params%i_2m, m1, num, lam
      call throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine get_slope_generic

  subroutine get_slope_generic_kf(ixy_inner, k, params, n0, lam, mu, lams, mass, Tk, num, &
                               m3)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GET_SLOPE_GENERIC_KF'

    integer, intent(in) :: ixy_inner
    integer, intent(in) :: k
    type(hydro_params), intent(inout) :: params
    real(wp), intent(out) :: n0, lam, mu
    real(wp), intent(out) :: lams(2)
    real(wp), intent(in) :: mass, Tk
    real(wp), intent(in), optional :: num, m3

    real(wp) :: m1, ms, p1, p2, p3
    real(wp) :: n_p, cficei(1), m_s(1), Ta(1), rhoa(1), qcf(1)
    real(wp) :: j1, j2
    real(wp) :: na, nb

    integer  :: points

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    m1=mass/params%c_x
    p1=params%p1
    p2=params%p2
    p3=params%p3

    if (l_kfsm) then

      j1=p1+params%b_x
      j2=p1+b_i
      ! modify p2 so it is equal to j1
      params%p2 = j1

      points = 1
      n_p = 0.0
      m_s(1)=0.0
      cficei(1) = 1.0
      rhoa(1) = rho(k,ixy_inner)
      qcf(1) = mass
      Ta(1) = Tk

    end if

      mu=params%fix_mu
      n0=params%fix_n0
!#if DEF_MODEL==UM
!AH-KID changes - is this needed?
!prf make n0 fixed param consistent with Field 2017
      if (params%id==3 .and. l_kfsm) then
        qcf(1)=max(qcf(1),1e-8) !stop lam going to zero
        call casim_moments(params,points,rhoa,Ta,qcf,cficei,0.0_wp,m_s) !get concentration
        n0=m_s(1)
        m_s(1)=0.0
      endif
!prf
!AH-KID changes 
!#endif
      call get_lam_n0(m1, params, lam, n0, params%id)

      if (l_kfsm) then
        ! Code for Kalli's single moment work
        if (params%id==2) then

          ! Rain: Abel and Boutle distribution
          na=0.22
          nb=2.2
          n0=na*params%gam_1_mu*lam**(nb-1.0-params%fix_mu)

        else if (params%id==5) then

          ! Graupel: Swann PSD
          ! JW: Note we could eventually look to using Field et al. (2019)
          !     PSD - JAMC.
          na=5.0e25
          nb=-4.0
          n0=na*params%gam_1_mu*lam**(nb-1.0-params%fix_mu)

        else if (params%id==3) then
!#if DEF_MODEL==UM
          ! Ice: need lsp_moments
          call casim_moments(params,points,rhoa,Ta,qcf,cficei,j1,m_s)

          ms=m_s(1)

          call casim_moments(params,points,rhoa,Ta,qcf,cficei,j2,m_s)

          if (m_s(1) < ms*params%a_x/a_i) then
            j1=j2
            ms=m_s(1)
            a_s(k)=a_i
            b_s(k)=b_i
          else
            a_s(k) = params%a_x
            b_s(k) = params%b_x
          end if ! m_s(1)

          ! lams(1) is the slope, while lams(2) is n0
          call get_lam_n0(m1, ms, params, lams(1), lams(2), params%id)
!#endif
        end if ! params%id

      end if ! l_kfsm

    if (lam <= 0) then
      write(std_msg, *) 'ERROR in lookup', params%id, params%i_2m, m1, num, lam, n0, qcf(1)
      call throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine get_slope_generic_kf

  subroutine get_mu(m1, m2, m3, p1, p2, p3, mu, mu_g_o, mu_i_o)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GET_MU'

    real(wp), intent(in) :: m1, m2, m3
    real(wp), intent(in) :: p1, p2, p3
    real(wp), intent(in), optional :: mu_g_o(:), mu_i_o(:)
    real(wp), intent(out) :: mu

    real(wp) :: G
    real(wp) :: muG(nmu), muI(nmu)
    integer  :: i
    real(wp) :: k1, k2, pos

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (present(mu_g_o) .and. present(mu_i_o)) then
      muG=mu_g_o
      muI=mu_i_o
    else
      muG=mu_g
      muI=mu_i
    end if

    k1=p3-p2
    k2=p1-p3

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

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine get_mu

  ! 3M version
  subroutine get_lam_n0_3M(m1, m2, m3, params, mu, lam, n0)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GET_LAM_N0_3M'

    real(wp), intent(in) :: m1, m2, m3, mu
                                !! m3 NEEDED for 3rd moment code
    type(hydro_params), intent(in) :: params
    !real(wp), intent(in) :: p1, p2, p3
    real(wp), intent(out) :: lam, n0

    real(wp) :: p, m
    real(wp) :: l2

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !    l2 = 1./(p2-p3)
    !
    !    lam = ((GammaFunc(1.+mu+p2)/GammaFunc(1.+mu+p3)) &
    !         *(m3/m2))**l2
    ! Make sure we use m1 and m2 to calculate lambda, in case we've gone
    ! out of bounds with mu and need to modify m3.
    l2=1.0/(params%p2-params%p1)

    lam=((GammaFunc(1.0+mu+params%p2)/GammaFunc(1.0+mu+params%p1))*(m1/m2))**l2

    m=m2
    p=params%p2

    n0=lam**(p)*m*GammaFunc(1.0+mu)/GammaFunc(1.0+mu+p)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine get_lam_n0_3M

  ! 2M version
  subroutine get_lam_n0_2M(m1, m2, params, lam, n0, id)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GET_LAM_N0_2M'

    type(hydro_params), intent(in) :: params
    real(wp), intent(in) :: m1, m2 !, mu
    !real(wp), intent(in) :: p1, p2
    real(wp), intent(out) :: lam, n0
    integer, intent(in) :: id

    real(wp) :: p, m

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    if (l_kfsm .and. id == 3) then 
       ! field formulation of ice psd, p2 can vary, so need to calc the gamma and inverse  
       ! each timestep. Note p2 is set outside the call - nasty!
       lam=((params%gam_1_mu_p1/GammaFunc(1.0+params%fix_mu+params%p2)) &
            *(m2/m1))**(1.0/(params%p1-params%p2))

       m=m2
       p=params%p2

       n0=lam**(p)*m*params%gam_1_mu/GammaFunc(1.0+params%fix_mu+p)
    else
       lam = ((params%gam_1_mu_p1/params%gam_1_mu_p2)* &
            (m2/m1))**(params%inv_p1_p2)
       
 
       n0=lam**(params%p2)*m2*params%gam_1_mu/params%gam_1_mu_p2
    endif
 
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine get_lam_n0_2M

  ! 1M version
  subroutine get_lam_n0_1M(m1, params, lam, n0)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GET_LAM_N0_1M'

    type(hydro_params), intent(in) :: params
    real(wp), intent(in) :: m1, n0 
    real(wp), intent(out) :: lam

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !
    ! Fixing Nx is equivalent to having n0=na*lam**(1+mu)
    ! (c.f. LEM formulation, na, nb)
    ! if we want to fix na and nb, we can do this and lambda is
    ! given by:
    ! lam=(na*gamma(1+mu+p1)/m1)**(1/(1+mu+p1-nb))
    !

    lam = (n0*params%gam_1_mu_p1/params%gam_1_mu*m1**(-1.0))**(params%inv_p1)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine get_lam_n0_1M

  subroutine get_lam_n0_1M_KF(m1, params, lam, n0, id)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GET_LAM_N0_1M_KF'

    type(hydro_params), intent(in) :: params
    real(wp), intent(in) :: m1, n0
    !real(wp), intent(in) :: p1
    real(wp), intent(out) :: lam
    integer, intent(in) :: id
    real(wp) :: na, nb

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !
    ! Fixing Nx is equivalent to having n0=na*lam**(1+mu)
    ! (c.f. LEM formulation, na, nb)
    ! if we want to fix na and nb, we can do this and lambda is
    ! given by:
    ! lam=(na*gamma(1+mu+p1)/m1)**(1/(1+mu+p1-nb))
    !

    if (l_kfsm .and. id == 2) then
      ! Rain: Abel and Boutle (2012)
      na  = 0.22
      nb  = 2.2
      lam = (na*params%gam_1_mu_p1*m1**(-1.0)) &
           **(1.0/(1+params%fix_mu+params%p1-nb))
    else if (l_kfsm .and. id == 5) then
      ! Graupel: Swann (LEM); Forbes and Halliwell (2003, internal report)
      na  = 5.0e25
      nb  = -4.0
      lam=(na*params%gam_1_mu_p1*m1**(-1.0)) &
           **(1.0/(1+params%fix_mu+params%p1-nb))
    else
      ! not l_kfsm or species other than rain or graupel
      lam = (n0*params%gam_1_mu_p1/params%gam_1_mu*m1**(-1.0))**(params%inv_p1)
    end if ! id / l_kfsm

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine get_lam_n0_1M_KF

  ! Get n0 given a moment and lamda and mu
  subroutine get_n0(m, p, mu, lam, n0)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GET_N0'

    real(wp), intent(in) :: m, mu, lam
    real(wp), intent(in) :: p
    real(wp), intent(out) :: n0

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    n0=m/(GammaFunc(1+mu+p)*lam**(-p)/GammaFunc(1+mu))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine get_n0

  function moment(n0,lam,mu,p)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='MOMENT'

    real(wp), intent(in) :: n0, lam, mu, p
    real(wp) :: moment

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    moment=n0*GammaFunc(1+mu+p)*lam**(-p)/GammaFunc(1+mu)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function moment
end module lookup
