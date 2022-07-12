module sweepout_rate
  use variable_precision, only: wp
  use special, only: pi, Gammafunc
  use mphys_constants, only: rho0
  use mphys_parameters, only: hydro_params
  
  implicit none
  private
  character(len=*), parameter, private :: ModuleName='SWEEPOUT_RATE'

  public sweepout, binary_collection, sweepout_1M2M, binary_collection_1M2M
contains

  ! Calculate the sweepout rate given the distribution
  ! and fallspeed parameters
  function sweepout(n0, lam, mu, params, rho, mass_weight)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SWEEPOUT'

    real(wp), intent(in) :: mu, n0, lam
    type(hydro_params), intent(in) :: params
    real(wp), intent(in) :: rho ! air density
    ! if present and true use mass-weighted sweepout
    logical, intent(in), optional :: mass_weight
    real(wp) :: sweepout

    ! local variables
    real(wp) :: G3_b_mu, G1_mu !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    real(wp) :: arg3 !< argument for gamma function
    real(wp) :: coef !< coefficient = c_x if mass-weighting

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    arg3=3.0+params%b_x+mu
    coef=1.0
    if (present(mass_weight)) then
      if (mass_weight) then
        arg3=arg3+params%d_x
        coef=params%c_x
      end if
    end if

    G3_b_mu=GammaFunc(arg3)
    G1_mu=GammaFunc(1.0+mu)
    sweepout=coef*(pi*n0*params%a_x/4.0)*G3_b_mu/G1_mu* (1.0 + params%f_x/lam)**(-arg3)&
         *lam**(1 + mu - arg3)* (rho0/rho)**(params%g_x)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function sweepout

 ! Calculate the sweepout rate given the distribution
  ! and fallspeed parameters
  function sweepout_1M2M(n0, lam, params, rho, mass_weight)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SWEEPOUT_1M2M'

    real(wp), intent(in) :: n0, lam
    type(hydro_params), intent(in) :: params
    real(wp), intent(in) :: rho ! air density
    ! if present and true use mass-weighted sweepout
    logical, intent(in), optional :: mass_weight
    real(wp) :: sweepout_1M2M

    ! local variables
    real(wp) :: G3_b_mu, G1_mu !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    real(wp) :: arg3 !< argument for gamma function
    real(wp) :: coef !< coefficient = c_x if mass-weighting

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    arg3=3.0+params%b_x+params%fix_mu
    G3_b_mu=params%gam_3_bx_mu
    coef=1.0

    if (present(mass_weight)) then
      if (mass_weight) then
         arg3=arg3+params%d_x
         G3_b_mu=params%gam_3_bx_mu_dx
         coef=params%c_x
      end if
    end if

    G1_mu=params%gam_1_mu

    sweepout_1M2M=coef*(pi*n0*params%a_x/4.0)*G3_b_mu/G1_mu* (1.0 + params%f_x/lam)**(-arg3)&
         *lam**(1 + params%fix_mu - arg3)* (rho0/rho)**(params%g_x)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function sweepout_1M2M

  ! Calculate the number of species Y collected by species X through
  ! binary collisions as both species sediment
  function binary_collection(n0_Xin, lam_X, mu_X, n0_Yin, lam_Y, mu_Y,   &
       params_X, params_Y, rho, mass_weight)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='BINARY_COLLECTION'

    real(wp), intent(in) :: mu_X, n0_Xin, lam_X
    real(wp), intent(in) :: mu_Y, n0_Yin, lam_Y
    real(wp) :: n0_X, n0_Y
    type(hydro_params), intent(in) :: params_X
    type(hydro_params), intent(in) :: params_Y
    real(wp), intent(in) :: rho ! air density
    ! if present and true use mass-weighted value (i.e. total mass accreted)
    logical, intent(in), optional :: mass_weight
    real(wp) :: binary_collection

    ! local variables
    real(wp) :: G1_X, G2_X, G3_X  !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    real(wp) :: G1_Y, G2_Y, G3_Y  !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    real(wp) :: arg1_X !< argument for gamma function
    real(wp) :: arg2_X !< argument for gamma function
    real(wp) :: arg3_X !< argument for gamma function
    real(wp) :: arg1_Y !< argument for gamma function
    real(wp) :: arg2_Y !< argument for gamma function
    real(wp) :: arg3_Y !< argument for gamma function
    real(wp) :: l_X1 !< lam_X^-arg1_X
    real(wp) :: l_X2 !< lam_X^-arg2_X
    real(wp) :: l_X3 !< lam_X^-arg3_X
    real(wp) :: l_Y1 !< lam_Y^-arg1_Y
    real(wp) :: l_Y2 !< lam_Y^-arg2_Y
    real(wp) :: l_Y3 !< lam_Y^-arg3_Y

    real(wp) :: coef !< coefficient = c_x if mass-weighting
    real(wp) :: V_X  !< mass-weighted fall velocity for X
    real(wp) :: V_Y  !< mass-weighted fall velocity for Y
    real(wp) :: delV !< bulk fall-speed differential

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    arg1_X=1.0+mu_X
    arg2_X=2.0+mu_X
    arg3_X=3.0+mu_X
    arg1_Y=1.0+mu_Y
    arg2_Y=2.0+mu_Y
    arg3_Y=3.0+mu_Y
    coef=1.0
    if (present(mass_weight)) then
      if (mass_weight) then
        arg1_Y=arg1_Y+params_Y%d_x
        arg2_Y=arg2_Y+params_Y%d_x
        arg3_Y=arg3_Y+params_Y%d_x
        coef=params_Y%c_x
      end if
    end if

    G1_X=GammaFunc(arg1_X)
    G2_X=GammaFunc(arg2_X)
    G3_X=GammaFunc(arg3_X)
    G1_Y=GammaFunc(arg1_Y)
    G2_Y=GammaFunc(arg2_Y)
    G3_Y=GammaFunc(arg3_Y)

    l_X1=lam_X**(-arg1_X)
    l_X2=lam_X**(-arg2_X)
    l_X3=lam_X**(-arg3_X)
    l_Y1=lam_Y**(-arg1_Y)
    l_Y2=lam_Y**(-arg2_Y)
    l_Y3=lam_Y**(-arg3_Y)

    n0_X=n0_Xin *lam_X**(mu_X+1.0)/GammaFunc(1.0+mu_X)
    n0_Y=n0_Yin *lam_Y**(mu_Y+1.0)/GammaFunc(1.0+mu_Y)
    
    V_X=params_X%a_x * lam_X**(-params_X%b_x)*(rho0/rho)**(params_X%g_x)     &
         *GammaFunc(1.0+mu_X+params_X%d_x+params_X%b_x)/GammaFunc(1.0+mu_X+params_X%d_x)

    V_Y=params_Y%a_x * lam_Y**(-params_Y%b_x)*(rho0/rho)**(params_Y%g_x)     &
         *GammaFunc(1.0+mu_Y+params_Y%d_x+params_Y%b_x)/GammaFunc(1.0+mu_Y+params_Y%d_x)

    delV=max(max(V_X,V_Y)/4.0, abs(V_X-V_Y))

    binary_collection=coef*0.25*pi*n0_X*n0_Y*delV*(l_X3*l_Y1*G1_Y*G3_X+2.0*l_X2*l_Y2*G2_Y*G2_X+l_X1*l_Y3*G3_Y*G1_X)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function binary_collection

  function binary_collection_1M2M(n0_Xin, lam_X, n0_Yin, lam_Y,   &
       params_X, params_Y, rho, mass_weight)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='BINARY_COLLECTION_1M2M'

    real(wp), intent(in) :: n0_Xin, lam_X
    real(wp), intent(in) :: n0_Yin, lam_Y
    real(wp) :: n0_X, n0_Y, mu_X, mu_Y
    
    type(hydro_params), intent(in) :: params_X
    type(hydro_params), intent(in) :: params_Y
    real(wp), intent(in) :: rho ! air density
    ! if present and true use mass-weighted value (i.e. total mass accreted)
    logical, intent(in), optional :: mass_weight
    real(wp) :: binary_collection_1M2M

    ! local variables
    real(wp) :: G1_X, G2_X, G3_X  !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    real(wp) :: G1_Y, G2_Y, G3_Y  !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    real(wp) :: arg1_X !< argument for gamma function
    real(wp) :: arg2_X !< argument for gamma function
    real(wp) :: arg3_X !< argument for gamma function
    real(wp) :: arg1_Y !< argument for gamma function
    real(wp) :: arg2_Y !< argument for gamma function
    real(wp) :: arg3_Y !< argument for gamma function
    real(wp) :: l_X1 !< lam_X^-arg1_X
    real(wp) :: l_X2 !< lam_X^-arg2_X
    real(wp) :: l_X3 !< lam_X^-arg3_X
    real(wp) :: l_Y1 !< lam_Y^-arg1_Y
    real(wp) :: l_Y2 !< lam_Y^-arg2_Y
    real(wp) :: l_Y3 !< lam_Y^-arg3_Y

    real(wp) :: coef !< coefficient = c_x if mass-weighting
    real(wp) :: V_X  !< mass-weighted fall velocity for X
    real(wp) :: V_Y  !< mass-weighted fall velocity for Y
    real(wp) :: delV !< bulk fall-speed differential

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    mu_X=params_X%fix_mu
    mu_Y=params_Y%fix_mu

    arg1_X=1.0+params_X%fix_mu
    arg2_X=2.0+params_X%fix_mu
    arg3_X=3.0+params_X%fix_mu
    arg1_Y=1.0+params_Y%fix_mu
    arg2_Y=2.0+params_Y%fix_mu
    arg3_Y=3.0+params_Y%fix_mu
    
    G1_X=params_X%gam_1_mu
    G2_X=params_X%gam_2_mu
    G3_X=params_X%gam_3_mu
    G1_Y=params_Y%gam_1_mu
    G2_Y=params_Y%gam_2_mu
    G3_Y=params_Y%gam_3_mu
    coef=1.0
    if (present(mass_weight)) then
       if (mass_weight) then
          arg1_Y=arg1_Y+params_Y%d_x
          arg2_Y=arg2_Y+params_Y%d_x
          arg3_Y=arg3_Y+params_Y%d_x
          G1_Y=params_Y%gam_1_mu_dx
          G2_Y=params_Y%gam_2_mu_dx
          G3_Y=params_Y%gam_3_mu_dx
          coef=params_Y%c_x
       end if
    end if
        
    l_X1=lam_X**(-arg1_X)
    l_X2=lam_X**(-arg2_X)
    l_X3=lam_X**(-arg3_X)
    l_Y1=lam_Y**(-arg1_Y)
    l_Y2=lam_Y**(-arg2_Y)
    l_Y3=lam_Y**(-arg3_Y)

    n0_X=n0_Xin *lam_X**(mu_X+1.0)/GammaFunc(1.0+mu_X)
    n0_Y=n0_Yin *lam_Y**(mu_Y+1.0)/GammaFunc(1.0+mu_Y)

    V_X=params_X%a_x * lam_X**(-params_X%b_x)*(rho0/rho)**(params_X%g_x)     &
         *params_X%gam_1_mu_dx_bx/params_X%gam_1_mu_dx
         !*GammaFunc(1.0+mu_X+params_X%d_x+params_X%b_x)/GammaFunc(1.0+mu_X+params_X%d_x)

    V_Y=params_Y%a_x * lam_Y**(-params_Y%b_x)*(rho0/rho)**(params_Y%g_x)     &
         *params_Y%gam_1_mu_dx_bx/params_Y%gam_1_mu_dx
         !*GammaFunc(1.0+mu_Y+params_Y%d_x+params_Y%b_x)/GammaFunc(1.0+mu_Y+params_Y%d_x)

    delV=max(max(V_X,V_Y)/4.0, abs(V_X-V_Y))

    binary_collection_1M2M=coef*0.25*pi*n0_X*n0_Y*delV*(l_X3*l_Y1*G1_Y*G3_X+2.0*l_X2*l_Y2*G2_Y*G2_X+l_X1*l_Y3*G3_Y*G1_X)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function binary_collection_1M2M

end module sweepout_rate
