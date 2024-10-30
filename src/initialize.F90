module initialize
  use mphys_die, only: throw_mphys_error, incorrect_opt, std_msg
  use variable_precision, only: wp
  use lookup, only: set_mu_lookup, mu_i, mu_g, mu_i_sed, mu_g_sed, nmu
  use derived_constants, only: set_constants
  use mphys_parameters, only: p1, p2, p3, sp1, sp2, sp3, snow_params
  use mphys_switches, only: aerosol_option, l_warm, cloud_params, &
       rain_params, ice_params, snow_params, graupel_params,&
       iopt_act, l_g, l_sg, l_override_checks, i_am10, i_an10, &
       isol, iinsol, active_rain, active_cloud, aero_index, active_number, process_level, iopt_act, l_process
  use gauss_casim_micro, only: gaussfunclookup
  use micro_main, only : initialise_micromain, finalise_micromain
  use sedimentation, only : initialise_sedr, finalise_sedr

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='INITIALIZE'

  public mphys_init, mphys_finalise

contains 

  subroutine mphys_init(il, iu, jl, ju, kl, ku,                 &       
       is_in, ie_in, js_in, je_in, ks_in, ke_in, l_tendency)


    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: il, iu ! upper and lower i levels
    integer, intent(in) :: jl, ju ! upper and lower j levels
    integer, intent(in) :: kl, ku ! upper and lower k levels

    integer, intent(in), optional :: is_in, ie_in ! upper and lower i levels which are to be used
    integer, intent(in), optional :: js_in, je_in ! upper and lower j levels
    integer, intent(in), optional :: ks_in, ke_in ! upper and lower k levels

    ! New optional l_tendency logical added...
    ! if true then a tendency is returned (i.e. units/s)
    ! if false then an increment is returned (i.e. units/timestep)
    logical, intent(in), optional :: l_tendency

    ! Local variables
    real(wp) :: tmp

    integer :: i_start, i_end ! upper and lower i levels which are to be used
    integer :: j_start, j_end ! upper and lower j levels
    integer :: k_start, k_end ! upper and lower k levels
    logical :: l_tendency_loc

    character(len=*), parameter :: RoutineName='MPHYS_INIT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP SINGLE

    call check_options()
    call set_constants()

!$OMP END SINGLE

    ! Set grid extents to operate on
    if (present(is_in)) i_start=is_in
    if (present(ie_in)) i_end=ie_in
    if (present(js_in)) j_start=js_in
    if (present(je_in)) j_end=je_in
    if (present(ks_in)) k_start=ks_in
    if (present(ke_in)) k_end=ke_in

    ! if not passed in, then default to full grid
    if (.not. present(is_in)) i_start=il
    if (.not. present(ie_in)) i_end=iu
    if (.not. present(js_in)) j_start=jl
    if (.not. present(je_in)) j_end=ju
    if (.not. present(ks_in)) k_start=kl
    if (.not. present(ke_in)) k_end=ku
    
    if (present(l_tendency)) then
      l_tendency_loc=l_tendency
    else
      l_tendency_loc=.true.
    end if

    call initialise_micromain(il, iu, jl, ju, kl, ku, i_start, i_end, j_start, j_end, k_start, k_end, l_tendency_loc)

    call initialise_sedr()

!$OMP SINGLE
    call initialise_lookup_tables()
!$OMP END SINGLE

    !Gaussfunc
    call gaussfunclookup(snow_params%id, tmp, a=snow_params%fix_mu, b=snow_params%b_x)

!$OMP SINGLE
    !Initialise the gamma function so not calced on every timestep
    call gamma_initialize()
!$OMP END SINGLE

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine mphys_init

  subroutine mphys_finalise()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Local variables
    character(len=*), parameter :: RoutineName='MPHYS_FINALISE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    call finalise_sedr()
    call finalise_micromain()

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine mphys_finalise

  subroutine initialise_lookup_tables()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Local variables
    character(len=*), parameter :: RoutineName='INITIALISE_LOOKUP_TABLES'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !------------------------------------------------------
    ! Look up tables... (These need to be extended)
    !------------------------------------------------------
    ! mu lookup

    allocate(mu_i(nmu))
    allocate(mu_g(nmu))
    call set_mu_lookup(p1, p2, p3, mu_i, mu_g)
    allocate(mu_i_sed(nmu))
    allocate(mu_g_sed(nmu))
    call set_mu_lookup(sp1, sp2, sp3, mu_i_sed, mu_g_sed)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine initialise_lookup_tables

  ! Check that the options that have been selected are consitent
  subroutine check_options()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Local variables
    character(len=*), parameter :: RoutineName='CHECK_OPTIONS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (.not. l_override_checks) then
      if (l_warm .and. l_process) then
        write(std_msg, '(A)') 'processing does not currently work with l_warm=.true.'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if

      if (.not. ice_params%l_1m .and. .not. l_warm) then
        write(std_msg, '(A)') 'l_warm must be true if not using ice'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if

      if (cloud_params%l_2m .and. aerosol_option == 0     &
           .and. (iopt_act== 1 .or. iopt_act==3)) then
        write(std_msg, '(A)') 'for double moment cloud you must have aerosol_option>0'// &
             'or else activation should be independent of aerosol'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
             
      end if

      if (.not. cloud_params%l_2m .and. aerosol_option > 0) then
        write(std_msg, '(A)') 'aerosol_option must be 0 if not using double moment microphysics'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if

      if (l_g .and. .not. l_sg) then
        write(std_msg, '(A)') 'Cannot run with graupel but not with snow'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
        l_g=.false.
      end if

      ! Some options are not yet working or well tested so don't let anyone use these

      if (i_am10 > 0 .or. i_an10 > 0) then
        write(std_msg, '(A)') 'Accumulation mode dust is not yet used.'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if

      ! Aerosol consistency
      if (active_rain(isol) .and. .not. active_cloud(isol)) then
        write(std_msg, '(A)') 'active_rain(isol) .and. .not. active_cloud(isol)'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if

      ! Aerosol consistency
      if (aero_index%i_accum == 0 .and. aero_index%i_coarse==0 .and. &
         (iopt_act==1 .or. iopt_act==3) )                            then
        write(std_msg, '(A)') 'Must have accumulation or coarse mode aerosol '//&
                              'for chosen activation option'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if

      ! Aerosol consistency
      if (aero_index%i_accum == 0 .and. aero_index%nccn /= 1 .and. &
          active_number(isol))                                     then
      
        write(std_msg, '(A)') 'Soluble modes must only be accumulation mode '//&
                              'with soluble active_number'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if

      ! Aerosol consistency
      if (aero_index%nin > 0) then 
         if (aero_index%i_accum_dust == 0 .and. aero_index%nin /= 1 .and. &
              active_number(iinsol))                                  then
            write(std_msg, '(A)') 'Dust modes must only be accumulation mode '//&
                 'with insoluble active_number'
            call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                 std_msg)
         end if
      endif

      ! Aerosol processing consistency
      if (process_level > 0 .and. iopt_act < 3) then
        write(std_msg, '(A)') 'If processing aerosol, must use higher '//&
                              'level activation code, i.e check iopt_act'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine check_options
  
  subroutine gamma_initialize()
    ! Purpose: When CASIM is used in 1M or 2M mode with a fixed shape parameter, 
    !          all the gamma functions can be calculated at the beginning of the job.
    !          This routine does this calculation so that results can be used in 
    !          sedimentation, lookup...
    
    !use mphys_parameters, only: cloud_params, rain_params, ice_params, snow_params, &
    !     graupel_params
    use special, only: Gammafunc

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='GAMMA_INITIALIZE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! moment gamma functions

    ! Cloud
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    cloud_params%gam_1_mu_p1 = GammaFunc(1.0_wp+cloud_params%fix_mu+cloud_params%p1)
    cloud_params%gam_1_mu_p2 = GammaFunc(1.0_wp+cloud_params%fix_mu+cloud_params%p2)
    cloud_params%gam_1_mu = GammaFunc(1.0_wp+cloud_params%fix_mu)
    ! Moment exponent for lamba
    cloud_params%inv_p1_p2 = (1.0_wp)/(cloud_params%p1 - cloud_params%p2)
    cloud_params%inv_p1 = (1.0_wp)/(cloud_params%p1)

    ! Rain
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    rain_params%gam_1_mu_p1 = GammaFunc(1.0_wp+rain_params%fix_mu+rain_params%p1)
    rain_params%gam_1_mu_p2 = GammaFunc(1.0_wp+rain_params%fix_mu+rain_params%p2)
    rain_params%gam_1_mu = GammaFunc(1.+rain_params%fix_mu)
    ! Moment exponent for lamba
    rain_params%inv_p1_p2 = (1.0_wp)/(rain_params%p1 - rain_params%p2)
    rain_params%inv_p1 = (1.0_wp)/(rain_params%p1)

    ! Ice
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    ice_params%gam_1_mu_p1 = GammaFunc(1.0_wp+ice_params%fix_mu+ice_params%p1)
    ice_params%gam_1_mu_p2 = GammaFunc(1.0_wp+ice_params%fix_mu+ice_params%p2)
    ice_params%gam_1_mu = GammaFunc(1.0_wp+ice_params%fix_mu)
    ! Moment exponent for lamba
    ice_params%inv_p1_p2 = (1.0_wp)/(ice_params%p1 - ice_params%p2)
    ice_params%inv_p1 = (1.0_wp)/(ice_params%p1)

    ! Snow
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    snow_params%gam_1_mu_p1 = GammaFunc(1.0_wp+snow_params%fix_mu+snow_params%p1)
    snow_params%gam_1_mu_p2 = GammaFunc(1.0_wp+snow_params%fix_mu+snow_params%p2)
    snow_params%gam_1_mu = GammaFunc(1.0_wp+snow_params%fix_mu)
    ! Moment exponent for lamba
    snow_params%inv_p1_p2 = (1.0_wp)/(snow_params%p1 - snow_params%p2)
    snow_params%inv_p1 = (1.0_wp)/(snow_params%p1)

    ! Graupel
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    graupel_params%gam_1_mu_p1 = GammaFunc(1.0_wp+graupel_params%fix_mu+graupel_params%p1)
    graupel_params%gam_1_mu_p2 = GammaFunc(1.0_wp+graupel_params%fix_mu+graupel_params%p2)
    graupel_params%gam_1_mu = GammaFunc(1.0_wp+graupel_params%fix_mu)
    ! Moment exponent for lamba
    graupel_params%inv_p1_p2 = (1.0_wp)/(graupel_params%p1 - graupel_params%p2)
    graupel_params%inv_p1 = (1.0_wp)/(graupel_params%p1)   

    ! gamma functions for sedimentation moments and ventilation
    
    ! Cloud
    ! gamma function constants (used in sedimentation)
    cloud_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+cloud_params%fix_mu+cloud_params%sp1)
    cloud_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+cloud_params%fix_mu+cloud_params%sp1+cloud_params%b_x)
    cloud_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+cloud_params%fix_mu+cloud_params%sp2)
    cloud_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+cloud_params%fix_mu+cloud_params%sp2+cloud_params%b_x)
    ! exponent of lambda for sedimentation
    cloud_params%exp_1_mu_sp1 = 1.0+cloud_params%fix_mu+cloud_params%sp1
    cloud_params%exp_1_mu_sp1_bx = 1.0 + cloud_params%fix_mu + cloud_params%sp1 + cloud_params%b_x
    cloud_params%exp_1_mu_sp2 = 1.0+cloud_params%fix_mu+cloud_params%sp2
    cloud_params%exp_1_mu_sp2_bx = 1.0 + cloud_params%fix_mu + cloud_params%sp2 + cloud_params%b_x
    ! gamma function used in ventilation (not really needed for cloud, just added for completeness)
    cloud_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*cloud_params%b_x+cloud_params%fix_mu+2.5)

    ! Rain
    ! gamma function constants (used in sedimentation)
    rain_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp1)
    rain_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp1+rain_params%b_x)
    rain_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp2)
    rain_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp2+rain_params%b_x)
    ! gamma function constants for Abel and shipway (used in sedimentation)
    rain_params%gam_1_mu_sp1_b2x = & 
         Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp1+rain_params%b2_x)
    rain_params%gam_1_mu_sp2_b2x = & 
         Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp2+rain_params%b2_x)
    ! exponent of lambda for sedimentation
    rain_params%exp_1_mu_sp1 = 1.0+rain_params%fix_mu+rain_params%sp1
    rain_params%exp_1_mu_sp1_bx = 1.0 + rain_params%fix_mu + rain_params%sp1 + rain_params%b_x
    rain_params%exp_1_mu_sp2 = 1.0+rain_params%fix_mu+rain_params%sp2
    rain_params%exp_1_mu_sp2_bx = 1.0 + rain_params%fix_mu + rain_params%sp2 + rain_params%b_x
    ! exponent of lambda for Abel and Shipway sedimentation
    rain_params%exp_1_mu_sp1_b2x = &
         1.0 + rain_params%fix_mu + rain_params%sp1 + rain_params%b2_x
    rain_params%exp_1_mu_sp2_b2x = &
         1.0 + rain_params%fix_mu + rain_params%sp2 + rain_params%b2_x
    ! gamma function used in ventilation
    rain_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*rain_params%b_x+rain_params%fix_mu+2.5)

    ! Ice
    ! gamma function constants (used in sedimentation)
    ice_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+ice_params%fix_mu+ice_params%sp1)
    ice_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+ice_params%fix_mu+ice_params%sp1+ice_params%b_x)
    ice_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+ice_params%fix_mu+ice_params%sp2)
    ice_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+ice_params%fix_mu+ice_params%sp2+ice_params%b_x)
    ! exponent of lambda for sedimentation
    ice_params%exp_1_mu_sp1 = 1.0+ice_params%fix_mu+ice_params%sp1
    ice_params%exp_1_mu_sp1_bx = 1.0 + ice_params%fix_mu + ice_params%sp1 + ice_params%b_x
    ice_params%exp_1_mu_sp2 = 1.0+ice_params%fix_mu+ice_params%sp2
    ice_params%exp_1_mu_sp2_bx = 1.0 + ice_params%fix_mu + ice_params%sp2 + ice_params%b_x
    ! gamma function used in ventilation
    ice_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*ice_params%b_x+ice_params%fix_mu+2.5)

    
    ! Snow
    ! gamma function constants (used in sedimentation)
    snow_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+snow_params%fix_mu+snow_params%sp1)
    snow_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+snow_params%fix_mu+snow_params%sp1+snow_params%b_x)
    snow_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+snow_params%fix_mu+snow_params%sp2)
    snow_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+snow_params%fix_mu+snow_params%sp2+snow_params%b_x)
    ! exponent of lambda for sedimentation
    snow_params%exp_1_mu_sp1 = 1.0+snow_params%fix_mu+snow_params%sp1
    snow_params%exp_1_mu_sp1_bx = 1.0 + snow_params%fix_mu + snow_params%sp1 + snow_params%b_x
    snow_params%exp_1_mu_sp2 = 1.0+snow_params%fix_mu+snow_params%sp2
    snow_params%exp_1_mu_sp2_bx = 1.0 + snow_params%fix_mu + snow_params%sp2 + snow_params%b_x
    ! gamma function used in ventilation
    snow_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*snow_params%b_x+snow_params%fix_mu+2.5)
    
    
    ! Graupel
    ! gamma function constants (used in sedimentation)
    graupel_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+graupel_params%fix_mu+graupel_params%sp1)
    graupel_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+graupel_params%fix_mu+graupel_params%sp1+graupel_params%b_x)
    graupel_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+graupel_params%fix_mu+graupel_params%sp2)
    graupel_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+graupel_params%fix_mu+graupel_params%sp2+graupel_params%b_x)
    ! exponent of lambda for sedimentation
    graupel_params%exp_1_mu_sp1 = 1.0+graupel_params%fix_mu+graupel_params%sp1
    graupel_params%exp_1_mu_sp1_bx = 1.0 + graupel_params%fix_mu + graupel_params%sp1 + graupel_params%b_x
    graupel_params%exp_1_mu_sp2 = 1.0+graupel_params%fix_mu+graupel_params%sp2
    graupel_params%exp_1_mu_sp2_bx = 1.0 + graupel_params%fix_mu + graupel_params%sp2 + graupel_params%b_x
    ! gamma function used in ventilation
    graupel_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*graupel_params%b_x+graupel_params%fix_mu+2.5)

! gamma functions for the calculation of the ice_accretion using the functions sweepout and binary collection

    ! Cloud 
    ! gamma function (used in sweepout)
    cloud_params%gam_3_bx_mu = GammaFunc(3.0+cloud_params%b_x+cloud_params%fix_mu)
    cloud_params%gam_3_bx_mu_dx = GammaFunc(3.0+cloud_params%b_x+cloud_params%fix_mu+cloud_params%d_x)
    ! gamma function constants (used in binary_collection)
    cloud_params%gam_2_mu = GammaFunc(2.0_wp+cloud_params%fix_mu)
    cloud_params%gam_3_mu = GammaFunc(3.0_wp+cloud_params%fix_mu)
    cloud_params%gam_2_mu_dx = GammaFunc(2.0_wp+cloud_params%fix_mu+cloud_params%d_x)
    cloud_params%gam_3_mu_dx = GammaFunc(3.0_wp+cloud_params%fix_mu+cloud_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    cloud_params%gam_1_mu_dx_bx = GammaFunc(1.0 + cloud_params%fix_mu + cloud_params%d_x + cloud_params%b_x)
    cloud_params%gam_1_mu_dx = GammaFunc(1.0 + cloud_params%fix_mu + cloud_params%d_x) 
    
    ! Rain
    ! gamma function (used in sweepout)
    rain_params%gam_3_bx_mu = GammaFunc(3.0+rain_params%b_x+rain_params%fix_mu)
    rain_params%gam_3_bx_mu_dx = GammaFunc(3.0+rain_params%b_x+rain_params%fix_mu+rain_params%d_x)
    ! gamma function constants (used in binary_collection)
    rain_params%gam_2_mu = GammaFunc(2.0_wp+rain_params%fix_mu)
    rain_params%gam_3_mu = GammaFunc(3.0_wp+rain_params%fix_mu)
    rain_params%gam_2_mu_dx = GammaFunc(2.0_wp+rain_params%fix_mu+rain_params%d_x)
    rain_params%gam_3_mu_dx = GammaFunc(3.0_wp+rain_params%fix_mu+rain_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    rain_params%gam_1_mu_dx_bx = GammaFunc(1.0 + rain_params%fix_mu + rain_params%d_x + rain_params%b_x)
    rain_params%gam_1_mu_dx = GammaFunc(1.0 + rain_params%fix_mu + rain_params%d_x) 
    
    ! Ice
    ! gamma function (used in sweepout)
    ice_params%gam_3_bx_mu = GammaFunc(3.0+ice_params%b_x+ice_params%fix_mu)
    ice_params%gam_3_bx_mu_dx = GammaFunc(3.0+ice_params%b_x+ice_params%fix_mu+ice_params%d_x)
    ! gamma function constants (used in binary_collection)
    ice_params%gam_2_mu = GammaFunc(2.0_wp+ice_params%fix_mu)
    ice_params%gam_3_mu = GammaFunc(3.0_wp+ice_params%fix_mu)
    ice_params%gam_2_mu_dx = GammaFunc(2.0_wp+ice_params%fix_mu+ice_params%d_x)
    ice_params%gam_3_mu_dx = GammaFunc(3.0_wp+ice_params%fix_mu+ice_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    ice_params%gam_1_mu_dx_bx = GammaFunc(1.0 + ice_params%fix_mu + ice_params%d_x + ice_params%b_x)
    ice_params%gam_1_mu_dx = GammaFunc(1.0 + ice_params%fix_mu + ice_params%d_x) 
    
    ! Snow
    ! gamma function (used in sweepout)
    snow_params%gam_3_bx_mu = GammaFunc(3.0+snow_params%b_x+snow_params%fix_mu)
    snow_params%gam_3_bx_mu_dx = GammaFunc(3.0+snow_params%b_x+snow_params%fix_mu+snow_params%d_x)
    ! gamma function constants (used in binary_collection)
    snow_params%gam_2_mu = GammaFunc(2.0_wp+snow_params%fix_mu)
    snow_params%gam_3_mu = GammaFunc(3.0_wp+snow_params%fix_mu)
    snow_params%gam_2_mu_dx = GammaFunc(2.0_wp+snow_params%fix_mu+snow_params%d_x)
    snow_params%gam_3_mu_dx = GammaFunc(3.0_wp+snow_params%fix_mu+snow_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    snow_params%gam_1_mu_dx_bx = GammaFunc(1.0 + snow_params%fix_mu + snow_params%d_x + snow_params%b_x)
    snow_params%gam_1_mu_dx = GammaFunc(1.0 + snow_params%fix_mu + snow_params%d_x) 

    ! Graupel
    ! gamma function (used in sweepout)
    graupel_params%gam_3_bx_mu = GammaFunc(3.0+graupel_params%b_x+graupel_params%fix_mu)
    graupel_params%gam_3_bx_mu_dx = GammaFunc(3.0+graupel_params%b_x+graupel_params%fix_mu+graupel_params%d_x)
    ! gamma function constants (used in binary_collection)
    graupel_params%gam_2_mu = GammaFunc(2.0_wp+graupel_params%fix_mu)
    graupel_params%gam_3_mu = GammaFunc(3.0_wp+graupel_params%fix_mu)
    graupel_params%gam_2_mu_dx = GammaFunc(2.0_wp+graupel_params%fix_mu+graupel_params%d_x)
    graupel_params%gam_3_mu_dx = GammaFunc(3.0_wp+graupel_params%fix_mu+graupel_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    graupel_params%gam_1_mu_dx_bx = GammaFunc(1.0 + graupel_params%fix_mu + graupel_params%d_x + graupel_params%b_x)
    graupel_params%gam_1_mu_dx = GammaFunc(1.0 + graupel_params%fix_mu + graupel_params%d_x) 
    

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    
  end subroutine gamma_initialize

end module initialize
