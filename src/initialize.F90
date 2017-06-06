module initialize
  use mphys_die, only: throw_mphys_error, incorrect_opt, std_msg
  use variable_precision, only: wp
  use lookup, only: set_mu_lookup, mu_i, mu_g, mu_i_sed, mu_g_sed, nmu
  use derived_constants, only: set_constants
  use mphys_parameters, only: p1, p2, p3, sp1, sp2, sp3, nprocs, naeroprocs, snow_params
  use mphys_switches, only: option, aerosol_option, l_warm, hydro_complexity, aero_complexity, cloud_params, &
       rain_params, ice_params, snow_params, graupel_params, iopt_act, l_g, l_sg, l_override_checks, i_am10, i_an10, &
       isol, iinsol, active_rain, active_cloud, aero_index, active_number, process_level, iopt_act, l_process
  use gauss_casim_micro, only: gaussfunclookup
  use micro_main, only : DTPUD, initialise_micromain, finalise_micromain
  use sedimentation, only : initialise_sedr, finalise_sedr

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='INITIALIZE'

  ! indices for process conversion rates
  integer :: n_cond_q, n_cond_n, n_aut_qr, n_aut_nr, n_aut_m3, n_aut_nl, n_acw_qr, n_acw_m3, n_acw_nl, n_evp_qr, &
       n_evp_nr, n_evp_m3, n_acr_nr, n_acr_m3, n_dep_qi, n_dep_ni, n_sub_qi, n_sub_ni, n_iacw_qi, n_iacw_ni, n_saut_qi, &
       n_saut_ni, n_raci_qi, n_raci_ni, n_gshd_qi, n_gshd_ni, n_hal_qi, n_hal_ni, n_iagg_ni, n_inuc_qi, n_inuc_ni, &
       n_homc_ni, n_homc_qi, n_gdep_qg, n_gsub_qg, n_gsub_ng, n_sdep_qs, n_ssub_qs, n_ssub_ns, n_sacw_ql, n_sacw_qs, &
       n_sacw_qg, n_sacw_ns, n_sacw_nl, n_sacw_ng, n_saci_qs, n_saci_qi, n_saci_ns, n_saci_ni, n_sacr_qr, n_sacr_qs, &
       n_sacr_qg, n_sacr_ns, n_sacr_nr, n_sacr_ng, n_gacr_qg, n_gacr_qr, n_gacr_nr, n_gacr_ng, n_gacw_qg, n_gacw_ql, &
       n_gacw_nl, n_gacw_ng, n_gaci_qg, n_gaci_qi, n_gaci_ni, n_gaci_ng, n_gacs_qs, n_gacs_qg, n_gacs_ns, n_gacs_ng, &
       n_sagg_ns, n_gagg_ng, n_sbrk_ns, n_gshd_qg, n_gshd_ng, n_gshd_qr, n_gshd_nr, n_gshd_qs, n_gshd_ns, n_hal_qs,  &
       n_hal_ns, n_hal_qg, n_hal_ng, n_smlt_ns, n_smlt_qs, n_gmlt_ng, n_gmlt_qg, n_homr_nr, n_homr_qr

  public mphys_init, mphys_finalise
contains 

  subroutine mphys_init(il, iu, jl, ju, kl, ku,                 &       
       is_in, ie_in, js_in, je_in, ks_in, ke_in, l_tendency)

    implicit none

    character(len=*), parameter :: RoutineName='MPHYS_INIT'

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

    integer :: iproc
    real(wp) :: tmp

    integer :: is, ie ! upper and lower i levels which are to be used
    integer :: js, je ! upper and lower j levels
    integer :: ks, ke ! upper and lower k levels    
    logical :: l_tendency_loc

    call check_options()
    call set_constants()


    ! Set grid extents to operate on
    if (present(is_in)) is=is_in
    if (present(ie_in)) ie=ie_in
    if (present(js_in)) js=js_in
    if (present(je_in)) je=je_in
    if (present(ks_in)) ks=ks_in
    if (present(ke_in)) ke=ke_in

    ! if not passed in, then default to full grid
    if (.not. present(is_in)) is=il
    if (.not. present(ie_in)) ie=iu
    if (.not. present(js_in)) js=jl
    if (.not. present(je_in)) je=ju
    if (.not. present(ks_in)) ks=kl
    if (.not. present(ke_in)) ke=ku
    
    if (present(l_tendency)) then
      l_tendency_loc=l_tendency
    else
      l_tendency_loc=.true.
    end if

    call initialise_micromain(il, iu, jl, ju, kl, ku, is, ie, js, je, ks, ke, l_tendency_loc)

    call initialise_sedr()
    
    call initialise_lookup_tables()
    !Gaussfunc
    call gaussfunclookup(snow_params%id, tmp, a=snow_params%fix_mu, b=snow_params%d_x, init=.true.)
  end subroutine mphys_init

  subroutine mphys_finalise()

    implicit none
    character(len=*), parameter :: RoutineName='MPHYS_FINALISE'

    call finalise_sedr()
    call finalise_micromain()
  end subroutine mphys_finalise

  subroutine initialise_lookup_tables()

    implicit none
    character(len=*), parameter :: RoutineName='INITIALISE_LOOKUP_TABLES'

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
  end subroutine initialise_lookup_tables

  ! Check that the options that have been selected are consitent
  subroutine check_options()

    implicit none
    character(len=*), parameter :: RoutineName='CHECK_OPTIONS'

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
      if (aerosol_option == 3) then
        write(std_msg, '(A)') 'This aerosol_option is not yet well tested.'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if

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
      if (aero_index%i_accum_dust == 0 .and. aero_index%nin /= 1 .and. &
          active_number(iinsol))                                  then
        write(std_msg, '(A)') 'Dust modes must only be accumulation mode '//&
                              'with insoluble active_number'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if

      ! Aerosol processing consistency
      if (process_level > 0 .and. iopt_act < 3) then
        write(std_msg, '(A)') 'If processing aerosol, must use higher '//&
                              'level activation code, i.e check iopt_act'
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      end if
    end if
  end subroutine check_options
end module initialize
