module initialize
  use mphys_die, only: throw_mphys_error
  use variable_precision, only: wp
  use lookup, only: set_mu_lookup, mu_i, mu_g, mu_i_sed, mu_g_sed, nmu
  use derived_constants, only: set_constants
  use mphys_parameters, only: p1, p2, p3, sp1, sp2, sp3, nprocs, naeroprocs, snow_params
  use mphys_switches, only: option, aerosol_option, l_warm, hydro_complexity, aero_complexity, cloud_params, &
       rain_params, ice_params, snow_params, graupel_params, iopt_act, l_g, l_sg, l_override_checks, i_am10, i_an10, &
       isol, iinsol, active_rain, active_cloud, aero_index, active_number, process_level, iopt_act, l_process
  use gauss_casim_micro, only: gaussfunclookup

  implicit none
  private

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

  real :: DTPUD ! Time step for puddle diagnostic

  public mphys_init
contains

  subroutine mphys_init()
    integer :: iproc
    real(wp) :: tmp

    call check_options()
    call set_constants()

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

    !Gaussfunc
    call gaussfunclookup(snow_params%id, tmp, a=snow_params%fix_mu, b=snow_params%d_x, init=.true.)

  end subroutine mphys_init

  ! Check that the options that have been selected are consitent
  subroutine check_options()
    if (.not. l_override_checks) then
      if (l_warm .and. l_process) then
        call throw_mphys_error(1, 'check_options', 'processing does not currently work with l_warm=.true.')
      end if

      if (.not. ice_params%l_1m .and. .not. l_warm) then
        call throw_mphys_error(1, 'check_options', 'l_warm must be true if not using ice')
      end if

      if (cloud_params%l_2m .and. aerosol_option == 0     &
           .and. (iopt_act== 1 .or. iopt_act==3)) then
        call throw_mphys_error(1, 'check_options', 'for double moment cloud you must have aerosol_option>0'// &
             'or else activation should be independent of aerosol')
      end if

      if (.not. cloud_params%l_2m .and. aerosol_option > 0) then
        call throw_mphys_error(1, 'check_options', 'aerosol_option must be 0 if not using double moment microphysics')
      end if

      if (l_g .and. .not. l_sg) then
        !      call throw_mphys_error(1, 'check_options', 'cannot run with graupel but not with snow')
        l_g=.false.
      end if

      ! Some options are not yet working or well tested so don't let anyone use these
      if (aerosol_option == 3) then
        call throw_mphys_error(1, 'check_options', 'This aerosol_option is not yet well tested.')
      end if

      if (i_am10 > 0 .or. i_an10 > 0) then
        call throw_mphys_error(1, 'check_options', 'Accumulation mode dust is not yet used.')
      end if

      ! Aerosol consistency
      if (active_rain(isol) .and. .not. active_cloud(isol)) then
        call throw_mphys_error(1, 'check_options', 'active_rain(isol) .and. .not. active_cloud(isol)')
      end if

      ! Aerosol consistency
      if (aero_index%i_accum == 0 .and. aero_index%i_coarse==0 .and. (iopt_act==1 .or. iopt_act==3) ) then
        call throw_mphys_error(1, 'check_options', 'Must have accumulation or coarse mode aerosol for chosen activation option')
      end if

      ! Aerosol consistency
      if (aero_index%i_accum == 0 .and. aero_index%nccn /= 1 .and. active_number(isol)) then
        call throw_mphys_error(1, 'check_options', 'Soluble modes must only be accumulation mode with soluble active_number')
      end if

      ! Aerosol consistency
      if (aero_index%i_accum_dust == 0 .and. aero_index%nin /= 1 .and. active_number(iinsol)) then
        call throw_mphys_error(1, 'check_options', 'Dust modes must only be accumulation mode with insoluble active_number')
      end if

      ! Aerosol processing consistency
      if (process_level > 0 .and. iopt_act < 3) then
        call throw_mphys_error(1, 'check_options', &
             'If processing aerosol, must use higher level activatio code, i.e check iopt_act')
      end if
    end if
  end subroutine check_options
end module initialize
