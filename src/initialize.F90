MODULE initialize

USE mphys_die, ONLY: throw_mphys_error
USE variable_precision, ONLY: wp
USE lookup, ONLY: set_mu_lookup, mu_i, mu_g, mu_i_sed, mu_g_sed, nmu
USE derived_constants, ONLY: set_constants
USE mphys_parameters, ONLY: p1, p2, p3, sp1, sp2, sp3   &
     , nprocs, naeroprocs, snow_params
USE mphys_switches, ONLY: option, aerosol_option, l_warm   &
     , hydro_complexity, aero_complexity,   &
     cloud_params, rain_params, ice_params, snow_params, graupel_params &
     , iopt_act, l_g, l_sg, l_override_checks &
     , i_am10, i_an10, isol, iinsol, active_rain, active_cloud, aero_index &
     , active_number, process_level, iopt_act, l_process

USE gauss_casim_micro, ONLY: gaussfunclookup

#if DEF_MODEL==MODEL_UM
USE mphys_casim_diagnostics, ONLY: ProcessRates, nProcessDiags, ProcessKeys, ProcessQs, &
       PhaseChanges
USE atm_fields_bounds_mod, ONLY: tdims
USE process_routines, ONLY: i_idep, i_sdep
#elif DEF_MODEL==MODEL_LEM
USE process_routines, ONLY:  process_rate   &
     , allocate_procs, deallocate_procs &
     , i_cond, i_praut, i_pracw, i_pracr, i_prevp, i_psedr, i_tidy, i_tidy2 &
     , i_aact, i_aaut, i_aacw, i_aevp, i_asedr, i_arevp, i_atidy, i_atidy2   &
     , i_inuc, i_homc
USE passive_fields, ONLY: dt
USE com_prametr, ONLY: ilatentdgp, nprc_prametr => nprc, iforscalp
USE com_dgstore, ONLY: procrate, fallq
USE com_dgchars, ONLY: procchar
USE com_rain, ONLY: puddle
#endif

IMPLICIT NONE

INTEGER ::   & ! indices for process conversion rates
          n_cond_q,  &
          n_cond_n,  &
          n_aut_qr,  &
          n_aut_nr,  &
          n_aut_m3,  &
          n_aut_nl,  &
          n_acw_qr,  &
          n_acw_m3,  &
          n_acw_nl,  &
          n_evp_qr,  &
          n_evp_nr,  &
          n_evp_m3,  &
          n_acr_nr,  &
          n_acr_m3,  &
          n_dep_qi,  &
          n_dep_ni,  &
          n_sub_qi,  &
          n_sub_ni,  &
          n_iacw_qi, &
          n_iacw_ni, &
          n_saut_qi, &
          n_saut_ni, &
          n_raci_qi, &
          n_raci_ni, &
          n_gshd_qi, &
          n_gshd_ni, &
          n_hal_qi,  &
          n_hal_ni,  &
          n_iagg_ni, &
          n_inuc_qi, &
          n_inuc_ni, &
          n_homc_ni, &
          n_homc_qi, &

          n_gdep_qg, &
          n_gsub_qg, &
          n_gsub_ng, &
          n_sdep_qs, &
          n_ssub_qs, &
          n_ssub_ns, &
          n_sacw_ql, &
          n_sacw_qs, &
          n_sacw_qg, &
          n_sacw_ns, &
          n_sacw_nl, &
          n_sacw_ng, &
          n_saci_qs, &
          n_saci_qi, &
          n_saci_ns, &
          n_saci_ni, &
          n_sacr_qr, &
          n_sacr_qs, &
          n_sacr_qg, &
          n_sacr_ns, &
          n_sacr_nr, &
          n_sacr_ng, &
          n_gacr_qg, &
          n_gacr_qr, &
          n_gacr_nr, &
          n_gacr_ng, &
          n_gacw_qg, &
          n_gacw_ql, &
          n_gacw_nl, &
          n_gacw_ng, &
          n_gaci_qg, &
          n_gaci_qi, &
          n_gaci_ni, &
          n_gaci_ng, &
          n_gacs_qs, &
          n_gacs_qg, &
          n_gacs_ns, &
          n_gacs_ng, &
          n_sagg_ns, &
          n_gagg_ng, &
          n_sbrk_ns, &
          n_gshd_qg, &
          n_gshd_ng, &
          n_gshd_qr, &
          n_gshd_nr, &
          n_gshd_qs, &
          n_gshd_ns, &
          n_hal_qs,  &
          n_hal_ns,  &
          n_hal_qg,  &
          n_hal_ng,  &
          n_smlt_ns,  &
          n_smlt_qs,  &
          n_gmlt_ng,  &
          n_gmlt_qg,  &
          n_homr_nr, &
          n_homr_qr

REAL :: DTPUD ! Time step for puddle diagnostic

CONTAINS

SUBROUTINE mphys_init

INTEGER :: iproc

REAL(wp) :: tmp

CALL check_options
CALL set_constants

    !------------------------------------------------------
    ! Look up tables... (These need to be extended)
    !------------------------------------------------------
    ! mu lookup
ALLOCATE(mu_i(nmu))
ALLOCATE(mu_g(nmu))
CALL set_mu_lookup(p1, p2, p3, mu_i, mu_g)
ALLOCATE(mu_i_sed(nmu))
ALLOCATE(mu_g_sed(nmu))
CALL set_mu_lookup(sp1, sp2, sp3, mu_i_sed, mu_g_sed)

    !Gaussfunc
CALL gaussfunclookup(snow_params%id, tmp, a=snow_params%fix_mu, b=snow_params%d_x, init=.TRUE.)


#if DEF_MODEL==MODEL_UM
!!$ALLOCATE(ProcessRates(tdims%i_start : tdims%i_end,                           &
!!$                       	  tdims%j_start : tdims%j_end,                       &
!!$                          1 : tdims%k_end,                                   &
!!$                          nProcessDiags  ))
!!$
!!$ALLOCATE(PhaseChanges(tdims%i_start : tdims%i_end,                           &
!!$                       	  tdims%j_start : tdims%j_end,                       &
!!$                          1 : tdims%k_end,                                   &
!!$                          6 ))

!    ProcessKeys(:)=(/i_idep%unique_id, i_sdep%unique_id/)
!    ProcessQs(:)  =(/ice_params%i_1m,  snow_params%i_1m/)

#elif DEF_MODEL==MODEL_LEM

    ! LEM Standard Diagnostics...
DO iproc=1,nprc_prametr
  SELECT CASE (procchar(iproc))
    CASE ('cond_q')
      n_cond_q = iproc
    CASE ('cond_n')
      n_cond_n = iproc
    CASE ('aut_qr')
      n_aut_qr = iproc
    CASE ('aut_nr')
      n_aut_nr = iproc
    CASE ('aut_m3')
      n_aut_m3 = iproc
    CASE ('aut_nl')
      n_aut_nl = iproc
    CASE ('acw_qr')
      n_acw_qr = iproc
    CASE ('acw_m3')
      n_acw_m3 = iproc
    CASE ('acw_nl')
      n_aut_nl = iproc
    CASE ('evp_qr')
      n_evp_qr = iproc
    CASE ('evp_nr')
      n_evp_nr = iproc
    CASE ('evp_m3')
      n_evp_m3 = iproc
    CASE ('acr_nr')
      n_acr_nr = iproc
    CASE ('acr_m3')
      n_acr_m3 = iproc
    CASE ('dep_qi')
      n_dep_qi = iproc
    CASE ('dep_ni')
      n_dep_ni = iproc
    CASE ('sub_qi')
      n_sub_qi = iproc
    CASE ('sub_ni')
      n_sub_ni = iproc
    CASE ('iacw_qi')
      n_iacw_qi = iproc
    CASE ('iacw_ni')
      n_iacw_ni = iproc
    CASE ('saut_qi')
      n_saut_qi = iproc
    CASE ('saut_ni')
      n_saut_ni = iproc
    CASE ('raci_qi')
      n_raci_qi = iproc
    CASE ('raci_ni')
      n_raci_ni = iproc
    CASE ('gshd_qi')
      n_gshd_qi = iproc
    CASE ('gshd_ni')
      n_gshd_ni = iproc
    CASE ('hal_qi')
      n_hal_qi = iproc
    CASE ('hal_ni')
      n_hal_ni = iproc
    CASE ('iagg_ni')
      n_iagg_ni = iproc
    CASE ('inuc_qi')
      n_inuc_qi = iproc
    CASE ('inuc_ni')
      n_inuc_ni = iproc
    CASE ('homc_ni')
      n_homc_ni = iproc
    CASE ('homc_qi')
      n_homc_qi = iproc

    CASE ('gdep_qg')
      n_gdep_qg = iproc
    CASE ('gsub_qg')
      n_gsub_qg = iproc
    CASE ('gsub_ng')
      n_gsub_ng = iproc
    CASE ('sdep_qs')
      n_sdep_qs = iproc
    CASE ('ssub_qs')
      n_ssub_qs = iproc
    CASE ('ssub_ns')
      n_ssub_ns = iproc
    CASE ('sacw_ql')
      n_sacw_ql = iproc
    CASE ('sacw_qs')
      n_sacw_qs = iproc
    CASE ('sacw_qg')
      n_sacw_qg = iproc
    CASE ('sacw_ns')
      n_sacw_ns = iproc
    CASE ('sacw_nl')
      n_sacw_nl = iproc
    CASE ('sacw_ng')
      n_sacw_ng = iproc
    CASE ('saci_qs')
      n_saci_qs = iproc
    CASE ('saci_qi')
      n_saci_qi = iproc
    CASE ('saci_ns')
      n_saci_ns = iproc
    CASE ('saci_ni')
      n_saci_ni = iproc
    CASE ('sacr_qr')
      n_sacr_qr = iproc
    CASE ('sacr_qs')
      n_sacr_qs = iproc
    CASE ('sacr_qg')
      n_sacr_qg = iproc
    CASE ('sacr_ns')
      n_sacr_ns = iproc
    CASE ('sacr_nr')
      n_sacr_nr = iproc
    CASE ('sacr_ng')
      n_sacr_ng = iproc
    CASE ('gacr_qg')
      n_gacr_qg = iproc
    CASE ('gacr_qr')
      n_gacr_qr = iproc
    CASE ('gacr_nr')
      n_gacr_nr = iproc
    CASE ('gacr_ng')
      n_gacr_ng = iproc
    CASE ('gacw_qg')
      n_gacw_qg = iproc
    CASE ('gacw_ql')
      n_gacw_ql = iproc
    CASE ('gacw_nl')
      n_gacw_nl = iproc
    CASE ('gacw_ng')
      n_gacw_ng = iproc
    CASE ('gaci_qg')
      n_gaci_qg = iproc
    CASE ('gaci_qi')
      n_gaci_qi = iproc
    CASE ('gaci_ni')
      n_gaci_ni = iproc
    CASE ('gaci_ng')
      n_gaci_ng = iproc
    CASE ('gacs_qs')
      n_gacs_qs = iproc
    CASE ('gacs_qg')
      n_gacs_qg = iproc
    CASE ('gacs_ns')
      n_gacs_ns = iproc
    CASE ('gacs_ng')
      n_gacs_ng = iproc
    CASE ('sagg_ns')
      n_sagg_ns = iproc
    CASE ('gagg_ng')
      n_gagg_ng = iproc
    CASE ('sbrk_ns')
      n_sbrk_ns = iproc
    CASE ('gshd_qg')
      n_gshd_qg = iproc
    CASE ('gshd_ng')
      n_gshd_ng = iproc
    CASE ('gshd_qr')
      n_gshd_qr = iproc
    CASE ('gshd_nr')
      n_gshd_nr = iproc
    CASE ('gshd_qs')
      n_gshd_qs = iproc
    CASE ('gshd_ns')
      n_gshd_ns = iproc
    CASE ('hal_qs')
      n_hal_qs = iproc
    CASE ('hal_ns')
      n_hal_ns = iproc
    CASE ('hal_qg')
      n_hal_qg = iproc
    CASE ('hal_ng')
      n_hal_ng = iproc
    CASE ('smlt_ns')
      n_smlt_ns = iproc
    CASE ('smlt_qs')
      n_smlt_qs = iproc
    CASE ('gmlt_ng')
      n_gmlt_ng = iproc
    CASE ('gmlt_qg')
      n_gmlt_qg = iproc
    CASE ('homr_nr')
      n_homr_nr = iproc
    CASE ('homr_qr')
      n_homr_qr = iproc
  END SELECT
END DO

IF (IFORSCALP == 0) THEN
  DTPUD = 0.5*DT
ELSE
  DTPUD = DT
END IF

#endif

END SUBROUTINE mphys_init


SUBROUTINE check_options
    ! Check that the options that have been selected are consitent


    ! if (.not. l_override_checks)then
    !   call throw_mphys_error(1, 'check_options', 'You havent yet sorted out the passive numbers: first thing is to look at thte tidy routines!')
    ! end if
IF (.NOT. l_override_checks) THEN

IF (l_warm .AND. l_process) THEN
  CALL throw_mphys_error(1, 'check_options', 'processing does not currently work with l_warm=.true.')
END IF

IF (.NOT. ice_params%l_1m .AND. .NOT. l_warm) THEN
  CALL throw_mphys_error(1, 'check_options', 'l_warm must be true if not using ice')
END IF

IF (cloud_params%l_2m .AND. aerosol_option == 0     &
       .AND. (iopt_act== 1 .OR. iopt_act==3)) THEN
  CALL throw_mphys_error(1, 'check_options', 'for double moment cloud you must have aerosol_option>0'// &
     'or else activation should be independent of aerosol')
END IF

IF (.NOT. cloud_params%l_2m .AND. aerosol_option > 0) THEN
  CALL throw_mphys_error(1, 'check_options', 'aerosol_option must be 0 if not using double moment microphysics')
END IF

IF (l_g .AND. .NOT. l_sg) THEN
!      call throw_mphys_error(1, 'check_options', 'cannot run with graupel but not with snow')
  l_g=.FALSE.
END IF


    ! Some options are not yet working or well tested so
    ! don't let anyone use these
  IF (aerosol_option==3) THEN
    CALL throw_mphys_error(1, 'check_options', 'This aerosol_option is not yet well tested.')
  END IF


  IF (i_am10 > 0 .OR. i_an10 > 0) THEN
    CALL throw_mphys_error(1, 'check_options', 'Accumulation mode dust is not yet used.')
  END IF


    ! Aerosol consistency
IF (active_rain(isol) .AND. .NOT. active_cloud(isol)) THEN
  CALL throw_mphys_error(1, 'check_options', 'active_rain(isol) .and. .not. active_cloud(isol)')
END IF


    ! Aerosol consistency
IF (aero_index%i_accum == 0 .AND. aero_index%i_coarse==0     &
       .AND. (iopt_act==1 .OR. iopt_act==3) ) THEN
  CALL throw_mphys_error(1, 'check_options', 'Must have accumulation or coarse mode aerosol for chosen activation option')
END IF

    ! Aerosol consistency
IF (aero_index%i_accum == 0 .AND. aero_index%nccn /= 1       &
       .AND. active_number(isol)) THEN
  CALL throw_mphys_error(1, 'check_options', 'Soluble modes must only be accumulation mode with soluble active_number')
END IF

    ! Aerosol consistency
IF (aero_index%i_accum_dust == 0 .AND. aero_index%nin /= 1       &
       .AND. active_number(iinsol)) THEN
  CALL throw_mphys_error(1, 'check_options', 'Dust modes must only be accumulation mode with insoluble active_number')
END IF

    ! Aerosol processing consistency
IF (process_level > 0 .AND. iopt_act < 3) THEN
  CALL throw_mphys_error(1, 'check_options', 'If processing aerosol, must use higher level activatio code, i.e check iopt_act')
END IF
END IF

END SUBROUTINE check_options


END MODULE initialize



