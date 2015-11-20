MODULE diags_lem
#if DEF_MODEL==MODEL_LEM_DIAG
USE type_process, ONLY: process_name, process_rate
USE mphys_switches, ONLY: aero_complexity, hydro_complexity  &
     , ntotalq, hydro_names, aero_names, ntotala
USE mphys_parameters, ONLY: nprocs,naeroprocs,   &
     ice_params, cloud_params, snow_params, &
     graupel_params, rain_params, nz
USE process_routines, ONLY:    &
       i_cond, i_praut, i_pracw, i_pracr, i_prevp, i_psedr &
     , i_tidy, i_tidy2, i_inuc                             &
     , i_aact, i_asedr, i_arevp, i_aevp, i_atidy, i_atidy2 &
     , i_idep, i_isub, i_smlt, i_gmlt, i_gacs &
     , i_saut, i_gaci, i_saci, i_raci, i_homc, i_homr, i_ihal &
     , i_iagg, i_gshd, i_iacw, i_gdep, i_gsub, i_sdep, i_ssub &
     , i_sacw, i_sacr, i_gacw, i_gacr, i_sagg, i_gagg, i_sbrk
USE com_dgstore, ONLY: l_dodgs

USE initialize, ONLY:   &
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

USE com_dgstore, ONLY: procrate
USE com_prametr, ONLY : jjp, kkp, iipep
USE com_me, ONLY: mype


IMPLICIT NONE

CONTAINS

SUBROUTINE lem_procdgs(n, kl, ku, procs)

USE diaghelp_lem
USE extra_dgs

IMPLICIT NONE

INTEGER, INTENT(IN) :: n ! substep index
INTEGER, INTENT(IN) :: kl, ku ! k bounds
TYPE(process_rate) :: procs(:,:)

TYPE(process_name) :: iproc
INTEGER :: iq
INTEGER :: k

IF (l_dodgs) THEN

  IF (j_here>=jstartdg .AND. j_here<=jenddg) THEN

    koff=kl - 1
    idgproc = req_dgproc('QDEP')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_idep%id)%source(ice_params%i_1m) &
           + procs(k-koff, i_isub%id)%source(ice_params%i_1m)
      END DO
    END IF


    idgproc = req_dgproc('QACC')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_pracw%id)%source(rain_params%i_1m)
      END DO
    END IF

    idgproc = req_dgproc('QAUT')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_praut%id)%source(rain_params%i_1m)
      END DO
    END IF

    idgproc = req_dgproc('QREVAP')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_prevp%id)%source(rain_params%i_1m)
      END DO
    END IF

    idgproc = req_dgproc('QMELT')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_smlt%id)%source(rain_params%i_1m)  &
           + procs(k-koff, i_gmlt%id)%source(rain_params%i_1m)
      END DO
    END IF


    idgproc = req_dgproc('NINUC')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_inuc%id)%source(ice_params%i_2m)
      END DO
    END IF


    idgproc = req_dgproc('NIHAL')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_ihal%id)%source(ice_params%i_2m)
      END DO
    END IF


    idgproc = req_dgproc('NIFRW')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_homc%id)%source(ice_params%i_2m)
      END DO
    END IF


    idgproc = req_dgproc('NIACC')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_gaci%id)%source(ice_params%i_2m)   &
           + procs(k-koff, i_saci%id)%source(ice_params%i_2m)  &
           + procs(k-koff, i_raci%id)%source(ice_params%i_2m)
      END DO
    END IF

    idgproc = req_dgproc('QSAUT')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_saut%id)%source(snow_params%i_1m)
      END DO
    END IF

    idgproc = req_dgproc('GDEP')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_gdep%id)%source(graupel_params%i_1m) &
           + procs(k-koff, i_gsub%id)%source(graupel_params%i_1m)
      END DO
    END IF

    idgproc = req_dgproc('GACC')
    IF (idgproc > 0) THEN
      DO k=kl,ku
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = &
           procs(k-koff, i_raci%id)%source(graupel_params%i_1m) &
           + procs(k-koff, i_sacr%id)%source(graupel_params%i_1m) &
           + procs(k-koff, i_gacw%id)%source(graupel_params%i_1m) &
           + procs(k-koff, i_gacr%id)%source(graupel_params%i_1m) &
           + procs(k-koff, i_gaci%id)%source(graupel_params%i_1m) &
           + procs(k-koff, i_gacs%id)%source(graupel_params%i_1m)
      END DO
    END IF

  END IF
      !-------------------------------
      ! For conditional diagnostics
      !-------------------------------
  DO k=1,nz
    iproc=i_cond
    iq=cloud_params%i_1m
    procrate(j_here,k,n_cond_q) = procs(k,iproc%id)%source(iq)

    iq=cloud_params%i_2m
    procrate(j_here,k,n_cond_n) = procs(k,iproc%id)%source(iq)

    iproc=i_praut
    iq=rain_params%i_1m
    procrate(j_here,k,n_aut_qr) = procs(k,iproc%id)%source(iq)

    iq=rain_params%i_2m
    procrate(j_here,k,n_aut_nr) = procs(k,iproc%id)%source(iq)

!      iq=rain_params%i_3m
!      procrate(j_here,k,n_aut_m3) = procs(k,iproc%id)%m3_source(iq)

    iq=cloud_params%i_1m
    procrate(j_here,k,n_aut_nl) = procs(k,iproc%id)%source(iq)

    iproc=i_pracw
    iq=rain_params%i_1m
    procrate(j_here,k,n_acw_qr) = procs(k,iproc%id)%source(iq)

!      iq=rain_params%i_3m
!      procrate(j_here,k,n_acw_m3) = procs(k,iproc%id)%m3_source(iq)

    iq=cloud_params%i_2m
    procrate(j_here,k,n_acw_nl) = procs(k,iproc%id)%source(iq)

    iproc=i_prevp
    iq=rain_params%i_1m
    procrate(j_here,k,n_evp_qr) = procs(k,iproc%id)%source(iq)

    iq=rain_params%i_2m
    procrate(j_here,k,n_evp_nr) = procs(k,iproc%id)%source(iq)

!      iq=rain_params%i_3m
!      procrate(j_here,k,n_evp_m3) = procs(k,iproc%id)%m3_source(iq)

    iproc=i_pracr
    iq=rain_params%i_2m
    procrate(j_here,k,n_acr_nr) = procs(k,iproc%id)%source(iq)

!      iq=rain_params%i_3m
!      procrate(j_here,k,n_acr_m3) = procs(k,iproc%id)%m3_source(iq)

    iproc=i_idep
    iq=ice_params%i_1m
    procrate(j_here,k,n_dep_qi) = procs(k,iproc%id)%source(iq)

    iproc=i_idep
    iq=ice_params%i_2m
    procrate(j_here,k,n_dep_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_isub
    iq=ice_params%i_1m
    procrate(j_here,k,n_sub_qi) = procs(k,iproc%id)%source(iq)

    iproc=i_isub
    iq=ice_params%i_2m
    procrate(j_here,k,n_sub_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_iacw
    iq=ice_params%i_1m
    procrate(j_here,k,n_iacw_qi) = procs(k,iproc%id)%source(iq)

    iq=ice_params%i_2m
    procrate(j_here,k,n_iacw_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_saut
    iq=ice_params%i_1m
    procrate(j_here,k,n_saut_qi) = procs(k,iproc%id)%source(iq)

    iq=ice_params%i_2m
    procrate(j_here,k,n_saut_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_raci
    iq=ice_params%i_1m
    procrate(j_here,k,n_raci_qi) = procs(k,iproc%id)%source(iq)

    iq=ice_params%i_2m
    procrate(j_here,k,n_raci_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_gshd
    iq=ice_params%i_1m
    procrate(j_here,k,n_gshd_qi) = procs(k,iproc%id)%source(iq)

    iq=ice_params%i_2m
    procrate(j_here,k,n_gshd_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_ihal
    iq=ice_params%i_1m
    procrate(j_here,k,n_hal_qi) = procs(k,iproc%id)%source(iq)

    iq=ice_params%i_2m
    procrate(j_here,k,n_hal_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_iagg
    iq=ice_params%i_2m
    procrate(j_here,k,n_iagg_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_inuc
    iq=ice_params%i_1m
    procrate(j_here,k,n_inuc_qi) = procs(k,iproc%id)%source(iq)

    iq=ice_params%i_2m
    procrate(j_here,k,n_inuc_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_homc
    iq=ice_params%i_1m
    procrate(j_here,k,n_homc_qi) = procs(k,iproc%id)%source(iq)

    iq=ice_params%i_2m
    procrate(j_here,k,n_homc_ni) = procs(k,iproc%id)%source(iq)

!----
    iproc=i_gdep
    iq=graupel_params%i_1m
    procrate(j_here,k,n_gdep_qg) = procs(k,iproc%id)%source(iq)

    iproc=i_gsub
    iq=graupel_params%i_1m
    procrate(j_here,k,n_gsub_qg) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_2m
    procrate(j_here,k,n_gsub_ng) = procs(k,iproc%id)%source(iq)

    iproc=i_sdep
    iq=snow_params%i_1m
    procrate(j_here,k,n_sdep_qs) = procs(k,iproc%id)%source(iq)

    iproc=i_ssub
    iq=snow_params%i_1m
    procrate(j_here,k,n_ssub_qs) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_2m
    procrate(j_here,k,n_ssub_ns) = procs(k,iproc%id)%source(iq)

    iproc=i_sacw
    iq=cloud_params%i_1m
    procrate(j_here,k,n_sacw_ql) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_1m
    procrate(j_here,k,n_sacw_qs) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_1m
    procrate(j_here,k,n_sacw_qg) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_2m
    procrate(j_here,k,n_sacw_ns) = procs(k,iproc%id)%source(iq)
    iq=cloud_params%i_2m
    procrate(j_here,k,n_sacw_nl) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_2m
    procrate(j_here,k,n_sacw_ng) = procs(k,iproc%id)%source(iq)

    iproc=i_saci
    iq=snow_params%i_1m
    procrate(j_here,k,n_saci_qs) = procs(k,iproc%id)%source(iq)
    iq=ice_params%i_1m
    procrate(j_here,k,n_saci_qi) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_2m
    procrate(j_here,k,n_saci_ns) = procs(k,iproc%id)%source(iq)
    iq=ice_params%i_2m
    procrate(j_here,k,n_saci_ni) = procs(k,iproc%id)%source(iq)

    iproc=i_sacr
    iq=rain_params%i_1m
    procrate(j_here,k,n_sacr_qr) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_1m
    procrate(j_here,k,n_sacr_qs) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_1m
    procrate(j_here,k,n_sacr_qg) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_2m
    procrate(j_here,k,n_sacr_ns) = procs(k,iproc%id)%source(iq)
    iq=rain_params%i_2m
    procrate(j_here,k,n_sacr_nr) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_2m
    procrate(j_here,k,n_sacr_ng) = procs(k,iproc%id)%source(iq)

    iproc=i_gacr
    iq=graupel_params%i_1m
    procrate(j_here,k,n_gacr_qg) = procs(k,iproc%id)%source(iq)
    iq=rain_params%i_1m
    procrate(j_here,k,n_gacr_qr) = procs(k,iproc%id)%source(iq)
    iq=rain_params%i_2m
    procrate(j_here,k,n_gacr_nr) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_2m
    procrate(j_here,k,n_gacr_ng) = procs(k,iproc%id)%source(iq)

    iproc=i_gacw
    iq=graupel_params%i_1m
    procrate(j_here,k,n_gacw_qg) = procs(k,iproc%id)%source(iq)
    iq=cloud_params%i_1m
    procrate(j_here,k,n_gacw_ql) = procs(k,iproc%id)%source(iq)
    iq=cloud_params%i_2m
    procrate(j_here,k,n_gacw_nl) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_2m
    procrate(j_here,k,n_gacw_ng) = procs(k,iproc%id)%source(iq)

    iproc=i_gaci
    iq=graupel_params%i_1m
    procrate(j_here,k,n_gaci_qg) = procs(k,iproc%id)%source(iq)
    iq=ice_params%i_1m
    procrate(j_here,k,n_gaci_qi) = procs(k,iproc%id)%source(iq)
    iq=ice_params%i_2m
    procrate(j_here,k,n_gaci_ni) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_2m
    procrate(j_here,k,n_gaci_ng) = procs(k,iproc%id)%source(iq)

    iproc=i_gacs
    iq=snow_params%i_1m
    procrate(j_here,k,n_gacs_qs) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_1m
    procrate(j_here,k,n_gacs_qg) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_2m
    procrate(j_here,k,n_gacs_ns) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_2m
    procrate(j_here,k,n_gacs_ng) = procs(k,iproc%id)%source(iq)

    iproc=i_sagg
    iq=snow_params%i_2m
    procrate(j_here,k,n_sagg_ns) = procs(k,iproc%id)%source(iq)

    iproc=i_gagg
    iq=graupel_params%i_2m
    procrate(j_here,k,n_gagg_ng) = procs(k,iproc%id)%source(iq)

    iproc=i_sbrk
    iq=snow_params%i_2m
    procrate(j_here,k,n_sbrk_ns) = procs(k,iproc%id)%source(iq)

    iproc=i_gshd
    iq=graupel_params%i_1m
    procrate(j_here,k,n_gshd_qg) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_2m
    procrate(j_here,k,n_gshd_ng) = procs(k,iproc%id)%source(iq)
    iq=rain_params%i_1m
    procrate(j_here,k,n_gshd_qr) = procs(k,iproc%id)%source(iq)
    iq=rain_params%i_2m
    procrate(j_here,k,n_gshd_nr) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_1m
    procrate(j_here,k,n_gshd_qs) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_2m
    procrate(j_here,k,n_gshd_ns) = procs(k,iproc%id)%source(iq)

    iproc=i_ihal
    iq=snow_params%i_1m
    procrate(j_here,k,n_hal_qs) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_2m
    procrate(j_here,k,n_hal_ns) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_1m
    procrate(j_here,k,n_hal_qg) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_2m
    procrate(j_here,k,n_hal_ng) = procs(k,iproc%id)%source(iq)

    iproc=i_smlt
    iq=snow_params%i_2m
    procrate(j_here,k,n_smlt_ns) = procs(k,iproc%id)%source(iq)
    iq=snow_params%i_1m
    procrate(j_here,k,n_smlt_qs) = procs(k,iproc%id)%source(iq)

    iproc=i_gmlt
    iq=graupel_params%i_2m
    procrate(j_here,k,n_gmlt_ng) = procs(k,iproc%id)%source(iq)
    iq=graupel_params%i_1m
    procrate(j_here,k,n_gmlt_qg) = procs(k,iproc%id)%source(iq)

    iproc=i_homr
    iq=rain_params%i_2m
    procrate(j_here,k,n_homr_nr) = procs(k,iproc%id)%source(iq)
    iq=rain_params%i_1m
    procrate(j_here,k,n_homr_qr) = procs(k,iproc%id)%source(iq)
  END DO

END IF

END SUBROUTINE lem_procdgs


#endif
END MODULE diags_lem

