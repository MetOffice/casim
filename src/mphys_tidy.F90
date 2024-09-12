module mphys_tidy
  use variable_precision, only: wp
  use process_routines, only: process_rate,  process_name
  use aerosol_routines, only: aerosol_active
  use thresholds, only: thresh_tidy, thresh_atidy
  use passive_fields, only: exner, pressure
  use mphys_switches, only:                       &
       i_qv, i_ql, i_nl, i_qr, i_nr, i_m3r, i_th, &
       i_qi, i_ni, i_qs, i_ns, i_m3s,             &
       i_qg, i_ng, i_m3g,                         &
!       l_3mr, l_3mg, l_3ms,                      &
       l_2mc, l_2mr,                              &
       l_2mi, l_2ms, l_2mg,                       &
       i_an2, i_am2, i_am4, i_am5, l_warm,        &
       i_an6, i_am6, i_am7, i_am8, i_am9,         &
       i_an11, i_an12,                            &
       l_process, ntotalq, ntotala,               &
       i_qstart, i_nstart, i_m3start,             &
       l_separate_rain, l_tidy_conserve_E, l_tidy_conserve_q, &
       l_passivenumbers, l_passivenumbers_ice
  use mphys_constants, only: Lv, Ls, cp
  use qsat_funs, only: qisaturation
  use mphys_parameters, only: hydro_params
  use mphys_die, only: throw_mphys_error, bad_values, warn, std_msg

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='MPHYS_TIDY'

  logical, parameter :: l_rescale_on_number = .false.
! logical :: l_tidym3 = .false.  ! Don't tidy based on m3 values

  real(wp), allocatable :: thresh(:), athresh(:), qin_thresh(:)
!$OMP THREADPRIVATE(thresh, athresh, qin_thresh)
  
  logical :: current_l_negonly, current_qin_l_negonly
!$OMP THREADPRIVATE(current_l_negonly, current_qin_l_negonly)

  public initialise_mphystidy, finalise_mphystidy, qtidy, ensure_positive, ensure_saturated, tidy_qin, &
       tidy_ain, ensure_positive_aerosol
contains

  subroutine initialise_mphystidy()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='INITIALISE_MPHYSTIDY'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    allocate(thresh(lbound(thresh_tidy,1):ubound(thresh_tidy,1)), qin_thresh(lbound(thresh_tidy,1):ubound(thresh_tidy,1)))

    if (l_process) then
      allocate(athresh(lbound(thresh_atidy,1):ubound(thresh_atidy,1)))
    end if

    current_l_negonly=.true.
    call recompute_constants(.true.)
    current_qin_l_negonly=.true.
    call recompute_qin_constants(.true.)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine initialise_mphystidy

  subroutine recompute_qin_constants(l_negonly)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='RECOMPUTE_QIN_CONSTANTS'

    logical, intent(in) :: l_negonly

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    qin_thresh=thresh_tidy
    if (l_negonly) then
      qin_thresh=0.0*qin_thresh
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine recompute_qin_constants

  subroutine recompute_constants(l_negonly)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='RECOMPUTE_CONSTANTS'

    logical, intent(in) :: l_negonly

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (l_negonly) then
      thresh=0.0
    else
      thresh=thresh_tidy
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine recompute_constants

  subroutine finalise_mphystidy()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='FINALISE_MPHYSTIDY'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (l_process) then
      deallocate(athresh)
    end if
    deallocate(thresh, qin_thresh)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine finalise_mphystidy

  subroutine qtidy(ixy_inner, dt, nz, qfields, procs, aerofields, aeroact, dustact, aeroice, dustliq, &
       aeroprocs, i_proc, i_aproc, l_negonly)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='QTIDY'

    integer, intent(in) :: ixy_inner
    integer, intent(in) :: nz
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: qfields(:,:), aerofields(:,:)
    type(aerosol_active), intent(in) :: aeroact(:), dustact(:), aeroice(:), dustliq(:)
    type(process_rate), intent(inout) :: procs(:,:)
    type(process_rate), intent(inout) :: aeroprocs(:,:)
    type(process_name), intent(in) :: i_proc, i_aproc
    logical, intent(in), optional :: l_negonly

    logical :: ql_reset, nl_reset, qr_reset, nr_reset, m3r_reset
    logical :: qi_reset, ni_reset, qs_reset, ns_reset, m3s_reset
    logical :: qg_reset, ng_reset, m3g_reset
    logical :: am4_reset, am5_reset, am7_reset, am8_reset, am9_reset
    logical :: an11_reset, an12_reset

    real(wp) :: dmass, dnumber

    logical :: l_qsig(0:ntotalq), l_qpos, l_qsmall, l_qsneg(0:ntotalq)
    !    l_qpos: q variable is positive
    !    l_qsmall:   q variable is positive, but below tidy threshold
    !    l_qsneg:    q variable is small or negative

    logical :: l_qice, l_qliquid

    logical :: l_apos, l_asmall, l_asneg(ntotala), l_asig(ntotala)

    integer :: iq, k

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_qsig(0)=.false.
    l_qsneg(0)=.false.
    
    do k = 1, nz
    nr_reset=.false.
    m3r_reset=.false.
    qi_reset=.false.
    qs_reset=.false.
    ns_reset=.false.
    m3s_reset=.false.
    qg_reset=.false.
    ng_reset=.false.
    m3g_reset=.false.

    am4_reset=.false.
    am5_reset=.false.
    am7_reset=.false.
    am8_reset=.false.
    am9_reset=.false.
    an11_reset=.false.
    an12_reset=.false.

    if (present(l_negonly)) then
      if (l_negonly .neqv. current_l_negonly) then
        current_l_negonly=l_negonly
        call recompute_constants(l_negonly)
      end if
    else if (current_l_negonly) then
      current_l_negonly=.false.
      call recompute_constants(.false.)
    end if

    do iq=1, ntotalq
      l_qsig(iq)=qfields(k, iq) > thresh(iq)

      l_qpos=qfields(k, iq) > 0.0

      l_qsmall = (.not. l_qsig(iq)) .and. l_qpos

      l_qsneg(iq)=qfields(k, iq) < 0.0 .or. l_qsmall
    end do

    if (l_process) then
      athresh=thresh_atidy
      if (present(l_negonly)) then
        if (l_negonly) athresh=0.0
      end if

      do iq=1, ntotala
        l_asig(iq)=aerofields(k, iq) > athresh(iq)

        l_apos=aerofields(k, iq) > 0.0

        l_asmall=(.not. l_asig(iq)) .and. l_apos

        l_asneg(iq)=aerofields(k, iq) < 0.0 .or. l_asmall
      end do
    end if

    l_qliquid=l_qsig(i_ql) .or. l_qsig(i_qr)
    l_qice=l_qsig(i_qi) .or. l_qsig(i_qs) .or. l_qsig(i_qg)

    ! Tidying of small and negative numbers and/or incompatible numbers (e.g.nl>0 and ql=0)
    ! - Mass and energy conserving...
    !==============================
    ! What should be reset?
    !==============================
    ql_reset=l_qsneg(i_ql)
    if (l_2mc) then
      nl_reset=l_qsneg(i_nl) .or. (l_qsig(i_nl) .and. ql_reset)
      ql_reset=ql_reset .or. (l_qsig(i_ql) .and. nl_reset)
    end if
    qr_reset=l_qsneg(i_qr)
    if (l_2mr) then
      nr_reset=l_qsneg(i_nr) .or. (l_qsig(i_nr) .and. qr_reset)
      qr_reset=qr_reset .or. (l_qsig(i_qr) .and. nr_reset)
    end if
    ! if (l_3mr .and. l_tidym3) then
    !   m3r_reset=l_qsneg(i_m3r) .or. (l_qsig(i_m3r) .and. (qr_reset .or. nr_reset))
    !   nr_reset=nr_reset .or. (l_qsig(i_nr) .and. m3r_reset)
    !   qr_reset=qr_reset .or. (l_qsig(i_qr) .and. m3r_reset)
    ! end if

    if (.not. l_warm) then
      qi_reset=l_qsneg(i_qi)
      if (l_2mi) then
        ni_reset=l_qsneg(i_ni) .or. (l_qsig(i_ni) .and. qi_reset)
        qi_reset=qi_reset .or. (l_qsig(i_qi) .and. ni_reset)
      end if

      qs_reset=l_qsneg(i_qs)
      if (l_2ms) then
        ns_reset=l_qsneg(i_ns) .or. (l_qsig(i_ns) .and. qs_reset)
        qs_reset=qs_reset .or. (l_qsig(i_qs) .and. ns_reset)
      end if
      ! if (l_3ms .and. l_tidym3) then
      !   m3s_reset=l_qsneg(i_m3s) .or. (l_qsig(i_m3s) .and. (qs_reset .or. ns_reset))
      !   ns_reset=ns_reset .or. (l_qsig(i_ns) .and. m3s_reset)
      !   qs_reset=qs_reset .or. (l_qsig(i_qs) .and. m3s_reset)
      ! end if

      qg_reset=l_qsneg(i_qg)
      if (l_2mg) then
        ng_reset=l_qsneg(i_ng) .or. (l_qsig(i_ng) .and. qg_reset)
        qg_reset=qg_reset .or. (l_qsig(i_qg) .and. ng_reset)
      end if
      ! if (l_3mg .and. l_tidym3) then
      !   m3g_reset=l_qsneg(i_m3g) .or. (l_qsig(i_m3g) .and. (qg_reset .or. ng_reset))
      !   ng_reset=ng_reset .or. (l_qsig(i_ng) .and. m3g_reset)
      !   qg_reset=qg_reset .or. (l_qsig(i_qg) .and. m3g_reset)
      ! end if
    end if

    !===========================================================
    ! Aerosol tests...
    !===========================================================
    if (l_process) then
      ! Aerosols in liquid water

      ! If small/neg values...
      if (l_asneg(i_am4))am4_reset=.true.
      if (l_passivenumbers) then
        if (l_asneg(i_an11))an11_reset=.true.  ! can be in liq or ice
      endif
      if (l_passivenumbers_ice) then
        if (l_asneg(i_an12))an12_reset=.true.  ! can be in liq or ice
      endif

      if (l_separate_rain) then
        if (l_asneg(i_am5))am5_reset=.true.
      end if
      if (.not. l_Warm) then
        if (l_asneg(i_am9))am9_reset=.true.
      end if
      ! If no hydrometeors...
      if ((ql_reset .and. qr_reset) .or. .not. l_qliquid) then
        if (l_asig(i_am4))am4_reset=.true.
        if (.not.l_warm) then
          if (l_asig(i_am9))am9_reset=.true.
        end if
        if (l_separate_rain) then
          if (l_asig(i_am5)) am5_reset=.true.
        end if
      end if

      ! If no active aerosol, then we shouldn't have any hydrometeor...(what about SIP?)
      ql_reset=ql_reset .or. (am4_reset .and. am9_reset .and. l_qsig(i_ql))
      qr_reset=qr_reset .or. (am4_reset .and. am9_reset .and. l_qsig(i_qr))
      qr_reset=qr_reset .or. (am5_reset .and. l_qsig(i_qr))

      ! Aerosols in ice
      ! If small/neg values...
      if (.not. l_Warm)then
        if (l_asneg(i_am7)) am7_reset=.true.
        if (l_asneg(i_am8)) then
          am8_reset=.true.
        end if
        ! If no hydrometeors...
        if ((qi_reset .and. qs_reset .and. qg_reset) .or. .not. l_qice) then
          if (l_asig(i_am7)) am7_reset=.true.
          if (l_asig(i_am8)) then
            am8_reset=.true.
          end if
        end if
      end if

      ! If no active aerosol, then we shouldn't have any hydrometeor...
      qi_reset=qi_reset .or. (am7_reset .and. am8_reset .and. l_qsig(i_qi))
      qs_reset=qs_reset .or. (am7_reset .and. am8_reset .and. l_qsig(i_qs))
      qg_reset=qg_reset .or. (am7_reset .and. am8_reset .and. l_qsig(i_qg))
    end if

    !===========================================================
    ! Consistency following aerosol
    !===========================================================
    nl_reset=ql_reset .and. l_qsig(i_nl)
    nr_reset=nr_reset .and. l_qsig(i_nr)
    m3r_reset=m3r_reset .and. l_qsig(i_m3r)
    ni_reset=qi_reset .and. l_qsig(i_ni)
    ns_reset=ns_reset .and. l_qsig(i_ns)
    m3s_reset=m3s_reset .and. l_qsig(i_m3s)
    ng_reset=ng_reset .and. l_qsig(i_ng)
    m3g_reset=m3g_reset .and. l_qsig(i_m3g)

    !==============================
    ! Now reset things...
    !==============================
    if (ql_reset .or. nl_reset) then
      dmass=qfields(k, i_ql)/dt
      procs(i_ql,i_proc%id)%column_data(k)=-dmass
      if (l_tidy_conserve_q)procs(i_qv,i_proc%id)%column_data(k)=dmass
      if (l_tidy_conserve_E)procs(i_th,i_proc%id)%column_data(k)=-Lv*dmass/cp/exner(k,ixy_inner)
      if (l_2mc) then
        dnumber=qfields(k, i_nl)/dt
        procs(i_nl,i_proc%id)%column_data(k)=-dnumber
      end if
      !--------------------------------------------------
      ! aerosol - not reset, but adjusted acordingly
      !--------------------------------------------------
      if (l_process) then
        if (.not. am4_reset) then
          dmass=aeroact(k)%mact1/dt
          dnumber=aeroact(k)%nact1/dt
          if (dmass > 0.0) then
            aeroprocs(i_am4,i_aproc%id)%column_data(k)=-dmass
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=dnumber
          end if
        end if
        if (.not. am9_reset) then
          dmass=dustliq(k)%mact1/dt
          dnumber=dustliq(k)%nact1/dt
          if (dmass > 0.0) then
            aeroprocs(i_am9,i_aproc%id)%column_data(k)=-dmass
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=dnumber
          end if
        end if
      end if
    end if

    if (qr_reset .or. nr_reset .or. m3r_reset) then
      dmass=qfields(k, i_qr)/dt
      procs(i_qr,i_proc%id)%column_data(k)=-dmass
      if (l_tidy_conserve_q) procs(i_qv,i_proc%id)%column_data(k)=             &
                                   procs(i_qv,i_proc%id)%column_data(k) + dmass
      if (l_tidy_conserve_E) procs(i_th,i_proc%id)%column_data(k)=             &
          procs(i_th,i_proc%id)%column_data(k) - Lv*dmass/cp/exner(k,ixy_inner)
      if (l_2mr) then
        dnumber=qfields(k, i_nr)/dt
        procs(i_nr,i_proc%id)%column_data(k)=-dnumber
      end if
      ! if (l_3mr) then
      !   procs(i_m3r,i_proc%id)%column_data(k)=-qfields(k, i_m3r)/dt
      ! end if
      !--------------------------------------------------
      ! aerosol
      !--------------------------------------------------
      if (l_process) then
        if (.not. am4_reset) then
          dmass=aeroact(k)%mact2/dt
          dnumber=aeroact(k)%nact2/dt
          aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
          aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
          aeroprocs(i_am4,i_aproc%id)%column_data(k)=aeroprocs(i_am4,i_aproc%id)%column_data(k)-dmass
        end if
        if (.not. am5_reset) then
          if (l_separate_rain) then
            dmass=aeroact(k)%mact2/dt
            dnumber=aeroact(k)%nact2/dt
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am5,i_aproc%id)%column_data(k)=aeroprocs(i_am5,i_aproc%id)%column_data(k)-dmass
          end if
        end if
        if (.not. am9_reset .and. .not. l_warm) then
          dmass=dustliq(k)%mact2/dt
          dnumber=dustliq(k)%nact2/dt
          if (dmass>0.0) then
            aeroprocs(i_am9,i_aproc%id)%column_data(k)=-dmass
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=dnumber
          end if
        end if
      end if
    end if

    if (.not. l_warm) then
      if (qi_reset .or. ni_reset) then
        dmass=qfields(k, i_qi)/dt
        procs(i_qi,i_proc%id)%column_data(k)=-dmass
        if (l_tidy_conserve_q) procs(i_qv,i_proc%id)%column_data(k)=           &
                                     procs(i_qv,i_proc%id)%column_data(k)+dmass
        if (l_tidy_conserve_E) procs(i_th,i_proc%id)%column_data(k)=           &
            procs(i_th,i_proc%id)%column_data(k)-Ls*dmass/cp/exner(k,ixy_inner)
        if (l_2mi) then
          dnumber=qfields(k, i_ni)/dt
          procs(i_ni,i_proc%id)%column_data(k)=-dnumber
        end if
        !--------------------------------------------------
        !aerosol
        !--------------------------------------------------
        if (l_process) then
          if (.not. am7_reset) then
            dmass=dustact(k)%mact1/dt
            dnumber=dustact(k)%nact1/dt
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=aeroprocs(i_an6,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am7,i_aproc%id)%column_data(k)=aeroprocs(i_am7,i_aproc%id)%column_data(k)-dmass
          end if
          if (.not. am8_reset) then
            dmass=aeroice(k)%mact1/dt
            dnumber=aeroice(k)%nact1/dt
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am8,i_aproc%id)%column_data(k)=aeroprocs(i_am8,i_aproc%id)%column_data(k)-dmass
          end if
        end if
      end if

      if (qs_reset .or. ns_reset .or. m3s_reset) then
        dmass=qfields(k, i_qs)/dt
        procs(i_qs,i_proc%id)%column_data(k)=-dmass
        if (l_tidy_conserve_q) procs(i_qv,i_proc%id)%column_data(k)=           &
                                     procs(i_qv,i_proc%id)%column_data(k)+dmass
        if (l_tidy_conserve_E) procs(i_th,i_proc%id)%column_data(k)=           &
            procs(i_th,i_proc%id)%column_data(k)-Ls*dmass/cp/exner(k,ixy_inner)
        if (l_2ms) then
          dnumber=qfields(k, i_ns)/dt
          procs(i_ns,i_proc%id)%column_data(k)=-dnumber
        end if
        ! if (l_3ms) then
        !   procs(i_m3s,i_proc%id)%column_data(k)=-qfields(k, i_m3s)/dt
        ! end if
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
        if (l_process) then
          if (.not. am7_reset) then
            dmass=dustact(k)%mact2/dt
            dnumber=dustact(k)%nact2/dt
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=aeroprocs(i_an6,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am7,i_aproc%id)%column_data(k)=aeroprocs(i_am7,i_aproc%id)%column_data(k)-dmass
          end if
          if (.not. am8_reset) then
            dmass=aeroice(k)%mact2/dt
            dnumber=aeroice(k)%nact2/dt
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am8,i_aproc%id)%column_data(k)=aeroprocs(i_am8,i_aproc%id)%column_data(k)-dmass
          end if
        end if
      end if

      if (qg_reset .or. ng_reset .or. m3g_reset) then
        dmass=qfields(k, i_qg)/dt
        procs(i_qg,i_proc%id)%column_data(k)=-dmass
        if (l_tidy_conserve_q)procs(i_qv,i_proc%id)%column_data(k)=            &
                                     procs(i_qv,i_proc%id)%column_data(k)+dmass
        if (l_tidy_conserve_E)procs(i_th,i_proc%id)%column_data(k)=            &
            procs(i_th,i_proc%id)%column_data(k)-Ls*dmass/cp/exner(k,ixy_inner)
        if (l_2mg) then
          dnumber=qfields(k, i_ng)/dt
          procs(i_ng,i_proc%id)%column_data(k)=-dnumber
        end if
        ! if (l_3mg) then
        !   procs(i_m3g,i_proc%id)%column_data(k)=-qfields(k, i_m3g)/dt
        ! end if
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
        if (l_process) then
          if (.not. am7_reset) then
            dmass=dustact(k)%mact3/dt
            dnumber=dustact(k)%nact3/dt
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=aeroprocs(i_an6,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am7,i_aproc%id)%column_data(k)=aeroprocs(i_am7,i_aproc%id)%column_data(k)-dmass
          end if
          if (.not. am8_reset) then
            dmass=aeroice(k)%mact3/dt
            dnumber=aeroice(k)%nact3/dt
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am8,i_aproc%id)%column_data(k)=aeroprocs(i_am8,i_aproc%id)%column_data(k)-dmass
          end if
        end if
      end if
    end if
    !==============================
    ! Now reset aerosol...
    !==============================

    if (l_process) then
      if (am4_reset) then
        dmass=aerofields(k, i_am4)/dt
        aeroprocs(i_am4,i_aproc%id)%column_data(k)=aeroprocs(i_am4,i_aproc%id)%column_data(k)-dmass
        aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
      end if

      if (am5_reset .and. l_separate_rain) then
        dmass=aerofields(k, i_am5)/dt
        aeroprocs(i_am5,i_aproc%id)%column_data(k)=aeroprocs(i_am5,i_aproc%id)%column_data(k)-dmass
        aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
      end if

      if (.not. l_warm) then
        if (am7_reset) then
          dmass=aerofields(k, i_am7)/dt
          aeroprocs(i_am7,i_aproc%id)%column_data(k)=aeroprocs(i_am7,i_aproc%id)%column_data(k)-dmass
          aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
        end if

        if (am8_reset) then
          dmass=aerofields(k, i_am8)/dt
          aeroprocs(i_am8,i_aproc%id)%column_data(k)=aeroprocs(i_am8,i_aproc%id)%column_data(k)-dmass
          aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
        end if

        if (am9_reset) then
          dmass=aerofields(k, i_am9)/dt
          aeroprocs(i_am9,i_aproc%id)%column_data(k)=aeroprocs(i_am9,i_aproc%id)%column_data(k)-dmass
          aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
        end if

        !only do this is passive numbers are used!
        if ( l_passivenumbers ) then        
          if (an11_reset) then  !just reset it and put number into accum sol
            aeroprocs(i_an11,i_aproc%id)%column_data(k)=-aerofields(k,i_an11)/dt
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aerofields(k,i_an11)/dt
          end if
        endif
        if ( l_passivenumbers_ice ) then        
          if (an12_reset) then  !just reset it and put number into coarse insol
            aeroprocs(i_an12,i_aproc%id)%column_data(k)=-aerofields(k,i_an12)/dt
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=aerofields(k,i_an12)/dt
          end if
        endif

      end if
    end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine qtidy

  subroutine tidy_qin(ixy_inner, qfields, l_negonly)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='TIDY_QIN'

    integer, intent(in) :: ixy_inner
    real(wp), intent(inout) :: qfields(:,:)
    logical, intent(in), optional :: l_negonly

    logical :: ql_reset, nl_reset, qr_reset, nr_reset, m3r_reset
    logical :: qi_reset, ni_reset, qs_reset, ns_reset, m3s_reset
    logical :: qg_reset, ng_reset, m3g_reset
    integer :: k

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (present(l_negonly)) then
      if (l_negonly .neqv. current_qin_l_negonly) then
        current_qin_l_negonly=l_negonly
        call recompute_qin_constants(l_negonly)
      end if
    else if (current_qin_l_negonly) then
      current_qin_l_negonly=.false.
      call recompute_qin_constants(.false.)
    end if

    do k=1, ubound(qfields,1)
      nl_reset=.false.
      nr_reset=.false.
      m3r_reset=.false.
      qi_reset=.false.
      ni_reset=.false.
      qs_reset=.false.
      ns_reset=.false.
      m3s_reset=.false.
      qg_reset=.false.
      ng_reset=.false.
      m3g_reset=.false.

      ql_reset=qfields(k, i_ql) < 0.0 .or. (qfields(k, i_ql) < qin_thresh(i_ql) .and. qfields(k, i_ql) >0)
      if (l_2mc) then
        nl_reset=qfields(k, i_nl) < 0.0 .or. (qfields(k, i_nl) < qin_thresh(i_nl) .and. qfields(k, i_nl) >0) .or. &
             (qfields(k, i_nl) > 0.0 .and. qfields(k, i_ql) <= 0.0)
        ql_reset=ql_reset .or. (qfields(k, i_ql) > 0.0 .and. qfields(k, i_nl) <= 0.0)
      end if
      qr_reset=qfields(k, i_qr) < 0.0 .or. (qfields(k, i_qr) < qin_thresh(i_qr) .and. qfields(k, i_qr) > 0.0)

      if (l_2mr) then
        nr_reset=qfields(k, i_nr) < 0.0 .or. (qfields(k, i_nr) < qin_thresh(i_nr) .and. qfields(k, i_nr) > 0.0) .or. &
             (qfields(k, i_nr) > 0.0 .and. qfields(k, i_qr) <= 0.0)
        qr_reset=qr_reset .or. (qfields(k, i_qr) > 0.0 .and. qfields(k, i_nr) <= 0.0)
      end if

      ! if (l_3mr .and. l_tidym3) then
      !   m3r_reset=qfields(k, i_m3r) < 0.0 .or. (qfields(k, i_m3r) < qin_thresh(i_m3r) .and. qfields(k, i_m3r) >0) .or. &
      !        (qfields(k, i_m3r) > 0.0 .and. (qfields(k, i_qr) <=0.0 .or. qfields(k, i_nr) <=0.0))
      !   qr_reset=qr_reset .or. (qfields(k, i_qr) > 0.0 .and. qfields(k, i_m3r) <= 0.0)
      !   nr_reset=nr_reset .or. (qfields(k, i_nr) > 0.0 .and. qfields(k, i_m3r) <= 0.0)
      ! end if

      if (.not. l_warm) then
        qi_reset=qfields(k, i_qi) < 0.0 .or. (qfields(k, i_qi) < qin_thresh(i_qi) .and. qfields(k, i_qi) > 0.0)
        if (l_2mi) then
          ni_reset=qfields(k, i_ni) < 0.0 .or. (qfields(k, i_ni) < qin_thresh(i_ni) .and. qfields(k, i_ni) > 0.0) .or. &
               (qfields(k, i_ni) > 0.0 .and. qfields(k, i_qi) <= 0.0)
          qi_reset=qi_reset .or. (qfields(k, i_qi) > 0.0 .and. qfields(k, i_ni) <= 0.0)
        end if

        qs_reset=qfields(k, i_qs) < 0.0 .or. (qfields(k, i_qs) < qin_thresh(i_qs) .and. qfields(k, i_qs) > 0.0)
        if (l_2ms) then
          ns_reset=qfields(k, i_ns) < 0.0 .or. (qfields(k, i_ns) < qin_thresh(i_ns) .and. qfields(k, i_ns) > 0.0) .or. &
               (qfields(k, i_ns) > 0.0 .and. qfields(k, i_qs) <= 0.0)
          qs_reset=qs_reset .or. (qfields(k, i_qs) > 0.0 .and. qfields(k, i_ns) <= 0.0)
        end if
        ! if (l_3ms .and. l_tidym3) then
        !   m3s_reset=qfields(k, i_m3s) < 0.0 .or. (qfields(k, i_m3s) < qin_thresh(i_m3s) .and. qfields(k, i_m3s) >0) .or.&
        !        (qfields(k, i_m3s) > 0.0 .and. (qfields(k, i_qs) <=0.0 .or. qfields(k, i_ns) <=0.0))
        !   qs_reset=qs_reset .or. (qfields(k, i_qs) > 0.0 .and. qfields(k, i_m3s) <= 0.0)
        !   ns_reset=ns_reset .or. (qfields(k, i_ns) > 0.0 .and. qfields(k, i_m3s) <= 0.0)
        ! end if

        qg_reset=qfields(k, i_qg) < 0.0 .or. (qfields(k, i_qg) < qin_thresh(i_qg) .and. qfields(k, i_qg) > 0.0)
        if (l_2mg) then
          ng_reset=qfields(k, i_ng) < 0.0 .or. (qfields(k, i_ng) < qin_thresh(i_ng) .and. qfields(k, i_ng) > 0.0) .or. &
               (qfields(k, i_ng) > 0.0 .and. qfields(k, i_qg) <= 0.0)
          qg_reset=qg_reset .or. (qfields(k, i_qg) > 0.0 .and. qfields(k, i_ng) <= 0.0)
        end if
        ! if (l_3mg .and. l_tidym3) then
        !   m3g_reset=qfields(k, i_m3g) < 0.0 .or. (qfields(k, i_m3g) < qin_thresh(i_m3g) .and. qfields(k, i_m3g) >0) .or.&
        !        (qfields(k, i_m3g) > 0.0 .and. (qfields(k, i_qg) <=0.0 .or. qfields(k, i_ng) <=0.0))
        !   qg_reset=qg_reset .or. (qfields(k, i_qg) > 0.0 .and. qfields(k, i_m3g) <= 0.0)
        !   ng_reset=ng_reset .or. (qfields(k, i_ng) > 0.0 .and. qfields(k, i_m3g) <= 0.0)
        ! end if
      end if

      !==============================
      ! Now reset things...
      !==============================
      if (ql_reset .or. nl_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Lv/cp*qfields(k,i_ql)/exner(k,ixy_inner)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_ql)
        qfields(k,i_ql)=0.0
        if (l_2mc) then
          qfields(k,i_nl)=0.0
        end if
      end if

      if (qr_reset .or. nr_reset .or. m3r_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Lv/cp*qfields(k,i_qr)/exner(k,ixy_inner)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qr)
        qfields(k,i_qr)=0.0
        if (l_2mr) then
          qfields(k,i_nr)=0.0
        end if
        ! if (l_3mr) then
        !   qfields(k,i_m3r)=0.0
        ! end if
      end if

      if (qi_reset .or. ni_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qi)/exner(k,ixy_inner)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qi)
        qfields(k,i_qi)=0.0
        if (l_2mi) then
          qfields(k,i_ni)=0.0
        end if
      end if

      if (qs_reset .or. ns_reset .or. m3s_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qs)/exner(k,ixy_inner)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qs)
        qfields(k,i_qs)=0.0
        if (l_2ms) then
          qfields(k,i_ns)=0.0
        end if
        ! if (l_3ms) then
        !   qfields(k,i_m3s)=0.0
        ! end if
      end if

      if (qg_reset .or. ng_reset .or. m3g_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qg)/exner(k,ixy_inner)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qg)
        qfields(k,i_qg)=0.0
        if (l_2mg) then
          qfields(k,i_ng)=0.0
        end if
        ! if (l_3mg) then
        !   qfields(k,i_m3g)=0.0
        ! end if
      end if

    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine tidy_qin

  subroutine tidy_ain(qfields, aerofields)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='TIDY_AIN'

    real(wp), intent(in) :: qfields(:,:)
    real(wp), intent(inout) :: aerofields(:,:)

    integer :: k

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k=1, ubound(qfields,1)
      if ((qfields(k, i_ql)+qfields(k,i_qr) <=0.0 .and. aerofields(k,i_am4)>0.0) .or. aerofields(k,i_am4) < 0.0) then

        aerofields(k,i_am4)=0.0
      end if
      if (i_am9 > 0)then
        if ((qfields(k, i_ql)+qfields(k,i_qr) <=0.0 .and. aerofields(k,i_am9)>0.0) .or. aerofields(k,i_am9) < 0.0) then
          aerofields(k,i_am9)=0.0
        end if
      end if

      if (i_am5 > 0) then
        if (((qfields(k,i_qr) <=0.0 .and. aerofields(k,i_am5)>0.0) .or. aerofields(k,i_am5) < 0.0)) then
          aerofields(k,i_am5)=0.0
        end if
      end if
      if (i_am7 > 0) then
        if (((qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) <=0.0 .and. aerofields(k,i_am7)>0.0) &
             .or. aerofields(k,i_am7) < 0.0)) then
          aerofields(k,i_am7)=0.0
        end if
      end if
      if (i_am8 > 0) then
        if (((qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) <=0.0 .and. aerofields(k,i_am8)>0.0) &
             .or. aerofields(k,i_am8) < 0.0)) then
          aerofields(k,i_am8)=0.0
        end if
      end if
      if (i_an12 > 0) then
        if (((qfields(k, i_ql) + qfields(k,i_qr) + qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) &
                                                                 <=0.0 .and. aerofields(k,i_an12)>0.0) &
                                                            .or. aerofields(k,i_an12) < 0.0)) then
          aerofields(k,i_an12)=0.0
        end if
      end if
      if (i_an11 > 0) then
        if (((qfields(k, i_ql) + qfields(k,i_qr) + qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) &
                                                                 <=0.0 .and. aerofields(k,i_an12)>0.0) &
                                                            .or. aerofields(k,i_an11) < 0.0)) then
          aerofields(k,i_an11)=0.0
        end if
      end if

    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine tidy_ain

  ! Subroutine to ensure parallel processes don't remove more
    ! mass than is available and then rescales all processes
    ! (including number and other terms)
  subroutine ensure_positive(nz, dt, qfields, procs, params, iprocs_scalable, &
                             iprocs_nonscalable, aeroprocs, iprocs_dependent, &
                             iprocs_dependent_ns)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ENSURE_POSITIVE'

    integer, intent(in) :: nz
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: qfields(:,:)
    type(process_rate), intent(inout) :: procs(:,:)         ! microphysical process rates

    type(hydro_params), intent(in) :: params        ! parameters from hydrometeor variable to test
    type(process_name), intent(in) :: iprocs_scalable(:)    ! list of processes to rescale
    type(process_name), intent(in), optional ::      &
         iprocs_nonscalable(:) ! list of other processes which
    ! provide source or sink, but
    ! which we don't want to rescale
    type(process_rate), intent(inout), optional ::   &
         aeroprocs(:,:)        ! associated aerosol process rates
    type(process_name), intent(in), optional ::      &
         iprocs_dependent(:)   ! list of aerosol processes which
    ! are dependent on rescaled processes and
    ! so should be rescaled themselves
    type(process_name), intent(in), optional ::      &
         iprocs_dependent_ns(:)   ! list of aerosol processes which
    ! are dependent on rescaled processes but
    ! we don't want to rescale
    integer :: iproc, id, iq, k
    real(wp) :: delta_scalable, delta_nonscalable, ratio, maxratio
    integer :: i_1m, i_2m, i_3m
    logical :: l_rescaled

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    do k=1,nz
    i_1m=params%i_1m
    if (params%l_2m) i_2m=params%i_2m
    if (params%l_3m) i_3m=params%i_3m
    l_rescaled=.false.

    if (qfields(k, i_1m) > 0.0) then
      ! First calculate the unscaled increments
      delta_scalable=0.0
      delta_nonscalable=0.0
      do iproc=1,size(iprocs_scalable)
        if (iprocs_scalable(iproc)%on) then
          id=iprocs_scalable(iproc)%id
          delta_scalable=delta_scalable+procs(i_1m, id)%column_data(k)
        end if
      end do
      delta_scalable=delta_scalable*dt
      if (present(iprocs_nonscalable)) then
        do iproc=1, size(iprocs_nonscalable)
          if (iprocs_nonscalable(iproc)%on) then
            id=iprocs_nonscalable(iproc)%id
            delta_nonscalable=delta_nonscalable+procs(i_1m, id)%column_data(k)
          end if
        end do
      end if
      delta_nonscalable=delta_nonscalable*dt

      ! Test to see if we need to do anything
      if (delta_scalable+delta_nonscalable+qfields(k, i_1m) < spacing(qfields(k, i_1m)) &
           .and. abs(delta_scalable) > epsilon(delta_scalable)) then
        ratio=-(qfields(k, i_1m)+delta_nonscalable)/delta_scalable
        if (ratio > 1.0) then
          if (present(iprocs_nonscalable)) then
            do iproc=1, size(iprocs_nonscalable)
              if (iprocs_nonscalable(iproc)%on) then
                id=iprocs_nonscalable(iproc)%id
              end if
            end do
          end if
          do iproc=1, size(iprocs_scalable)
            if (iprocs_scalable(iproc)%on) then
              id=iprocs_scalable(iproc)%id
            end if
          end do

          do iproc=1, size(iprocs_scalable)
            if (iprocs_scalable(iproc)%on) then
              id=iprocs_scalable(iproc)%id
            end if
          end do
          do iproc=1, size(iprocs_nonscalable)
            if (iprocs_nonscalable(iproc)%on) then
              id=iprocs_nonscalable(iproc)%id
            end if
          end do

          write(std_msg, '(A, F7.4)') 'Problem with ratio > 1.0: ratio = ', ratio
          ! Tell mphys_error that this is due to bad values
          call throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)

        end if
        ! Now rescale the scalable processes
        do iproc=1, size(iprocs_scalable)
          if (iprocs_scalable(iproc)%on) then
            id=iprocs_scalable(iproc)%id
            do iq = i_qstart, i_nstart-1
               procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
            enddo
          end if
        end do
        ! Set flag to indicate a rescaling was performed
        l_rescaled=.true.
      end if

      ! Now we need to rescale additional moments
      if (l_rescaled) then
        if (params%l_2m) then ! second moment
          ! First calculate the unscaled increments
          delta_scalable=0.0
          delta_nonscalable=0.0
          do iproc=1,size(iprocs_scalable)
            if (iprocs_scalable(iproc)%on) then
              id=iprocs_scalable(iproc)%id
              delta_scalable=delta_scalable+procs(i_2m, id)%column_data(k)
            end if
          end do
          delta_scalable=delta_scalable*dt
          if (abs(delta_scalable) > epsilon(delta_scalable)) then
            if (present(iprocs_nonscalable)) then
              do iproc=1, size(iprocs_nonscalable)
                if (iprocs_nonscalable(iproc)%on) then
                  id=iprocs_nonscalable(iproc)%id
                  delta_nonscalable=delta_nonscalable+procs(i_2m, id)%column_data(k)
                end if
              end do
            end if
            delta_nonscalable=delta_nonscalable*dt
            ! ratio may now be greater than 1
            ratio=-(qfields(k, i_2m)+delta_nonscalable)/delta_scalable
            ! Now rescale the scalable processes
            do iproc=1, size(iprocs_scalable)
              if (iprocs_scalable(iproc)%on) then
                id=iprocs_scalable(iproc)%id
                do iq = i_nstart, i_m3start-1
                   procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
                enddo
              end if
            end do
          end if
        end if
        ! ! if (params%l_3m) then ! third moment
        ! !   ! First calculate the unscaled increments
        ! !   delta_scalable=0.0
        ! !   delta_nonscalable=0.0
        ! !   do iproc=1,size(iprocs_scalable)
        ! !     if (iprocs_scalable(iproc)%on) then
        ! !       id=iprocs_scalable(iproc)%id
        ! !       delta_scalable=delta_scalable+procs(i_3m, id)%column_data(k)
        ! !     end if
        ! !   end do
        ! !   delta_scalable=delta_scalable*dt
        ! !   if (abs(delta_scalable) > epsilon(delta_scalable)) then
        ! !     if (present(iprocs_nonscalable)) then
        ! !       do iproc=1, size(iprocs_nonscalable)
        ! !         if (iprocs_nonscalable(iproc)%on) then
        ! !           id=iprocs_nonscalable(iproc)%id
        ! !           delta_nonscalable=delta_nonscalable+procs(i_3m, id)%column_data(k)
        ! !         end if
        ! !       end do
        ! !     end if
        ! !     delta_nonscalable=delta_nonscalable*dt
        ! !     ! ratio may now be greater than 1
        ! !     ratio=-(qfields(k, i_3m)+delta_nonscalable)/(delta_scalable + epsilon(1.0))
        ! !     ! Now rescale the scalable processes
        ! !     do iproc=1, size(iprocs_scalable)
        ! !       if (iprocs_scalable(iproc)%on) then
        ! !         id=iprocs_scalable(iproc)%id
        ! !         do iq = i_m3start, ubound(procs,1)
        ! !            procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
        ! !         enddo
        ! !       end if
        ! !     end do
        ! !   end if
        ! ! end if

        !Now rescale the increments to aerosol
        ! How do we do this?????
        !        print*, 'WARNING: Should be rescaling aerosol?, but not done!'

      else ! What if we haven't rescaled mass, but number is now not conserved?
        if (l_rescale_on_number) then
          if (params%l_2m) then ! second moment
            ! First calculate the unscaled increments
            delta_scalable=0.0
            delta_nonscalable=0.0
            do iproc=1, size(iprocs_scalable)
              if (iprocs_scalable(iproc)%on) then
                id=iprocs_scalable(iproc)%id
                delta_scalable=delta_scalable+procs(i_2m, id)%column_data(k)
              end if
            end do
            delta_scalable=delta_scalable*dt
            if (present(iprocs_nonscalable)) then
              do iproc=1, size(iprocs_nonscalable)
                if (iprocs_nonscalable(iproc)%on) then
                  id=iprocs_nonscalable(iproc)%id
                  delta_nonscalable=delta_nonscalable+procs(i_2m, id)%column_data(k)
                end if
              end do
            end if
            delta_nonscalable=delta_nonscalable*dt
            ! Test to see if we need to do anything
            if (delta_scalable+delta_nonscalable+qfields(k, i_2m) < spacing(qfields(k, i_2m)) &
                 .and. abs(delta_scalable) > epsilon(delta_scalable) ) then
              maxratio=(1.0-spacing(delta_scalable))
              if (abs(delta_scalable) < spacing(delta_scalable)) then
                if (present(iprocs_nonscalable)) then
                  do iproc=1, size(iprocs_nonscalable)
                    if (iprocs_nonscalable(iproc)%on) then
                      id=iprocs_nonscalable(iproc)%id
                    end if
                  end do
                end if
                do iproc=1, size(iprocs_scalable)
                  if (iprocs_scalable(iproc)%on) then
                    id=iprocs_scalable(iproc)%id
                  end if
                end do
              end if
              ratio=-(maxratio*qfields(k, i_2m)+delta_nonscalable)/delta_scalable
              ratio=max(ratio, 0.0_wp)
              if (ratio==0.0_wp) then
                if (present(iprocs_nonscalable)) then
                  do iproc=1, size(iprocs_nonscalable)
                    if (iprocs_nonscalable(iproc)%on) id=iprocs_nonscalable(iproc)%id
                  end do
                end if
                do iproc=1, size(iprocs_scalable)
                  if (iprocs_scalable(iproc)%on) id=iprocs_scalable(iproc)%id
                end do
              end if

              if (ratio<0.95) then
                ! Some warnings for testing

                write(std_msg, *) 'WARNING: Significantly rescaled number, but not sure ' // &
                                   'what to do with other moments. id, ratio, bad',           &
                                    params%id, ratio, delta_scalable + delta_nonscalable +    &
                                    qfields(k, i_2m), spacing(qfields(k, i_2m)),              &
                                    (qfields(k, i_2m) + delta_nonscalable), delta_scalable,   &
                                    'qfields', qfields(k, i_1m), qfields(k, i_2m)

                call throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)


                if (present(iprocs_nonscalable)) then
                  do iproc=1, size(iprocs_nonscalable)
                    if (iprocs_nonscalable(iproc)%on) then
                      id=iprocs_nonscalable(iproc)%id
                    end if
                  end do
                end if
                do iproc=1, size(iprocs_scalable)
                  if (iprocs_scalable(iproc)%on) then
                    id=iprocs_scalable(iproc)%id
                  end if
                end do
              end if

              ! Now rescale the scalable processes
              do iproc=1, size(iprocs_scalable)
                if (iprocs_scalable(iproc)%on) then
                  id=iprocs_scalable(iproc)%id
                  if (i_2m > 0) then
                     do iq = i_nstart, i_m3start-1
                        procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
                     enddo
                  end if
                end if
              end do

              ! Set flag to indicate a rescaling was performed
              l_rescaled=.true.
            end if
          end if
        end if
      end if
    end if
    end do
  
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ensure_positive

  ! Subroutine to ensure parallel aerosol processes don't remove more
  ! mass than is available and then rescales all processes
  ! (NB this follows any rescaling due to the parent microphysical processes and
  ! we might lose consistency between number and mass here)
  subroutine ensure_positive_aerosol(nz, dt, aerofields, aerosol_procs, iprocs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ENSURE_POSITIVE_AEROSOL'

    integer, intent(in) :: nz
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: aerofields(:,:)
    type(process_rate), intent(inout) :: aerosol_procs(:,:)  ! aerosol process rates
    type(process_name), intent(in) :: iprocs(:)    ! list of processes to rescale

    integer :: iq, iproc, id, k
    real(wp) :: ratio, delta_scalable

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k=1,nz
      do iq=1, ntotala
        delta_scalable=0.0
        do iproc=1, size(iprocs)
          if (iprocs(iproc)%on) then
            id=iprocs(iproc)%id
            delta_scalable=delta_scalable + aerosol_procs(iq, id)%column_data(k)
          end if
        end do
        delta_scalable=delta_scalable*dt
        if (delta_scalable + aerofields(k, iq) < spacing(aerofields(k, iq))    &
             .and. abs(delta_scalable) > spacing(aerofields(k, iq))) then
          ratio=(spacing(aerofields(k, iq))-aerofields(k, iq))/(delta_scalable)

          do iproc=1 ,size(iprocs)
            if (iprocs(iproc)%on) then
              id=iprocs(iproc)%id
              aerosol_procs(iq,id)%column_data(k)=aerosol_procs(iq,id)%column_data(k)*ratio
            end if
          end do
        end if
      end do
    end do
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ensure_positive_aerosol

  ! Subroutine to ensure parallel ice processes don't remove more
    ! vapour than is available (i.e. so become subsaturated)
    ! and then rescales processes

    ! Modified so it can also prevent sublimation processes from putting
    ! back too much vapour and so become supersaturated
  subroutine ensure_saturated(ixy_inner, nz, l_Tcold, dt, qfields, procs, iprocs_scalable)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ENSURE_SATURATED'

    integer, intent(in) :: ixy_inner
    integer, intent(in) :: nz
    logical, intent(in) :: l_Tcold(:) 
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: qfields(:,:)
    type(process_rate), intent(inout) :: procs(:,:)         ! microphysical process rates
    type(process_name), intent(in) :: iprocs_scalable(:)    ! list of processes to rescale

    integer :: iproc, id, iq, k
    real(wp) :: delta_scalable, ratio, delta_sat
    real(wp) :: th, qis

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k=1, nz
    if (l_Tcold(k)) then
   
    delta_scalable=0.0
    do iproc=1, size(iprocs_scalable)
      if (iprocs_scalable(iproc)%on) then
        id=iprocs_scalable(iproc)%id
        delta_scalable=delta_scalable+procs(i_qv, id)%column_data(k)*dt
      end if
    end do

    delta_scalable=abs(delta_scalable)

    if (delta_scalable > spacing(delta_scalable)) then
      th=qfields(k, i_th)

      qis=qisaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)

      delta_sat=abs(qis-qfields(k, i_qv))

      if (delta_scalable > delta_sat) then
        ratio=delta_sat/delta_scalable
        do iproc=1, size(iprocs_scalable)
          if (iprocs_scalable(iproc)%on) then
            id=iprocs_scalable(iproc)%id
            do iq = i_qstart, i_nstart-1
               procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
            enddo
          end if
        end do
      end if
    end if
    end if
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ensure_saturated

end module mphys_tidy
