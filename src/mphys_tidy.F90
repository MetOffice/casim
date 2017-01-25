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
       l_2mc, l_2mr, l_3mr,                       &
       l_2mi, l_2ms, l_3ms, l_2mg, l_3mg,         &
       i_an2, i_am2, i_am4, i_am5, l_warm,        &
       i_an6, i_am6, i_am7, i_am8, i_am9,         &
       l_process, ntotalq, ntotala,               &
       i_qstart, i_nstart, i_m3start,             &
       active_rain, active_cloud, isol, iinsol,   &
       l_separate_rain, l_tidy_conserve_E, l_tidy_conserve_q
  use mphys_constants, only: Lv, Ls, cp
  use qsat_funs, only: qisaturation
  use mphys_parameters, only: hydro_params
  use mphys_die, only: throw_mphys_error

#if DEF_MODEL==MODEL_KiD
  ! Kid modules
  use diagnostics, only: save_dg, i_dgtime, i_here, j_here, k_here, nx
  use runtime, only: time
#elif DEF_MODEL==MODEL_UM

  use UM_ParCore, only: mype
  use diaghelp_um, only: i_here, j_here, k_here
  use timestep_mod, only: time => timestep_number
#elif DEF_MODEL==MODEL_LEM
  use diaghelp_lem, only: i_here, j_here, k_here
  use com_params, only: time
#elif  DEF_MODEL==MODEL_MONC
  use diaghelp_monc, only: i_here, j_here, time
#endif
  implicit none
  private

  logical :: l_rescale_on_number
  logical :: l_tidym3 = .false.  ! Don't tidy based on m3 values

  real(wp), allocatable :: thresh(:), athresh(:), qin_thresh(:)
  logical,  allocatable :: l_qpos(:), l_qpresent(:), l_qsmall(:), l_qsneg(:), l_qsig(:)
  !    l_qpos: q variable is positive
  !    l_qpresent: q variable is not zero
  !    l_qsmall:   q variable is positive, but below tidy threshold
  !    l_qsneg:    q variable is small or negative
  logical,  allocatable :: l_apos(:), l_apresent(:), l_asmall(:), l_asneg(:), l_asig(:)
  logical :: l_qice, l_qliquid, current_l_negonly, current_qin_l_negonly

  public initialise_mphystidy, finalise_mphystidy, qtidy, ensure_positive, ensure_saturated, tidy_qin, &
       tidy_ain, ensure_positive_aerosol
contains

  subroutine initialise_mphystidy()
    allocate(thresh(lbound(thresh_tidy,1):ubound(thresh_tidy,1)), qin_thresh(lbound(thresh_tidy,1):ubound(thresh_tidy,1)), &
         l_qsig(0:ntotalq), l_qpos(0:ntotalq), l_qsmall(0:ntotalq), l_qsneg(0:ntotalq), l_qpresent(0:ntotalq))

    if (l_process) then
      allocate(athresh(lbound(thresh_atidy,1):ubound(thresh_atidy,1)), l_asig(ntotala), l_apos(ntotala), &
           l_asmall(ntotala), l_asneg(ntotala), l_apresent(ntotala))
    end if

    current_l_negonly=.true.
    call recompute_constants(.true.)
    current_qin_l_negonly=.true.
    call recompute_qin_constants(.true.)

    l_qsig(0)=.false.
    l_qpos(0)=.false.
    l_qsmall(0)=.false.
    l_qsneg(0)=.false.
    l_qpresent(0)=.false.
  end subroutine initialise_mphystidy

  subroutine recompute_qin_constants(l_negonly)
    logical, intent(in) :: l_negonly

    qin_thresh=thresh_tidy
    if (l_negonly) then
      qin_thresh=0.0*qin_thresh
    end if   
  end subroutine recompute_qin_constants  

  subroutine recompute_constants(l_negonly)
    logical, intent(in) :: l_negonly
   
    if (l_negonly) then
      thresh=0.0
    else
      thresh=thresh_tidy
    end if   
  end subroutine recompute_constants

  subroutine finalise_mphystidy()
    if (l_process) then
      deallocate(l_asig, l_apresent, l_asneg, l_apos, l_asmall, athresh)
    end if
    deallocate(l_qsig, l_qpresent, l_qsneg, l_qpos, l_qsmall, thresh, qin_thresh)
  end subroutine finalise_mphystidy  

  subroutine qtidy(dt, k, qfields, procs, aerofields, aeroact, dustact, aeroice, dustliq, &
       aeroprocs, i_proc, i_aproc, l_negonly)
    integer, intent(in) :: k
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

    integer :: iq

    ql_reset=.false.
    nl_reset=.false.
    qr_reset=.false.
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

    am4_reset=.false.
    am5_reset=.false.
    am7_reset=.false.
    am8_reset=.false.
    am9_reset=.false.

    an11_reset=.false. !? What to do with these???
    an12_reset=.false. !? What to do with these???

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
    end do
    
    do iq=1, ntotalq
      l_qpos(iq)=qfields(k, iq) > 0.0
    end do
    
    do iq=1,ntotalq
      l_qsmall(iq) = (.not. l_qsig(iq)) .and. l_qpos(iq)
    end do
    
    do iq=1, ntotalq
      l_qsneg(iq)=qfields(k, iq) < 0.0 .or. l_qsmall(iq)
    end do
    
    do iq=1, ntotalq
      l_qpresent(iq)=l_qpos(iq) .or. l_qsneg(iq)
    end do

    if (l_process) then      
      athresh=thresh_atidy
      if (present(l_negonly)) then
        if (l_negonly) athresh=0.0
      end if

      do iq=1, ntotala
        l_asig(iq)=aerofields(k, iq) > athresh(iq)
      end do

      do iq=1, ntotala
        l_apos(iq)=aerofields(k, iq) > 0.0
      end do

      do iq=1, ntotala
        l_asmall(iq)=(.not. l_asig(iq)) .and. l_apos(iq)
      end do

      do iq=1, ntotala
        l_asneg(iq)=aerofields(k, iq) < 0.0 .or. l_asmall(iq)
      end do

      do iq=1, ntotala
        l_apresent(iq)=l_apos(iq) .or. l_asneg(iq)
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
    if (l_3mr .and. l_tidym3) then
      m3r_reset=l_qsneg(i_m3r) .or. (l_qsig(i_m3r) .and. (qr_reset .or. nr_reset))
      nr_reset=nr_reset .or. (l_qsig(i_nr) .and. m3r_reset)
      qr_reset=qr_reset .or. (l_qsig(i_qr) .and. m3r_reset)
    end if

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
      if (l_3ms .and. l_tidym3) then
        m3s_reset=l_qsneg(i_m3s) .or. (l_qsig(i_m3s) .and. (qs_reset .or. ns_reset))
        ns_reset=ns_reset .or. (l_qsig(i_ns) .and. m3s_reset)
        qs_reset=qs_reset .or. (l_qsig(i_qs) .and. m3s_reset)
      end if

      qg_reset=l_qsneg(i_qg)
      if (l_2mg) then
        ng_reset=l_qsneg(i_ng) .or. (l_qsig(i_ng) .and. qg_reset)
        qg_reset=qg_reset .or. (l_qsig(i_qg) .and. ng_reset)
      end if
      if (l_3mg .and. l_tidym3) then
        m3g_reset=l_qsneg(i_m3g) .or. (l_qsig(i_m3g) .and. (qg_reset .or. ng_reset))
        ng_reset=ng_reset .or. (l_qsig(i_ng) .and. m3g_reset)
        qg_reset=qg_reset .or. (l_qsig(i_qg) .and. m3g_reset)
      end if
    end if

    !===========================================================
    ! Aerosol tests...
    !===========================================================
    if (l_process) then
      ! Aerosols in liquid water

      ! If small/neg values...
      if (l_asneg(i_am4))am4_reset=.true.
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

      ! If no active aerosol, then we shouldn't have any hydrometeor...
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
      procs(k,i_proc%id)%source(i_ql)=-dmass
      if (l_tidy_conserve_q)procs(k,i_proc%id)%source(i_qv)=dmass
      if (l_tidy_conserve_E)procs(k,i_proc%id)%source(i_th)=-Lv*dmass/cp
      if (l_2mc) then
        dnumber=qfields(k, i_nl)/dt
        procs(k,i_proc%id)%source(i_nl)=-dnumber
      end if
      !--------------------------------------------------
      ! aerosol - not reset, but adjusted acordingly
      !--------------------------------------------------
      if (l_process) then
        if (.not. am4_reset) then
          dmass=aeroact(k)%mact1/dt
          aeroprocs(k,i_aproc%id)%source(i_am4)=-dmass
          aeroprocs(k,i_aproc%id)%source(i_am2)=dmass
          aeroprocs(k,i_aproc%id)%source(i_an2)=dnumber
        end if
        if (.not. am9_reset) then
          dmass=dustliq(k)%mact1/dt
          if (dmass > 0.0) then
            aeroprocs(k,i_aproc%id)%source(i_am9)=-dmass
            aeroprocs(k,i_aproc%id)%source(i_am6)=dmass
            aeroprocs(k,i_aproc%id)%source(i_an6)=dnumber
          end if
        end if
      end if
    end if

    if (qr_reset .or. nr_reset .or. m3r_reset) then
      dmass=qfields(k, i_qr)/dt
      procs(k,i_proc%id)%source(i_qr)=-dmass
      if (l_tidy_conserve_q) procs(k,i_proc%id)%source(i_qv)=procs(k,i_proc%id)%source(i_qv) + dmass
      if (l_tidy_conserve_E) procs(k,i_proc%id)%source(i_th)=procs(k,i_proc%id)%source(i_th) - Lv*dmass/cp
      if (l_2mr) then
        dnumber=qfields(k, i_nr)/dt
        procs(k,i_proc%id)%source(i_nr)=-dnumber
      end if
      if (l_3mr) then
        procs(k,i_proc%id)%source(i_m3r)=-qfields(k, i_m3r)/dt
      end if
      !--------------------------------------------------
      ! aerosol
      !--------------------------------------------------
      if (l_process) then
        if (.not. am4_reset) then
          dmass=aeroact(k)%mact2/dt
          aeroprocs(k,i_aproc%id)%source(i_am2)=aeroprocs(k,i_aproc%id)%source(i_am2)+dmass
          aeroprocs(k,i_aproc%id)%source(i_an2)=aeroprocs(k,i_aproc%id)%source(i_an2)+dnumber
          aeroprocs(k,i_aproc%id)%source(i_am4)=aeroprocs(k,i_aproc%id)%source(i_am4)-dmass
        end if
        if (.not. am5_reset) then
          if (l_separate_rain) then
            dmass=aeroact(k)%mact2/dt
            aeroprocs(k,i_aproc%id)%source(i_am2)=aeroprocs(k,i_aproc%id)%source(i_am2)+dmass
            aeroprocs(k,i_aproc%id)%source(i_an2)=aeroprocs(k,i_aproc%id)%source(i_an2)+dnumber
            aeroprocs(k,i_aproc%id)%source(i_am5)=aeroprocs(k,i_aproc%id)%source(i_am5)-dmass
          end if
        end if
        if (.not. am9_reset .and. .not. l_warm) then
          dmass=dustliq(k)%mact2/dt
          if (dmass>0.0) then
            aeroprocs(k,i_aproc%id)%source(i_am9)=-dmass
            aeroprocs(k,i_aproc%id)%source(i_am6)=dmass
            aeroprocs(k,i_aproc%id)%source(i_an6)=dnumber
          end if
        end if
      end if
    end if

    if (.not. l_warm) then
      if (qi_reset .or. ni_reset) then
        dmass=qfields(k, i_qi)/dt
        procs(k,i_proc%id)%source(i_qi)=-dmass
        if (l_tidy_conserve_q) procs(k,i_proc%id)%source(i_qv)=procs(k,i_proc%id)%source(i_qv)+dmass
        if (l_tidy_conserve_E) procs(k,i_proc%id)%source(i_th)=procs(k,i_proc%id)%source(i_th)-Ls*dmass/cp
        if (l_2mi) then
          dnumber=qfields(k, i_ni)/dt
          procs(k,i_proc%id)%source(i_ni)=-dnumber
        end if
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
        if (l_process) then
          if (.not. am7_reset) then
            dmass=dustact(k)%mact1/dt
            aeroprocs(k,i_aproc%id)%source(i_am6)=aeroprocs(k,i_aproc%id)%source(i_am6)+dmass
            aeroprocs(k,i_aproc%id)%source(i_an6)=aeroprocs(k,i_aproc%id)%source(i_an6)+dnumber
            aeroprocs(k,i_aproc%id)%source(i_am7)=aeroprocs(k,i_aproc%id)%source(i_am7)-dmass
          end if
          if (.not. am8_reset) then
            dmass=aeroice(k)%mact1/dt
            aeroprocs(k,i_aproc%id)%source(i_am2)=aeroprocs(k,i_aproc%id)%source(i_am2)+dmass
            aeroprocs(k,i_aproc%id)%source(i_an2)=aeroprocs(k,i_aproc%id)%source(i_an2)+dnumber
            aeroprocs(k,i_aproc%id)%source(i_am8)=aeroprocs(k,i_aproc%id)%source(i_am8)-dmass
          end if
        end if
      end if

      if (qs_reset .or. ns_reset .or. m3s_reset) then
        dmass=qfields(k, i_qs)/dt
        procs(k,i_proc%id)%source(i_qs)=-dmass
        if (l_tidy_conserve_q) procs(k,i_proc%id)%source(i_qv)=procs(k,i_proc%id)%source(i_qv)+dmass
        if (l_tidy_conserve_E) procs(k,i_proc%id)%source(i_th)=procs(k,i_proc%id)%source(i_th)-Ls*dmass/cp
        if (l_2ms) then
          dnumber=qfields(k, i_ns)/dt
          procs(k,i_proc%id)%source(i_ns)=-dnumber
        end if
        if (l_3ms) then
          procs(k,i_proc%id)%source(i_m3s)=-qfields(k, i_m3s)/dt
        end if
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
        if (l_process) then
          if (.not. am7_reset) then
            dmass=dustact(k)%mact2/dt
            aeroprocs(k,i_aproc%id)%source(i_am6)=aeroprocs(k,i_aproc%id)%source(i_am6)+dmass
            aeroprocs(k,i_aproc%id)%source(i_an6)=aeroprocs(k,i_aproc%id)%source(i_an6)+dnumber
            aeroprocs(k,i_aproc%id)%source(i_am7)=aeroprocs(k,i_aproc%id)%source(i_am7)-dmass
          end if
          if (.not. am8_reset) then
            dmass=aeroice(k)%mact2/dt
            aeroprocs(k,i_aproc%id)%source(i_am2)=aeroprocs(k,i_aproc%id)%source(i_am2)+dmass
            aeroprocs(k,i_aproc%id)%source(i_an2)=aeroprocs(k,i_aproc%id)%source(i_an2)+dnumber
            aeroprocs(k,i_aproc%id)%source(i_am8)=aeroprocs(k,i_aproc%id)%source(i_am8)-dmass
          end if
        end if
      end if

      if (qg_reset .or. ng_reset .or. m3g_reset) then
        dmass=qfields(k, i_qg)/dt
        procs(k,i_proc%id)%source(i_qg)=-dmass
        if (l_tidy_conserve_q)procs(k,i_proc%id)%source(i_qv)=procs(k,i_proc%id)%source(i_qv)+dmass
        if (l_tidy_conserve_E)procs(k,i_proc%id)%source(i_th)=procs(k,i_proc%id)%source(i_th)-Ls*dmass/cp
        if (l_2mg) then
          dnumber=qfields(k, i_ng)/dt
          procs(k,i_proc%id)%source(i_ng)=-dnumber
        end if
        if (l_3mg) then
          procs(k,i_proc%id)%source(i_m3g)=-qfields(k, i_m3g)/dt
        end if
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
        if (l_process) then
          if (.not. am7_reset) then
            dmass=dustact(k)%mact3/dt
            aeroprocs(k,i_aproc%id)%source(i_am6)=aeroprocs(k,i_aproc%id)%source(i_am6)+dmass
            aeroprocs(k,i_aproc%id)%source(i_an6)=aeroprocs(k,i_aproc%id)%source(i_an6)+dnumber
            aeroprocs(k,i_aproc%id)%source(i_am7)=aeroprocs(k,i_aproc%id)%source(i_am7)-dmass
          end if
          if (.not. am8_reset) then
            dmass=aeroice(k)%mact3/dt
            aeroprocs(k,i_aproc%id)%source(i_am2)=aeroprocs(k,i_aproc%id)%source(i_am2)+dmass
            aeroprocs(k,i_aproc%id)%source(i_an2)=aeroprocs(k,i_aproc%id)%source(i_an2)+dnumber
            aeroprocs(k,i_aproc%id)%source(i_am8)=aeroprocs(k,i_aproc%id)%source(i_am8)-dmass
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
        aeroprocs(k,i_aproc%id)%source(i_am4)=aeroprocs(k,i_aproc%id)%source(i_am4)-dmass
        aeroprocs(k,i_aproc%id)%source(i_am2)=aeroprocs(k,i_aproc%id)%source(i_am2)+dmass
      end if

      if (am5_reset .and. l_separate_rain) then
        dmass=aerofields(k, i_am5)/dt
        aeroprocs(k,i_aproc%id)%source(i_am5)=aeroprocs(k,i_aproc%id)%source(i_am5)-dmass
        aeroprocs(k,i_aproc%id)%source(i_am2)=aeroprocs(k,i_aproc%id)%source(i_am2)+dmass
      end if

      if (.not. l_warm) then
        if (am7_reset) then
          dmass=aerofields(k, i_am7)/dt
          aeroprocs(k,i_aproc%id)%source(i_am7)=aeroprocs(k,i_aproc%id)%source(i_am7)-dmass
          aeroprocs(k,i_aproc%id)%source(i_am6)=aeroprocs(k,i_aproc%id)%source(i_am6)+dmass
        end if

        if (am8_reset) then
          dmass=aerofields(k, i_am8)/dt
          aeroprocs(k,i_aproc%id)%source(i_am8)=aeroprocs(k,i_aproc%id)%source(i_am8)-dmass
          aeroprocs(k,i_aproc%id)%source(i_am2)=aeroprocs(k,i_aproc%id)%source(i_am2)+dmass
        end if

        if (am9_reset) then
          dmass=aerofields(k, i_am9)/dt
          aeroprocs(k,i_aproc%id)%source(i_am9)=aeroprocs(k,i_aproc%id)%source(i_am9)-dmass
          aeroprocs(k,i_aproc%id)%source(i_am6)=aeroprocs(k,i_aproc%id)%source(i_am6)+dmass
        end if
      end if
    end if
  end subroutine qtidy

  subroutine tidy_qin(qfields, l_negonly)
    real(wp), intent(inout) :: qfields(:,:)
    logical, intent(in), optional :: l_negonly

    real(wp), allocatable :: thresh(:)
    logical :: ql_reset, nl_reset, qr_reset, nr_reset, m3r_reset
    logical :: qi_reset, ni_reset, qs_reset, ns_reset, m3s_reset
    logical :: qg_reset, ng_reset, m3g_reset
    integer :: k    

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
      ql_reset=.false.
      nl_reset=.false.
      qr_reset=.false.
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

      if (l_3mr .and. l_tidym3) then
        m3r_reset=qfields(k, i_m3r) < 0.0 .or. (qfields(k, i_m3r) < qin_thresh(i_m3r) .and. qfields(k, i_m3r) >0) .or. &
             (qfields(k, i_m3r) > 0.0 .and. (qfields(k, i_qr) <=0.0 .or. qfields(k, i_nr) <=0.0))
        qr_reset=qr_reset .or. (qfields(k, i_qr) > 0.0 .and. qfields(k, i_m3r) <= 0.0)
        nr_reset=nr_reset .or. (qfields(k, i_nr) > 0.0 .and. qfields(k, i_m3r) <= 0.0)
      end if

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
        if (l_3ms .and. l_tidym3) then
          m3s_reset=qfields(k, i_m3s) < 0.0 .or. (qfields(k, i_m3s) < qin_thresh(i_m3s) .and. qfields(k, i_m3s) >0) .or.&
               (qfields(k, i_m3s) > 0.0 .and. (qfields(k, i_qs) <=0.0 .or. qfields(k, i_ns) <=0.0))
          qs_reset=qs_reset .or. (qfields(k, i_qs) > 0.0 .and. qfields(k, i_m3s) <= 0.0)
          ns_reset=ns_reset .or. (qfields(k, i_ns) > 0.0 .and. qfields(k, i_m3s) <= 0.0)
        end if

        qg_reset=qfields(k, i_qg) < 0.0 .or. (qfields(k, i_qg) < qin_thresh(i_qg) .and. qfields(k, i_qg) > 0.0)
        if (l_2mg) then
          ng_reset=qfields(k, i_ng) < 0.0 .or. (qfields(k, i_ng) < qin_thresh(i_ng) .and. qfields(k, i_ng) > 0.0) .or. &
               (qfields(k, i_ng) > 0.0 .and. qfields(k, i_qg) <= 0.0)
          qg_reset=qg_reset .or. (qfields(k, i_qg) > 0.0 .and. qfields(k, i_ng) <= 0.0)
        end if
        if (l_3mg .and. l_tidym3) then
          m3g_reset=qfields(k, i_m3g) < 0.0 .or. (qfields(k, i_m3g) < qin_thresh(i_m3g) .and. qfields(k, i_m3g) >0) .or.&
               (qfields(k, i_m3g) > 0.0 .and. (qfields(k, i_qg) <=0.0 .or. qfields(k, i_ng) <=0.0))
          qg_reset=qg_reset .or. (qfields(k, i_qg) > 0.0 .and. qfields(k, i_m3g) <= 0.0)
          ng_reset=ng_reset .or. (qfields(k, i_ng) > 0.0 .and. qfields(k, i_m3g) <= 0.0)
        end if
      end if

      !==============================
      ! Now reset things...
      !==============================
      if (ql_reset .or. nl_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Lv/cp*qfields(k,i_ql)/exner(k)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_ql)
        qfields(k,i_ql)=0.0
        if (l_2mc) then
          qfields(k,i_nl)=0.0
        end if
      end if

      if (qr_reset .or. nr_reset .or. m3r_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Lv/cp*qfields(k,i_qr)/exner(k)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qr)
        qfields(k,i_qr)=0.0
        if (l_2mr) then
          qfields(k,i_nr)=0.0
        end if
        if (l_3mr) then
          qfields(k,i_m3r)=0.0
        end if
      end if

      if (qi_reset .or. ni_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qi)/exner(k)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qi)
        qfields(k,i_qi)=0.0
        if (l_2mi) then
          qfields(k,i_ni)=0.0
        end if
      end if

      if (qs_reset .or. ns_reset .or. m3s_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qs)/exner(k)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qs)
        qfields(k,i_qs)=0.0
        if (l_2ms) then
          qfields(k,i_ns)=0.0
        end if
        if (l_3ms) then
          qfields(k,i_m3s)=0.0
        end if
      end if

      if (qg_reset .or. ng_reset .or. m3g_reset) then
        if (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qg)/exner(k)
        if (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qg)
        qfields(k,i_qg)=0.0
        if (l_2mg) then
          qfields(k,i_ng)=0.0
        end if
        if (l_3mg) then
          qfields(k,i_m3g)=0.0
        end if
      end if

#if DEF_MODEL==MODEL_KiD
      if (qr_reset) call save_dg(k,1.0_wp, 'qin_reset_qr', i_dgtime)
      if (nr_reset) call save_dg(k,1.0_wp, 'qin_reset_nr', i_dgtime)
      if (m3r_reset) call save_dg(k,1.0_wp, 'qin_reset_m3r', i_dgtime)
      if (qs_reset) call save_dg(k,1.0_wp, 'qin_reset_qs', i_dgtime)
      if (ns_reset) call save_dg(k,1.0_wp, 'qin_reset_ns', i_dgtime)
      if (m3s_reset) call save_dg(k,1.0_wp, 'qin_reset_m3s', i_dgtime)
      if (qg_reset) call save_dg(k,1.0_wp, 'qin_reset_qg', i_dgtime)
      if (ng_reset) call save_dg(k,1.0_wp, 'qin_reset_ng', i_dgtime)
      if (m3g_reset) call save_dg(k,1.0_wp, 'qin_reset_m3g', i_dgtime)
#endif
    end do
  end subroutine tidy_qin

  subroutine tidy_ain(qfields, aerofields)
    real(wp), intent(in) :: qfields(:,:)
    real(wp), intent(inout) :: aerofields(:,:)

    integer :: k

    do k=1, ubound(qfields,1)
      if ((qfields(k, i_ql)+qfields(k,i_qr) <=0.0 .and. aerofields(k,i_am4)>0.0) .or. aerofields(k,i_am4) < 0.0) then
        !      print*, 'Correcting active liquid aerofields am4...', aerofields(k,i_am4), aerofields(k, i_am9)
#if DEF_MODEL==MODEL_KiD
        call save_dg(k, i_here, aerofields(k, i_am4), 'corr_am4', i_dgtime)
#endif
        aerofields(k,i_am4)=0.0
      end if
      if (i_am9 > 0)then
        if ((qfields(k, i_ql)+qfields(k,i_qr) <=0.0 .and. aerofields(k,i_am9)>0.0) .or. aerofields(k,i_am9) < 0.0) then
          !      print*, 'Correcting active liquid aerofields am9...', aerofields(k,i_am9), aerofields(k, i_am9)
#if DEF_MODEL==MODEL_KiD
          call save_dg(k, i_here, aerofields(k, i_am9), 'corr_am9', i_dgtime)
#endif
          aerofields(k,i_am9)=0.0
        end if
      end if

      if (i_am5 > 0) then
        if (((qfields(k,i_qr) <=0.0 .and. aerofields(k,i_am5)>0.0) .or. aerofields(k,i_am5) < 0.0)) then
          !      print*, 'Correcting active liquid aerofields am5...', aerofields(k,i_am5), aerofields(k, i_am5)
#if DEF_MODEL==MODEL_KiD
          call save_dg(k, i_here, aerofields(k, i_am5), 'corr_am5', i_dgtime)
#endif
          aerofields(k,i_am5)=0.0
        end if
      end if
      if (i_am7 > 0) then
        if (((qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) <=0.0 .and. aerofields(k,i_am7)>0.0) &
             .or. aerofields(k,i_am7) < 0.0)) then
          !      print*, 'Correcting active liquid aerofields am7...', aerofields(k,i_am7), aerofields(k, i_am7)
#if DEF_MODEL==MODEL_KiD
          call save_dg(k, i_here, aerofields(k, i_am7), 'corr_am7', i_dgtime)
#endif
          aerofields(k,i_am7)=0.0
        end if
      end if
      if (i_am8 > 0) then
        if (((qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) <=0.0 .and. aerofields(k,i_am8)>0.0) &
             .or. aerofields(k,i_am8) < 0.0)) then
          !     print*, 'Correcting active liquid aerofields am8...', aerofields(k,i_am8), aerofields(k, i_am8)
#if DEF_MODEL==MODEL_KiD
          call save_dg(k, i_here, aerofields(k, i_am8), 'corr_am8', i_dgtime)
#endif
          aerofields(k,i_am8)=0.0
        end if
      end if
    end do
  end subroutine tidy_ain

  ! Subroutine to ensure parallel processes don't remove more
    ! mass than is available and then rescales all processes
    ! (including number and other terms)
  subroutine ensure_positive(k, dt, qfields, procs, params, iprocs_scalable, iprocs_nonscalable &
       , aeroprocs, iprocs_dependent, iprocs_dependent_ns)
    integer, intent(in) :: k
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
    integer :: iproc, id
    real(wp) :: delta_scalable, delta_nonscalable, ratio, maxratio
    integer :: i_1m, i_2m, i_3m
    logical :: l_rescaled

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
          delta_scalable=delta_scalable+procs(k, id)%source(i_1m)
        end if
      end do
      delta_scalable=delta_scalable*dt
      if (present(iprocs_nonscalable)) then
        do iproc=1, size(iprocs_nonscalable)
          if (iprocs_nonscalable(iproc)%on) then
            id=iprocs_nonscalable(iproc)%id
            delta_nonscalable=delta_nonscalable+procs(k, id)%source(i_1m)
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
          stop
        end if
        ! Now rescale the scalable processes
        do iproc=1, size(iprocs_scalable)
          if (iprocs_scalable(iproc)%on) then
            id=iprocs_scalable(iproc)%id
            procs(k, id)%source(i_qstart:i_nstart-1)=procs(k, id)%source(i_qstart:i_nstart-1)*ratio
#if DEF_MODEL==MODEL_KiD
            call save_dg(k, i_here, ratio, 'rescale_mass_ratio_'//trim(iprocs_scalable(iproc)%name), i_dgtime)
#endif
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
              delta_scalable=delta_scalable+procs(k, id)%source(i_2m)
            end if
          end do
          delta_scalable=delta_scalable*dt
          if (abs(delta_scalable) > epsilon(delta_scalable)) then
            if (present(iprocs_nonscalable)) then
              do iproc=1, size(iprocs_nonscalable)
                if (iprocs_nonscalable(iproc)%on) then
                  id=iprocs_nonscalable(iproc)%id
                  delta_nonscalable=delta_nonscalable+procs(k, id)%source(i_2m)
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
                procs(k, id)%source(i_nstart:i_m3start-1)=procs(k, id)%source(i_nstart:i_m3start-1)*ratio
              end if
            end do
          end if
        end if
        if (params%l_3m) then ! third moment
          ! First calculate the unscaled increments
          delta_scalable=0.0
          delta_nonscalable=0.0
          do iproc=1,size(iprocs_scalable)
            if (iprocs_scalable(iproc)%on) then
              id=iprocs_scalable(iproc)%id
              delta_scalable=delta_scalable+procs(k, id)%source(i_3m)
            end if
          end do
          delta_scalable=delta_scalable*dt
          if (abs(delta_scalable) > epsilon(delta_scalable)) then
            if (present(iprocs_nonscalable)) then
              do iproc=1, size(iprocs_nonscalable)
                if (iprocs_nonscalable(iproc)%on) then
                  id=iprocs_nonscalable(iproc)%id
                  delta_nonscalable=delta_nonscalable+procs(k, id)%source(i_3m)
                end if
              end do
            end if
            delta_nonscalable=delta_nonscalable*dt
            ! ratio may now be greater than 1
            ratio=-(qfields(k, i_3m)+delta_nonscalable)/(delta_scalable + epsilon(1.0))
            ! Now rescale the scalable processes
            do iproc=1, size(iprocs_scalable)
              if (iprocs_scalable(iproc)%on) then
                id=iprocs_scalable(iproc)%id
                procs(k, id)%source(i_m3start:)=procs(k, id)%source(i_m3start:)*ratio
              end if
            end do
          end if
        end if

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
                delta_scalable=delta_scalable+procs(k, id)%source(i_2m)
              end if
            end do
            delta_scalable=delta_scalable*dt
            if (present(iprocs_nonscalable)) then
              do iproc=1, size(iprocs_nonscalable)
                if (iprocs_nonscalable(iproc)%on) then
                  id=iprocs_nonscalable(iproc)%id
                  delta_nonscalable=delta_nonscalable+procs(k, id)%source(i_2m)
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
                print*, 'WARNING: Significantly rescaled number, but not sure what to do with other moments'
                print*, 'id, ratio, bad', params%id, ratio, delta_scalable + delta_nonscalable + qfields(k, i_2m)
                print*, spacing(qfields(k, i_2m)), (qfields(k, i_2m)+delta_nonscalable), delta_scalable
                print*, 'qfields', qfields(k, i_1m), qfields(k, i_2m)
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
                    procs(k, id)%source(i_nstart:i_m3start-1)=procs(k, id)%source(i_nstart:i_m3start-1)*ratio
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

#if DEF_MODEL==MODEL_KiD
    if (l_rescaled) then
      if (nx>1) then
        call save_dg(k, i_here, 1.0_wp, 'rescaled', i_dgtime)
      else
        call save_dg(k, 1.0_wp, 'rescaled', i_dgtime)
      end if
    end if
#endif
  end subroutine ensure_positive

  ! Subroutine to ensure parallel aerosol processes don't remove more
  ! mass than is available and then rescales all processes
  ! (NB this follows any rescaling due to the parent microphysical processes and
  ! we might lose consistency between number and mass here)
  subroutine ensure_positive_aerosol(k, dt, aerofields, aerosol_procs, iprocs)
    integer, intent(in) :: k
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: aerofields(:,:)
    type(process_rate), intent(inout) :: aerosol_procs(:,:)  ! aerosol process rates
    type(process_name), intent(in) :: iprocs(:)    ! list of processes to rescale

    integer :: iq, iproc, id
    real(wp) :: ratio, delta_scalable

    do iq=1, ntotala
      delta_scalable=0.0
      do iproc=1, size(iprocs)
        if (iprocs(iproc)%on) then
          id=iprocs(iproc)%id
          delta_scalable=delta_scalable + aerosol_procs(k, id)%source(iq)
        end if
      end do
      delta_scalable=delta_scalable*dt
      if (delta_scalable + aerofields(k, iq) < spacing(aerofields(k, iq))    &
           .and. abs(delta_scalable) > spacing(aerofields(k, iq))) then
        ratio=(spacing(aerofields(k, iq))-aerofields(k, iq))/(delta_scalable)

        do iproc=1 ,size(iprocs)
          if (iprocs(iproc)%on) then
            id=iprocs(iproc)%id
            aerosol_procs(k,id)%source(iq)=aerosol_procs(k,id)%source(iq)*ratio
          end if
        end do
      end if
    end do
  end subroutine ensure_positive_aerosol

  ! Subroutine to ensure parallel ice processes don't remove more
    ! vapour than is available (i.e. so become subsaturated)
    ! and then rescales processes

    ! Modified so it can also prevent sublimation processes from putting
    ! back too much vapour and so become supersaturated
  subroutine ensure_saturated(k, dt, qfields, procs, iprocs_scalable)
    integer, intent(in) :: k
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: qfields(:,:)
    type(process_rate), intent(inout) :: procs(:,:)         ! microphysical process rates
    type(process_name), intent(in) :: iprocs_scalable(:)    ! list of processes to rescale

    integer :: iproc, id
    real(wp) :: delta_scalable, ratio, delta_sat
    real(wp) :: th, qv, qis

    delta_scalable=0.0
    do iproc=1, size(iprocs_scalable)
      if (iprocs_scalable(iproc)%on) then
        id=iprocs_scalable(iproc)%id
        delta_scalable=delta_scalable+procs(k, id)%source(i_qv)*dt
      end if
    end do

    delta_scalable=abs(delta_scalable)

    if (delta_scalable > spacing(delta_scalable)) then
      qv=qfields(k, i_qv)
      th=qfields(k, i_th)

      qis=qisaturation(th*exner(k), pressure(k)/100.0)

      delta_sat=abs(qis-qfields(k, i_qv))

      if (delta_scalable > delta_sat) then
        ratio=delta_sat/delta_scalable
        do iproc=1, size(iprocs_scalable)
          if (iprocs_scalable(iproc)%on) then
            id=iprocs_scalable(iproc)%id
            procs(k, id)%source(i_qstart:i_nstart-1)=procs(k, id)%source(i_qstart:i_nstart-1)*ratio
          end if
        end do
      end if
    end if
  end subroutine ensure_saturated

  ! Once tidied, there may be some negative numbers due to
  ! rounding error, so should be safe to lose these
  subroutine cleanup_q(qfields)
    real(wp), intent(inout) :: qfields(:,:)

    integer :: i1,i2 ! use loops instead of where construct
    integer :: lbounds(2), ubounds(2)

    lbounds=lbound(qfields)
    ubounds=ubound(qfields)
    do i2=lbounds(2), ubounds(2)
      do i1=lbounds(1), ubounds(1)
        if (qfields(i1,i2) < 0.0) then
          qfields(i1,i2)=0.0
        end if
      end do
    end do
  end subroutine cleanup_q
end module mphys_tidy
