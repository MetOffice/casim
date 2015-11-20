MODULE mphys_tidy

USE variable_precision, ONLY: wp
USE process_routines, ONLY: process_rate,  process_name
USE aerosol_routines, ONLY: aerosol_active
USE thresholds, ONLY: thresh_tidy, thresh_atidy
USE passive_fields, ONLY: exner, pressure
USE mphys_switches, ONLY:                       &
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
USE mphys_constants, ONLY: Lv, Ls, cp
USE qsat_funs, ONLY: qisaturation
USE mphys_parameters, ONLY: hydro_params
USE mphys_die, ONLY: throw_mphys_error

#if DEF_MODEL==MODEL_KiD
  ! Kid modules
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, j_here, k_here, nx
USE runtime, ONLY: time
#elif DEF_MODEL==MODEL_UM

USE UM_ParCore, ONLY: mype
USE diaghelp_um, ONLY: i_here, j_here, k_here
USE timestep_mod, ONLY: time => timestep_number
#elif DEF_MODEL==MODEL_LEM
USE diaghelp_lem, ONLY: i_here, j_here, k_here
USE com_params, ONLY: time
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here, time
#endif
IMPLICIT NONE

LOGICAL :: l_rescale_on_number
LOGICAL :: l_tidym3 = .FALSE.  ! Don't tidy based on m3 values

CONTAINS

SUBROUTINE qtidy(dt, k, qfields, procs, aerofields, aeroact, dustact, aeroice, dustliq, &
     aeroprocs, i_proc, i_aproc, l_negonly)

INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN) :: dt
REAL(wp), INTENT(IN) :: qfields(:,:), aerofields(:,:)
TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustact(:), aeroice(:), dustliq(:)
TYPE(process_rate), INTENT(INOUT) :: procs(:,:)
TYPE(process_rate), INTENT(INOUT) :: aeroprocs(:,:)
TYPE(process_name), INTENT(IN) :: i_proc, i_aproc
LOGICAL, INTENT(IN), OPTIONAL :: l_negonly

LOGICAL :: ql_reset, nl_reset, qr_reset, nr_reset, m3r_reset
LOGICAL :: qi_reset, ni_reset, qs_reset, ns_reset, m3s_reset
LOGICAL :: qg_reset, ng_reset, m3g_reset
LOGICAL :: am4_reset, am5_reset, am7_reset, am8_reset, am9_reset
LOGICAL :: an11_reset, an12_reset

REAL(wp) :: dmass, dnumber
REAL(wp), ALLOCATABLE :: thresh(:), athresh(:)
LOGICAL,  ALLOCATABLE :: l_qpos(:), l_qpresent(:), l_qsmall(:), l_qsneg(:), l_qsig(:)
    !    l_qpos: q variable is positive
    !    l_qpresent: q variable is not zero
    !    l_qsmall:   q variable is positive, but below tidy threshold
    !    l_qsneg:    q variable is small or negative
LOGICAL,  ALLOCATABLE :: l_apos(:), l_apresent(:), l_asmall(:), l_asneg(:), l_asig(:)

LOGICAL :: l_qice, l_qliquid

INTEGER :: iq

ql_reset=.FALSE.
nl_reset=.FALSE.
qr_reset=.FALSE.
nr_reset=.FALSE.
m3r_reset=.FALSE.
qi_reset=.FALSE.
ni_reset=.FALSE.
qs_reset=.FALSE.
ns_reset=.FALSE.
m3s_reset=.FALSE.
qg_reset=.FALSE.
ng_reset=.FALSE.
m3g_reset=.FALSE.

am4_reset=.FALSE.
am5_reset=.FALSE.
am7_reset=.FALSE.
am8_reset=.FALSE.
am9_reset=.FALSE.

an11_reset=.FALSE. !? What to do with these???
an12_reset=.FALSE. !? What to do with these???

ALLOCATE(thresh(LBOUND(thresh_tidy,1):UBOUND(thresh_tidy,1)))
thresh=thresh_tidy
IF (PRESENT(l_negonly)) THEN
  IF (l_negonly) thresh=0.0
END IF

ALLOCATE(l_qsig(0:ntotalq))
l_qsig(0)=.FALSE.
DO iq=1,ntotalq
  l_qsig(iq) = qfields(k, iq) > thresh(iq)
END DO

ALLOCATE(l_qpos(0:ntotalq))
l_qpos(0)=.FALSE.
DO iq=1,ntotalq
  l_qpos(iq) = qfields(k, iq) > 0.0
END DO

ALLOCATE(l_qsmall(0:ntotalq))
l_qsmall(0)=.FALSE.
DO iq=1,ntotalq
  l_qsmall(iq) = (.NOT. l_qsig(iq)) .AND. l_qpos(iq)
END DO

ALLOCATE(l_qsneg(0:ntotalq))
l_qsneg(0)=.FALSE.
DO iq=1,ntotalq
  l_qsneg(iq) = qfields(k, iq) < 0.0 .OR. l_qsmall(iq)
END DO

ALLOCATE(l_qpresent(0:ntotalq))
l_qpresent(0)=.FALSE.
DO iq=1,ntotalq
  l_qpresent(iq) = l_qpos(iq) .OR. l_qsneg(iq)
END DO

IF (l_process) THEN

  ALLOCATE(athresh(LBOUND(thresh_atidy,1):UBOUND(thresh_atidy,1)))
  athresh=thresh_atidy
  IF (PRESENT(l_negonly)) THEN
    IF (l_negonly) athresh=0.0
  END IF

  ALLOCATE(l_asig(ntotala))
  DO iq=1,ntotala
    l_asig(iq) = aerofields(k, iq) > athresh(iq)
  END DO

  ALLOCATE(l_apos(ntotala))
  DO iq=1,ntotala
    l_apos(iq) = aerofields(k, iq) > 0.0
  END DO

  ALLOCATE(l_asmall(ntotala))
  DO iq=1,ntotala
    l_asmall(iq) = (.NOT. l_asig(iq)) .AND. l_apos(iq)
  END DO

  ALLOCATE(l_asneg(ntotala))
  DO iq=1,ntotala
    l_asneg(iq) = aerofields(k, iq) < 0.0 .OR. l_asmall(iq)
  END DO

  ALLOCATE(l_apresent(ntotala))
  DO iq=1,ntotala
    l_apresent(iq) = l_apos(iq) .OR. l_asneg(iq)
  END DO

END IF

l_qliquid = l_qsig(i_ql) .OR. l_qsig(i_qr)
l_qice    = l_qsig(i_qi) .OR. l_qsig(i_qs) .OR. l_qsig(i_qg)

    ! Tidying of small and negative numbers and/or incompatible numbers (e.g.nl>0 and ql=0)
    ! - Mass and energy conserving...
    !==============================
    ! What should be reset?
    !==============================
ql_reset=l_qsneg(i_ql)
IF (l_2mc) THEN
  nl_reset=l_qsneg(i_nl) .OR. (l_qsig(i_nl) .AND. ql_reset)
  ql_reset=ql_reset .OR. (l_qsig(i_ql) .AND. nl_reset)
END IF
qr_reset=l_qsneg(i_qr)
IF (l_2mr) THEN
  nr_reset=l_qsneg(i_nr) .OR. (l_qsig(i_nr) .AND. qr_reset)
  qr_reset=qr_reset .OR. (l_qsig(i_qr) .AND. nr_reset)
END IF
IF (l_3mr .AND. l_tidym3) THEN
  m3r_reset=l_qsneg(i_m3r) .OR. (l_qsig(i_m3r) .AND. (qr_reset .OR. nr_reset))
  nr_reset=nr_reset .OR. (l_qsig(i_nr) .AND. m3r_reset)
  qr_reset=qr_reset .OR. (l_qsig(i_qr) .AND. m3r_reset)
END IF

IF (.NOT. l_warm) THEN
  qi_reset=l_qsneg(i_qi)
  IF (l_2mi) THEN
    ni_reset=l_qsneg(i_ni) .OR. (l_qsig(i_ni) .AND. qi_reset)
    qi_reset=qi_reset .OR. (l_qsig(i_qi) .AND. ni_reset)
  END IF

  qs_reset=l_qsneg(i_qs)
  IF (l_2ms) THEN
    ns_reset=l_qsneg(i_ns) .OR. (l_qsig(i_ns) .AND. qs_reset)
    qs_reset=qs_reset .OR. (l_qsig(i_qs) .AND. ns_reset)
  END IF
  IF (l_3ms .AND. l_tidym3) THEN
    m3s_reset=l_qsneg(i_m3s) .OR. (l_qsig(i_m3s) .AND. (qs_reset .OR. ns_reset))
    ns_reset=ns_reset .OR. (l_qsig(i_ns) .AND. m3s_reset)
    qs_reset=qs_reset .OR. (l_qsig(i_qs) .AND. m3s_reset)
  END IF

  qg_reset=l_qsneg(i_qg)
  IF (l_2mg) THEN
    ng_reset=l_qsneg(i_ng) .OR. (l_qsig(i_ng) .AND. qg_reset)
    qg_reset=qg_reset .OR. (l_qsig(i_qg) .AND. ng_reset)
  END IF
  IF (l_3mg .AND. l_tidym3) THEN
    m3g_reset=l_qsneg(i_m3g) .OR. (l_qsig(i_m3g) .AND. (qg_reset .OR. ng_reset))
    ng_reset=ng_reset .OR. (l_qsig(i_ng) .AND. m3g_reset)
    qg_reset=qg_reset .OR. (l_qsig(i_qg) .AND. m3g_reset)
  END IF

END IF

    !===========================================================
    ! Aerosol tests...
    !===========================================================
IF (l_process) THEN
      ! Aerosols in liquid water

      ! If small/neg values...
  IF (l_asneg(i_am4))am4_reset=.TRUE.
  IF (l_separate_rain) THEN
    IF (l_asneg(i_am5))am5_reset=.TRUE.
  END IF
  if (.not. l_Warm)then
     IF (l_asneg(i_am9))am9_reset=.TRUE.
  end if
      ! If no hydrometeors...
  IF ((ql_reset .AND. qr_reset) .OR. .NOT. l_qliquid) THEN
    IF (l_asig(i_am4))am4_reset=.TRUE.
    if (.not.l_warm)then
       IF (l_asig(i_am9))am9_reset=.TRUE.
    end if
    IF (l_separate_rain) THEN
      IF (l_asig(i_am5)) am5_reset=.TRUE.
    END IF
  END IF

      ! If no active aerosol, then we shouldn't have any hydrometeor...
  ql_reset=ql_reset .OR. (am4_reset .AND. am9_reset .AND. l_qsig(i_ql))
  qr_reset=qr_reset .OR. (am4_reset .AND. am9_reset .AND. l_qsig(i_qr))
  qr_reset=qr_reset .OR. (am5_reset .AND. l_qsig(i_qr))

      ! Aerosols in ice
      ! If small/neg values...
  if (.not. l_Warm)then
     IF (l_asneg(i_am7))am7_reset=.TRUE.
     IF (l_asneg(i_am8)) THEN
        am8_reset=.TRUE.
     END IF
      ! If no hydrometeors...
  IF ((qi_reset .AND. qs_reset .AND. qg_reset) .OR. .NOT. l_qice) THEN
    IF (l_asig(i_am7))am7_reset=.TRUE.
    IF (l_asig(i_am8)) THEN
      am8_reset=.TRUE.
    END IF
  END IF
  end if

      ! If no active aerosol, then we shouldn't have any hydrometeor...
  qi_reset=qi_reset .OR. (am7_reset .AND. am8_reset .AND. l_qsig(i_qi))
  qs_reset=qs_reset .OR. (am7_reset .AND. am8_reset .AND. l_qsig(i_qs))
  qg_reset=qg_reset .OR. (am7_reset .AND. am8_reset .AND. l_qsig(i_qg))
END IF

    !===========================================================
    ! Consistency following aerosol
    !===========================================================
nl_reset=ql_reset .AND. l_qsig(i_nl)
nr_reset=nr_reset .AND. l_qsig(i_nr)
m3r_reset=m3r_reset .AND. l_qsig(i_m3r)
ni_reset=qi_reset .AND. l_qsig(i_ni)
ns_reset=ns_reset .AND. l_qsig(i_ns)
m3s_reset=m3s_reset .AND. l_qsig(i_m3s)
ng_reset=ng_reset .AND. l_qsig(i_ng)
m3g_reset=m3g_reset .AND. l_qsig(i_m3g)

    !==============================
    ! Now reset things...
    !==============================
IF (ql_reset .OR. nl_reset) THEN
  dmass=qfields(k, i_ql)/dt
  procs(k,i_proc%id)%source(i_ql)=-dmass
  if (l_tidy_conserve_q)procs(k,i_proc%id)%source(i_qv)=dmass
  if (l_tidy_conserve_E)procs(k,i_proc%id)%source(i_th)= -Lv*dmass/cp
  IF (l_2mc) THEN
    dnumber=qfields(k, i_nl)/dt
    procs(k,i_proc%id)%source(i_nl)=-dnumber
  END IF
      !--------------------------------------------------
      ! aerosol - not reset, but adjusted acordingly
      !--------------------------------------------------
  IF (l_process) THEN
    IF (.NOT. am4_reset) THEN
      dmass=aeroact(k)%mact1/dt
      aeroprocs(k,i_aproc%id)%source(i_am4)=-dmass
      aeroprocs(k,i_aproc%id)%source(i_am2)=dmass
      aeroprocs(k,i_aproc%id)%source(i_an2)=dnumber
    END IF
    IF (.NOT. am9_reset) THEN
      dmass=dustliq(k)%mact1/dt
      IF (dmass > 0.0) THEN
        aeroprocs(k,i_aproc%id)%source(i_am9)=-dmass
        aeroprocs(k,i_aproc%id)%source(i_am6)=dmass
        aeroprocs(k,i_aproc%id)%source(i_an6)=dnumber
      END IF
    END IF
  END IF
END IF

IF (qr_reset .OR. nr_reset .OR. m3r_reset) THEN
  dmass=qfields(k, i_qr)/dt
  procs(k,i_proc%id)%source(i_qr)=-dmass
  if (l_tidy_conserve_q)procs(k,i_proc%id)%source(i_qv)=     &
         procs(k,i_proc%id)%source(i_qv) + dmass
  if (l_tidy_conserve_E)procs(k,i_proc%id)%source(i_th)=     &
         procs(k,i_proc%id)%source(i_th) - Lv*dmass/cp
  IF (l_2mr) THEN
    dnumber=qfields(k, i_nr)/dt
    procs(k,i_proc%id)%source(i_nr)=-dnumber
  END IF
  IF (l_3mr) THEN
    procs(k,i_proc%id)%source(i_m3r)=-qfields(k, i_m3r)/dt
  END IF
      !--------------------------------------------------
      ! aerosol
      !--------------------------------------------------
  IF (l_process) THEN
    IF (.NOT. am4_reset) THEN
      dmass=aeroact(k)%mact2/dt
      aeroprocs(k,i_aproc%id)%source(i_am2)=     &
             aeroprocs(k,i_aproc%id)%source(i_am2) + dmass
      aeroprocs(k,i_aproc%id)%source(i_an2)=     &
             aeroprocs(k,i_aproc%id)%source(i_an2) + dnumber
      aeroprocs(k,i_aproc%id)%source(i_am4)=    &
             aeroprocs(k,i_aproc%id)%source(i_am4) - dmass
    END IF
    IF (.NOT. am5_reset) THEN
      IF (l_separate_rain) THEN
        dmass=aeroact(k)%mact2/dt
        aeroprocs(k,i_aproc%id)%source(i_am2)=     &
               aeroprocs(k,i_aproc%id)%source(i_am2) + dmass
        aeroprocs(k,i_aproc%id)%source(i_an2)=     &
               aeroprocs(k,i_aproc%id)%source(i_an2) + dnumber
        aeroprocs(k,i_aproc%id)%source(i_am5)=    &
               aeroprocs(k,i_aproc%id)%source(i_am5) - dmass
      END IF
    END IF
    IF (.NOT. am9_reset .AND. .NOT. l_warm) THEN
      dmass=dustliq(k)%mact2/dt
      IF (dmass>0.0) THEN
        aeroprocs(k,i_aproc%id)%source(i_am9)=-dmass
        aeroprocs(k,i_aproc%id)%source(i_am6)=dmass
        aeroprocs(k,i_aproc%id)%source(i_an6)=dnumber
      END IF
    END IF

  END IF

END IF

IF (.NOT. l_warm) THEN
  IF (qi_reset .OR. ni_reset) THEN
    dmass=qfields(k, i_qi)/dt
    procs(k,i_proc%id)%source(i_qi)=-dmass
    if (l_tidy_conserve_q)procs(k,i_proc%id)%source(i_qv)=     &
           procs(k,i_proc%id)%source(i_qv) + dmass
    if (l_tidy_conserve_E)procs(k,i_proc%id)%source(i_th)=     &
           procs(k,i_proc%id)%source(i_th) - Ls*dmass/cp
    IF (l_2mi) THEN
      dnumber=qfields(k, i_ni)/dt
      procs(k,i_proc%id)%source(i_ni)=-dnumber
    END IF
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
    IF (l_process) THEN
      IF (.NOT. am7_reset) THEN
        dmass=dustact(k)%mact1/dt
        aeroprocs(k,i_aproc%id)%source(i_am6)=     &
               aeroprocs(k,i_aproc%id)%source(i_am6) + dmass
        aeroprocs(k,i_aproc%id)%source(i_an6)=     &
               aeroprocs(k,i_aproc%id)%source(i_an6) + dnumber
        aeroprocs(k,i_aproc%id)%source(i_am7)=    &
               aeroprocs(k,i_aproc%id)%source(i_am7) - dmass
      END IF
      IF (.NOT. am8_reset) THEN
        dmass=aeroice(k)%mact1/dt
        aeroprocs(k,i_aproc%id)%source(i_am2)=     &
               aeroprocs(k,i_aproc%id)%source(i_am2) + dmass
        aeroprocs(k,i_aproc%id)%source(i_an2)=     &
               aeroprocs(k,i_aproc%id)%source(i_an2) + dnumber
        aeroprocs(k,i_aproc%id)%source(i_am8)=    &
               aeroprocs(k,i_aproc%id)%source(i_am8) - dmass
      END IF
    END IF
  END IF

  IF (qs_reset .OR. ns_reset .OR. m3s_reset) THEN
    dmass=qfields(k, i_qs)/dt
    procs(k,i_proc%id)%source(i_qs)=-dmass
    if (l_tidy_conserve_q)procs(k,i_proc%id)%source(i_qv)=     &
           procs(k,i_proc%id)%source(i_qv) + dmass
    if (l_tidy_conserve_E)procs(k,i_proc%id)%source(i_th)=     &
           procs(k,i_proc%id)%source(i_th) - Ls*dmass/cp
    IF (l_2ms) THEN
      dnumber=qfields(k, i_ns)/dt
      procs(k,i_proc%id)%source(i_ns)=-dnumber
    END IF
    IF (l_3ms) THEN
      procs(k,i_proc%id)%source(i_m3s)=-qfields(k, i_m3s)/dt
    END IF
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
    IF (l_process) THEN
      IF (.NOT. am7_reset) THEN
        dmass=dustact(k)%mact2/dt
        aeroprocs(k,i_aproc%id)%source(i_am6)=     &
               aeroprocs(k,i_aproc%id)%source(i_am6) + dmass
        aeroprocs(k,i_aproc%id)%source(i_an6)=     &
               aeroprocs(k,i_aproc%id)%source(i_an6) + dnumber
        aeroprocs(k,i_aproc%id)%source(i_am7)=    &
               aeroprocs(k,i_aproc%id)%source(i_am7) - dmass
      END IF
      IF (.NOT. am8_reset) THEN
        dmass=aeroice(k)%mact2/dt
        aeroprocs(k,i_aproc%id)%source(i_am2)=     &
               aeroprocs(k,i_aproc%id)%source(i_am2) + dmass
        aeroprocs(k,i_aproc%id)%source(i_an2)=     &
               aeroprocs(k,i_aproc%id)%source(i_an2) + dnumber
        aeroprocs(k,i_aproc%id)%source(i_am8)=    &
               aeroprocs(k,i_aproc%id)%source(i_am8) - dmass
      END IF
    END IF
  END IF


  IF (qg_reset .OR. ng_reset .OR. m3g_reset) THEN
    dmass=qfields(k, i_qg)/dt
    procs(k,i_proc%id)%source(i_qg)=-dmass
    if (l_tidy_conserve_q)procs(k,i_proc%id)%source(i_qv)=     &
           procs(k,i_proc%id)%source(i_qv) + dmass
    if (l_tidy_conserve_E)procs(k,i_proc%id)%source(i_th)=     &
           procs(k,i_proc%id)%source(i_th) - Ls*dmass/cp
    IF (l_2mg) THEN
      dnumber=qfields(k, i_ng)/dt
      procs(k,i_proc%id)%source(i_ng)=-dnumber
    END IF
    IF (l_3mg) THEN
      procs(k,i_proc%id)%source(i_m3g)=-qfields(k, i_m3g)/dt
    END IF
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
    IF (l_process) THEN
      IF (.NOT. am7_reset) THEN
        dmass=dustact(k)%mact3/dt
        aeroprocs(k,i_aproc%id)%source(i_am6)=     &
               aeroprocs(k,i_aproc%id)%source(i_am6) + dmass
        aeroprocs(k,i_aproc%id)%source(i_an6)=     &
               aeroprocs(k,i_aproc%id)%source(i_an6) + dnumber
        aeroprocs(k,i_aproc%id)%source(i_am7)=    &
               aeroprocs(k,i_aproc%id)%source(i_am7) - dmass
      END IF
      IF (.NOT. am8_reset) THEN
        dmass=aeroice(k)%mact3/dt
        aeroprocs(k,i_aproc%id)%source(i_am2)=     &
               aeroprocs(k,i_aproc%id)%source(i_am2) + dmass
        aeroprocs(k,i_aproc%id)%source(i_an2)=     &
               aeroprocs(k,i_aproc%id)%source(i_an2) + dnumber
        aeroprocs(k,i_aproc%id)%source(i_am8)=    &
               aeroprocs(k,i_aproc%id)%source(i_am8) - dmass
      END IF
    END IF
  END IF

END IF
    !==============================
    ! Now reset aerosol...
    !==============================

IF (l_process) THEN
  IF (am4_reset) THEN
    dmass=aerofields(k, i_am4)/dt
    aeroprocs(k,i_aproc%id)%source(i_am4)=     &
           aeroprocs(k,i_aproc%id)%source(i_am4) - dmass
    aeroprocs(k,i_aproc%id)%source(i_am2)=     &
           aeroprocs(k,i_aproc%id)%source(i_am2) + dmass
  END IF

  IF (am5_reset .AND. l_separate_rain) THEN
    dmass=aerofields(k, i_am5)/dt
    aeroprocs(k,i_aproc%id)%source(i_am5)=     &
           aeroprocs(k,i_aproc%id)%source(i_am5) - dmass
    aeroprocs(k,i_aproc%id)%source(i_am2)=     &
           aeroprocs(k,i_aproc%id)%source(i_am2) + dmass
  END IF

  IF (.NOT. l_warm) THEN
    IF (am7_reset) THEN
      dmass=aerofields(k, i_am7)/dt
      aeroprocs(k,i_aproc%id)%source(i_am7)=     &
             aeroprocs(k,i_aproc%id)%source(i_am7) - dmass
      aeroprocs(k,i_aproc%id)%source(i_am6)=     &
             aeroprocs(k,i_aproc%id)%source(i_am6) + dmass
    END IF

    IF (am8_reset) THEN
      dmass=aerofields(k, i_am8)/dt
      aeroprocs(k,i_aproc%id)%source(i_am8)=     &
             aeroprocs(k,i_aproc%id)%source(i_am8) - dmass
      aeroprocs(k,i_aproc%id)%source(i_am2)=     &
             aeroprocs(k,i_aproc%id)%source(i_am2) + dmass
    END IF

    IF (am9_reset) THEN
      dmass=aerofields(k, i_am9)/dt
      aeroprocs(k,i_aproc%id)%source(i_am9)=     &
             aeroprocs(k,i_aproc%id)%source(i_am9) - dmass
      aeroprocs(k,i_aproc%id)%source(i_am6)=     &
             aeroprocs(k,i_aproc%id)%source(i_am6) + dmass
    END IF
  END IF

END IF

IF (l_process) THEN
  DEALLOCATE(l_asig)
  DEALLOCATE(l_apresent)
  DEALLOCATE(l_asneg)
  DEALLOCATE(l_apos)
  DEALLOCATE(l_asmall)
  DEALLOCATE(athresh)
END IF
DEALLOCATE(l_qsig)
DEALLOCATE(l_qpresent)
DEALLOCATE(l_qsneg)
DEALLOCATE(l_qpos)
DEALLOCATE(l_qsmall)

DEALLOCATE(thresh)

END SUBROUTINE qtidy

SUBROUTINE tidy_qin(qfields, l_negonly)

REAL(wp), INTENT(INOUT) :: qfields(:,:)
LOGICAL, INTENT(IN), OPTIONAL :: l_negonly

REAL(wp), ALLOCATABLE :: thresh(:)

LOGICAL :: ql_reset, nl_reset, qr_reset, nr_reset, m3r_reset
LOGICAL :: qi_reset, ni_reset, qs_reset, ns_reset, m3s_reset
LOGICAL :: qg_reset, ng_reset, m3g_reset

INTEGER :: k


ALLOCATE(thresh(LBOUND(thresh_tidy,1):UBOUND(thresh_tidy,1)))
thresh=thresh_tidy
IF (PRESENT(l_negonly)) THEN
  IF (l_negonly)thresh=0.0*thresh
END IF

DO k=1,UBOUND(qfields,1)

  ql_reset=.FALSE.
  nl_reset=.FALSE.
  qr_reset=.FALSE.
  nr_reset=.FALSE.
  m3r_reset=.FALSE.
  qi_reset=.FALSE.
  ni_reset=.FALSE.
  qs_reset=.FALSE.
  ns_reset=.FALSE.
  m3s_reset=.FALSE.
  qg_reset=.FALSE.
  ng_reset=.FALSE.
  m3g_reset=.FALSE.

  ql_reset=qfields(k, i_ql) < 0.0 .OR.  &
       (qfields(k, i_ql) < thresh(i_ql) .AND. qfields(k, i_ql) >0)
  IF (l_2mc) THEN
    nl_reset=qfields(k, i_nl) < 0.0 .OR.  &
         (qfields(k, i_nl) < thresh(i_nl) .AND. qfields(k, i_nl) >0) .OR. &
         (qfields(k, i_nl) > 0.0 .AND. qfields(k, i_ql) <= 0.0)
    ql_reset=ql_reset .OR. (qfields(k, i_ql) > 0.0 .AND. qfields(k, i_nl) <= 0.0)
  END IF
  qr_reset=qfields(k, i_qr) < 0.0 .OR.  &
       (qfields(k, i_qr) < thresh(i_qr) .AND. qfields(k, i_qr) > 0.0)

  IF (l_2mr) THEN
    nr_reset=qfields(k, i_nr) < 0.0 .OR.  &
         (qfields(k, i_nr) < thresh(i_nr) .AND. qfields(k, i_nr) > 0.0) .OR. &
         (qfields(k, i_nr) > 0.0 .AND. qfields(k, i_qr) <= 0.0)
    qr_reset=qr_reset .OR. (qfields(k, i_qr) > 0.0 .AND. qfields(k, i_nr) <= 0.0)
  END IF

  IF (l_3mr .AND. l_tidym3) THEN
    m3r_reset=qfields(k, i_m3r) < 0.0 .OR.  &
         (qfields(k, i_m3r) < thresh(i_m3r) .AND. qfields(k, i_m3r) >0) .OR.&
         (qfields(k, i_m3r) > 0.0 .AND. (qfields(k, i_qr) <=0.0 .OR. qfields(k, i_nr) <=0.0))
    qr_reset=qr_reset .OR. (qfields(k, i_qr) > 0.0 .AND. qfields(k, i_m3r) <= 0.0)
    nr_reset=nr_reset .OR. (qfields(k, i_nr) > 0.0 .AND. qfields(k, i_m3r) <= 0.0)
  END IF
  
  IF (.NOT. l_warm) THEN
    qi_reset=qfields(k, i_qi) < 0.0 .OR.  &
         (qfields(k, i_qi) < thresh(i_qi) .AND. qfields(k, i_qi) > 0.0)
    IF (l_2mi) THEN
      ni_reset=qfields(k, i_ni) < 0.0 .OR.  &
           (qfields(k, i_ni) < thresh(i_ni) .AND. qfields(k, i_ni) > 0.0) .OR. &
           (qfields(k, i_ni) > 0.0 .AND. qfields(k, i_qi) <= 0.0)
      qi_reset=qi_reset .OR. (qfields(k, i_qi) > 0.0 .AND. qfields(k, i_ni) <= 0.0)
    END IF

    qs_reset=qfields(k, i_qs) < 0.0 .OR.  &
         (qfields(k, i_qs) < thresh(i_qs) .AND. qfields(k, i_qs) > 0.0)
    IF (l_2ms) THEN
      ns_reset=qfields(k, i_ns) < 0.0 .OR.  &
           (qfields(k, i_ns) < thresh(i_ns) .AND. qfields(k, i_ns) > 0.0) .OR. &
           (qfields(k, i_ns) > 0.0 .AND. qfields(k, i_qs) <= 0.0)
      qs_reset=qs_reset .OR. (qfields(k, i_qs) > 0.0 .AND. qfields(k, i_ns) <= 0.0)
    END IF
    IF (l_3ms .AND. l_tidym3) THEN
      m3s_reset=qfields(k, i_m3s) < 0.0 .OR.  &
           (qfields(k, i_m3s) < thresh(i_m3s) .AND. qfields(k, i_m3s) >0) .OR.&
           (qfields(k, i_m3s) > 0.0 .AND. (qfields(k, i_qs) <=0.0 .OR. qfields(k, i_ns) <=0.0))
      qs_reset=qs_reset .OR. (qfields(k, i_qs) > 0.0 .AND. qfields(k, i_m3s) <= 0.0)
      ns_reset=ns_reset .OR. (qfields(k, i_ns) > 0.0 .AND. qfields(k, i_m3s) <= 0.0)
    END IF


    qg_reset=qfields(k, i_qg) < 0.0 .OR.  &
         (qfields(k, i_qg) < thresh(i_qg) .AND. qfields(k, i_qg) > 0.0)
    IF (l_2mg) THEN
      ng_reset=qfields(k, i_ng) < 0.0 .OR.  &
           (qfields(k, i_ng) < thresh(i_ng) .AND. qfields(k, i_ng) > 0.0) .OR. &
           (qfields(k, i_ng) > 0.0 .AND. qfields(k, i_qg) <= 0.0)
      qg_reset=qg_reset .OR. (qfields(k, i_qg) > 0.0 .AND. qfields(k, i_ng) <= 0.0)
    END IF
    IF (l_3mg .AND. l_tidym3) THEN
      m3g_reset=qfields(k, i_m3g) < 0.0 .OR.  &
           (qfields(k, i_m3g) < thresh(i_m3g) .AND. qfields(k, i_m3g) >0) .OR.&
           (qfields(k, i_m3g) > 0.0 .AND. (qfields(k, i_qg) <=0.0 .OR. qfields(k, i_ng) <=0.0))
      qg_reset=qg_reset .OR. (qfields(k, i_qg) > 0.0 .AND. qfields(k, i_m3g) <= 0.0)
      ng_reset=ng_reset .OR. (qfields(k, i_ng) > 0.0 .AND. qfields(k, i_m3g) <= 0.0)
    END IF
  END IF


    !==============================
    ! Now reset things...
    !==============================
  IF (ql_reset .OR. nl_reset) THEN
    if (l_tidy_conserve_E)qfields(k,i_th)=qfields(k,i_th)-Lv/cp*qfields(k,i_ql)/exner(k)
    if (l_tidy_conserve_q)qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_ql)
    qfields(k,i_ql)=0.0
    IF (l_2mc) THEN
      qfields(k,i_nl)=0.0
    END IF
  END IF


  IF (qr_reset .OR. nr_reset .OR. m3r_reset) THEN
    if (l_tidy_conserve_E)qfields(k,i_th)=qfields(k,i_th)-Lv/cp*qfields(k,i_qr)/exner(k)
    if (l_tidy_conserve_q)qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qr)
    qfields(k,i_qr)=0.0
    IF (l_2mr) THEN
      qfields(k,i_nr)=0.0
    END IF
    IF (l_3mr) THEN
      qfields(k,i_m3r)=0.0
    END IF
  END IF

  IF (qi_reset .OR. ni_reset) THEN
    if (l_tidy_conserve_E)qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qi)/exner(k)
    if (l_tidy_conserve_q)qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qi)
    qfields(k,i_qi)=0.0
    IF (l_2mi) THEN
      qfields(k,i_ni)=0.0
    END IF
  END IF

  IF (qs_reset .OR. ns_reset .OR. m3s_reset) THEN
    if (l_tidy_conserve_E)qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qs)/exner(k)
    if (l_tidy_conserve_q)qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qs)
    qfields(k,i_qs)=0.0
    IF (l_2ms) THEN
      qfields(k,i_ns)=0.0
    END IF
    IF (l_3ms) THEN
      qfields(k,i_m3s)=0.0
    END IF
  END IF

  IF (qg_reset .OR. ng_reset .OR. m3g_reset) THEN
    if (l_tidy_conserve_E)qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qg)/exner(k)
    if (l_tidy_conserve_q)qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qg)
    qfields(k,i_qg)=0.0
    IF (l_2mg) THEN
      qfields(k,i_ng)=0.0
    END IF
    IF (l_3mg) THEN
      qfields(k,i_m3g)=0.0
    END IF
  END IF

#if DEF_MODEL==MODEL_KiD
  IF (qr_reset) CALL save_dg(k,1.0_wp, 'qin_reset_qr', i_dgtime)
  IF (nr_reset) CALL save_dg(k,1.0_wp, 'qin_reset_nr', i_dgtime)
  IF (m3r_reset) CALL save_dg(k,1.0_wp, 'qin_reset_m3r', i_dgtime)
  IF (qs_reset) CALL save_dg(k,1.0_wp, 'qin_reset_qs', i_dgtime)
  IF (ns_reset) CALL save_dg(k,1.0_wp, 'qin_reset_ns', i_dgtime)
  IF (m3s_reset) CALL save_dg(k,1.0_wp, 'qin_reset_m3s', i_dgtime)
  IF (qg_reset) CALL save_dg(k,1.0_wp, 'qin_reset_qg', i_dgtime)
  IF (ng_reset) CALL save_dg(k,1.0_wp, 'qin_reset_ng', i_dgtime)
  IF (m3g_reset) CALL save_dg(k,1.0_wp, 'qin_reset_m3g', i_dgtime)
#endif


END DO

DEALLOCATE(thresh)

END SUBROUTINE tidy_qin

SUBROUTINE tidy_ain(qfields, aerofields)

REAL(wp), INTENT(IN) :: qfields(:,:)
REAL(wp), INTENT(INOUT) :: aerofields(:,:)

INTEGER :: k

DO k=1,UBOUND(qfields,1)
  IF ((qfields(k, i_ql)+qfields(k,i_qr) <=0.0 .AND. aerofields(k,i_am4)>0.0)   &
       .OR. aerofields(k,i_am4) < 0.0) THEN
!      print*, 'Correcting active liquid aerofields am4...', aerofields(k,i_am4), aerofields(k, i_am9)
#if DEF_MODEL==MODEL_KiD
    CALL save_dg(k, i_here, aerofields(k, i_am4), 'corr_am4', i_dgtime)
#endif
    aerofields(k,i_am4) = 0.0
  END IF
  if (i_am9 > 0)then
  IF ((qfields(k, i_ql)+qfields(k,i_qr) <=0.0 .AND. aerofields(k,i_am9)>0.0)   &
       .OR. aerofields(k,i_am9) < 0.0) THEN
!      print*, 'Correcting active liquid aerofields am9...', aerofields(k,i_am9), aerofields(k, i_am9)
#if DEF_MODEL==MODEL_KiD
    CALL save_dg(k, i_here, aerofields(k, i_am9), 'corr_am9', i_dgtime)
#endif
    aerofields(k,i_am9) = 0.0
  END IF
  end if

  IF (i_am5 > 0) THEN
    IF (((qfields(k,i_qr) <=0.0 .AND. aerofields(k,i_am5)>0.0)  &
       .OR. aerofields(k,i_am5) < 0.0)) THEN
!      print*, 'Correcting active liquid aerofields am5...', aerofields(k,i_am5), aerofields(k, i_am5)
#if DEF_MODEL==MODEL_KiD
      CALL save_dg(k, i_here, aerofields(k, i_am5), 'corr_am5', i_dgtime)
#endif
      aerofields(k,i_am5) = 0.0
    END IF
  END IF
  IF (i_am7 > 0) THEN
    IF (((qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) <=0.0 .AND. aerofields(k,i_am7)>0.0) &
       .OR. aerofields(k,i_am7) < 0.0)) THEN
!      print*, 'Correcting active liquid aerofields am7...', aerofields(k,i_am7), aerofields(k, i_am7)
#if DEF_MODEL==MODEL_KiD
      CALL save_dg(k, i_here, aerofields(k, i_am7), 'corr_am7', i_dgtime)
#endif
      aerofields(k,i_am7) = 0.0
    END IF
  END IF
  IF (i_am8 > 0) THEN
    IF (((qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) <=0.0 .AND. aerofields(k,i_am8)>0.0) &
       .OR. aerofields(k,i_am8) < 0.0)) THEN
 !     print*, 'Correcting active liquid aerofields am8...', aerofields(k,i_am8), aerofields(k, i_am8)
#if DEF_MODEL==MODEL_KiD
      CALL save_dg(k, i_here, aerofields(k, i_am8), 'corr_am8', i_dgtime)
#endif
      aerofields(k,i_am8) = 0.0
    END IF
  END IF

END DO

END SUBROUTINE tidy_ain


SUBROUTINE ensure_positive(k, dt, qfields, procs, params, iprocs_scalable, iprocs_nonscalable &
     , aeroprocs, iprocs_dependent, iprocs_dependent_ns)
  ! Subroutine to ensure parallel processes don't remove more
  ! mass than is available and then rescales all processes
  ! (including number and other terms)

INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN) :: dt
REAL(wp), INTENT(IN) :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT) :: procs(:,:)         ! microphysical process rates

TYPE(hydro_params), INTENT(IN) :: params        ! parameters from hydrometeor variable to test
TYPE(process_name), INTENT(IN) :: iprocs_scalable(:)    ! list of processes to rescale
TYPE(process_name), INTENT(IN), OPTIONAL ::      &
     iprocs_nonscalable(:) ! list of other processes which
  ! provide source or sink, but
  ! which we don't want to rescale
TYPE(process_rate), INTENT(INOUT), OPTIONAL ::   &
     aeroprocs(:,:)        ! associated aerosol process rates
TYPE(process_name), INTENT(IN), OPTIONAL ::      &
     iprocs_dependent(:)   ! list of aerosol processes which
  ! are dependent on rescaled processes and
  ! so should be rescaled themselves
TYPE(process_name), INTENT(IN), OPTIONAL ::      &
     iprocs_dependent_ns(:)   ! list of aerosol processes which
  ! are dependent on rescaled processes but
  ! we don't want to rescale
INTEGER :: iproc, id
REAL(wp) :: delta_scalable, delta_nonscalable, ratio, maxratio

INTEGER :: i_1m, i_2m, i_3m
LOGICAL :: l_rescaled

i_1m = params%i_1m
IF (params%l_2m)i_2m = params%i_2m
IF (params%l_3m)i_3m = params%i_3m
l_rescaled=.FALSE.

IF (qfields(k, i_1m) > 0.0) THEN
    ! First calculate the unscaled increments
  delta_scalable=0.0
  delta_nonscalable=0.0
  DO iproc=1,SIZE(iprocs_scalable)
    IF (iprocs_scalable(iproc)%on) THEN
      id= iprocs_scalable(iproc)%id
      delta_scalable = delta_scalable + procs(k, id)%source(i_1m)
    END IF
  END DO
  delta_scalable = delta_scalable*dt
  IF (PRESENT(iprocs_nonscalable)) THEN
    DO iproc=1,SIZE(iprocs_nonscalable)
      IF (iprocs_nonscalable(iproc)%on) THEN
        id= iprocs_nonscalable(iproc)%id
        delta_nonscalable = delta_nonscalable + procs(k, id)%source(i_1m)
      END IF
    END DO
  END IF
  delta_nonscalable = delta_nonscalable*dt

    ! Test to see if we need to do anything
  IF (delta_scalable + delta_nonscalable + qfields(k, i_1m) < SPACING(qfields(k, i_1m)) &
       .AND. ABS(delta_scalable) > EPSILON(delta_scalable)) THEN
    ratio = -(qfields(k, i_1m)+delta_nonscalable)/delta_scalable
    IF (ratio > 1.0) THEN
      PRINT*, 'DEBUG - something wrong in ensure_positive', ratio, k, i_here, j_here, params%id
      PRINT*, delta_scalable, delta_nonscalable
      IF (PRESENT(iprocs_nonscalable)) THEN
        DO iproc=1,SIZE(iprocs_nonscalable)
          IF (iprocs_nonscalable(iproc)%on) THEN
            id= iprocs_nonscalable(iproc)%id
            IF (procs(k, id)%source(i_1m) /= 0.0_wp) THEN
              PRINT*, 'NS:', iprocs_nonscalable(iproc)%name, '(', procs(k, id)%source(i_1m) ,')'
            END IF
          END IF
        END DO
      END IF
      DO iproc=1,SIZE(iprocs_scalable)
        IF (iprocs_scalable(iproc)%on) THEN
          id= iprocs_scalable(iproc)%id
          IF (procs(k, id)%source(i_1m) /= 0.0_wp) THEN
            PRINT*, 'S:', iprocs_scalable(iproc)%name, '(', procs(k, id)%source(i_1m) ,')'
          END IF
        END IF
      END DO

      DO iproc=1,SIZE(iprocs_scalable)
        IF (iprocs_scalable(iproc)%on) THEN
          id= iprocs_scalable(iproc)%id
          PRINT*, 'ALL S', iprocs_scalable(iproc)%name, procs(k, id)%source
        END IF
      END DO
      DO iproc=1,SIZE(iprocs_nonscalable)
        IF (iprocs_nonscalable(iproc)%on) THEN
          id= iprocs_nonscalable(iproc)%id
          PRINT*, 'ALL NS', iprocs_nonscalable(iproc)%name, procs(k, id)%source
        END IF
      END DO
      STOP
    END IF
      ! Now rescale the scalable processes
    DO iproc=1,SIZE(iprocs_scalable)
      IF (iprocs_scalable(iproc)%on) THEN
        id= iprocs_scalable(iproc)%id
        procs(k, id)%source(i_qstart:i_nstart-1) =     &
                 procs(k, id)%source(i_qstart:i_nstart-1) * ratio
#if DEF_MODEL==MODEL_KiD
        CALL save_dg(k, i_here, ratio, 'rescale_mass_ratio_'//TRIM(iprocs_scalable(iproc)%name), i_dgtime)
#endif
      END IF
    END DO
      ! Set flag to indicate a rescaling was performed
    l_rescaled=.TRUE.
  END IF

    ! Now we need to rescale additional moments
  IF (l_rescaled) THEN
    IF (params%l_2m) THEN ! second moment
        ! First calculate the unscaled increments
      delta_scalable=0.0
      delta_nonscalable=0.0
      DO iproc=1,SIZE(iprocs_scalable)
        IF (iprocs_scalable(iproc)%on) THEN
          id= iprocs_scalable(iproc)%id
          delta_scalable = delta_scalable + procs(k, id)%source(i_2m)
        END IF
      END DO
      delta_scalable = delta_scalable*dt
      IF (ABS(delta_scalable) > EPSILON(delta_scalable)) THEN
        IF (PRESENT(iprocs_nonscalable)) THEN
          DO iproc=1,SIZE(iprocs_nonscalable)
            IF (iprocs_nonscalable(iproc)%on) THEN
              id= iprocs_nonscalable(iproc)%id
              delta_nonscalable = delta_nonscalable + procs(k, id)%source(i_2m)
            END IF
          END DO
        END IF
        delta_nonscalable = delta_nonscalable*dt
          ! ratio may now be greater than 1
        ratio = -(qfields(k, i_2m)+delta_nonscalable)/delta_scalable
          ! Now rescale the scalable processes
        DO iproc=1,SIZE(iprocs_scalable)
          IF (iprocs_scalable(iproc)%on) THEN
            id= iprocs_scalable(iproc)%id
            procs(k, id)%source(i_nstart:i_m3start-1) =     &
                     procs(k, id)%source(i_nstart:i_m3start-1) * ratio
          END IF
        END DO
      END IF
    END IF
    IF (params%l_3m) THEN ! third moment
        ! First calculate the unscaled increments
      delta_scalable=0.0
      delta_nonscalable=0.0
      DO iproc=1,SIZE(iprocs_scalable)
        IF (iprocs_scalable(iproc)%on) THEN
          id= iprocs_scalable(iproc)%id
          delta_scalable = delta_scalable + procs(k, id)%source(i_3m)
        END IF
      END DO
      delta_scalable = delta_scalable*dt
      IF (ABS(delta_scalable) > EPSILON(delta_scalable)) THEN
        IF (PRESENT(iprocs_nonscalable)) THEN
          DO iproc=1,SIZE(iprocs_nonscalable)
            IF (iprocs_nonscalable(iproc)%on) THEN
              id= iprocs_nonscalable(iproc)%id
              delta_nonscalable = delta_nonscalable + procs(k, id)%source(i_3m)
            END IF
          END DO
        END IF
        delta_nonscalable = delta_nonscalable*dt
          ! ratio may now be greater than 1
        ratio = -(qfields(k, i_3m)+delta_nonscalable)/(delta_scalable + EPSILON(1.0))
          ! Now rescale the scalable processes
        DO iproc=1,SIZE(iprocs_scalable)
          IF (iprocs_scalable(iproc)%on) THEN
            id= iprocs_scalable(iproc)%id
            procs(k, id)%source(i_m3start:) =     &
                     procs(k, id)%source(i_m3start:) * ratio
          END IF
        END DO
      END IF
    END IF

      !Now rescale the increments to aerosol
      ! How do we do this?????
      !        print*, 'WARNING: Should be rescaling aerosol?, but not done!'

  ELSE ! What if we haven't rescaled mass, but number is now not conserved?
    IF (l_rescale_on_number) THEN
      IF (params%l_2m) THEN ! second moment
          ! First calculate the unscaled increments
        delta_scalable=0.0
        delta_nonscalable=0.0
        DO iproc=1,SIZE(iprocs_scalable)
          IF (iprocs_scalable(iproc)%on) THEN
            id= iprocs_scalable(iproc)%id
            delta_scalable = delta_scalable + procs(k, id)%source(i_2m)
          END IF
        END DO
        delta_scalable = delta_scalable*dt
        IF (PRESENT(iprocs_nonscalable)) THEN
          DO iproc=1,SIZE(iprocs_nonscalable)
            IF (iprocs_nonscalable(iproc)%on) THEN
              id= iprocs_nonscalable(iproc)%id
              delta_nonscalable = delta_nonscalable + procs(k, id)%source(i_2m)
            END IF
          END DO
        END IF
        delta_nonscalable = delta_nonscalable*dt
          ! Test to see if we need to do anything
        IF (delta_scalable + delta_nonscalable + qfields(k, i_2m) < SPACING(qfields(k, i_2m)) &
             .AND. ABS(delta_scalable) > EPSILON(delta_scalable) ) THEN
          maxratio=(1.0-SPACING(delta_scalable))
          IF (ABS(delta_scalable) < SPACING(delta_scalable)) THEN
            PRINT*, 'Rescaling problem encountered...', time
            PRINT*, delta_scalable, delta_nonscalable
            IF (PRESENT(iprocs_nonscalable)) THEN
              DO iproc=1,SIZE(iprocs_nonscalable)
                IF (iprocs_nonscalable(iproc)%on) THEN
                  id= iprocs_nonscalable(iproc)%id
                  IF (procs(k, id)%source(i_2m) /= 0.0_wp) THEN
                    PRINT*, 'NS:', iprocs_nonscalable(iproc)%name, procs(k, id)%source(i_2m) &
                              , '(', procs(k, id)%source(i_1m) ,')'
                  END IF
                END IF
              END DO
            END IF
            DO iproc=1,SIZE(iprocs_scalable)
              IF (iprocs_scalable(iproc)%on) THEN
                id= iprocs_scalable(iproc)%id
                IF (i_2m > 0) THEN
                  PRINT*, 'S:', iprocs_scalable(iproc)%name, procs(k, id)%source(i_2m) &
                            , '(', procs(k, id)%source(i_1m) ,')'
                END IF
              END IF
            END DO
          END IF
          ratio = -(maxratio*qfields(k, i_2m)+delta_nonscalable)/delta_scalable
          ratio = MAX(ratio, 0.0_wp)
          IF (ratio==0.0_wp) THEN
            PRINT*, 'Rescaling problem encountered - zero ratio...', time
            PRINT*, delta_scalable, delta_nonscalable
            IF (PRESENT(iprocs_nonscalable)) THEN
              DO iproc=1,SIZE(iprocs_nonscalable)
                IF (iprocs_nonscalable(iproc)%on) THEN
                  id= iprocs_nonscalable(iproc)%id
                  IF (procs(k, id)%source(i_2m) /= 0.0_wp)     &
                           PRINT*, 'NS:', iprocs_nonscalable(iproc)%name, procs(k, id)%source(i_2m)
                END IF
              END DO
            END IF
            DO iproc=1,SIZE(iprocs_scalable)
              IF (iprocs_scalable(iproc)%on) THEN
                id= iprocs_scalable(iproc)%id
                IF (procs(k, id)%source(i_2m) /= 0.0_wp)     &
                         PRINT*, 'S:', iprocs_scalable(iproc)%name, procs(k, id)%source(i_2m)
              END IF
            END DO
          END IF


          IF (ratio<0.95) THEN
              ! Some warnings for testing
            PRINT*, 'WARNING: Significantly rescaled number, but not sure what to do with other moments'
            PRINT*, 'id, ratio, bad', params%id, ratio, delta_scalable + delta_nonscalable + qfields(k, i_2m)
            PRINT*, SPACING(qfields(k, i_2m)), (qfields(k, i_2m)+delta_nonscalable), delta_scalable
            PRINT*, 'qfields', qfields(k, i_1m), qfields(k, i_2m)
            IF (PRESENT(iprocs_nonscalable)) THEN
              DO iproc=1,SIZE(iprocs_nonscalable)
                IF (iprocs_nonscalable(iproc)%on) THEN
                  id= iprocs_nonscalable(iproc)%id
                  IF (procs(k, id)%source(i_2m) /= 0.0_wp .OR.     &
                           procs(k, id)%source(i_1m) /= 0.0_wp) THEN
                    PRINT*, 'NS:', iprocs_nonscalable(iproc)%name, procs(k, id)%source(i_2m) &
                              , '(', procs(k, id)%source(i_1m) ,')'
                  END IF
                END IF
              END DO
            END IF
            DO iproc=1,SIZE(iprocs_scalable)
              IF (iprocs_scalable(iproc)%on) THEN
                id= iprocs_scalable(iproc)%id
                IF (procs(k, id)%source(i_2m) /= 0.0_wp .OR.     &
                         procs(k, id)%source(i_1m) /= 0.0_wp) THEN
                  PRINT*, 'S:', iprocs_scalable(iproc)%name, procs(k, id)%source(i_2m) &
                            , '(', procs(k, id)%source(i_1m) ,')'
                END IF
              END IF
            END DO
          END IF

            ! Now rescale the scalable processes
          DO iproc=1,SIZE(iprocs_scalable)
            IF (iprocs_scalable(iproc)%on) THEN
              id= iprocs_scalable(iproc)%id
              IF (i_2m > 0) THEN
                procs(k, id)%source(i_nstart:i_m3start-1) =      &
                          procs(k, id)%source(i_nstart:i_m3start-1) * ratio
              END IF
            END IF
          END DO

            ! Set flag to indicate a rescaling was performed
          l_rescaled=.TRUE.
        END IF
      END IF
    END IF


  END IF

END IF

#if DEF_MODEL==MODEL_KiD
IF (l_rescaled) THEN
  IF (nx>1) THEN
    CALL save_dg(k, i_here, 1.0_wp, 'rescaled', i_dgtime)
  ELSE
    CALL save_dg(k, 1.0_wp, 'rescaled', i_dgtime)
  END IF
END IF
#endif


END SUBROUTINE ensure_positive


SUBROUTINE ensure_positive_aerosol(k, dt, aerofields, aerosol_procs, iprocs)
  ! Subroutine to ensure parallel aerosol processes don't remove more
  ! mass than is available and then rescales all processes
  ! (NB this follows any rescaling due to the parent microphysical processes and
  ! we might lose consistency between number and mass here)

INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN) :: dt
REAL(wp), INTENT(IN) :: aerofields(:,:)
TYPE(process_rate), INTENT(INOUT) :: aerosol_procs(:,:)  ! aerosol process rates

TYPE(process_name), INTENT(IN) :: iprocs(:)    ! list of processes to rescale


INTEGER :: iq, iproc, id
REAL(wp) :: ratio, delta_scalable

DO iq=1,ntotala
  delta_scalable=0.0
  DO iproc=1,SIZE(iprocs)
    IF (iprocs(iproc)%on) THEN
      id= iprocs(iproc)%id
      delta_scalable = delta_scalable + aerosol_procs(k, id)%source(iq)
    END IF
  END DO
  delta_scalable=delta_scalable*dt
  IF (delta_scalable + aerofields(k, iq) < SPACING(aerofields(k, iq))    &
          .AND. ABS(delta_scalable) > SPACING(aerofields(k, iq))) THEN
    ratio = (SPACING(aerofields(k, iq))-aerofields(k, iq))/(delta_scalable)

    DO iproc=1,SIZE(iprocs)
      IF (iprocs(iproc)%on) THEN
        id= iprocs(iproc)%id
        aerosol_procs(k,id)%source(iq) = aerosol_procs(k,id)%source(iq)*ratio
      END IF
    END DO
  END IF
END DO

END SUBROUTINE ensure_positive_aerosol


SUBROUTINE ensure_saturated(k, dt, qfields, procs, iprocs_scalable)
  ! Subroutine to ensure parallel ice processes don't remove more
  ! vapour than is available (i.e. so become subsaturated)
  ! and then rescales processes

  ! Modified so it can also prevent sublimation processes from putting
  ! back too much vapour and so become supersaturated

INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN) :: dt
REAL(wp), INTENT(IN) :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT) :: procs(:,:)         ! microphysical process rates
TYPE(process_name), INTENT(IN) :: iprocs_scalable(:)    ! list of processes to rescale


INTEGER :: iproc, id
REAL(wp) :: delta_scalable, ratio, delta_sat
REAL(wp) :: th, qv, qis

delta_scalable=0.0
DO iproc=1,SIZE(iprocs_scalable)
  IF (iprocs_scalable(iproc)%on) THEN
    id= iprocs_scalable(iproc)%id
    delta_scalable = delta_scalable + procs(k, id)%source(i_qv)*dt
  END IF
END DO

delta_scalable = ABS(delta_scalable)

IF (delta_scalable > SPACING(delta_scalable)) THEN
  qv = qfields(k, i_qv)
  th = qfields(k, i_th)

  qis = qisaturation(th*exner(k), pressure(k)/100.0)

  delta_sat = ABS(qis - qfields(k, i_qv))

  IF (delta_scalable > delta_sat) THEN
    ratio = delta_sat/delta_scalable
    DO iproc=1,SIZE(iprocs_scalable)
      IF (iprocs_scalable(iproc)%on) THEN
        id= iprocs_scalable(iproc)%id
        procs(k, id)%source(i_qstart:i_nstart-1) =     &
                 procs(k, id)%source(i_qstart:i_nstart-1) * ratio
      END IF
    END DO
  END IF
END IF

END SUBROUTINE ensure_saturated

SUBROUTINE cleanup_q(qfields)
  !
  ! Once tidied, there may be some negative numbers due to
  ! rounding error, so should be safe to lose these
  !
REAL(wp), INTENT(INOUT) :: qfields(:,:)

INTEGER :: i1,i2 ! use loops instead of where construct
INTEGER :: lbounds(2), ubounds(2)

lbounds=LBOUND(qfields)
ubounds=UBOUND(qfields)
DO i2=lbounds(2),ubounds(2)
  DO i1=lbounds(1),ubounds(1)
    IF (qfields(i1,i2) < 0.0) THEN
        !          write(6,*) 'Warning zeroing negative element:', i1,i2,qfields(i1,i2)
      qfields(i1,i2) = 0.0
    END IF
  END DO
END DO

END SUBROUTINE cleanup_q

END MODULE mphys_tidy

