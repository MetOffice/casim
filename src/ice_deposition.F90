MODULE ice_deposition

USE variable_precision, ONLY: wp, iwp
USE passive_fields, ONLY: rho, pressure, w, exner, TdegK
USE mphys_switches, ONLY: i_qv, i_qi, i_ni, i_th   &
     , hydro_complexity, i_am6, i_an2, l_2mi, l_2ms, l_2mg &
     , i_ns, i_ng, iopt_inuc, i_am7, i_an6                   &
     , l_process, l_passivenumbers_ice, l_passivenumbers, active_number &
     , active_ice, iinsol, i_an12 &
     , i_qr, i_ql, i_qs, i_qg, i_am8, i_am2, i_an11
USE type_process, ONLY: process_name
USE process_routines, ONLY: process_rate, i_idep,    &
     i_dsub, i_sdep, i_gdep, i_dssub, i_dgsub &
     , i_isub, i_ssub, i_gsub, i_iacw, i_raci, i_sacw, i_sacr  &
     , i_gacw, i_gacr
USE mphys_parameters, ONLY: hydro_params, ice_params, rain_params, cloud_params
USE mphys_constants, ONLY: Ls, cp,  Lv, Lf, ka, Dv, Rv
USE qsat_funs, ONLY: qsaturation, qisaturation
USE thresholds, ONLY: thresh_small
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active

USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
USE ventfac, ONLY: ventilation
USE special, ONLY: pi, Gammafunc
USE m3_incs, ONLY: m3_inc_type2


#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, k_here
USE runtime, ONLY: time
#elif DEF_MODEL==MODEL_UM
USE timestep_mod, ONLY: timestep_number
USE diaghelp_um, ONLY: k_here, i_here, j_here
USE UM_ParCore, ONLY: mype
#endif
IMPLICIT NONE

LOGICAL :: l_latenteffects = .FALSE.

CONTAINS

SUBROUTINE idep(dt, k, params, qfields, procs,    &
     dustact, aeroice, aerosol_procs)

    !< Subroutine to determine the deposition/sublimation onto/from
    !< ice, snow and graupel.  There is no source/sink for number
    !< when undergoing deposition, but there is a sink when sublimating.

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
TYPE(hydro_params), INTENT(IN) :: params
REAL(wp), INTENT(IN) :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! aerosol fields
TYPE(aerosol_active), INTENT(IN) :: dustact(:), aeroice(:)

    ! optional aerosol fields to be processed
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

TYPE(process_name) :: iproc, iaproc  ! processes selected depending on
                                         ! which species we're depositing on.

TYPE(process_name) :: i_acw, i_acr ! collection processes with cloud and rain
REAL(wp) :: dmass, dnumber, dmad, dnumber_a, dnumber_d

REAL(wp) :: th
REAL(wp) :: qv, qh
REAL(wp) :: number, mass, m1,m2, m3, dm1,dm2,dm3
TYPE(process_rate), POINTER :: this_proc
TYPE(process_rate), POINTER :: aero_proc

REAL(wp) :: qs, qis, Si, Sw, limit
REAL(wp) :: n0, lam, mu
REAL(wp) :: V_x, AB
REAL(wp) :: frac ! fraction of ice below a threshold

LOGICAL :: l_suball

l_suball=.FALSE. ! do we want to sublimate everything?

mass = qfields(k, params%i_1m)

qv = qfields(k, i_qv)
th = qfields(k, i_th)

qis = qisaturation(th*exner(k), pressure(k)/100.0)

IF (qv>qis) THEN
  SELECT CASE (params%id)
    CASE (3_iwp) !ice
      iproc=i_idep
      iaproc=i_dsub
      i_acw=i_iacw
      i_acr=i_raci
    CASE (4_iwp) !snow
      iproc=i_sdep
      iaproc=i_dssub
      i_acw=i_sacw
      i_acr=i_sacr
    CASE (5_iwp) !graupel
      iproc=i_gdep
      iaproc=i_dgsub
      i_acw=i_gacw
      i_acr=i_gacr
  END SELECT
ELSE
  SELECT CASE (params%id)
    CASE (3_iwp) !ice
      iproc=i_isub
      iaproc=i_dsub
      i_acw=i_iacw
      i_acr=i_raci
    CASE (4_iwp) !snow
      iproc=i_ssub
      iaproc=i_dssub
      i_acw=i_sacw
      i_acr=i_sacr
    CASE (5_iwp) !graupel
      iproc=i_gsub
      iaproc=i_dgsub
      i_acw=i_gacw
      i_acr=i_gacr
  END SELECT
END IF

IF (mass > thresh_small(params%i_1m)) THEN ! if no existing ice, we don't grow/deplete it.

#if DEF_MODEL==MODEL_KiD
  IF (iproc%id==i_sdep%id) THEN
    CALL save_dg(k, i_here, mass, 'depthresh', i_dgtime)
  END IF
#endif

  this_proc => procs(k, iproc%id)
  IF (params%l_2m) number = qfields(k, params%i_2m)
  IF (params%l_3m) m3 = qfields(k, params%i_3m)

  n0 = dist_n0(k,params%id)
  mu = dist_mu(k,params%id)
  lam = dist_lambda(k,params%id)

  CALL ventilation(k, V_x, n0, lam, mu, params)

  AB = 1.0/(Ls*Ls/(Rv*ka*TdegK(k)*TdegK(k))*rho(k)    &
         + 1.0/(Dv*qis))
  dmass= (qv/qis-1.0)*V_x*AB


#if DEF_MODEL==MODEL_KiD
  IF (iproc%id==i_sdep%id) THEN
    CALL save_dg(k, i_here, V_x, 'Vs', i_dgtime)
  END IF
#endif

      ! Include latent heat effects of collection of rain and cloud
      ! as done in Milbrandt & Yau (2005)
  IF (l_latenteffects) THEN
    dmass = dmass - Lf*Ls/(Rv*ka*TdegK(k)*TdegK(k))         &
           *(procs(k, i_acw%id)%source(cloud_params%i_1m) &
           + procs(k, i_acr%id)%source(rain_params%i_1m))
  END IF

      ! Check we haven't become subsaturated and limit if we have (dep only)
      ! NB doesn't account for simultaneous ice/snow growth - checked elsewhere
  IF (dmass > 0.0)dmass = MIN((qv-qis)/dt,dmass)
      ! Check we don't remove too much (sub only)
  IF (dmass < 0.0)dmass = MAX(-mass/dt,dmass)

  this_proc%source(i_qv) = -dmass
  this_proc%source(params%i_1m) = dmass

  IF (params%l_2m) THEN
    dnumber=0.0
    IF (dmass < 0.0) dnumber=dmass * number / mass
  END IF

  IF (-dmass*dt >0.98*mass .OR. (params%l_2m .AND. -dnumber*dt > 0.98*number)) THEN
    l_suball=.TRUE.
    dmass=-mass/dt
    dnumber=-number/dt
  END IF

  IF (params%l_2m)this_proc%source(params%i_2m) = dnumber


  IF (dmass < 0.0 .AND. l_process) THEN ! Only process aerosol if sublimating

    aero_proc => aerosol_procs(k, iaproc%id)

    IF (iaproc%id==i_dsub%id) THEN
      dmad = dnumber*dustact(k)%mact1_mean*dustact(k)%nratio1
      dnumber_d =  dnumber*dustact(k)%nratio1
    ELSE IF (iaproc%id==i_dssub%id) THEN
      dmad = dnumber*dustact(k)%mact2_mean*dustact(k)%nratio2
      dnumber_d =  dnumber*dustact(k)%nratio2
    ELSE IF (iaproc%id==i_dgsub%id) THEN
      dmad = dnumber*dustact(k)%mact3_mean*dustact(k)%nratio3
      dnumber_d =  dnumber*dustact(k)%nratio3
    END IF

    aero_proc%source(i_am7) = dmad
    aero_proc%source(i_am6) = -dmad       ! <WARNING: putting back in coarse mode

    IF (iaproc%id==i_dsub%id) THEN
      dmad = dnumber*aeroice(k)%mact1_mean*aeroice(k)%nratio1
      dnumber_a =  dnumber*aeroice(k)%nratio1
    ELSE IF (iaproc%id==i_dssub%id) THEN
      dmad = dnumber*aeroice(k)%mact2_mean*aeroice(k)%nratio2
      dnumber_a =  dnumber*aeroice(k)%nratio2
    ELSE IF (iaproc%id==i_dgsub%id) THEN
      dmad = dnumber*aeroice(k)%mact3_mean*aeroice(k)%nratio3
      dnumber_a =  dnumber*aeroice(k)%nratio3
    END IF

    aero_proc%source(i_am8) = dmad
    aero_proc%source(i_am2) = -dmad    ! <WARNING: putting back in accumulation mode

    IF (l_passivenumbers_ice) THEN
      aero_proc%source(i_an12) = dnumber_d
    END IF
    aero_proc%source(i_an6) = -dnumber_d  ! <WARNING: putting back in coarse mode

    IF (l_passivenumbers) THEN
      aero_proc%source(i_an11) = dnumber_a
    END IF
    aero_proc%source(i_an2) = -dnumber_a  ! <WARNING: putting back in accumulation mode

    NULLIFY(aero_proc)
  END IF

  NULLIFY(this_proc)

END IF

END SUBROUTINE idep

END MODULE ice_deposition
