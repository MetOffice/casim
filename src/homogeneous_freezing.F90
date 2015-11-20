MODULE homogeneous

USE variable_precision, ONLY: wp
USE passive_fields, ONLY: rho, pressure, w, exner
USE mphys_switches, ONLY: i_qv, i_ql, i_qi, i_ni, i_th   &
     , hydro_complexity, i_am6, i_an2, l_2mi, l_2ms, l_2mg &
     , i_ns, i_ng, iopt_inuc, i_am7, i_an6, i_am9                 &
     , i_m3r, i_m3g, i_qr, i_qg, i_nr, i_ng, i_nl          &
     , isol, i_am4, i_am8, active_ice, l_process
USE process_routines, ONLY: process_rate, i_homr, i_homc,    &
     i_dhomc, i_dhomr
USE mphys_parameters, ONLY: rain_params, graupel_params, cloud_params,   &
     ice_params, T_hom_freeze
USE mphys_constants, ONLY: Lf, cp
USE qsat_funs, ONLY: qsaturation, qisaturation
USE thresholds, ONLY: thresh_small, thresh_tidy
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
USE special, ONLY: GammaFunc, pi
USE m3_incs, ONLY: m3_inc_type2, m3_inc_type4

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, k_here
#endif
IMPLICIT NONE

CONTAINS

SUBROUTINE ihom_rain(dt, k, qfields, aeroact, dustliq, procs, aerosol_procs)

    ! Calculates immersion freezing of rain drops
    ! See Bigg 1953

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustliq(:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! optional aerosol fields to be processed
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

REAL(wp) :: dmass, dnumber, dmac, coef, dmadl
REAL(wp) :: dm1, dm2, dm3, dm3_g, m1, m2, m3
REAL(wp) :: n0, lam, mu

REAL(wp) :: th
REAL(wp) :: qv, qr, nr
TYPE(process_rate), POINTER :: this_proc
TYPE(process_rate), POINTER :: aero_proc

REAL(wp) :: Tc

REAL(wp), PARAMETER :: A_bigg = 0.66, B_bigg = 100.0

LOGICAL :: l_condition, l_freezeall
LOGICAL :: l_ziegler=.TRUE.

th = qfields(k, i_th)
Tc = th*exner(k) - 273.15
qr = qfields(k, i_qr)
nr = qfields(k, i_nr)

l_condition=(Tc < -4.0 .AND. qr > thresh_small(i_qr))
l_freezeall=.FALSE.
IF (l_condition) THEN
  this_proc => procs(k, i_homr%id)

  n0 = dist_n0(k,rain_params%id)
  mu = dist_mu(k,rain_params%id)
  lam = dist_lambda(k,rain_params%id)

  IF (l_ziegler) THEN
    dnumber = B_bigg*(EXP(-A_bigg*Tc)-1.0)*rho(k)*qr/rain_params%density
    dnumber = MIN(dnumber, nr/dt)
    dmass = (qr/nr)*dnumber
  ELSE
    coef = B_bigg*(pi/6.0)*(EXP(-A_bigg*Tc)-1.0)/rho(k)               &
           * n0 /(lam*lam*lam)/GammaFunc(1.0 + mu)

    dmass = coef * rain_params%c_x * lam**(-rain_params%d_x)          &
           * GammaFunc(4.0 + mu + rain_params%d_x)

    dnumber = coef                                                    &
           * GammaFunc(4.0 + mu)

  END IF

      ! PRAGMATIC HACK - FIX ME
      ! If most of the drops are frozen, do all of them
  IF (dmass*dt >0.95*qr .OR. dnumber*dt > 0.95*nr) THEN
    dmass=qr/dt
    dnumber=nr/dt
    l_freezeall=.TRUE.
  END IF

  this_proc%source(i_qr) = -dmass
  this_proc%source(i_qg) = dmass

  IF (rain_params%l_2m) THEN
    this_proc%source(i_nr) = -dnumber
  END IF
  IF (graupel_params%l_2m) THEN
    this_proc%source(i_ng) = dnumber
  END IF


  IF (rain_params%l_3m) THEN
    IF (l_freezeall) THEN
      dm3=-qfields(k,i_m3r)/dt
    ELSE
      m1=qr/rain_params%c_x
      m2=qfields(k,i_nr)
      m3=qfields(k,i_m3r)

      dm1=-dt*dmass/rain_params%c_x
      dm2=-dt*dnumber

      CALL m3_inc_type2(m1, m2, m3, rain_params%p1, rain_params%p2, rain_params%p3, dm1, dm2, dm3)
      dm3=dm3/dt

    END IF

    this_proc%source(i_m3r) = dm3

  END IF

  IF (graupel_params%l_3m) THEN

    IF (rain_params%l_3m) THEN
      CALL m3_inc_type4(dm3, graupel_params%c_x, rain_params%c_x, rain_params%p3, dm3_g)
    ELSE
      m1=qfields(k,i_qg)/graupel_params%c_x
      m2=qfields(k,i_ng)
      m3=qfields(k,i_m3g)

      dm1=-dt*dmass/graupel_params%c_x
      dm2=-dt*dnumber

      CALL m3_inc_type2(m1, m2, m3, graupel_params%p1, graupel_params%p2, graupel_params%p3, dm1, dm2, dm3_g)
      dm3_g=dm3_g/dt
    END IF

    this_proc%source(i_m3g) = dm3_g

  END IF

  IF (l_process) THEN
    aero_proc => aerosol_procs(k, i_dhomr%id)
    dmac = dnumber*aeroact(k)%mact2_mean

    aero_proc%source(i_am8) = dmac
    aero_proc%source(i_am4) = -dmac

        ! Dust already in the liquid phase
    dmadl = dnumber*dustliq(k)%mact2_mean*dustliq(k)%nratio2
    IF (dmadl /=0.0) THEN
      aero_proc%source(i_am9) = -dmadl
      aero_proc%source(i_am7) = dmadl
    END IF
    NULLIFY(aero_proc)
  END IF

  NULLIFY(this_proc)
END IF

END SUBROUTINE ihom_rain

SUBROUTINE ihom_droplets(dt, k, qfields, aeroact, dustliq, procs, aerosol_procs)

    ! Calculates homogeneous freezing of cloud drops
    ! See Wisener 1972

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    ! aerosol fields
TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustliq(:)

TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

REAL(wp) :: dmass, dnumber, dmac, coef, dmadl
REAL(wp) :: dm1, dm2, dm3, m1, m2, m3
REAL(wp) :: n0, lam, mu

REAL(wp) :: th
REAL(wp) :: qv, ql
TYPE(process_rate), POINTER :: this_proc
TYPE(process_rate), POINTER :: aero_proc

REAL(wp) :: Tc

LOGICAL :: l_condition

th = qfields(k, i_th)
Tc = th*exner(k) - 273.15
ql = qfields(k, i_ql)

l_condition=(Tc < T_hom_freeze .AND. ql > thresh_tidy(i_ql))
IF (l_condition) THEN
  this_proc => procs(k, i_homc%id)

  dmass = MIN(ql, cp*(T_hom_freeze - Tc)/Lf)/dt

  dnumber = dmass*qfields(k, i_nl)/ql

  this_proc%source(i_ql) = -dmass
  this_proc%source(i_qi) = dmass

  IF (cloud_params%l_2m) THEN
    this_proc%source(i_nl) = -dnumber
  END IF
  IF (ice_params%l_2m) THEN
    this_proc%source(i_ni) = dnumber
  END IF

  IF (l_process) THEN
    aero_proc => aerosol_procs(k, i_dhomc%id)
    dmac = dnumber*aeroact(k)%mact1_mean

    aero_proc%source(i_am8) = dmac
    aero_proc%source(i_am4) = -dmac

        ! Dust already in the liquid phase
    dmadl = dnumber*dustliq(k)%mact1_mean*dustliq(k)%nratio1
    IF (dmadl /=0.0) THEN
      aero_proc%source(i_am9) = -dmadl
      aero_proc%source(i_am7) = dmadl
    END IF
    NULLIFY(aero_proc)
  END IF

  NULLIFY(this_proc)
END IF

END SUBROUTINE ihom_droplets
END MODULE homogeneous
