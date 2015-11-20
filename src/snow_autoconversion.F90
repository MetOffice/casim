MODULE snow_autoconversion

USE variable_precision, ONLY: wp, iwp
USE passive_fields, ONLY: rho, exner, TdegK, pressure
USE mphys_switches, ONLY: iopt_auto,   &
     i_qi, i_qs, i_ni, i_ns, i_m3s, hydro_complexity, &
     l_dsaut, i_qv, i_th, &
     l_harrington
USE mphys_constants, ONLY: rhow, fixed_ice_number,   &
     Ls, cp,  Lv,ka, Dv, Rv
USE mphys_parameters, ONLY: mu_saut, snow_params, ice_params,   &
     DImax, tau_saut, DI2S
USE process_routines, ONLY: process_rate, i_saut
USE qsat_funs, ONLY: qsaturation, qisaturation
USE thresholds, ONLY: thresh_small
USE special, ONLY: pi
USE m3_incs, ONLY: m3_inc_type2, m3_inc_type3

USE lookup, ONLY: gfunc
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, i_here
#endif

IMPLICIT NONE


CONTAINS

SUBROUTINE saut(dt, k, qfields, aerofields, procs, aerosol_procs)

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN) :: qfields(:,:), aerofields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

REAL(wp) :: dmass, dnumber
REAL(wp) :: m1, m2, m3, dm1,dm2,dm3
REAL(wp) :: ice_n0, ice_lam, ice_mu

REAL(wp) :: ice_mass
REAL(wp) :: ice_number
REAL(wp) :: snow_mass
REAL(wp) :: snow_number
REAL(wp) :: snow_m3
REAL(wp) :: th
REAL(wp) :: qv
TYPE(process_rate), POINTER :: this_proc
REAL(wp) :: p1, p2, p3
REAL(wp) :: lami_min , AB, qis
INTEGER :: i

LOGICAL :: l_condition ! logical condition to switch on process

l_condition=.TRUE.
lami_min=0.0

ice_mass   = qfields(k, i_qi)
IF (ice_params%l_2m) THEN
  ice_number = qfields(k, i_ni)
ELSE
  ice_number = fixed_ice_number
END IF

IF (ice_mass < thresh_small(i_qi) .OR. ice_number < thresh_small(i_ni)) THEN
  l_condition = .FALSE.

ELSE

  ice_n0 = dist_n0(k,ice_params%id)
  ice_mu = dist_mu(k,ice_params%id)
  ice_lam = dist_lambda(k,ice_params%id)

  lami_min = (1.0 + ice_mu) / DImax

#if DEF_MODEL==MODEL_KiD
!      if (nx==1)then
  CALL save_dg(k, (1000.0*(ice_mass/ice_number)/ice_params%c_x)**(1.0/ice_params%d_x), &
       'Di_mean', i_dgtime)
  CALL save_dg(k, ice_mass, 'ice_mass_saut', i_dgtime)
  CALL save_dg(k, ice_number, 'ice_number_saut', i_dgtime)
!      else
!        call save_dg(k,i_here,((ice_mass/ice_number)/ice_params%c_x)**ice_params%d_x, &
!           'Di_mean', i_dgtime)
!      end if

#endif

END IF

IF (l_condition) THEN
  IF (l_harrington) THEN
    qv = qfields(k, i_qv)
    th = qfields(k, i_th)
    qis = qisaturation(th*exner(k), pressure(k)/100.0)
    l_condition=qv > qis
  ELSE
    l_condition=ice_lam < lami_min ! LEM autconversion
  END IF
END IF


IF (l_condition) THEN

  this_proc => procs(k, i_saut%id)

  p1=snow_params%p1
  p2=snow_params%p2
  p3=snow_params%p3

  snow_mass = qfields(k, i_qs)
  IF (snow_params%l_2m) snow_number = qfields(k, i_ns)
  IF (snow_params%l_3m) snow_m3 = qfields(k, i_m3s)


  IF (l_harrington) THEN

        !< AB This is used elsewhere, so we should do it more efficiently.
    AB = 1.0/(Lv*Lv/(Rv*ka*TdegK(k)*TdegK(k))*rho(k)    &
           + 1.0/(Dv*qis))
    dnumber = 4.0/DImax/ice_params%density*(qv-qis)    &
           *rho(k)*ice_number*EXP(-ice_lam*DImax)*Dv/AB
    dnumber=MIN(dnumber,0.9*ice_number/dt)
    dmass = pi/6.0*ice_params%density*DImax**3*dnumber
    dmass=MIN(dmass, 0.7*ice_mass/dt)

  ELSE ! LEM version

    dmass = ((lami_min/ice_lam)**ice_params%d_x - 1.0) * ice_mass/tau_saut

    IF (ice_mass < dmass*dt) THEN
#if DEF_MODEL==MODEL_KiD
          !print*, 'DEBUG', k, i_here, ice_mass, 'saut_m', i_dgtime
      CALL save_dg(k, i_here, ice_mass, 'saut_m', i_dgtime)
#endif
      dmass=.5*ice_mass/dt
    END IF

    dmass=MIN(dmass, .5*ice_mass/dt)

    dnumber = dmass/(ice_params%c_x*DI2S**ice_params%d_x)

  END IF

  IF (ice_number < dnumber*dt) THEN
#if DEF_MODEL==MODEL_KiD
    PRINT*, 'DEBUG', k, i_here, ice_number, 'saut_n', i_dgtime
    CALL save_dg(k, i_here, ice_number, 'saut_n', i_dgtime)
#endif
    dnumber = dmass/ice_mass * ice_number
  END IF


  IF (dmass*dt > thresh_small(i_qs)) THEN
    IF (snow_params%l_3m) THEN
      dm1 = dt*dmass/snow_params%c_x
      dm2 = dt*dnumber
      CALL m3_inc_type3(p1, p2, p3,     &
           dm1, dm2, dm3, mu_saut)
      dm3=dm3/dt
    END IF

    this_proc%source(i_qi) = -dmass
    this_proc%source(i_qs) = dmass

    IF (ice_params%l_2m) THEN
      this_proc%source(i_ni) = -dnumber
    END IF
    IF (snow_params%l_2m) THEN
      this_proc%source(i_ns) = dnumber
    END IF
    IF (snow_params%l_3m) THEN
      this_proc%source(i_m3s) = dm3
    END IF
  END IF
      !==============================
      ! No aerosol processing needed
      !==============================

  NULLIFY(this_proc)



END IF

END SUBROUTINE saut

END MODULE snow_autoconversion
