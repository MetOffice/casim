MODULE aggregation

USE variable_precision, ONLY: wp, iwp
USE passive_fields, ONLY: rho, TdegC
USE mphys_switches, ONLY:   &
     i_qr, i_nr, i_m3r, hydro_complexity, l_2mr, l_3mr
USE process_routines, ONLY: process_rate,  process_name,   &
     i_pracr, i_iagg, i_sagg, i_gagg
USE mphys_parameters, ONLY: c_r, d_r, p1, p2, p3, hydro_params,   &
     rain_params
USE mphys_constants, ONLY: rhow, rho0
USE thresholds, ONLY: qr_small, thresh_small

USE sweepout_rate, ONLY: sweepout, binary_collection
USE m3_incs, ONLY: m3_inc_type2
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
USE special, ONLY: pi
USE gauss_casim_micro, ONLY: gauss_casim_func, gaussfunclookup, gaussfunclookup_2d

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE racr(dt, k, qfields, procs)

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

REAL(wp) :: dmass, dnumber, Dr, Eff
REAL(wp) :: m1, m2, m3, dm1, dm2, dm3

REAL(wp) :: rain_mass
REAL(wp) :: rain_number
REAL(wp) :: rain_m3
TYPE(process_rate), POINTER :: this_proc

LOGICAL :: l_beheng=.TRUE.

    !debug
REAL(wp) :: mu, t1, t2 ,H, n0, lam, k1,k2,k3

IF (l_2mr) THEN

  rain_mass = qfields(k, i_qr)
  rain_number = qfields(k, i_nr)
  IF (l_3mr)rain_m3 = qfields(k, i_m3r)

  IF (rain_mass > qr_small .AND. rain_number>0) THEN

    this_proc => procs(k, i_pracr%id)

    IF (l_beheng) THEN

      Dr = (.75/pi)*(rain_mass/rain_number/rhow)**(1.0/d_r)
      IF ( Dr < 600.0e-6) THEN
            ! Modified from original
        Eff=.5
      ELSE
        Eff=0.0
      END IF

      dnumber = Eff*8.0*rain_number*rain_mass*rho(k)

    ELSE

          !???

    END IF

    this_proc%source(i_nr) = -dnumber

    IF (l_3mr) THEN
      m1=rain_mass/rain_params%c_x
      m2=rain_number
      m3=rain_m3

      dm1=0.0
      dm2=-dt*dnumber

      CALL m3_inc_type2(m1, m2, m3, p1, p2, p3, dm1, dm2, dm3)

      dm3=dm3/dt

      this_proc%source(i_m3r) = dm3

    END IF

    NULLIFY(this_proc)

  END IF

END IF

END SUBROUTINE racr

SUBROUTINE ice_aggregation(dt, k, params, qfields, procs)

    !< Subroutine to determine the aggregation of
    !< ice, snow and graupel.  Aggregation causes a sink in number
    !< and for triple moment species a corresponding change in the
    !< 3rd moment assuming shape parameter is not changed
    !< NB: Aerosol mass is not modified by this process
    !
    !< CODE TIDYING: Move efficiencies into parameters

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
TYPE(hydro_params), INTENT(IN) :: params
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)


TYPE(process_name) :: iproc ! processes selected depending on
                                ! which species we're modifying

REAL(wp) :: dnumber, dm1, dm2, dm3

REAL(wp) :: Eff ! collection efficiencies need to re-evaluate these
                    ! and put them in properly to mphys_parameters

REAL(wp) :: number, mass, m1, m2, m3, gaussterm
TYPE(process_rate), POINTER :: this_proc

REAL(wp) :: n0, lam, mu


Eff = MIN(1.0_wp, 0.2*EXP(0.08*TdegC(k)))
SELECT CASE (params%id)
  CASE (3_iwp) !ice
    iproc=i_iagg
  CASE (4_iwp) !snow
    iproc=i_sagg
  CASE (5_iwp) !graupel
    iproc=i_gagg
END SELECT

mass = qfields(k, params%i_1m)

IF (mass > thresh_small(params%i_1m) .AND. params%l_2m) THEN ! if no significant ice, we don't bother

  this_proc => procs(k, iproc%id)
  IF (params%l_2m) number = qfields(k, params%i_2m)
  IF (params%l_3m) m3 = qfields(k, params%i_3m)

  n0 = dist_n0(k,params%id)
  mu = dist_mu(k,params%id)
  lam = dist_lambda(k,params%id)

  IF (params%l_3m) THEN
        !gaussterm = gauss_casim_func(mu, params%b_x)
    CALL gaussfunclookup_2d(params%id, gaussterm, mu, params%b_x)
  ELSE
    CALL gaussfunclookup(params%id, gaussterm)
  END IF

  dnumber = -1.0*gaussterm * pi*0.125*(rho0/rho(k))**params%g_x     &
         * params%a_x*n0*n0*Eff*lam**(-(4.0 + 2.0*mu + params%b_x))

  this_proc%source(params%i_2m) = dnumber

  NULLIFY(this_proc)

END IF

END SUBROUTINE ice_aggregation

END MODULE aggregation
