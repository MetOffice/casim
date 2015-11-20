MODULE graupel_wetgrowth

USE variable_precision, ONLY: wp
USE process_routines, ONLY: process_rate, process_name,   &
     i_gshd, i_gacw, i_gacr, i_gaci, i_gacs
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
USE passive_fields, ONLY: TdegC,  qws0, rho

USE mphys_parameters, ONLY: ice_params, snow_params, graupel_params,   &
     rain_params, T_hom_freeze, DR_melt
USE mphys_switches, ONLY: i_qv
USE mphys_constants, ONLY: Lv, Lf, Ka, Cwater, Cice
USE thresholds, ONLY: thresh_small

USE m3_incs, ONLY: m3_inc_type2
USE ventfac, ONLY: ventilation
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE wetgrowth(dt, k, qfields, procs, aerophys, aerochem, aeroact   &
     , aerosol_procs)

    !< Subroutine to determine if all liquid accreted by graupel can be
    !< frozen or if there will be some shedding
    !<
    !< CODE TIDYING: Should move efficiencies into parameters


REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! aerosol fields
TYPE(aerosol_phys), INTENT(IN) :: aerophys(:)
TYPE(aerosol_chem), INTENT(IN) :: aerochem(:)
TYPE(aerosol_active), INTENT(IN) :: aeroact(:)

    ! optional aerosol fields to be processed
TYPE(process_rate), INTENT(INOUT), OPTIONAL :: aerosol_procs(:,:)

TYPE(process_name) :: iproc ! processes selected depending on
                                ! which species we're modifying

REAL(wp) :: qv

REAL(wp) :: dmass, dnumber, dm1, dm2, dm3

REAL(wp) :: number, mass, m1, m2, m3

REAL(wp) :: gacw, gacr, gaci, gacs  ! graupel accretion process rates
REAL(wp) :: dnumber_s, dnumber_g  ! number conversion rate from snow/graupel
REAL(wp) :: dmass_s, dmass_g      ! mass conversion rate from snow/graupel

TYPE(process_rate), POINTER :: this_proc

REAL(wp) :: n0, lam, mu, V_x

REAL(wp) :: pgacw, pgacr, pgaci, pgacs
REAL(wp) :: Eff, sdryfac, idryfac
REAL(wp) :: pgaci_dry, pgacs_dry, pgdry

REAL(wp) :: pgwet       !< Amount of liquid that graupel can freeze withouth shedding

REAL(wp) :: pgacsum     !< Sum of all graupel accretion terms


mass = qfields(k, graupel_params%i_1m)

IF (mass > thresh_small(graupel_params%i_1m)) THEN

  pgacw = procs(k, i_gacw%id)%source(graupel_params%i_1m)
  pgacr = procs(k, i_gacr%id)%source(graupel_params%i_1m)
  pgaci = procs(k, i_gaci%id)%source(graupel_params%i_1m)
  pgacs = procs(k, i_gacs%id)%source(graupel_params%i_1m)

      ! Factors for converting wet collection efficiencies to dry ones
  Eff = 1.0
  sdryfac = MIN(1.0_wp, 0.2*EXP(0.08*TdegC(k)) / Eff )
  Eff = 1.0
  idryfac = MIN(1.0_wp, 0.2*EXP(0.08*TdegC(k)) / Eff )

  pgaci_dry = idryfac*pgaci
  pgacs_dry = sdryfac*pgacs

  pgdry = pgacw + pgacr + pgaci_dry + pgacs_dry

  pgacsum = pgacw + pgacr + pgaci + pgacs

  IF (pgacsum > thresh_small(graupel_params%i_1m)) THEN

    qv = qfields(k, i_qv)
    mass = qfields(k, graupel_params%i_1m)
    IF (graupel_params%l_2m) number = qfields(k, graupel_params%i_2m)

    n0 = dist_n0(k,graupel_params%id)
    mu = dist_mu(k,graupel_params%id)
    lam = dist_lambda(k,graupel_params%id)

    CALL ventilation(k, V_x, n0, lam, mu, graupel_params)

    pgwet = (910.0/graupel_params%density)**0.625            &
           *(Lv*(qws0(k) - qv) - Ka*TdegC(k)/rho(k))       &
           /(Lf + Cwater*TdegC(k)) * V_x
    pgwet = pgwet + (pgaci + pgacs)*(1.0 - Cice*TdegC(k)/(Lf + Cwater*TdegC(k)))

    IF (pgdry < pgwet .OR. TdegC(k) < T_hom_freeze) THEN ! Dry growth, so use recalculated gaci, gacs

      IF (pgaci + pgacs > 0.0) THEN
        this_proc => procs(k, i_gacs%id)
        dmass = pgacs_dry
        this_proc%source(graupel_params%i_1m) = dmass
        this_proc%source(snow_params%i_1m) = -dmass
        IF (snow_params%l_2m) THEN
          dnumber = this_proc%source(snow_params%i_2m)*sdryfac
          this_proc%source(snow_params%i_2m) = dnumber
        END IF

        this_proc => procs(k, i_gaci%id)
        dmass = pgaci_dry
        this_proc%source(graupel_params%i_1m) = dmass
        this_proc%source(ice_params%i_1m) = -dmass
        IF (ice_params%l_2m) THEN
          dnumber = this_proc%source(ice_params%i_2m)*idryfac
          this_proc%source(ice_params%i_2m) = dnumber
        END IF
      END IF

    ELSE ! Wet growth mode so recalculate gacr, gshd

      this_proc => procs(k, i_gacr%id)
      IF (ABS(procs(k, i_gacr%id)%source(graupel_params%i_1m)) > thresh_small(graupel_params%i_1m)) THEN
        dmass = pgacr + pgwet - pgacsum
        IF (dmass > 0.0) THEN
          this_proc%source(graupel_params%i_1m) = dmass
          this_proc%source(rain_params%i_1m) = -dmass
        ELSE
          iproc = i_gshd
          this_proc => procs(k, iproc%id)
          dmass = MIN(pgacw, -1.0*dmass)
          this_proc%source(rain_params%i_1m) = dmass
          IF (rain_params%l_2m) THEN
            dnumber = dmass/(rain_params%c_x*DR_melt**3)
            this_proc%source(rain_params%i_2m)     &
                   = dnumber
          END IF

        END IF
      END IF

    END IF

  END IF

      !----------------------
      ! Aerosol processing...
      !----------------------

  NULLIFY(this_proc)

END IF

END SUBROUTINE wetgrowth

END MODULE graupel_wetgrowth
