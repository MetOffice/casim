MODULE ice_accretion

USE variable_precision, ONLY: wp, iwp
USE passive_fields, ONLY: rho, pressure, w, exner, TdegC
USE mphys_switches, ONLY: hydro_complexity,  i_ql, l_process,   &
     cloud_params, rain_params, i_am4, i_am7, i_am8, i_am9
USE type_process, ONLY: process_name
USE process_routines, ONLY: process_rate,   &
     i_iacw, i_sacw, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs, &
     i_diacw, i_dsacw, i_dgacw, i_dsacr, i_dgacr, i_draci
USE mphys_parameters, ONLY: hydro_params
USE mphys_constants, ONLY: Ls, cp
USE qsat_funs, ONLY: qsaturation, qisaturation
USE thresholds, ONLY: thresh_small
USE activation, ONLY: activate
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active

USE m3_incs, ONLY: m3_inc_type2, m3_inc_type3
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
USE sweepout_rate, ONLY: sweepout, binary_collection

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, k_here
USE runtime, ONLY: time
#elif DEF_MODEL==MODEL_LEM
USE com_params, ONLY: time
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here
#else
USE diaghelp_um, ONLY: i_here,j_here
USE UM_ParCore, ONLY: mype
#endif
USE special, ONLY: pi, Gammafunc

IMPLICIT NONE

CONTAINS

SUBROUTINE iacc(dt, k, params_X, params_Y, params_Z, qfields, procs,   &
     aeroact, dustliq, aerosol_procs)

    !
    !< CODE TIDYING: Move efficiencies into parameters

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
TYPE(hydro_params), INTENT(IN) :: params_X !< parameters for species which does the collecting
TYPE(hydro_params), INTENT(IN) :: params_Y !< parameters for species which is collected
TYPE(hydro_params), INTENT(IN) :: params_Z !< parameters for species to which resulting amalgamation is sent
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! aerosol fields
TYPE(aerosol_active), INTENT(IN) :: aeroact(:)
TYPE(aerosol_active), INTENT(IN) :: dustliq(:)

    ! optional aerosol fields to be processed
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

TYPE(process_name) :: iproc, iaproc  ! processes selected depending on
                                         ! which species we're depositing on.

REAL(wp) :: dmass_Y, dnumber_Y, dmac, dmad
REAL(wp) :: dmass_X, dmass_Z, dnumber_X, dnumber_Z

REAL(wp) :: number_X, mass_X, m3_X, m1, m2, m3
REAL(wp) :: mass_Y, number_Y, m3_Y
TYPE(process_rate), POINTER :: this_proc, aero_proc

REAL(wp) :: n0_X, lam_X, mu_X
REAL(wp) :: n0_Y, lam_Y, mu_Y
REAL(wp) :: dm1, dm2

REAL(wp) :: Eff ! collection efficiencies need to re-evaluate these
                    ! and put them in properly to mphys_parameters

LOGICAL :: l_condition, l_alternate_Z, l_freezeall_X, l_freezeall_y
LOGICAL :: l_slow ! fall speed is slow compared to interactive species

LOGICAL :: l_aero ! If true then this process will modify aerosol

mass_Y = qfields(k, params_Y%i_1m)
mass_X = qfields(k, params_X%i_1m)

l_condition = mass_X > thresh_small(params_X%i_1m) .AND.     &
       mass_Y > thresh_small(params_Y%i_1m)
l_freezeall_X=.FALSE.
l_freezeall_Y=.FALSE.

IF (l_condition) THEN

      ! initialize variables which may not have been set
  number_X=0.0
  m3_X=0.0
  dnumber_Y=0.0

  l_alternate_Z = params_X%id /= params_Z%id
  l_aero=.FALSE.

  SELECT CASE (params_Y%id)
    CASE(1_iwp) ! cloud liquid is collected
      l_slow=.TRUE. ! Y catagory falls slowly compared to X
      SELECT CASE (params_X%id)
        CASE (3_iwp) !ice collects
          iproc=i_iacw
          Eff=1.0
          iaproc=i_diacw
          l_aero=.TRUE.
        CASE (4_iwp) !snow collects
          iproc=i_sacw
          Eff=1.0
          iaproc=i_dsacw
          l_aero=.TRUE.
        CASE (5_iwp) !graupel collects
          iproc=i_gacw
          Eff=1.0
          iaproc=i_dgacw
          l_aero=.TRUE.
      END SELECT
    CASE(2_iwp) ! rain is collected
      l_slow=.FALSE. ! Y catagory falls slowly compared to X
      SELECT CASE (params_X%id)
        CASE (4_iwp) !snow collects
          iproc=i_sacr
          Eff=1.0
          iaproc=i_dsacr
          l_aero=.TRUE.
        CASE (5_iwp) !graupel collects
          iproc=i_gacr
          Eff=1.0
          iaproc=i_dgacr
          l_aero=.TRUE.
      END SELECT
    CASE(3_iwp)! cloud ice is collected
      l_slow=.TRUE. ! Y catagory falls slowly compared to X
      SELECT CASE (params_X%id)
        CASE (2_iwp) !rain collects
          iproc=i_raci
          Eff=1.0
          iaproc=i_draci
          l_aero=.TRUE.
        CASE (4_iwp) !snow collects
          iproc=i_saci
          !Eff=1.0
          Eff = MIN(1.0_wp, 0.2*EXP(0.08*TdegC(k)))
        CASE (5_iwp) !graupel collects
          iproc=i_gaci
          !Eff=1.0
          Eff = MIN(1.0_wp, 0.2*EXP(0.08*TdegC(k)))
      END SELECT
    CASE(4_iwp) ! snow is collected
      l_slow=.FALSE. ! Y catagory falls slowly compared to X
      SELECT CASE (params_X%id)
        CASE (5_iwp) !graupel collects
          iproc=i_gacs
          !Eff=1.0
          Eff = MIN(1.0_wp, 0.2*EXP(0.08*TdegC(k)))
      END SELECT
  END SELECT

  this_proc => procs(k, iproc%id)

  IF (params_X%l_2m) number_X = qfields(k, params_X%i_2m)
  IF (params_X%l_3m) m3_X = qfields(k, params_X%i_3m)

  n0_X = dist_n0(k,params_X%id)
  mu_X = dist_mu(k,params_X%id)
  lam_X = dist_lambda(k,params_X%id)

  IF (l_slow) THEN  ! collected species is approximated to have zero fallspeed

    dmass_Y = -Eff*sweepout(n0_X, lam_X, mu_X, params_X, rho(k)) * mass_Y
    dmass_Y = MAX(dmass_Y, -mass_Y/dt)

    IF (params_Y%l_2m) THEN
      number_Y = qfields(k, params_Y%i_2m)
      dnumber_Y = dmass_Y * number_Y / mass_Y
    END IF

    IF (l_alternate_Z) THEN ! We move resulting collision to another species.
      IF (params_Y%l_2m) THEN
        number_Y = qfields(k, params_Y%i_2m)
      ELSE
            ! we need an else in here
      END IF
      dmass_X =  -Eff*sweepout(n0_X, lam_X, mu_X, params_X, rho(k), mass_weight=.TRUE.) * number_Y
      dmass_X = MAX(dmass_X, -mass_X/dt)
      dnumber_X = dnumber_Y !dmass_Y*number_X/mass_Y
      dmass_Z = -1.0*( dmass_X + dmass_Y )
      dnumber_Z = -1.0*dnumber_X
    END IF

  ELSE  ! both species have significant fall velocity

    IF (params_Y%l_2m) THEN
      number_Y = qfields(k, params_Y%i_2m)
    ELSE
          ! we need an else in here
    END IF
    IF (params_Y%l_3m) m3_Y = qfields(k, params_Y%i_3m)

    n0_Y = dist_n0(k,params_Y%id)
    mu_Y = dist_mu(k,params_Y%id)
    lam_Y = dist_lambda(k,params_Y%id)

    dmass_Y   = -Eff*binary_collection(n0_X, lam_X, mu_X, n0_Y, lam_Y, mu_Y, params_X, params_Y, rho(k), mass_weight=.TRUE.)
    dmass_Y = MAX(dmass_Y, -mass_Y/dt)

    IF (params_Y%l_2m) THEN
      dnumber_Y = -Eff*binary_collection(n0_X, lam_X, mu_X, n0_Y, lam_Y, mu_Y, params_X, params_Y, rho(k))
    END IF

    IF (l_alternate_Z) THEN ! We move resulting collision to another species.
      dmass_X =  -Eff*binary_collection(n0_Y, lam_Y, mu_Y, n0_X, lam_X, mu_X, params_Y, params_X, rho(k), mass_weight=.TRUE.)
      dmass_X = MAX(dmass_X, -mass_X/dt)
      dnumber_X = dnumber_Y
      dmass_Z = -1.0*( dmass_X + dmass_Y )
      dnumber_Z = -1.0*dnumber_X
    END IF

  END IF


      ! PRAGMATIC HACK - FIX ME (ALTHOUGH I QUITE LIKE IT)
      ! If most of the collected particles are to be removed then remove all of them
  IF (-dmass_Y*dt >0.95*mass_Y .OR. (params_Y%l_2m .AND. -dnumber_Y*dt > 0.95*number_Y)) THEN
    dmass_Y=-mass_Y/dt
    dnumber_Y=-number_Y/dt
    l_freezeall_Y=.TRUE.
  END IF
      ! If most of the collecting particles are to be removed then remove all of them
  IF (l_alternate_Z .AND. (-dmass_X*dt >0.95*mass_X .OR.     &
         (params_X%l_2m .AND. -dnumber_X*dt > 0.95*number_X))) THEN
    dmass_X=-mass_X/dt
    dnumber_X=-number_X/dt
    l_freezeall_X=.TRUE.
    dmass_Z = -1.0*( dmass_X + dmass_Y )
    dnumber_Z = -1.0*dnumber_X
  END IF

  IF (.NOT. l_alternate_Z) THEN
    dmass_X=-dmass_Y
    dnumber_X=0.0
    dnumber_Z=0.0
    dmass_Z=0.0
  END IF


  this_proc%source(params_Y%i_1m) = dmass_Y
  IF (params_Y%l_2m) THEN
    this_proc%source(params_Y%i_2m) = dnumber_Y
  END IF

  this_proc%source(params_X%i_1m) = dmass_X
  IF (params_X%l_2m) THEN
    this_proc%source(params_X%i_2m) = dnumber_X
  END IF

  IF (l_alternate_Z) THEN
    this_proc%source(params_Z%i_1m) = dmass_Z
    IF (params_Z%l_2m) THEN
      this_proc%source(params_Z%i_2m) = dnumber_Z
    END IF
  END IF

      !----------------------
      ! Aerosol processing...
      !----------------------

  IF (l_process .AND. l_aero) THEN
        ! We note that all processes result in a source of frozen water
        ! i.e. params_Z is either ice, snow or graupel

    aero_proc => aerosol_procs(k, iaproc%id)

    IF (params_Y%id == cloud_params%id) THEN
      dmac = ABS(dnumber_Y)*aeroact(k)%mact1_mean*aeroact(k)%nratio1
      dmad = ABS(dnumber_Y)*dustliq(k)%mact1_mean*dustliq(k)%nratio1
    ELSE IF (params_Y%id == rain_params%id) THEN
      dmac = ABS(dnumber_Y)*aeroact(k)%mact2_mean*aeroact(k)%nratio2
      dmad = ABS(dnumber_Y)*dustliq(k)%mact2_mean*dustliq(k)%nratio2
    ELSE IF (params_X%id == cloud_params%id) THEN ! This is never the case !
      dmac = ABS(dnumber_X)*aeroact(k)%mact1_mean*aeroact(k)%nratio1
      dmad = ABS(dnumber_X)*dustliq(k)%mact1_mean*dustliq(k)%nratio1
    ELSE IF (params_X%id == rain_params%id) THEN
      dmac = ABS(dnumber_X)*aeroact(k)%mact2_mean*aeroact(k)%nratio2
      dmad = ABS(dnumber_X)*dustliq(k)%mact2_mean*dustliq(k)%nratio2
    END IF

    aero_proc%source(i_am8) = dmac
    aero_proc%source(i_am4) = -dmac

    aero_proc%source(i_am7) = dmad
    aero_proc%source(i_am9) = -dmad


    NULLIFY(aero_proc)

  END IF


  NULLIFY(this_proc)

END IF

END SUBROUTINE iacc

END MODULE ice_accretion
