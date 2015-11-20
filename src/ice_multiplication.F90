MODULE ice_multiplication

USE variable_precision, ONLY: wp, iwp
USE process_routines, ONLY: process_rate, process_name,   &
     i_gacw, i_sacw, i_ihal
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
USE passive_fields, ONLY: TdegC

USE mphys_parameters, ONLY: ice_params, snow_params, graupel_params,   &
     dN_hallet_mossop, M0_hallet_mossop
USE thresholds, ONLY: thresh_small

USE m3_incs, ONLY: m3_inc_type2
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, k_here
USE runtime, ONLY: time
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE hallet_mossop(dt, k, qfields, procs, aerophys, aerochem, aeroact   &
     , aerosol_procs)

    !< Subroutine to determine the ice splintering by Hallet-Mossop
    !< This effect requires prior calculation of the accretion rate of
    !< graupel and snow.
    !< This is a source of ice number and mass and a sink of liquid
    !< (but this is done via the accretion processes already so is
    !< represented here as a sink of snow/graupel)
    !< For triple moment species there is a corresponding change in the
    !< 3rd moment assuming shape parameter is not changed
    !<
    !< AEROSOL: All aerosol sinks/sources are assumed to come from soluble modes
    !
    !< OPTIMISATION POSSIBILITIES:

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

REAL(wp) :: dnumber, dm1, dm2, dm3

REAL(wp) :: number, mass, m1, m2, m3

REAL(wp) :: gacw, sacw  ! accretion process rates
REAL(wp) :: dnumber_s, dnumber_g  ! number conversion rate from snow/graupel
REAL(wp) :: dmass_s, dmass_g      ! mass conversion rate from snow/graupel

TYPE(process_rate), POINTER :: this_proc

TYPE(process_name) :: iproc ! processes selected depending on
                                ! which species we're modifying

REAL(wp) :: n0, lam, mu

REAL(wp) :: Eff  !< splintering efficiency

Eff = 1.0 - ABS(TdegC(k) + 5.0)/2.5 ! linear increase between -2.5/-7.5 and -5C

IF (Eff > 0.0) THEN

  sacw = 0.0
  gacw = 0.0
  IF (snow_params%i_1m > 0)   sacw =     &
         procs(k, i_sacw%id)%source(snow_params%i_1m)
  IF (graupel_params%i_1m > 0) gacw =     &
         procs(k, i_gacw%id)%source(graupel_params%i_1m)

  IF ((sacw + gacw)*dt > thresh_small(snow_params%i_1m)) THEN

    iproc = i_ihal
    this_proc => procs(k, iproc%id)

    dnumber_g = dN_hallet_mossop * Eff * (gacw) ! Number of splinters from graupel
    dnumber_s = dN_hallet_mossop * Eff * (sacw) ! Number of splinters from snow

    dnumber_g = MIN(dnumber_g, 0.5*gacw/M0_hallet_mossop) ! don't remove more than 50% of rimed liquid
    dnumber_s = MIN(dnumber_s, 0.5*sacw/M0_hallet_mossop) ! don't remove more than 50% of rimed liquid

    dmass_g = dnumber_g * M0_hallet_mossop
    dmass_s = dnumber_s * M0_hallet_mossop

        !-------------------
        ! Sources for ice...
        !-------------------
    this_proc%source(ice_params%i_1m) = dmass_g + dmass_s
    this_proc%source(ice_params%i_2m) = dnumber_g + dnumber_s

        !-------------------
        ! Sinks for snow...
        !-------------------
    IF (sacw > 0.0) THEN
      this_proc%source(snow_params%i_1m) = -dmass_s
      this_proc%source(snow_params%i_2m) = 0.0

    END IF

        !---------------------
        ! Sinks for graupel...
        !---------------------
    IF (gacw > 0.0) THEN
      this_proc%source(graupel_params%i_1m) = -dmass_g
      this_proc%source(graupel_params%i_2m) = 0.0

    END IF

        !----------------------
        ! Aerosol processing...
        !----------------------

  END IF

  NULLIFY(this_proc)

END IF

END SUBROUTINE hallet_mossop

END MODULE ice_multiplication
