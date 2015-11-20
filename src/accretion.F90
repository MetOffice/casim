MODULE accretion

USE variable_precision, ONLY: wp
USE passive_fields, ONLY: rho
USE mphys_switches, ONLY:                            &
     i_ql, i_qr, i_nl, i_nr, i_m3r                   &
     , l_2mc, l_2mr                                  &
     , l_aacc, i_am4, i_am5, l_process               &
     , active_rain, isol
USE mphys_constants, ONLY: fixed_cloud_number
USE mphys_parameters, ONLY: hydro_params
USE process_routines, ONLY: process_rate, i_pracw, i_aacw
USE thresholds, ONLY: ql_small, qr_small
USE sweepout_rate, ONLY: sweepout

USE distributions, ONLY: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, k_here, i_here
#endif

IMPLICIT NONE

CONTAINS
  !---------------------------------------------------------------------------
  !> @author
  !> Ben Shipway
  !
  !> @brief
  !> This subroutine calculates increments due to the accretion of
  !> cloud water by rain
  !--------------------------------------------------------------------------- !

SUBROUTINE racw(dt, k, qfields, aerofields, procs, params, aerosol_procs)

REAL(wp), INTENT(IN) :: dt !< microphysics time increment (s)
INTEGER,  INTENT(IN)  :: k  !< k index
REAL(wp), INTENT(IN) :: qfields(:,:)     !< hydrometeor fields
REAL(wp), INTENT(IN) :: aerofields(:,:)  !< aerosol fields
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)         !< hydrometeor process rates
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:) !< aerosol process rates
TYPE(hydro_params), INTENT(IN) :: params !< parameters describing hydrometor size distribution/fallspeeds etc.

REAL(wp) :: dmass, dnumber, damass
REAL(wp) :: m1, m2, dm1, dm2

REAL(wp) :: cloud_mass
REAL(wp) :: cloud_number
REAL(wp) :: rain_mass
REAL(wp) :: rain_number

TYPE(process_rate), POINTER :: this_proc
TYPE(process_rate), POINTER :: aero_proc

REAL(wp) :: mu, n0, lam
LOGICAL :: l_kk_acw=.TRUE.

cloud_mass   = qfields(k, i_ql)
IF (l_2mc) THEN
  cloud_number = qfields(k, i_nl)
ELSE
  cloud_number = fixed_cloud_number
END IF
rain_mass = qfields(k, i_qr)
IF (l_2mr)rain_number = qfields(k, i_nr)

IF (cloud_mass > ql_small .AND. rain_mass > qr_small) THEN

  this_proc => procs(k, i_pracw%id)

  IF (l_kk_acw) THEN
    dmass = MIN(0.9*cloud_mass, 67.0*(cloud_mass*rain_mass)**1.15)
  ELSE
    n0 = dist_n0(k,params%id)
    mu = dist_mu(k,params%id)
    lam = dist_lambda(k,params%id)
    dmass = sweepout(n0, lam, mu, params, rho(k)) * cloud_mass
  END IF

  IF (l_2mc)dnumber = dmass/(cloud_mass/cloud_number)

  this_proc%source(i_ql) = -dmass
  this_proc%source(i_qr) = dmass
  IF (l_2mc) THEN
    this_proc%source(i_nl) = -dnumber
  END IF

  IF (l_aacc .AND. l_process) THEN
    aero_proc => aerosol_procs(k, i_aacw%id)

    IF (active_rain(isol)) THEN
      damass = dmass/cloud_mass*aerofields(k,i_am4)
      aero_proc%source(i_am4) = -damass
      aero_proc%source(i_am5) = damass
    END IF
    NULLIFY(aero_proc)
  END IF

  NULLIFY(this_proc)

END IF

END SUBROUTINE racw

END MODULE accretion
