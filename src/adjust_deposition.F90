MODULE adjust_deposition

  ! As in Harrington et al. 1995 move some of the depositional
  ! growth on ice into the snow category

USE variable_precision, ONLY: wp, iwp
USE passive_fields, ONLY: rho
USE mphys_switches, ONLY: i_qi, i_qs, i_ns, i_m3s
USE mphys_parameters, ONLY: DImax, snow_params, ice_params
USE distributions, ONLY: dist_lambda, dist_mu
USE process_routines, ONLY: process_rate, i_idep, i_sdep, i_saut
USE m3_incs, ONLY: m3_inc_type2

IMPLICIT NONE


CONTAINS

SUBROUTINE adjust_dep(dt, k, procs, qfields)

    ! only grow ice which is not autoconverted to
    ! snow, c.f. Harrington et al (1995)
    ! This assumes that mu_ice==0, so the fraction becomes
    ! P(mu+2, lambda*DImax) (see Abramowitz & Stegun 6.5.13)

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
REAL(wp), INTENT(IN) :: qfields(:,:)

TYPE(process_rate), POINTER :: ice_dep, snow_dep, ice_aut
REAL(wp) :: lam, frac, dmass
INTEGER :: i_pqi, i_pqai, i_pqs, i_pns, i_pm3s
REAL(wp) :: m1,m2,m3,dm1,dm2,dm3

ice_dep => procs(k, i_idep%id)
snow_dep => procs(k, i_sdep%id)
ice_aut => procs(k, i_saut%id)

IF (ice_aut%source(ice_params%i_1m) > 0) THEN

  lam = dist_lambda(k,ice_params%id)
  frac = 1.0-EXP(-lam*DImax)*(1.0+lam*DImax)
  dmass = frac*ice_dep%source(ice_params%i_1m)

  ice_dep%source(ice_params%i_1m) = ice_dep%source(ice_params%i_1m) - dmass
  snow_dep%source(snow_params%i_1m) = snow_dep%source(snow_params%i_1m) + dmass

END IF

NULLIFY(ice_dep)
NULLIFY(ice_aut)
NULLIFY(snow_dep)

END SUBROUTINE adjust_dep

END MODULE adjust_deposition
