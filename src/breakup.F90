MODULE breakup

USE variable_precision, ONLY: wp, iwp
USE process_routines, ONLY: process_rate, process_name,   &
     i_sbrk

USE mphys_parameters, ONLY: hydro_params, DSbrk, tau_sbrk
USE thresholds, ONLY: qr_small, thresh_small

USE m3_incs, ONLY: m3_inc_type2
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE ice_breakup(dt, k, params, qfields, procs)

    !< Subroutine to determine the breakup of large particles
    !< This code is specified for just snow, but could be used for
    !< other species (e.g. rain)
    !< For triple moment species there is a corresponding change in the
    !< 3rd moment assuming shape parameter is not changed
    !< NB: Aerosol mass is not modified by this process
    !
    !< OPTIMISATION POSSIBILITIES: strip out shape parameters

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
TYPE(hydro_params), INTENT(IN) :: params
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)


TYPE(process_name) :: iproc ! processes selected depending on
                                ! which species we're modifying

REAL(wp) :: dnumber, dm1, dm2, dm3

REAL(wp) :: number, mass, m1, m2, m3
TYPE(process_rate), POINTER :: this_proc

REAL(wp) :: n0, lam, mu

REAL(wp) :: Dm ! Mass-weighted mean diameter

SELECT CASE (params%id)
  CASE (4_iwp) !snow
    iproc=i_sbrk
END SELECT

mass = qfields(k, params%i_1m)

IF (mass > thresh_small(params%i_1m) .AND. params%l_2m) THEN ! if no existing ice, we don't bother

  this_proc => procs(k, iproc%id)
  number = qfields(k, params%i_2m)
  IF (params%l_3m) m3 = qfields(k, params%i_3m)

  n0 = dist_n0(k,params%id)
  mu = dist_mu(k,params%id)
  lam = dist_lambda(k,params%id)

  Dm = (1.0 + params%d_x + mu)/lam

  IF (Dm > DSbrk) THEN ! Mean size exceeds threshold

    dnumber = (Dm/DSbrk - 1.0)**params%d_x * number / tau_sbrk

    this_proc%source(params%i_2m) = dnumber

    NULLIFY(this_proc)

  END IF

END IF

END SUBROUTINE ice_breakup

END MODULE breakup
