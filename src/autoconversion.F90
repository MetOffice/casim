MODULE autoconversion

USE variable_precision, ONLY: wp
USE passive_fields, ONLY: rho
USE mphys_switches, ONLY:    &
     i_ql, i_qr, i_nl, i_nr, i_m3r, hydro_complexity &
     , l_2mc, l_2mr, l_3mr                    &
     , l_aaut, i_am4, i_am5, &
     cloud_params, rain_params, l_process        &
     , l_separate_rain
USE mphys_constants, ONLY: rhow, fixed_cloud_number
USE mphys_parameters, ONLY: mu_aut, rain_params
USE process_routines, ONLY: process_rate   &
     , i_praut, i_aaut
USE thresholds, ONLY: ql_small, nl_small, qr_small, qr_sig
USE special, ONLY: pi

USE m3_incs, ONLY: m3_inc_type2, m3_inc_type3

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, i_here
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE raut(dt, k, qfields, aerofields, procs, aerosol_procs)

REAL(wp), INTENT(IN) :: dt
INTEGER,  INTENT(IN) :: k
REAL(wp), INTENT(IN) :: qfields(:,:), aerofields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

REAL(wp) :: dmass, dnumber1, dnumber2, damass
REAL(wp) :: m1, m2, m3, dm1,dm2,dm3
REAL(wp) :: n0, lam, mu

REAL(wp) :: cloud_mass
REAL(wp) :: cloud_number
REAL(wp) :: rain_mass
REAL(wp) :: rain_number
REAL(wp) :: rain_m3
TYPE(process_rate), POINTER :: this_proc
TYPE(process_rate), POINTER :: aero_proc
REAL(wp) :: p1, p2, p3
REAL(wp) :: k1, k2, k3

REAL(wp) :: mu_qc ! < cloud shape parameter (currently only used diagnostically here)

INTEGER :: i

cloud_mass   = qfields(k, i_ql)
IF (l_2mc) THEN
  cloud_number = qfields(k, i_nl)
ELSE
  cloud_number = fixed_cloud_number
END IF

rain_mass = qfields(k, i_qr)
IF (l_2mr) rain_number = qfields(k, i_nr)
IF (l_3mr) rain_m3 = qfields(k, i_m3r)

IF (cloud_mass > ql_small .AND. cloud_number > nl_small) THEN

  this_proc => procs(k, i_praut%id)
  dmass = 1350.0*cloud_mass**2.47*     &
         (cloud_number/1.0e6*rho(k))**(-1.79)
  dmass = MIN(.25*cloud_mass/dt, dmass)
  IF (l_2mc)dnumber1 = dmass/(cloud_mass/cloud_number)

  mu_qc = MIN(15.0_wp, (1000.0E6/cloud_number + 2.0))

  IF (l_2mr)dnumber2 = dmass/(rain_params%c_x*(mu_qc/3.0)*(50.0E-6)**3)

  IF (l_3mr) THEN
    dm1 = dt*dmass/rain_params%c_x
    dm2 = dt*dnumber2
    p1=rain_params%p1
    p2=rain_params%p2
    p3=rain_params%p3
    m1=rain_mass/rain_params%c_x
    m2=rain_number
    m3=rain_m3
    CALL m3_inc_type3(p1, p2, p3,     &
         dm1, dm2, dm3, mu_aut)
    dm3=dm3/dt
  END IF
  this_proc%source(i_ql) = -dmass
  this_proc%source(i_qr) = dmass

  IF (cloud_params%l_2m) THEN
    this_proc%source(i_nl) = -dnumber1
  END IF
  IF (rain_params%l_2m) THEN
    this_proc%source(i_nr) = dnumber2
  END IF
  IF (rain_params%l_3m) THEN
    this_proc%source(i_m3r) = dm3
  END IF


  IF (l_aaut .AND. l_process) THEN
    aero_proc => aerosol_procs(k, i_aaut%id)

        ! Standard Single soluble mode, 2 activated species
    IF (l_separate_rain) THEN
      damass = dmass/cloud_mass*aerofields(k,i_am4)
      aero_proc%source(i_am4) = -damass
      aero_proc%source(i_am5) = damass
    END IF

    NULLIFY(aero_proc)

  END IF


  NULLIFY(this_proc)

END IF

END SUBROUTINE raut

END MODULE autoconversion
