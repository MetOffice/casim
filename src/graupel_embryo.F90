MODULE graupel_embryo

USE variable_precision, ONLY: wp
USE process_routines, ONLY: process_rate, process_name,   &
     i_sacw
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
USE passive_fields, ONLY: TdegC, TdegK, qws0, rho

USE mphys_parameters, ONLY: ice_params, snow_params, graupel_params,   &
     rain_params, cloud_params
USE mphys_constants, ONLY: rho0
USE thresholds, ONLY: thresh_sig

USE special, ONLY: pi, Gammafunc
USE m3_incs, ONLY: m3_inc_type2, m3_inc_type3
USE ventfac, ONLY: ventilation
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, i_here
USE runtime, ONLY: time
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE graupel_embryos(dt, k, qfields, procs, aerophys, aerochem, aeroact  &
     , aerosol_procs)

    !< Subroutine to convert some of the small rimed snow to graupel
    !< (Ikawa & Saito 1991)
    !<    !
    !< OPTIMISATION POSSIBILITIES: See gamma functions and distribution calculations
    !<
    !< AEROSOL: NOT DONE YET - internal category transfer

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

REAL(wp) :: cloud_mass, snow_mass, snow_number
REAL(wp) :: dmass, dnumber
REAL(wp) :: m1, m2, m3, dm1, dm2, dm3
REAL(wp) :: snow_n0, snow_lam, snow_mu ! distribution parameters

REAL(wp) :: dnembryo ! rate of embryo creation
REAL(wp) :: embryo_mass = 1.6e-10 ! mass(kg) of a new graupel embryo (should be in parameters)

TYPE(process_rate), POINTER :: this_proc

REAL(wp) :: pgsacw  ! Rate of mass transfer to graupel
REAL(wp) :: Eff

REAL(wp) :: alpha = 4.0 ! Tunable parameter see Reisner 1998, Ikawa+Saito 1991
REAL(wp) :: rhogms      ! Difference between graupel density and snow density

REAL(wp) :: Garg ! Argument for Gamma function
INTEGER  :: pid  ! Process id

snow_mass = qfields(k, snow_params%i_1m)
cloud_mass = qfields(k, cloud_params%i_1m)

IF (snow_mass > thresh_sig(snow_params%i_1m) .AND. cloud_mass > thresh_sig(cloud_params%i_1m)) THEN

  IF (snow_params%l_2m)snow_number=qfields(k, snow_params%i_2m)

  snow_n0 = dist_n0(k,snow_params%id)
  snow_mu = dist_mu(k,snow_params%id)
  snow_lam = dist_lambda(k,snow_params%id)

      !< This efficiency should be the  same as used in sacw calculation, i.e.
      !< they should be defined consistently in the parameters
  Eff = 1.0
  rhogms = graupel_params%density - snow_params%density

  Garg=2.0+2*snow_params%b_x + snow_mu
  pgsacw = (0.75*alpha*dt*pi/rhogms)*Eff*rho(k)*rho(k)*cloud_mass*cloud_mass   &
         *snow_params%a_x*snow_params%a_x*snow_n0 &
         *GammaFunc(Garg) &
         *(2*snow_params%f_x + 2*snow_lam)**(-Garg) &
         *(rho(k)/rho0)**(2*snow_params%g_x)

  dnembryo = MAX(snow_params%density*pgsacw/rhogms/embryo_mass/rho(k), 0.0_wp)
  dnumber = MIN(dnembryo, 0.95*snow_number/dt)

  pid=i_sacw%id
  this_proc => procs(k, pid)
  dmass = this_proc%source(snow_params%i_1m) - pgsacw
  this_proc%source(snow_params%i_1m) = dmass
  this_proc%source(graupel_params%i_1m) = pgsacw
  IF (graupel_params%l_2m) THEN
    this_proc%source(graupel_params%i_2m) = dnumber
  END IF
  IF (snow_params%l_2m) THEN
    this_proc%source(snow_params%i_2m) = -dnumber
  END IF

      !----------------------
      ! Aerosol processing...
      !----------------------

  NULLIFY(this_proc)

END IF

END SUBROUTINE graupel_embryos

END MODULE graupel_embryo
