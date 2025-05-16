MODULE casim_reflec_mod

! This file belongs in section: Large Scale Precipitation
! Adapted for use with CASIM microphysics by Annette K. Miltenberger
! in March 2016
! Updated Summer 2017 by Jonathan Wilkinson to fit with latest code
! structure

! AKM_TODO: currently no cloud fractions used!!!
! (Cloud fractions are now available)

! Description:
! Calculates radar reflectivity in dBZ for all available
! hydrometer species.

USE precision, ONLY: wp
! special module is a CASIM module, the pi value is the same as that in the UM
USE special, ONLY : pi

IMPLICIT NONE
PRIVATE
! Constants for reflectivity
REAL(wp), PARAMETER :: kliq  = 0.93
REAL(wp), PARAMETER :: kice  = 0.174
REAL(wp), PARAMETER :: mm6m3 = 1.0e18

! Define reflectivity limit
REAL(wp), PARAMETER :: ref_lim = -35.0 ! dBZ
! Convert this to linear units (mm6 m-3) using 10.0**p
! Where p = -35 dBZ / 10.0 .
REAL(wp), PARAMETER :: ref_lim_lin = 3.1623e-4

! Mixing ratio limit (below which we ignore the species)
! Set this to be the same as the absolute value used in the
! rest of the microphysics.
REAL(wp), PARAMETER :: mr_lim = 1.0E-8

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CASIM_REFLEC_MOD'

REAL(wp) :: kgraup       ! Reflectivity prefactor due to graupel
REAL(wp) :: krain        ! Reflectivity prefactor due to rain
REAL(wp) :: kclw         ! Reflectivity prefactor due to cloud liquid water
REAL(wp) :: kice_a       ! Reflectivity prefactor due to ice aggregates
REAL(wp) :: kice_c       ! Reflectivity prefactor due to ice crystals

LOGICAL :: l_reflec_warn = .TRUE.

PUBLIC setup_reflec_constants, casim_reflec, kclw, krain, kgraup, kice_a, &
       kice_c, ref_lim

CONTAINS

SUBROUTINE setup_reflec_constants()

USE mphys_parameters, ONLY: cloud_params, rain_params, graupel_params

! Dr Hook Modules
USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim

IMPLICIT NONE

REAL(wp), PARAMETER :: s_ice_density = 900.0 ! [kg m-3]

! Declarations for Dr Hook:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETUP_REFLEC_CONSTANTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

kclw   = kliq / 0.93 * (6.0 / pi / cloud_params%density)**2
krain  = kliq / 0.93 * (6.0 / pi / rain_params%density)**2
kgraup = kice / 0.93 * (6.0 / pi / graupel_params%density)**2
kice_a = kice / 0.93 * (6.0 / pi / s_ice_density)**2
kice_c = kice / 0.93 * (6.0 / pi / s_ice_density)**2

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE setup_reflec_constants

SUBROUTINE casim_reflec( ixy_inner, points, nq, rho, qfields, cffields, &
                         dbz_tot_c, dbz_g_c, dbz_i_c, &
                         dbz_i2_c, dbz_l_c, dbz_r_c )

USE mphys_radar_mod,      ONLY: kliq, kice, mm6m3, ref_lim, ref_lim_lin, &
                                mr_lim

USE distributions,        ONLY: dist_lambda, dist_mu, dist_n0

USE mphys_parameters,     ONLY: cloud_params, rain_params, ice_params,   &
                                snow_params, graupel_params

USE casim_stph,           ONLY: l_rp2_casim, snow_a_x_rp, ice_a_x_rp
USE mphys_switches,       ONLY: l_g, l_warm, l_cfrac_casim_diag_scheme, &
                                i_cfl, i_cfr, i_cfi,       &
                                i_cfs, i_cfg
USE special,              ONLY: Gammafunc
USE mphys_die,            ONLY: throw_mphys_error, warn, std_msg

!!use thresholds, only: cfliq_small ! 
!! not used now removed due to comment in #374

USE distributions,        ONLY: query_distributions, dist_lambda, dist_mu, dist_n0


! Dr Hook Modules
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

!------------------------------------------------------------------------------
! Purpose:
!   Calculates any radar reflectivity required for diagnostics
!   Documentation: UMDP 26.
!------------------------------------------------------------------------------
! Subroutine Arguments

integer, intent(in) :: ixy_inner

INTEGER, INTENT(IN) :: points      ! Number of points calculation is done over
INTEGER, INTENT(IN) :: nq

REAL(wp), INTENT(IN) :: rho(points)    ! Air density [kg m-3]
REAL(wp), INTENT(INOUT) :: qfields(points,nq) ! Graupel mixing ratio [kg kg-1]

REAL(wp), INTENT(OUT) :: dbz_tot_c(points) ! Total reflectivity [dBZ]
REAL(wp), INTENT(OUT) :: dbz_g_c(points)   ! Reflectivity due to graupel [dBZ]
REAL(wp), INTENT(OUT) :: dbz_i_c(points)   ! Reflectivity due to ice
                                           ! aggregates [dBZ]
REAL(wp), INTENT(OUT) :: dbz_r_c(points)   ! Reflectivity due to rain [dBZ]
REAL(wp), INTENT(OUT) :: dbz_i2_c(points)  ! Reflectivity due to ice crystals
                                           ! [dBZ]
REAL(wp), INTENT(OUT) :: dbz_l_c(points)   ! Reflectivity due to liquid
                                           ! cloud [dBZ]

real(wp), intent(in) :: cffields(points,5)

!------------------------------------------------------------------------------
! Local Variables

INTEGER :: i

REAL(wp) :: kgraup       ! Reflectivity prefactor due to graupel
REAL(wp) :: krain        ! Reflectivity prefactor due to rain
REAL(wp) :: kclw         ! Reflectivity prefactor due to cloud liquid water
REAL(wp) :: kice_a       ! Reflectivity prefactor due to ice aggregates
REAL(wp) :: kice_c       ! Reflectivity prefactor due to ice crystals

REAL(wp) :: ze_g(points)   ! Linear reflectivity due to graupel [mm6 m-3]
REAL(wp) :: ze_i(points)   ! Linear reflectivity due to ice agg [mm6 m-3]
REAL(wp) :: ze_i2(points)  ! Linear reflectivity due to ice cry [mm6 m-3]
REAL(wp) :: ze_r(points)   ! Linear reflectivity due to rain    [mm6 m-3]
REAL(wp) :: ze_l(points)   ! Linear reflectivity due to liq cld [mm6 m-3]
REAL(wp) :: ze_tot(points) ! Total linear reflectivity [mm6 m-3]

REAL(wp) :: arg3, arg2     ! arguments for gamma-function

! Declarations for Dr Hook:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CASIM_REFLEC'

!==============================================================================
! Start of calculations
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Apply RP scheme
IF ( l_rp2_casim ) THEN
  snow_params%a_x = snow_a_x_rp
  ice_params%a_x = ice_a_x_rp
END IF

IF (l_cfrac_casim_diag_scheme .AND. l_reflec_warn ) THEN
  WRITE(std_msg, '(A)') 'Using radar reflectivity with cloud fraction scheme.' &
                     //' Assumed cloud fraction is 1'
  CALL throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)

  l_reflec_warn = .FALSE. ! Switch off the warning message
END IF ! l_cfrac_casim_diag_scheme


! Initialise all output dBZ variables. These are initialised to the
! minimum reflectivity as this is usually less than zero.
dbz_tot_c(:) = ref_lim
dbz_g_c(:)   = ref_lim
dbz_i_c(:)   = ref_lim
dbz_r_c(:)   = ref_lim
dbz_i2_c(:)  = ref_lim
dbz_l_c(:)   = ref_lim

! Initialise local variables:
ze_g(:)   = 0.0
ze_i(:)   = 0.0
ze_i2(:)  = 0.0
ze_r(:)   = 0.0
ze_l(:)   = 0.0


kclw   = kliq / 0.93 * (6.0 / pi / cloud_params%density)**2
krain  = kliq / 0.93 * (6.0 / pi / rain_params%density)**2
kgraup = kice / 0.93 * (6.0 / pi / graupel_params%density)**2
kice_a = kice / 0.93 * (6.0 / pi / 900.0)**2!snow_params%density)**2
kice_c = kice / 0.93 * (6.0 / pi / 900.0)**2!ice_params%density)**2


!dist_lambda(:,:) = 0.0 !prob  dont need this as set to zero in query dist
!dist_mu(:,:) = 0.0
!dist_n0(:,:) = 0.0

! determine n0, lambda, mu

CALL query_distributions(ixy_inner, cloud_params, qfields, cffields)
CALL query_distributions(ixy_inner, rain_params, qfields, cffields)
IF (.NOT. l_warm) THEN
  CALL query_distributions(ixy_inner, ice_params, qfields, cffields)
  CALL query_distributions(ixy_inner, snow_params, qfields, cffields)
  IF (l_g) THEN
     CALL query_distributions(ixy_inner, graupel_params, qfields, cffields)
  END IF
END IF

! calculate z_lin
DO i = 1,points
  ! Reflectivity for cloud water
  ! ------------------------------------------------------------------
  IF (qfields(i,cloud_params%i_1m)>mr_lim) THEN
     arg2 = 1.0+dist_mu(i,cloud_params%id)+cloud_params%p2
     arg3 = 2.0*cloud_params%d_x+dist_mu(i,cloud_params%id)+1.0
     ze_l(i) = rho(i)*mm6m3*kclw*cloud_params%c_x**2*(dist_n0(i,cloud_params%id)* &
         dist_lambda(i,cloud_params%id)**(arg2)/Gammafunc(arg2)) / &
         dist_lambda(i,cloud_params%id)**(arg3)*Gammafunc(arg3) * cffields(i,i_cfl) 
     ! Convert from linear (mm^6 m^-3) to dBZ
     IF (ze_l(i)   > ref_lim_lin) dbz_l_c(i)   = 10.0 * LOG10(ze_l(i))
  END IF

  ! Reflectivity for rain water
  ! ------------------------------------------------------------------
  IF (qfields(i,rain_params%i_1m)>mr_lim) THEN
     arg2 = 1.0+dist_mu(i,rain_params%id)+rain_params%p2
     arg3 = 2.0*rain_params%d_x+dist_mu(i,rain_params%id)+1.0
     ze_r(i) = rho(i)*mm6m3*kclw*rain_params%c_x**2*(dist_n0(i,rain_params%id)* &
         dist_lambda(i,rain_params%id)**(arg2)/Gammafunc(arg2)) / &
         dist_lambda(i,rain_params%id)**(arg3)*Gammafunc(arg3) * cffields(i,i_cfr) 
     IF (ze_r(i)   > ref_lim_lin) dbz_r_c(i)   = 10.0 * LOG10(ze_r(i))
  END IF

  IF (.NOT. l_warm) THEN
     ! Reflectivity for ice
     ! ------------------------------------------------------------------
     IF (qfields(i,ice_params%i_1m)>mr_lim) THEN
        arg2 = 1.0+dist_mu(i,ice_params%id)+ice_params%p2
        arg3 = 2.0*ice_params%d_x+dist_mu(i,ice_params%id)+1.0
        ze_i(i) = rho(i)*mm6m3*kice_c*ice_params%c_x**2*(dist_n0(i,ice_params%id)* &
            dist_lambda(i,ice_params%id)**(arg2)/Gammafunc(arg2)) / &
            dist_lambda(i,ice_params%id)**(arg3)*Gammafunc(arg3) * cffields(i,i_cfi)
        IF (ze_i(i)   > ref_lim_lin) dbz_i_c(i)   = 10.0 * LOG10(ze_i(i))
     END IF

     ! Reflectivity for snow
     ! ------------------------------------------------------------------
     IF (qfields(i,snow_params%i_1m)>mr_lim) THEN
        arg2 = 1.0+dist_mu(i,snow_params%id)+snow_params%p2
        arg3 = 2.0*snow_params%d_x+dist_mu(i,snow_params%id)+1.0
        ze_i2(i) = rho(i)*mm6m3*kice_a*snow_params%c_x**2*(dist_n0(i,snow_params%id)* &
            dist_lambda(i,snow_params%id)**(arg2)/Gammafunc(arg2)) / &
            dist_lambda(i,snow_params%id)**(arg3)*Gammafunc(arg3) * cffields(i,i_cfs)
        IF (ze_i2(i)  > ref_lim_lin) dbz_i2_c(i)  = 10.0 * LOG10(ze_i2(i))
     END IF

     IF (l_g) THEN
        ! Reflectivity for  graupel
        ! ------------------------------------------------------------------
        IF (qfields(i,graupel_params%i_1m)>mr_lim) THEN
           arg2 = 1.0+dist_mu(i,graupel_params%id)+graupel_params%p2
           arg3 = 2.0*graupel_params%d_x+dist_mu(i,graupel_params%id)+1.0
           ze_g(i) = rho(i)*mm6m3*kgraup*graupel_params%c_x**2*(dist_n0(i,graupel_params%id)* &
               dist_lambda(i,graupel_params%id)**(arg2)/Gammafunc(arg2)) / &
               dist_lambda(i,graupel_params%id)**(arg3)*Gammafunc(arg3) * cffields(i,i_cfg) 
           IF (ze_g(i)   > ref_lim_lin) dbz_g_c(i)   = 10.0 * LOG10(ze_g(i))
        END IF
     END IF
  END IF

  ! Compute total reflectivity
  ze_tot(i) = ze_g(i) + ze_i(i) + ze_i2(i) + ze_r(i) + ze_l(i)

  ! Convert from linear (mm^6 m^-3) to dBZ
  IF (ze_tot(i) > ref_lim_lin) dbz_tot_c(i) = 10.0 * LOG10(ze_tot(i))

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE casim_reflec
END MODULE casim_reflec_mod
