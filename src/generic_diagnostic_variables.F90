MODULE generic_diagnostic_variables

USE casim_reflec_mod, ONLY: ref_lim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
ModuleName = 'GENERIC_DIAGNOSTIC_VARIABLES'

SAVE

! Method:
! To add a new diagnostic item, include a logical flag for that
! diagnostic and a real allocatable array of the appropriate
! dimensions to the list below.

! Each parent model sets the logical flag (outside of CASIM)
! and then CASIM uses that flag to allocate the appropriate space
! as required.

!  
! 
!--- see Field et al. https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.4414 
! 
! p = process rate
! n = number tendancy 
! 
!--- hydrometeors
! c = cloud
! w = cloud liquid 
! r = rain
! i = cloud ice
! s = snow
! g = grauple
! 
!--- processes
! con = condensation
! aut = autoconversion # form rain, snow
! ac = acreation  also refered to as capture, riming (mixed phase)
!- for acretion collector hydrometeor appears first, captured hydrometeor last
! sed = sedimentation
! evp = evaporation
! sub = sublimation
! mlt = melt
! dep = deposition, 
! inu = ice nucleation 
! hom = homogeneous 
! hal = Hallet Mossop secondary ice production
! dps = droplet shatter secondary ice production
! ic = ice collision secondary ice production
! 
!--- vertically integrated quantities
! lwp = liquid water path
! rwp = rain water path
! iwp = ice water path
! swp = snow water path
! gwp = grauple water path
! 
!--- rates of change
! dqv = delta q vapour (q = humidity)
! dqc = delta q cloud
! dqr = delta q rain
! dqi = delta q ice
! dqs = delta q snow
! dqg = delta q grauple
! 
!--- radar reflectivity
! dbz = radar reflectivity
! 
!--- process rates
! gacs	Graupel-snow accretion rate
! gacw	Graupel-cloud water accretion rate
! gdep	Deposition rate for graupel
! gmlt	Graupel melting rate
! gshd	Graupel shedding of rain rate
! gsub	Graupel evaporation (sublimation) rate
! homc	Homogeneous nucleation rate of cloud
! homr	Homogeneous freezing of rain
! iacw	Ice-water accretion rate; that is, riming rate of ice crystals
! idep	Deposition rate of ice crystals
! imlt	Melting rate of ice crystals
! inuc	Heterogeneous nucleation rate
! isub	Evaporation (sublimation) of ice crystals
! racw	Rain-cloud water accretion rate
! raut	Rain autoconversion rate
! revp	Rain evaporation rate
! saci	Snow-ice accretion rate
! sacr	Snow-rain accretion rate
! sacw	Snow-cloud water accretion rate; that is, riming rate of snow aggregates
! saut	Snow autoconversion rate (from ice crystals)
! sdep	Deposition rate of snow aggregates
! sedg	Graupel sedimentation rate
! sedi	Ice crystal sedimentation rate
! sedl	Cloud liquid sedimentation rate
! sedr	Rain sedimentation rate
! seds	Snow aggregate sedimentation rate
! smlt	Melting rate of snow aggregates
! ssub	Evaporation (sublimation) of snow aggregates.
! condensation	Condensation/evaporation rate from the Unified Model cloud-fraction scheme
! activation	Droplet number rate from activation



TYPE diaglist

  !---------------------------------
  ! 2D variable logical flags
  !---------------------------------
  LOGICAL :: l_precip         = .FALSE.
  LOGICAL :: l_surface_cloud  = .FALSE.
  LOGICAL :: l_surface_rain   = .FALSE.
  LOGICAL :: l_surface_snow   = .FALSE.
  LOGICAl :: l_surface_graup  = .FALSE.

  !---------------------------------
  ! 3D variable logical flags
  !---------------------------------
  LOGICAL :: l_radar          = .FALSE.
  LOGICAL :: l_process_rates  = .FALSE.
  LOGICAL :: l_phomc          = .FALSE.
  LOGICAL :: l_pinuc          = .FALSE.
  LOGICAL :: l_pidep          = .FALSE.
  LOGICAL :: l_psdep          = .FALSE.
  LOGICAL :: l_piacw          = .FALSE.
  LOGICAL :: l_psacw          = .FALSE.
  LOGICAL :: l_psacr          = .FALSE.
  LOGICAL :: l_pisub          = .FALSE.
  LOGICAL :: l_pssub          = .FALSE.
  LOGICAL :: l_pimlt          = .FALSE.
  LOGICAL :: l_psmlt          = .FALSE.
  LOGICAL :: l_psaut          = .FALSE.
  LOGICAL :: l_psaci          = .FALSE.
  LOGICAL :: l_praut          = .FALSE.
  LOGICAL :: l_pracw          = .FALSE.
  LOGICAL :: l_pracr          = .FALSE.
  LOGICAL :: l_prevp          = .FALSE.
  LOGICAL :: l_pgacw          = .FALSE.
  LOGICAL :: l_pgacs          = .FALSE.
  LOGICAL :: l_pgmlt          = .FALSE.
  LOGICAL :: l_pgsub          = .FALSE.
  LOGICAL :: l_psedi          = .FALSE.
  LOGICAL :: l_pseds          = .FALSE.
  LOGICAL :: l_psedr          = .FALSE.
  LOGICAL :: l_psedg          = .FALSE.
  LOGICAL :: l_psedl          = .FALSE.
  LOGICAL :: l_pcond          = .FALSE.
  LOGICAL :: l_phomr          = .FALSE.
  LOGICAL :: l_nhomc          = .FALSE.
  LOGICAL :: l_nhomr          = .FALSE.
  LOGICAL :: l_nihal          = .FALSE.
  LOGICAL :: l_ninuc          = .FALSE.
  LOGICAL :: l_nsedi          = .FALSE.
  LOGICAL :: l_nseds          = .FALSE.
  LOGICAL :: l_nsedg          = .FALSE.
  LOGICAL :: l_nraut          = .FALSE.
  LOGICAL :: l_nsedl          = .FALSE.
  LOGICAL :: l_nracw          = .FALSE.
  LOGICAL :: l_nracr          = .FALSE.
  LOGICAL :: l_nsedr          = .FALSE.
  LOGICAL :: l_nrevp          = .FALSE.
  LOGICAL :: l_nisub          = .FALSE.
  LOGICAL :: l_nssub          = .FALSE.
  LOGICAL :: l_nsaut          = .FALSE.
  LOGICAL :: l_nsaci          = .FALSE.
  LOGICAL :: l_ngacs          = .FALSE.
  LOGICAL :: l_ngsub          = .FALSE.
  LOGICAL :: l_niacw          = .FALSE.
  LOGICAL :: l_nsacw          = .FALSE.
  LOGICAL :: l_nsacr          = .FALSE.
  LOGICAL :: l_nimlt          = .FALSE.
  LOGICAL :: l_nsmlt          = .FALSE.
  LOGICAL :: l_ngacw          = .FALSE.
  LOGICAL :: l_ngmlt          = .FALSE.
  LOGICAL :: l_pihal          = .FALSE.
  LOGICAL :: l_praci_g        = .FALSE.
  LOGICAL :: l_praci_r        = .FALSE.
  LOGICAL :: l_praci_i        = .FALSE.
  LOGICAL :: l_nraci_g        = .FALSE.
  LOGICAL :: l_nraci_r        = .FALSE.
  LOGICAL :: l_nraci_i        = .FALSE. 
  LOGICAL :: l_pidps          = .FALSE.
  LOGICAL :: l_nidps          = .FALSE.
  LOGICAL :: l_pgaci          = .FALSE.
  LOGICAL :: l_ngaci          = .FALSE.
  LOGICAL :: l_niics_s        = .FALSE.
  LOGICAL :: l_niics_i        = .FALSE. 
  LOGICAL :: l_rainfall_3d    = .FALSE.
  LOGICAL :: l_snowfall_3d    = .FALSE.
  LOGICAL :: l_snowonly_3d    = .FALSE.
  LOGICAL :: l_graupfall_3d   = .FALSE.
  LOGICAL :: l_mphys_pts      = .FALSE.


!PRF water path
  LOGICAL :: l_lwp          = .FALSE.
  LOGICAL :: l_rwp          = .FALSE.
  LOGICAL :: l_iwp          = .FALSE.
  LOGICAL :: l_swp          = .FALSE.
  LOGICAL :: l_gwp          = .FALSE.
!PRF



  !---------------------------------
  ! logical flags for theta tendencies
  ! (based on LEM and MONC, should
  ! work with UM)
  !--------------------------------
  LOGICAL :: l_tendency_dg    = .FALSE.
  LOGICAL :: l_dth            = .FALSE.
   
  !---------------------------------
  ! logical flags for
  ! mass tendencies (based on LEM and
  ! MONC, should work with UM)
  !--------------------------------
  LOGICAL :: l_dqv          = .FALSE.
  LOGICAL :: l_dqc          = .FALSE.
  LOGICAL :: l_dqr          = .FALSE.
  LOGICAL :: l_dqi          = .FALSE.
  LOGICAL :: l_dqs          = .FALSE.
  LOGICAL :: l_dqg          = .FALSE.

  
  !--------------------------------
  ! 2D variable arrays
  !--------------------------------
  ! Surface Precipitation rates
  REAL, ALLOCATABLE :: precip(:,:)
  REAL, ALLOCATABLE :: SurfaceRainR(:,:)
  REAL, ALLOCATABLE :: SurfaceCloudR(:,:)
  REAL, ALLOCATABLE :: SurfaceSnowR(:,:)
  REAL, ALLOCATABLE :: SurfaceGraupR(:,:)

!PRF
  REAL, ALLOCATABLE :: lwp(:,:)
  REAL, ALLOCATABLE :: rwp(:,:)
  REAL, ALLOCATABLE :: iwp(:,:)
  REAL, ALLOCATABLE :: swp(:,:)
  REAL, ALLOCATABLE :: gwp(:,:)



  !--------------------------------
  ! 3D variable arrays
  !--------------------------------
  ! Radar Reflectivity arrays
  REAL, ALLOCATABLE :: dbz_tot(:,:,:)
  REAL, ALLOCATABLE :: dbz_g(:,:,:)
  REAL, ALLOCATABLE :: dbz_i(:,:,:)
  REAL, ALLOCATABLE :: dbz_s(:,:,:)
  REAL, ALLOCATABLE :: dbz_l(:,:,:)
  REAL, ALLOCATABLE :: dbz_r(:,:,:)

  ! Process rate diagnostics
  REAL, ALLOCATABLE :: phomc(:,:,:)
  REAL, ALLOCATABLE :: pinuc(:,:,:)
  REAL, ALLOCATABLE :: pidep(:,:,:)
  REAL, ALLOCATABLE :: psdep(:,:,:)
  REAL, ALLOCATABLE :: piacw(:,:,:)
  REAL, ALLOCATABLE :: psacw(:,:,:)
  REAL, ALLOCATABLE :: psacr(:,:,:)
  REAL, ALLOCATABLE :: pisub(:,:,:)
  REAL, ALLOCATABLE :: pssub(:,:,:)
  REAL, ALLOCATABLE :: pimlt(:,:,:)
  REAL, ALLOCATABLE :: psmlt(:,:,:)
  REAL, ALLOCATABLE :: psaut(:,:,:)
  REAL, ALLOCATABLE :: psaci(:,:,:)
  REAL, ALLOCATABLE :: praut(:,:,:)
  REAL, ALLOCATABLE :: pracw(:,:,:)
  REAL, ALLOCATABLE :: pracr(:,:,:)
  REAL, ALLOCATABLE :: prevp(:,:,:)
  REAL, ALLOCATABLE :: pgacw(:,:,:)
  REAL, ALLOCATABLE :: pgacs(:,:,:)
  REAL, ALLOCATABLE :: pgmlt(:,:,:)
  REAL, ALLOCATABLE :: pgsub(:,:,:)
  REAL, ALLOCATABLE :: psedi(:,:,:)
  REAL, ALLOCATABLE :: pseds(:,:,:)
  REAL, ALLOCATABLE :: psedr(:,:,:)
  REAL, ALLOCATABLE :: psedg(:,:,:)
  REAL, ALLOCATABLE :: psedl(:,:,:)
  REAL, ALLOCATABLE :: pcond(:,:,:)
  REAL, ALLOCATABLE :: phomr(:,:,:)
  REAL, ALLOCATABLE :: nhomc(:,:,:)
  REAL, ALLOCATABLE :: nhomr(:,:,:)
  REAL, ALLOCATABLE :: nihal(:,:,:)
  REAL, ALLOCATABLE :: ninuc(:,:,:)
  REAL, ALLOCATABLE :: nsedi(:,:,:)
  REAL, ALLOCATABLE :: nseds(:,:,:)
  REAL, ALLOCATABLE :: nsedg(:,:,:)
  REAL, ALLOCATABLE :: nraut(:,:,:)   
  REAL, ALLOCATABLE :: nsedl(:,:,:)   
  REAL, ALLOCATABLE :: nracw(:,:,:)   
  REAL, ALLOCATABLE :: nracr(:,:,:)   
  REAL, ALLOCATABLE :: nsedr(:,:,:)   
  REAL, ALLOCATABLE :: nrevp(:,:,:)   
  REAL, ALLOCATABLE :: nisub(:,:,:)   
  REAL, ALLOCATABLE :: nssub(:,:,:)   
  REAL, ALLOCATABLE :: nsaut(:,:,:)   
  REAL, ALLOCATABLE :: nsaci(:,:,:)   
  REAL, ALLOCATABLE :: ngacs(:,:,:)   
  REAL, ALLOCATABLE :: ngsub(:,:,:)   
  REAL, ALLOCATABLE :: niacw(:,:,:)
  REAL, ALLOCATABLE :: nsacw(:,:,:)
  REAL, ALLOCATABLE :: nsacr(:,:,:)   
  REAL, ALLOCATABLE :: nimlt(:,:,:)  
  REAL, ALLOCATABLE :: nsmlt(:,:,:)   
  REAL, ALLOCATABLE :: ngacw(:,:,:)   
  REAL, ALLOCATABLE :: ngmlt(:,:,:)   
  REAL, ALLOCATABLE :: pihal(:,:,:)   
  REAL, ALLOCATABLE :: praci_g(:,:,:) 
  REAL, ALLOCATABLE :: praci_r(:,:,:) 
  REAL, ALLOCATABLE :: praci_i(:,:,:) 
  REAL, ALLOCATABLE :: nraci_g(:,:,:) 
  REAL, ALLOCATABLE :: nraci_r(:,:,:) 
  REAL, ALLOCATABLE :: nraci_i(:,:,:) 
  REAL, ALLOCATABLE :: pidps(:,:,:)   
  REAL, ALLOCATABLE :: nidps(:,:,:)   
  REAL, ALLOCATABLE :: pgaci(:,:,:)   
  REAL, ALLOCATABLE :: ngaci(:,:,:)   
  REAL, ALLOCATABLE :: niics_s(:,:,:) 
  REAL, ALLOCATABLE :: niics_i(:,:,:) 

  !---------------------------------
  ! 3D variables logical for
  ! theta tendencies (based on LEM and
  ! MONC, should work with UM)
  !--------------------------------
  REAL, ALLOCATABLE :: dth_total(:,:,:)
  REAL, ALLOCATABLE :: dth_cond_evap(:,:,:)

  !---------------------------------
  ! 3D variables for
  ! mass tendencies (based on LEM and
  ! MONC, should work with UM)
  !--------------------------------
  REAL, ALLOCATABLE :: dqv_total(:,:,:)
  REAL, ALLOCATABLE :: dqv_cond_evap(:,:,:)
  REAL, ALLOCATABLE :: dqc(:,:,:)
  REAL, ALLOCATABLE :: dqr(:,:,:)
  REAL, ALLOCATABLE :: dqi(:,:,:)
  REAL, ALLOCATABLE :: dqs(:,:,:)
  REAL, ALLOCATABLE :: dqg(:,:,:)

  !---------------------------------
  ! 3D rainfall and snowfall rates
  !---------------------------------
  REAL, ALLOCATABLE :: rainfall_3d(:,:,:)
  REAL, ALLOCATABLE :: snowfall_3d(:,:,:)
  REAL, ALLOCATABLE :: snowonly_3d(:,:,:)
  REAL, ALLOCATABLE :: graupfall_3d(:,:,:)

  !---------------------------------
  ! Points on which we do microphysics
  !---------------------------------
  LOGICAL, ALLOCATABLE :: mphys_pts(:,:,:)

END TYPE diaglist

TYPE (diaglist) :: casdiags

CONTAINS

SUBROUTINE allocate_diagnostic_space(i_start, i_end, j_start, j_end, k_start, k_end)

USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
USE mphys_parameters, ONLY: zero_real_wp

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_start ! Start of i array
INTEGER, INTENT(IN) :: i_end ! End of i array
INTEGER, INTENT(IN) :: j_start ! Start of j array
INTEGER, INTENT(IN) :: j_end ! End of j array
INTEGER, INTENT(IN) :: k_start ! Start of k array
INTEGER, INTENT(IN) :: k_end ! End of k array

! Local variables
CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_DIAGNOSTIC_SPACE'
INTEGER :: k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( casdiags % l_precip) THEN

  ALLOCATE ( casdiags % precip(i_start:i_end, j_start:j_end) )
  casdiags % precip(:,:) = zero_real_wp

END IF ! casdiags % l_precip

IF ( casdiags % l_surface_cloud ) THEN

  ALLOCATE ( casdiags % SurfaceCloudR(i_start:i_end, j_start:j_end) )
  casdiags % SurfaceCloudR(:,:) = zero_real_wp

END IF

IF ( casdiags % l_surface_rain ) THEN

  ALLOCATE ( casdiags % SurfaceRainR(i_start:i_end, j_start:j_end) )
  casdiags % SurfaceRainR(:,:) = zero_real_wp

END IF

IF ( casdiags % l_surface_snow ) THEN

  ALLOCATE ( casdiags % SurfaceSnowR(i_start:i_end, j_start:j_end) )
  casdiags % SurfaceSnowR(:,:) = zero_real_wp

END IF

IF ( casdiags % l_surface_graup ) THEN

  ALLOCATE ( casdiags % SurfaceGraupR(i_start:i_end, j_start:j_end) )
  casdiags % SurfaceGraupR(:,:) = zero_real_wp

END IF

IF ( casdiags % l_radar ) THEN

  ! A single logical allocates all variables
  ALLOCATE ( casdiags % dbz_tot(i_start:i_end, j_start:j_end, k_start:k_end) )
  ALLOCATE ( casdiags % dbz_g(i_start:i_end, j_start:j_end, k_start:k_end) )
  ALLOCATE ( casdiags % dbz_i(i_start:i_end, j_start:j_end, k_start:k_end) )
  ALLOCATE ( casdiags % dbz_s(i_start:i_end, j_start:j_end, k_start:k_end) )
  ALLOCATE ( casdiags % dbz_l(i_start:i_end, j_start:j_end, k_start:k_end) )
  ALLOCATE ( casdiags % dbz_r(i_start:i_end, j_start:j_end, k_start:k_end) )

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % dbz_tot(:,:,k) = ref_lim
    casdiags % dbz_g(:,:,k)   = ref_lim
    casdiags % dbz_i(:,:,k)   = ref_lim
    casdiags % dbz_s(:,:,k)   = ref_lim
    casdiags % dbz_l(:,:,k)   = ref_lim
    casdiags % dbz_r(:,:,k)   = ref_lim
  END DO
!$OMP END PARALLEL DO

END IF ! casdiags % l_radar

IF (casdiags % l_phomc) THEN
  ALLOCATE ( casdiags % phomc(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % phomc(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nhomc) THEN
  ALLOCATE ( casdiags % nhomc(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % nhomc(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pinuc) THEN
  ALLOCATE ( casdiags % pinuc(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pinuc(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_ninuc) THEN
  ALLOCATE ( casdiags % ninuc(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % ninuc(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pidep) THEN
  ALLOCATE ( casdiags % pidep(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pidep(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psdep) THEN
  ALLOCATE ( casdiags % psdep(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psdep(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_piacw) THEN
  ALLOCATE ( casdiags % piacw(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % piacw(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psacw) THEN
  ALLOCATE ( casdiags % psacw(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psacw(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psacr) THEN
  ALLOCATE ( casdiags % psacr(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psacr(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pisub) THEN
  ALLOCATE ( casdiags % pisub(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pisub(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pssub) THEN
  ALLOCATE ( casdiags % pssub(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pssub(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pimlt) THEN
  ALLOCATE ( casdiags % pimlt(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pimlt(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psmlt) THEN
  ALLOCATE ( casdiags % psmlt(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psmlt(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psaut) THEN
  ALLOCATE ( casdiags % psaut(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psaut(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psaci) THEN
  ALLOCATE ( casdiags % psaci(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psaci(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_praut) THEN
  ALLOCATE ( casdiags % praut(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % praut(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pracw) THEN
  ALLOCATE ( casdiags % pracw(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pracw(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pracr) THEN
  ALLOCATE ( casdiags % pracr(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pracr(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_prevp) THEN
  ALLOCATE ( casdiags % prevp(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % prevp(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pgacw) THEN
  ALLOCATE ( casdiags % pgacw(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pgacw(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pgacs) THEN
  ALLOCATE ( casdiags % pgacs(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pgacs(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pgmlt) THEN
  ALLOCATE ( casdiags % pgmlt(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pgmlt(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pgsub) THEN
  ALLOCATE ( casdiags % pgsub(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pgsub(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psedi) THEN
  ALLOCATE ( casdiags % psedi(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psedi(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsedi) THEN
  ALLOCATE ( casdiags % nsedi(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % nsedi(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pseds) THEN
  ALLOCATE ( casdiags % pseds(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % pseds(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nseds) THEN
  ALLOCATE ( casdiags % nseds(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % nseds(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psedr) THEN
  ALLOCATE ( casdiags % psedr(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psedr(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psedg) THEN
  ALLOCATE ( casdiags % psedg(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psedg(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsedg) THEN
  ALLOCATE ( casdiags % nsedg(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % nsedg(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psedl) THEN
  ALLOCATE ( casdiags % psedl(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % psedl(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pcond) THEN
  ALLOCATE ( casdiags % pcond(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
  casdiags % pcond(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_phomr) THEN
  ALLOCATE ( casdiags % phomr(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
  casdiags % phomr(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nhomr) THEN
  ALLOCATE ( casdiags % nhomr(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % nhomr(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nihal) THEN
  ALLOCATE ( casdiags % nihal(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nihal(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nraut) THEN
  ALLOCATE ( casdiags % nraut(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nraut(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsedl) THEN
  ALLOCATE ( casdiags % nsedl(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nsedl(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nracw) THEN
  ALLOCATE ( casdiags % nracw(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nracw(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nracr) THEN
  ALLOCATE ( casdiags % nracr(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nracr(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsedr) THEN
  ALLOCATE ( casdiags % nsedr(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nsedr(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nrevp) THEN
  ALLOCATE ( casdiags % nrevp(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nrevp(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nisub) THEN
  ALLOCATE ( casdiags % nisub(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nisub(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nssub) THEN
  ALLOCATE ( casdiags % nssub(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nssub(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsaut) THEN
  ALLOCATE ( casdiags % nsaut(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nsaut(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsaci) THEN
  ALLOCATE ( casdiags % nsaci(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nsaci(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_ngacs) THEN
  ALLOCATE ( casdiags % ngacs(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % ngacs(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_ngsub) THEN
  ALLOCATE ( casdiags % ngsub(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % ngsub(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_niacw) THEN
  ALLOCATE ( casdiags % niacw(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % niacw(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsacw) THEN
  ALLOCATE ( casdiags % nsacw(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nsacw(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsacr) THEN
  ALLOCATE ( casdiags % nsacr(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nsacr(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nimlt) THEN
  ALLOCATE ( casdiags % nimlt(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nimlt(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsmlt) THEN
  ALLOCATE ( casdiags % nsmlt(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nsmlt(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_ngacw) THEN
  ALLOCATE ( casdiags % ngacw(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % ngacw(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_ngmlt) THEN
  ALLOCATE ( casdiags % ngmlt(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % ngmlt(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pihal) THEN
  ALLOCATE ( casdiags % pihal(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % pihal(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_praci_g) THEN
  ALLOCATE ( casdiags % praci_g(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % praci_g(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_praci_r) THEN
  ALLOCATE ( casdiags % praci_r(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % praci_r(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_praci_i) THEN
  ALLOCATE ( casdiags % praci_i(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % praci_i(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nraci_g) THEN
  ALLOCATE ( casdiags % nraci_g(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nraci_g(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nraci_r) THEN
  ALLOCATE ( casdiags % nraci_r(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nraci_r(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nraci_i) THEN
  ALLOCATE ( casdiags % nraci_i(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nraci_i(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pidps) THEN
  ALLOCATE ( casdiags % pidps(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % pidps(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nidps) THEN
  ALLOCATE ( casdiags % nidps(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % nidps(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pgaci) THEN
  ALLOCATE ( casdiags % pgaci(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % pgaci(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_ngaci) THEN
  ALLOCATE ( casdiags % ngaci(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % ngaci(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_niics_s) THEN
  ALLOCATE ( casdiags % niics_s(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % niics_s(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_niics_i) THEN
  ALLOCATE ( casdiags % niics_i(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % niics_i(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

! potential temp and mass tendencies
IF (casdiags % l_dth) THEN
  ALLOCATE ( casdiags % dth_total(i_start:i_end, j_start:j_end, k_start:k_end) )
  ALLOCATE ( casdiags % dth_cond_evap(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % dth_total(:,:,k) = zero_real_wp
    casdiags % dth_cond_evap(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqv) THEN
  ALLOCATE ( casdiags % dqv_total(i_start:i_end, j_start:j_end, k_start:k_end) )
  ALLOCATE ( casdiags % dqv_cond_evap(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % dqv_total(:,:,k) = zero_real_wp
    casdiags % dqv_cond_evap(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqc) THEN
  ALLOCATE ( casdiags % dqc(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % dqc(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqr) THEN
  ALLOCATE ( casdiags % dqr(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % dqr(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqi) THEN
  ALLOCATE ( casdiags % dqi(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % dqi(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqs) THEN
  ALLOCATE ( casdiags % dqs(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % dqs(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqg) THEN
  ALLOCATE ( casdiags % dqg(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % dqg(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_rainfall_3d) THEN
  ALLOCATE ( casdiags % rainfall_3d(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % rainfall_3d(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
END IF

IF (casdiags % l_snowfall_3d) THEN
  ALLOCATE ( casdiags % snowfall_3d(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % snowfall_3d(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
END IF

IF (casdiags % l_snowonly_3d) THEN
  ALLOCATE ( casdiags % snowonly_3d(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % snowonly_3d(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
END IF

IF (casdiags % l_graupfall_3d) THEN
  ALLOCATE ( casdiags % graupfall_3d(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % graupfall_3d(:,:,k) = zero_real_wp
  END DO
!$OMP END PARALLEL DO
END IF

IF (casdiags % l_mphys_pts) THEN
  ALLOCATE ( casdiags % mphys_pts(i_start:i_end, j_start:j_end, k_start:k_end) )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)                    &
!$OMP SHARED(k_start, k_end, casdiags)
  DO k = k_start, k_end
    casdiags % mphys_pts(:,:,k) = .FALSE.
  END DO
!$OMP END PARALLEL DO
END IF

!PRF water paths
IF (casdiags % l_lwp) THEN
  ALLOCATE ( casdiags % lwp(i_start:i_end, j_start:j_end) )
  casdiags % lwp(:,:) = zero_real_wp
END IF
IF (casdiags % l_rwp) THEN
  ALLOCATE ( casdiags % rwp(i_start:i_end, j_start:j_end) )
  casdiags % rwp(:,:) = zero_real_wp
END IF
IF (casdiags % l_iwp) THEN
  ALLOCATE ( casdiags % iwp(i_start:i_end, j_start:j_end) )
  casdiags % iwp(:,:) = zero_real_wp
END IF
IF (casdiags % l_swp) THEN
  ALLOCATE ( casdiags % swp(i_start:i_end, j_start:j_end) )
  casdiags % swp(:,:) = zero_real_wp
END IF
IF (casdiags % l_gwp) THEN
  ALLOCATE ( casdiags % gwp(i_start:i_end, j_start:j_end) )
  casdiags % gwp(:,:) = zero_real_wp
END IF
!!!!





IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE allocate_diagnostic_space

SUBROUTINE deallocate_diagnostic_space()

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Local variables
CHARACTER(LEN=*), PARAMETER :: RoutineName='DEALLOCATE_DIAGNOSTIC_SPACE'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Remember to deallocate all memory in the reverse order that it was
! allocated.


!PRF water paths
IF ( ALLOCATED ( casdiags % gwp ) ) THEN
  DEALLOCATE ( casdiags % gwp )
END IF

IF ( ALLOCATED ( casdiags % swp ) ) THEN
  DEALLOCATE ( casdiags % swp )
END IF

IF ( ALLOCATED ( casdiags % iwp ) ) THEN
  DEALLOCATE ( casdiags % iwp )
END IF

IF ( ALLOCATED ( casdiags % rwp ) ) THEN
  DEALLOCATE ( casdiags % rwp )
END IF

IF ( ALLOCATED ( casdiags % lwp ) ) THEN
  DEALLOCATE ( casdiags % lwp )
END IF

IF ( ALLOCATED ( casdiags % mphys_pts ) ) THEN
  DEALLOCATE ( casdiags % mphys_pts )
END IF

IF ( ALLOCATED ( casdiags % graupfall_3d ) ) THEN
  DEALLOCATE ( casdiags % graupfall_3d )
END IF

IF ( ALLOCATED ( casdiags % snowonly_3d ) ) THEN
  DEALLOCATE ( casdiags % snowonly_3d )
END IF

IF ( ALLOCATED ( casdiags % snowfall_3d ) ) THEN
  DEALLOCATE ( casdiags % snowfall_3d )
END IF

IF ( ALLOCATED ( casdiags % rainfall_3d ) ) THEN
  DEALLOCATE ( casdiags % rainfall_3d )
END IF

IF ( ALLOCATED ( casdiags %  dqg )) THEN
  DEALLOCATE ( casdiags % dqg )
END IF

IF ( ALLOCATED ( casdiags %  dqs )) THEN
  DEALLOCATE ( casdiags % dqs )
END IF

IF ( ALLOCATED ( casdiags %  dqi )) THEN
  DEALLOCATE ( casdiags % dqi )
END IF

IF ( ALLOCATED ( casdiags %  dqr )) THEN
  DEALLOCATE ( casdiags % dqr )
END IF

IF ( ALLOCATED ( casdiags %  dqc )) THEN
  DEALLOCATE ( casdiags % dqc )
END IF

IF (casdiags % l_dqv) THEN
   IF ( ALLOCATED ( casdiags %  dqv_cond_evap )) THEN
      DEALLOCATE ( casdiags % dqv_cond_evap )
   END IF

   IF ( ALLOCATED ( casdiags %  dqv_total )) THEN
      DEALLOCATE ( casdiags % dqv_total )
   END IF
ENDIF

IF (casdiags % l_dth) THEN
   ! potential temp and mass tendencies
   IF ( ALLOCATED ( casdiags %  dth_cond_evap )) THEN
      DEALLOCATE ( casdiags % dth_cond_evap )
   END IF
   IF ( ALLOCATED ( casdiags %  dth_total )) THEN
      DEALLOCATE ( casdiags % dth_total )
   END IF
ENDIF

IF ( ALLOCATED ( casdiags % niics_i ) ) THEN
  DEALLOCATE ( casdiags % niics_i )
END IF

IF ( ALLOCATED ( casdiags % niics_s ) ) THEN
  DEALLOCATE ( casdiags % niics_s )
END IF

IF ( ALLOCATED ( casdiags % ngaci ) ) THEN
  DEALLOCATE ( casdiags % ngaci )
END IF

IF ( ALLOCATED ( casdiags % pgaci ) ) THEN
  DEALLOCATE ( casdiags % pgaci )
END IF

IF ( ALLOCATED ( casdiags % nidps ) ) THEN
  DEALLOCATE ( casdiags % nidps )
END IF

IF ( ALLOCATED ( casdiags % pidps ) ) THEN
  DEALLOCATE ( casdiags % pidps )
END IF

IF ( ALLOCATED ( casdiags % nraci_i ) ) THEN
  DEALLOCATE ( casdiags % nraci_i )
END IF

IF ( ALLOCATED ( casdiags % nraci_r ) ) THEN
  DEALLOCATE ( casdiags % nraci_r )
END IF

IF ( ALLOCATED ( casdiags % nraci_g ) ) THEN
  DEALLOCATE ( casdiags % nraci_g )
END IF

IF ( ALLOCATED ( casdiags % praci_i ) ) THEN
  DEALLOCATE ( casdiags % praci_i )
END IF

IF ( ALLOCATED ( casdiags % praci_r ) ) THEN
  DEALLOCATE ( casdiags % praci_r )
END IF

IF ( ALLOCATED ( casdiags % praci_g ) ) THEN
  DEALLOCATE ( casdiags % praci_g )
END IF

IF ( ALLOCATED ( casdiags % pihal ) ) THEN
  DEALLOCATE ( casdiags % pihal )
END IF

IF ( ALLOCATED ( casdiags % ngmlt ) ) THEN
  DEALLOCATE ( casdiags % ngmlt )
END IF

IF ( ALLOCATED ( casdiags % ngacw ) ) THEN
  DEALLOCATE ( casdiags % ngacw )
END IF

IF ( ALLOCATED ( casdiags % nsmlt ) ) THEN
  DEALLOCATE ( casdiags % nsmlt )
END IF

IF ( ALLOCATED ( casdiags % nimlt ) ) THEN
  DEALLOCATE ( casdiags % nimlt )
END IF

IF ( ALLOCATED ( casdiags % nsacr ) ) THEN
  DEALLOCATE ( casdiags % nsacr )
END IF

IF ( ALLOCATED ( casdiags % nsacw ) ) THEN
  DEALLOCATE ( casdiags % nsacw )
END IF

IF ( ALLOCATED ( casdiags % niacw ) ) THEN
  DEALLOCATE ( casdiags % niacw )
END IF

IF ( ALLOCATED ( casdiags % ngsub ) ) THEN
  DEALLOCATE ( casdiags % ngsub )
END IF

IF ( ALLOCATED ( casdiags % ngacs ) ) THEN
  DEALLOCATE ( casdiags % ngacs )
END IF

IF ( ALLOCATED ( casdiags % nsaci ) ) THEN
  DEALLOCATE ( casdiags % nsaci )
END IF

IF ( ALLOCATED ( casdiags % nsaut ) ) THEN
  DEALLOCATE ( casdiags % nsaut )
END IF

IF ( ALLOCATED ( casdiags % nssub ) ) THEN
  DEALLOCATE ( casdiags % nssub )
END IF

IF ( ALLOCATED ( casdiags % nisub ) ) THEN
  DEALLOCATE ( casdiags % nisub )
END IF

IF ( ALLOCATED ( casdiags % nrevp ) ) THEN
  DEALLOCATE ( casdiags % nrevp )
END IF

IF ( ALLOCATED ( casdiags % nsedr ) ) THEN
  DEALLOCATE ( casdiags % nsedr )
END IF

IF ( ALLOCATED ( casdiags % nracr ) ) THEN
  DEALLOCATE ( casdiags % nracr )
END IF

IF ( ALLOCATED ( casdiags % nracw ) ) THEN
  DEALLOCATE ( casdiags % nracw )
END IF

IF ( ALLOCATED ( casdiags % nsedl ) ) THEN
  DEALLOCATE ( casdiags % nsedl )
END IF

IF ( ALLOCATED ( casdiags % nraut ) ) THEN
  DEALLOCATE ( casdiags % nraut )
END IF

IF ( ALLOCATED ( casdiags % nihal ) ) THEN
  DEALLOCATE ( casdiags % nihal )
END IF

IF ( ALLOCATED ( casdiags % nhomr ) ) THEN
  DEALLOCATE ( casdiags % nhomr )
END IF

IF ( ALLOCATED ( casdiags % phomr ) ) THEN
  DEALLOCATE ( casdiags % phomr )
END IF

IF ( ALLOCATED ( casdiags % pcond ) ) THEN
  DEALLOCATE ( casdiags % pcond )
END IF

IF ( ALLOCATED ( casdiags % psedl ) ) THEN
  DEALLOCATE ( casdiags % psedl )
END IF

IF ( ALLOCATED ( casdiags % nsedg ) ) THEN
  DEALLOCATE ( casdiags % nsedg )
END IF

IF ( ALLOCATED ( casdiags % psedg ) ) THEN
  DEALLOCATE ( casdiags % psedg )
END IF

IF ( ALLOCATED ( casdiags % psedr ) ) THEN
  DEALLOCATE ( casdiags % psedr )
END IF

IF ( ALLOCATED ( casdiags % nseds ) ) THEN
  DEALLOCATE ( casdiags % nseds )
END IF

IF ( ALLOCATED ( casdiags % pseds ) ) THEN
  DEALLOCATE ( casdiags % pseds )
END IF

IF ( ALLOCATED ( casdiags % nsedi ) ) THEN
  DEALLOCATE ( casdiags % nsedi )
END IF

IF ( ALLOCATED ( casdiags % psedi ) ) THEN
  DEALLOCATE ( casdiags % psedi )
END IF

IF ( ALLOCATED ( casdiags % pgsub ) ) THEN
  DEALLOCATE ( casdiags % pgsub )
END IF

IF ( ALLOCATED ( casdiags % pgmlt ) ) THEN
  DEALLOCATE( casdiags % pgmlt )
END IF

IF ( ALLOCATED ( casdiags % pgacs ) ) THEN
  DEALLOCATE ( casdiags % pgacs )
END IF

IF ( ALLOCATED ( casdiags % pgacw ) ) THEN
  DEALLOCATE ( casdiags % pgacw )
END IF

IF ( ALLOCATED ( casdiags % prevp ) ) THEN
  DEALLOCATE ( casdiags % prevp )
END IF

IF ( ALLOCATED ( casdiags % pracr ) ) THEN
  DEALLOCATE (casdiags % pracr )
END IF

IF ( ALLOCATED ( casdiags % pracw ) ) THEN
  DEALLOCATE (casdiags % pracw )
END IF

IF ( ALLOCATED ( casdiags % praut ) ) THEN
  DEALLOCATE ( casdiags % praut )
END IF

IF ( ALLOCATED ( casdiags % psaci ) ) THEN
  DEALLOCATE( casdiags % psaci )
END IF

IF ( ALLOCATED ( casdiags % psaut ) ) THEN
  DEALLOCATE ( casdiags % psaut )
END IF

IF ( ALLOCATED ( casdiags % psmlt ) ) THEN
  DEALLOCATE ( casdiags % psmlt )
END IF

IF ( ALLOCATED ( casdiags % pimlt ) ) THEN
  DEALLOCATE ( casdiags % pimlt )
END IF

IF ( ALLOCATED ( casdiags % pssub ) ) THEN
  DEALLOCATE ( casdiags % pssub )
END IF

IF ( ALLOCATED ( casdiags % pisub ) ) THEN
  DEALLOCATE ( casdiags % pisub )
END IF

IF ( ALLOCATED ( casdiags % psacr ) ) THEN
  DEALLOCATE ( casdiags % psacr )
END IF

IF ( ALLOCATED ( casdiags % psacw ) ) THEN
  DEALLOCATE ( casdiags % psacw )
END IF

IF ( ALLOCATED ( casdiags % piacw ) ) THEN
  DEALLOCATE ( casdiags % piacw )
END IF

IF ( ALLOCATED ( casdiags % psdep ) ) THEN
  DEALLOCATE ( casdiags % psdep )
END IF

IF ( ALLOCATED ( casdiags % pidep ) ) THEN
  DEALLOCATE ( casdiags % pidep )
END IF

IF ( ALLOCATED ( casdiags % ninuc ) ) THEN
  DEALLOCATE ( casdiags % ninuc )
END IF

IF ( ALLOCATED ( casdiags % pinuc ) ) THEN
  DEALLOCATE ( casdiags % pinuc )
  casdiags % l_pinuc = .FALSE.
END IF

IF ( ALLOCATED ( casdiags % nhomc ) ) THEN
  DEALLOCATE ( casdiags % nhomc )
END IF

IF ( ALLOCATED ( casdiags % phomc ) ) THEN
  DEALLOCATE ( casdiags % phomc )
END IF

IF (casdiags % l_radar) THEN
   IF ( ALLOCATED ( casdiags % dbz_r ) ) THEN
      DEALLOCATE( casdiags % dbz_r )
   END IF
   
   IF ( ALLOCATED ( casdiags % dbz_l ) ) THEN
      DEALLOCATE ( casdiags % dbz_l )
   END IF
   
   IF ( ALLOCATED ( casdiags % dbz_s ) ) THEN
      DEALLOCATE ( casdiags % dbz_s )
   END IF
   
   IF ( ALLOCATED ( casdiags % dbz_i ) ) THEN
      DEALLOCATE ( casdiags % dbz_i )
   END IF
   
   IF ( ALLOCATED ( casdiags % dbz_g ) ) THEN
      DEALLOCATE ( casdiags % dbz_g )
   END IF
   
   IF ( ALLOCATED ( casdiags % dbz_tot ) ) THEN
      DEALLOCATE ( casdiags % dbz_tot )
   END IF
ENDIF

IF ( ALLOCATED ( casdiags %  SurfaceGraupR )) THEN
   DEALLOCATE ( casdiags % SurfaceGraupR )
END IF

IF ( ALLOCATED ( casdiags %  SurfaceSnowR )) THEN
   DEALLOCATE ( casdiags % SurfaceSnowR )
END IF

IF ( ALLOCATED ( casdiags %  SurfaceRainR )) THEN
   DEALLOCATE ( casdiags % SurfaceRainR )
END IF

IF ( ALLOCATED ( casdiags %  SurfaceCloudR )) THEN
   DEALLOCATE ( casdiags % SurfaceCloudR )
END IF

IF ( ALLOCATED ( casdiags %  precip )) THEN
   DEALLOCATE ( casdiags % precip )
END IF
!!

! Set to False all switches which affect groups of more than one diagnostic
casdiags % l_process_rates = .FALSE.
casdiags % l_tendency_dg   = .FALSE.
casdiags % l_radar         = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE deallocate_diagnostic_space

END MODULE generic_diagnostic_variables
