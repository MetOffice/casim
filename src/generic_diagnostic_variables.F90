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

SUBROUTINE allocate_diagnostic_space(is, ie, js, je, ks, ke)

USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
USE mphys_parameters, ONLY: zero_real_wp

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: is ! Start of i array
INTEGER, INTENT(IN) :: ie ! End of i array
INTEGER, INTENT(IN) :: js ! Start of j array
INTEGER, INTENT(IN) :: je ! End of j array
INTEGER, INTENT(IN) :: ks ! Start of k array
INTEGER, INTENT(IN) :: ke ! End of k array

! Local variables
CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_DIAGNOSTIC_SPACE'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( casdiags % l_precip) THEN

  ALLOCATE ( casdiags % precip(is:ie, js:je) )
  casdiags % precip(:,:) = zero_real_wp

END IF ! casdiags % l_precip

IF ( casdiags % l_surface_cloud ) THEN

  ALLOCATE ( casdiags % SurfaceCloudR(is:ie, js:je) )
  casdiags % SurfaceCloudR(:,:) = zero_real_wp

END IF

IF ( casdiags % l_surface_rain ) THEN

  ALLOCATE ( casdiags % SurfaceRainR(is:ie, js:je) )
  casdiags % SurfaceRainR(:,:) = zero_real_wp

END IF

IF ( casdiags % l_surface_snow ) THEN

  ALLOCATE ( casdiags % SurfaceSnowR(is:ie, js:je) )
  casdiags % SurfaceSnowR(:,:) = zero_real_wp

END IF

IF ( casdiags % l_surface_graup ) THEN

  ALLOCATE ( casdiags % SurfaceGraupR(is:ie, js:je) )
  casdiags % SurfaceGraupR(:,:) = zero_real_wp

END IF

IF ( casdiags % l_radar ) THEN

  ! A single logical allocates all variables
  ALLOCATE ( casdiags % dbz_tot(is:ie, js:je, ks:ke) )
  ALLOCATE ( casdiags % dbz_g(is:ie, js:je, ks:ke) )
  ALLOCATE ( casdiags % dbz_i(is:ie, js:je, ks:ke) )
  ALLOCATE ( casdiags % dbz_s(is:ie, js:je, ks:ke) )
  ALLOCATE ( casdiags % dbz_l(is:ie, js:je, ks:ke) )
  ALLOCATE ( casdiags % dbz_r(is:ie, js:je, ks:ke) )
  
  casdiags % dbz_tot(:,:,:) = ref_lim
  casdiags % dbz_g(:,:,:)   = ref_lim
  casdiags % dbz_i(:,:,:)   = ref_lim
  casdiags % dbz_s(:,:,:)   = ref_lim
  casdiags % dbz_l(:,:,:)   = ref_lim
  casdiags % dbz_r(:,:,:)   = ref_lim

END IF ! casdiags % l_radar

IF (casdiags % l_phomc) THEN
  ALLOCATE ( casdiags % phomc(is:ie, js:je, ks:ke) )
  casdiags % phomc(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nhomc) THEN
  ALLOCATE ( casdiags % nhomc(is:ie, js:je, ks:ke) )
  casdiags % nhomc(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pinuc) THEN
  ALLOCATE ( casdiags % pinuc(is:ie, js:je, ks:ke) )
  casdiags % pinuc(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_ninuc) THEN
  ALLOCATE ( casdiags % ninuc(is:ie, js:je, ks:ke) )
  casdiags % ninuc(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pidep) THEN
  ALLOCATE ( casdiags % pidep(is:ie, js:je, ks:ke) )
  casdiags % pidep(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psdep) THEN
  ALLOCATE ( casdiags % psdep(is:ie, js:je, ks:ke) )
  casdiags % psdep(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_piacw) THEN
  ALLOCATE ( casdiags % piacw(is:ie, js:je, ks:ke) )
  casdiags % piacw(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psacw) THEN
  ALLOCATE ( casdiags % psacw(is:ie, js:je, ks:ke) )
  casdiags % psacw(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psacr) THEN
  ALLOCATE ( casdiags % psacr(is:ie, js:je, ks:ke) )
  casdiags % psacr(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pisub) THEN
  ALLOCATE ( casdiags % pisub(is:ie, js:je, ks:ke) )
  casdiags % pisub(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pssub) THEN
  ALLOCATE ( casdiags % pssub(is:ie, js:je, ks:ke) )
  casdiags % pssub(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pimlt) THEN
  ALLOCATE ( casdiags % pimlt(is:ie, js:je, ks:ke) )
  casdiags % pimlt(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psmlt) THEN
  ALLOCATE ( casdiags % psmlt(is:ie, js:je, ks:ke) )
  casdiags % psmlt(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psaut) THEN
  ALLOCATE ( casdiags % psaut(is:ie, js:je, ks:ke) )
  casdiags % psaut(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psaci) THEN
  ALLOCATE ( casdiags % psaci(is:ie, js:je, ks:ke) )
  casdiags % psaci(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_praut) THEN
  ALLOCATE ( casdiags % praut(is:ie, js:je, ks:ke) )
  casdiags % praut(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pracw) THEN
  ALLOCATE ( casdiags % pracw(is:ie, js:je, ks:ke) )
  casdiags % pracw(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pracr) THEN
  ALLOCATE ( casdiags % pracr(is:ie, js:je, ks:ke) )
  casdiags % pracr(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_prevp) THEN
  ALLOCATE ( casdiags % prevp(is:ie, js:je, ks:ke) )
  casdiags % prevp(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pgacw) THEN
  ALLOCATE ( casdiags % pgacw(is:ie, js:je, ks:ke) )
  casdiags % pgacw(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pgacs) THEN
  ALLOCATE ( casdiags % pgacs(is:ie, js:je, ks:ke) )
  casdiags % pgacs(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pgmlt) THEN
  ALLOCATE ( casdiags % pgmlt(is:ie, js:je, ks:ke) )
  casdiags % pgmlt(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pgsub) THEN
  ALLOCATE ( casdiags % pgsub(is:ie, js:je, ks:ke) )
  casdiags % pgsub(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psedi) THEN
  ALLOCATE ( casdiags % psedi(is:ie, js:je, ks:ke) )
  casdiags % psedi(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsedi) THEN
  ALLOCATE ( casdiags % nsedi(is:ie, js:je, ks:ke) )
  casdiags % nsedi(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pseds) THEN
  ALLOCATE ( casdiags % pseds(is:ie, js:je, ks:ke) )
  casdiags % pseds(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nseds) THEN
  ALLOCATE ( casdiags % nseds(is:ie, js:je, ks:ke) )
  casdiags % nseds(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psedr) THEN
  ALLOCATE ( casdiags % psedr(is:ie, js:je, ks:ke) )
  casdiags % psedr(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psedg) THEN
  ALLOCATE ( casdiags % psedg(is:ie, js:je, ks:ke) )
  casdiags % psedg(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nsedg) THEN
  ALLOCATE ( casdiags % nsedg(is:ie, js:je, ks:ke) )
  casdiags % nsedg(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_psedl) THEN
  ALLOCATE ( casdiags % psedl(is:ie, js:je, ks:ke) )
  casdiags % psedl(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_pcond) THEN
  ALLOCATE ( casdiags % pcond(is:ie, js:je, ks:ke) )
  casdiags % pcond(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_phomr) THEN
  ALLOCATE ( casdiags % phomr(is:ie, js:je, ks:ke) )
  casdiags % phomr(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nhomr) THEN
  ALLOCATE ( casdiags % nhomr(is:ie, js:je, ks:ke) )
  casdiags % nhomr(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

IF (casdiags % l_nihal) THEN
  ALLOCATE ( casdiags % nihal(is:ie, js:je, ks:ke) )
  casdiags % nihal(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

! potential temp and mass tendencies
IF (casdiags % l_dth) THEN
  ALLOCATE ( casdiags % dth_total(is:ie, js:je, ks:ke) )
  casdiags % dth_total(:,:,:) = zero_real_wp
  ALLOCATE ( casdiags % dth_cond_evap(is:ie, js:je, ks:ke) )
  casdiags % dth_cond_evap(:,:,:) = zero_real_wp
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqv) THEN
  ALLOCATE ( casdiags % dqv_total(is:ie, js:je, ks:ke) )
  casdiags % dqv_total(:,:,:) = zero_real_wp
  ALLOCATE ( casdiags % dqv_cond_evap(is:ie, js:je, ks:ke) )
  casdiags % dqv_cond_evap(:,:,:) = zero_real_wp
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqc) THEN
  ALLOCATE ( casdiags % dqc(is:ie, js:je, ks:ke) )
  casdiags % dqc(:,:,:) = zero_real_wp
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqr) THEN
  ALLOCATE ( casdiags % dqr(is:ie, js:je, ks:ke) )
  casdiags % dqr(:,:,:) = zero_real_wp
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqi) THEN
  ALLOCATE ( casdiags % dqi(is:ie, js:je, ks:ke) )
  casdiags % dqi(:,:,:) = zero_real_wp
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqs) THEN
  ALLOCATE ( casdiags % dqs(is:ie, js:je, ks:ke) )
  casdiags % dqs(:,:,:) = zero_real_wp
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_dqg) THEN
  ALLOCATE ( casdiags % dqg(is:ie, js:je, ks:ke) )
  casdiags % dqg(:,:,:) = zero_real_wp
  casdiags % l_tendency_dg = .TRUE.
END IF

IF (casdiags % l_rainfall_3d) THEN
  ALLOCATE ( casdiags % rainfall_3d(is:ie, js:je, ks:ke) )
  casdiags % rainfall_3d(:,:,:) = zero_real_wp
END IF

IF (casdiags % l_snowfall_3d) THEN
  ALLOCATE ( casdiags % snowfall_3d(is:ie, js:je, ks:ke) )
  casdiags % snowfall_3d(:,:,:) = zero_real_wp
END IF

IF (casdiags % l_snowonly_3d) THEN
  ALLOCATE ( casdiags % snowonly_3d(is:ie, js:je, ks:ke) )
  casdiags % snowonly_3d(:,:,:) = zero_real_wp
END IF

IF (casdiags % l_graupfall_3d) THEN
  ALLOCATE ( casdiags % graupfall_3d(is:ie, js:je, ks:ke) )
  casdiags % graupfall_3d(:,:,:) = zero_real_wp
END IF

IF (casdiags % l_mphys_pts) THEN
  ALLOCATE ( casdiags % mphys_pts(is:ie, js:je, ks:ke) )
  casdiags % mphys_pts(:,:,:) = .FALSE.
END IF

!PRF water paths
IF (casdiags % l_lwp) THEN
  ALLOCATE ( casdiags % lwp(is:ie, js:je) )
  casdiags % lwp(:,:) = zero_real_wp
END IF
IF (casdiags % l_rwp) THEN
  ALLOCATE ( casdiags % rwp(is:ie, js:je) )
  casdiags % rwp(:,:) = zero_real_wp
END IF
IF (casdiags % l_iwp) THEN
  ALLOCATE ( casdiags % iwp(is:ie, js:je) )
  casdiags % iwp(:,:) = zero_real_wp
END IF
IF (casdiags % l_swp) THEN
  ALLOCATE ( casdiags % swp(is:ie, js:je) )
  casdiags % swp(:,:) = zero_real_wp
END IF
IF (casdiags % l_gwp) THEN
  ALLOCATE ( casdiags % gwp(is:ie, js:je) )
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

IF ( ALLOCATED ( casdiags %  precip )) THEN
   DEALLOCATE ( casdiags % precip )
END IF

IF (casdiags % l_dth) THEN
   ! potential temp and mass tendencies
   IF ( ALLOCATED ( casdiags %  dth_total )) THEN
      DEALLOCATE ( casdiags % dth_total )
   END IF

   IF ( ALLOCATED ( casdiags %  dth_cond_evap )) THEN
      DEALLOCATE ( casdiags % dth_cond_evap )
   END IF
ENDIF

IF (casdiags % l_dqv) THEN
   IF ( ALLOCATED ( casdiags %  dqv_total )) THEN
      DEALLOCATE ( casdiags % dqv_total )
   END IF

   IF ( ALLOCATED ( casdiags %  dqv_cond_evap )) THEN
      DEALLOCATE ( casdiags % dqv_cond_evap )
   END IF
ENDIF

IF ( ALLOCATED ( casdiags %  dqc )) THEN
  DEALLOCATE ( casdiags % dqc )
END IF

IF ( ALLOCATED ( casdiags %  dqr )) THEN
  DEALLOCATE ( casdiags % dqr )
END IF

IF ( ALLOCATED ( casdiags %  dqi )) THEN
  DEALLOCATE ( casdiags % dqi )
END IF

IF ( ALLOCATED ( casdiags %  dqs )) THEN
  DEALLOCATE ( casdiags % dqs )
END IF

IF ( ALLOCATED ( casdiags %  dqg )) THEN
  DEALLOCATE ( casdiags % dqg )
END IF


!PRF water paths
IF ( ALLOCATED ( casdiags % lwp ) ) THEN
  DEALLOCATE ( casdiags % lwp )
END IF
IF ( ALLOCATED ( casdiags % rwp ) ) THEN
  DEALLOCATE ( casdiags % rwp )
END IF
IF ( ALLOCATED ( casdiags % iwp ) ) THEN
  DEALLOCATE ( casdiags % iwp )
END IF
IF ( ALLOCATED ( casdiags % swp ) ) THEN
  DEALLOCATE ( casdiags % swp )
END IF
IF ( ALLOCATED ( casdiags % gwp ) ) THEN
  DEALLOCATE ( casdiags % gwp )
END IF

!!



! Set to False all switches which affect groups of more than one diagnostic
casdiags % l_process_rates = .FALSE.
casdiags % l_tendency_dg   = .FALSE.
casdiags % l_radar         = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE deallocate_diagnostic_space

END MODULE generic_diagnostic_variables
