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

  !-----------------------------------
  ! aerosol STASH test
  !-----------------------------------
  ! Aerosol codes
  !  am1     Aitken Sol Mass
  !  an1     Aitken Sol Number
  !  am2     Accum Sol Mass
  !  an2     Accum Sol Number
  !  am3     Coarse Sol Mass
  !  an3     Coarse Sol Number
  !  am4     Act Sol Liq in casim
  !  am5     Act Sol Rain in casim  <-this is not currently available
  !  am6     Coarse Dust Mass
  !  an6     Coarse Dust Number
  !  am7     Act Insol Ice in casim
  !  am8     Act Sol Ice in casim
  !  am9     Act Inso lLiq in casim
  !  am10    Accum Dust Mass
  !  an10    Accum Dust Number
  !  an11    Act Sol Number in casim
  !  an12    Act Insol Number in casim <- dust in hydrometeors ice and liquid
  !  ak1     Aitken Sol Bk  <-for ukca coupling
  !  ak2     Accum Sol Bk   <-for ukca coupling
  !  ak3     Coarse Sol Bk  <-for ukca coupling
  !-----------------------------------
  LOGICAL :: l_aact_am1       = .FALSE. !601
  LOGICAL :: l_aact_an1       = .FALSE. !602
  LOGICAL :: l_aact_am2       = .FALSE. !603
  LOGICAL :: l_aact_an2       = .FALSE. !604
  LOGICAL :: l_aact_am3       = .FALSE. !605
  LOGICAL :: l_aact_an3       = .FALSE. !606
  LOGICAL :: l_aact_am9       = .FALSE. !607
  LOGICAL :: l_aact_an6       = .FALSE. !608
  LOGICAL :: l_aaut           = .FALSE. !609
  LOGICAL :: l_aacw           = .FALSE. !610
  LOGICAL :: l_asedr_am       = .FALSE. !611
  LOGICAL :: l_asedr_an11     = .FALSE. !612
  LOGICAL :: l_asedr_an12     = .FALSE. !613
  LOGICAL :: l_arevp_am2      = .FALSE. !614
  LOGICAL :: l_arevp_an2      = .FALSE. !615
  LOGICAL :: l_arevp_am3      = .FALSE. !616
  LOGICAL :: l_arevp_an3      = .FALSE. !617
  LOGICAL :: l_arevp_am4      = .FALSE. !618
  LOGICAL :: l_arevp_am5      = .FALSE. !619
  LOGICAL :: l_arevp_am6      = .FALSE. !620
  LOGICAL :: l_arevp_an6      = .FALSE. !621
  LOGICAL :: l_asedl_am4      = .FALSE. !622
  LOGICAL :: l_asedl_an11     = .FALSE. !623
  LOGICAL :: l_asedl_an12     = .FALSE. !624
  LOGICAL :: l_dnuc_am8       = .FALSE. !625
  LOGICAL :: l_dnuc_am6       = .FALSE. !626
  LOGICAL :: l_dnuc_am9       = .FALSE. !627
  LOGICAL :: l_dnuc_an6       = .FALSE. !628
  LOGICAL :: l_dsub_am2       = .FALSE. !629
  LOGICAL :: l_dsub_an2       = .FALSE. !630
  LOGICAL :: l_dsub_am6       = .FALSE. !631
  LOGICAL :: l_dsub_an6       = .FALSE. !632
  LOGICAL :: l_dsedi_am7      = .FALSE. !633
  LOGICAL :: l_dsedi_am8      = .FALSE. !634
  LOGICAL :: l_dsedi_an11     = .FALSE. !635
  LOGICAL :: l_dsedi_an12     = .FALSE. !636
  LOGICAL :: l_dseds_am7      = .FALSE. !637
  LOGICAL :: l_dseds_am8      = .FALSE. !638
  LOGICAL :: l_dseds_an11     = .FALSE. !639
  LOGICAL :: l_dseds_an12     = .FALSE. !640
  LOGICAL :: l_dsedg_am7      = .FALSE. !641
  LOGICAL :: l_dsedg_am8      = .FALSE. !642
  LOGICAL :: l_dsedg_an11     = .FALSE. !643
  LOGICAL :: l_dsedg_an12     = .FALSE. !644
  LOGICAL :: l_dssub_am2      = .FALSE. !645
  LOGICAL :: l_dssub_an2      = .FALSE. !646
  LOGICAL :: l_dssub_am6      = .FALSE. !647
  LOGICAL :: l_dssub_an6      = .FALSE. !648
  LOGICAL :: l_dgsub_am2      = .FALSE. !649
  LOGICAL :: l_dgsub_an2      = .FALSE. !650
  LOGICAL :: l_dgsub_am6      = .FALSE. !651
  LOGICAL :: l_dgsub_an6      = .FALSE. !652
  LOGICAL :: l_dhomc_am8      = .FALSE. !653
  LOGICAL :: l_dhomc_am7      = .FALSE. !654
  LOGICAL :: l_dhomr_am8      = .FALSE. !655
  LOGICAL :: l_dhomr_am7      = .FALSE. !656
  LOGICAL :: l_dimlt_am4      = .FALSE. !657
  LOGICAL :: l_dimlt_am9      = .FALSE. !658
  LOGICAL :: l_dsmlt_am4      = .FALSE. !659
  LOGICAL :: l_dsmlt_am9      = .FALSE. !660
  LOGICAL :: l_dgmlt_am4      = .FALSE. !661
  LOGICAL :: l_dgmlt_am9      = .FALSE. !662
  LOGICAL :: l_diacw_am8      = .FALSE. !663
  LOGICAL :: l_diacw_am7      = .FALSE. !664
  LOGICAL :: l_dsacw_am8      = .FALSE. !665
  LOGICAL :: l_dsacw_am7      = .FALSE. !666
  LOGICAL :: l_dgacw_am8      = .FALSE. !667
  LOGICAL :: l_dgacw_am7      = .FALSE. !668
  LOGICAL :: l_dsacr_am8      = .FALSE. !669
  LOGICAL :: l_dsacr_am7      = .FALSE. !670
  LOGICAL :: l_dgacr_am8      = .FALSE. !671
  LOGICAL :: l_dgacr_am7      = .FALSE. !672
  LOGICAL :: l_draci_am8      = .FALSE. !673
  LOGICAL :: l_draci_am7      = .FALSE. !674
  LOGICAL :: l_asedl_am9      = .FALSE. !675
  LOGICAL :: l_asedr_am9      = .FALSE. !676


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

  !-------------------------------------------
  ! aerosol
  !-------------------------------------------
  REAL, ALLOCATABLE :: aact_am1(:,:,:)   !601
  REAL, ALLOCATABLE :: aact_an1(:,:,:)   !602
  REAL, ALLOCATABLE :: aact_am2(:,:,:)   !603
  REAL, ALLOCATABLE :: aact_an2(:,:,:)   !604
  REAL, ALLOCATABLE :: aact_am3(:,:,:)   !605
  REAL, ALLOCATABLE :: aact_an3(:,:,:)   !606
  REAL, ALLOCATABLE :: aact_am9(:,:,:)   !607
  REAL, ALLOCATABLE :: aact_an6(:,:,:)   !608
  REAL, ALLOCATABLE :: aaut(:,:,:)       !609
  REAL, ALLOCATABLE :: aacw(:,:,:)       !610
  REAL, ALLOCATABLE :: asedr_am(:,:,:)   !611
  REAL, ALLOCATABLE :: asedr_an11(:,:,:) !612
  REAL, ALLOCATABLE :: asedr_an12(:,:,:) !613
  REAL, ALLOCATABLE :: arevp_am2(:,:,:)  !614
  REAL, ALLOCATABLE :: arevp_an2(:,:,:)  !615
  REAL, ALLOCATABLE :: arevp_am3(:,:,:)  !616
  REAL, ALLOCATABLE :: arevp_an3(:,:,:)  !617
  REAL, ALLOCATABLE :: arevp_am4(:,:,:)  !618
  REAL, ALLOCATABLE :: arevp_am5(:,:,:)  !619
  REAL, ALLOCATABLE :: arevp_am6(:,:,:)  !620
  REAL, ALLOCATABLE :: arevp_an6(:,:,:)  !621
  REAL, ALLOCATABLE :: asedl_am4(:,:,:)  !622
  REAL, ALLOCATABLE :: asedl_an11(:,:,:) !623
  REAL, ALLOCATABLE :: asedl_an12(:,:,:) !624
  REAL, ALLOCATABLE :: dnuc_am8(:,:,:)   !625
  REAL, ALLOCATABLE :: dnuc_am6(:,:,:)   !626
  REAL, ALLOCATABLE :: dnuc_am9(:,:,:)   !627
  REAL, ALLOCATABLE :: dnuc_an6(:,:,:)   !628
  REAL, ALLOCATABLE :: dsub_am2(:,:,:)   !629
  REAL, ALLOCATABLE :: dsub_an2(:,:,:)   !630
  REAL, ALLOCATABLE :: dsub_am6(:,:,:)   !631
  REAL, ALLOCATABLE :: dsub_an6(:,:,:)   !632
  REAL, ALLOCATABLE :: dsedi_am7(:,:,:)  !633
  REAL, ALLOCATABLE :: dsedi_am8(:,:,:)  !634
  REAL, ALLOCATABLE :: dsedi_an11(:,:,:) !635
  REAL, ALLOCATABLE :: dsedi_an12(:,:,:) !636
  REAL, ALLOCATABLE :: dseds_am7(:,:,:)  !637
  REAL, ALLOCATABLE :: dseds_am8(:,:,:)  !638
  REAL, ALLOCATABLE :: dseds_an11(:,:,:) !639
  REAL, ALLOCATABLE :: dseds_an12(:,:,:) !640
  REAL, ALLOCATABLE :: dsedg_am7(:,:,:)  !641
  REAL, ALLOCATABLE :: dsedg_am8(:,:,:)  !642
  REAL, ALLOCATABLE :: dsedg_an11(:,:,:) !643
  REAL, ALLOCATABLE :: dsedg_an12(:,:,:) !644
  REAL, ALLOCATABLE :: dssub_am2(:,:,:)  !645
  REAL, ALLOCATABLE :: dssub_an2(:,:,:)  !646
  REAL, ALLOCATABLE :: dssub_am6(:,:,:)  !647
  REAL, ALLOCATABLE :: dssub_an6(:,:,:)  !648
  REAL, ALLOCATABLE :: dgsub_am2(:,:,:)  !649
  REAL, ALLOCATABLE :: dgsub_an2(:,:,:)  !650
  REAL, ALLOCATABLE :: dgsub_am6(:,:,:)  !651
  REAL, ALLOCATABLE :: dgsub_an6(:,:,:)  !652
  REAL, ALLOCATABLE :: dhomc_am8(:,:,:)  !653
  REAL, ALLOCATABLE :: dhomc_am7(:,:,:)  !654
  REAL, ALLOCATABLE :: dhomr_am8(:,:,:)  !655
  REAL, ALLOCATABLE :: dhomr_am7(:,:,:)  !656
  REAL, ALLOCATABLE :: dimlt_am4(:,:,:)  !657
  REAL, ALLOCATABLE :: dimlt_am9(:,:,:)  !658
  REAL, ALLOCATABLE :: dsmlt_am4(:,:,:)  !659
  REAL, ALLOCATABLE :: dsmlt_am9(:,:,:)  !660
  REAL, ALLOCATABLE :: dgmlt_am4(:,:,:)  !661
  REAL, ALLOCATABLE :: dgmlt_am9(:,:,:)  !662
  REAL, ALLOCATABLE :: diacw_am8(:,:,:)  !663
  REAL, ALLOCATABLE :: diacw_am7(:,:,:)  !664
  REAL, ALLOCATABLE :: dsacw_am8(:,:,:)  !665
  REAL, ALLOCATABLE :: dsacw_am7(:,:,:)  !666
  REAL, ALLOCATABLE :: dgacw_am8(:,:,:)  !667
  REAL, ALLOCATABLE :: dgacw_am7(:,:,:)  !668
  REAL, ALLOCATABLE :: dsacr_am8(:,:,:)  !669
  REAL, ALLOCATABLE :: dsacr_am7(:,:,:)  !670
  REAL, ALLOCATABLE :: dgacr_am8(:,:,:)  !671
  REAL, ALLOCATABLE :: dgacr_am7(:,:,:)  !672
  REAL, ALLOCATABLE :: draci_am8(:,:,:)  !673
  REAL, ALLOCATABLE :: draci_am7(:,:,:)  !674
  REAL, ALLOCATABLE :: asedl_am9(:,:,:)  !675
  REAL, ALLOCATABLE :: asedr_am9(:,:,:)  !676

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

!-----------------------------------------------------
! aerosol testing
!-----------------------------------------------------

!601
IF (casdiags % l_aact_am1) THEN
  ALLOCATE ( casdiags % aact_am1(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aact_am1(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!602
IF (casdiags % l_aact_an1) THEN
  ALLOCATE ( casdiags % aact_an1(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aact_an1(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!603
IF (casdiags % l_aact_am2) THEN
  ALLOCATE ( casdiags % aact_am2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aact_am2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!604
IF (casdiags % l_aact_an2) THEN
  ALLOCATE ( casdiags % aact_an2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aact_an2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!605
IF (casdiags % l_aact_am3) THEN
  ALLOCATE ( casdiags % aact_am3(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aact_am3(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!606
IF (casdiags % l_aact_an3) THEN
  ALLOCATE ( casdiags % aact_an3(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aact_an3(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!607
IF (casdiags % l_aact_am9) THEN
  ALLOCATE ( casdiags % aact_am9(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aact_am9(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!608
IF (casdiags % l_aact_an6) THEN
  ALLOCATE ( casdiags % aact_an6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aact_an6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!609
IF (casdiags % l_aaut) THEN
  ALLOCATE ( casdiags % aaut(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aaut(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!610
IF (casdiags % l_aacw) THEN
  ALLOCATE ( casdiags % aacw(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % aacw(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!611
IF (casdiags % l_asedr_am) THEN
  ALLOCATE ( casdiags % asedr_am(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % asedr_am(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!612
IF (casdiags % l_asedr_an11) THEN
  ALLOCATE ( casdiags % asedr_an11(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % asedr_an11(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!613
IF (casdiags % l_asedr_an12) THEN
  ALLOCATE ( casdiags % asedr_an12(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % asedr_an12(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!614
IF (casdiags % l_arevp_am2) THEN
  ALLOCATE ( casdiags % arevp_am2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % arevp_am2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!615
IF (casdiags % l_arevp_an2) THEN
  ALLOCATE ( casdiags % arevp_an2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % arevp_an2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!616
IF (casdiags % l_arevp_am3) THEN
  ALLOCATE ( casdiags % arevp_am3(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % arevp_am3(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!617
IF (casdiags % l_arevp_an3) THEN
  ALLOCATE ( casdiags % arevp_an3(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % arevp_an3(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!618
IF (casdiags % l_arevp_am4) THEN
  ALLOCATE ( casdiags % arevp_am4(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % arevp_am4(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!619
IF (casdiags % l_arevp_am5) THEN
  ALLOCATE ( casdiags % arevp_am5(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % arevp_am5(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!620
IF (casdiags % l_arevp_am6) THEN
  ALLOCATE ( casdiags % arevp_am6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % arevp_am6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!621
IF (casdiags % l_arevp_an6) THEN
  ALLOCATE ( casdiags % arevp_an6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % arevp_an6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!622
IF (casdiags % l_asedl_am4) THEN
  ALLOCATE ( casdiags % asedl_am4(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % asedl_am4(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!623
IF (casdiags % l_asedl_an11) THEN
  ALLOCATE ( casdiags % asedl_an11(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % asedl_an11(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!624
IF (casdiags % l_asedl_an12) THEN
  ALLOCATE ( casdiags % asedl_an12(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % asedl_an12(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!625
IF (casdiags % l_dnuc_am8) THEN
  ALLOCATE ( casdiags % dnuc_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dnuc_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!626
IF (casdiags % l_dnuc_am6) THEN
  ALLOCATE ( casdiags % dnuc_am6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dnuc_am6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!627
IF (casdiags % l_dnuc_am9) THEN
  ALLOCATE ( casdiags % dnuc_am9(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dnuc_am9(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!628
IF (casdiags % l_dnuc_an6) THEN
  ALLOCATE ( casdiags % dnuc_an6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dnuc_an6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!629
IF (casdiags % l_dsub_am2) THEN
  ALLOCATE ( casdiags % dsub_am2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsub_am2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!630
IF (casdiags % l_dsub_an2) THEN
  ALLOCATE ( casdiags % dsub_an2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsub_an2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!631
IF (casdiags % l_dsub_am6) THEN
  ALLOCATE ( casdiags % dsub_am6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsub_am6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!632
IF (casdiags % l_dsub_an6) THEN
  ALLOCATE ( casdiags % dsub_an6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsub_an6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!633
IF (casdiags % l_dsedi_am7) THEN
  ALLOCATE ( casdiags % dsedi_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsedi_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!634
IF (casdiags % l_dsedi_am8) THEN
  ALLOCATE ( casdiags % dsedi_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsedi_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!635
IF (casdiags % l_dsedi_an11) THEN
  ALLOCATE ( casdiags % dsedi_an11(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsedi_an11(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!636
IF (casdiags % l_dsedi_an12) THEN
  ALLOCATE ( casdiags % dsedi_an12(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsedi_an12(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!637
IF (casdiags % l_dseds_am7) THEN
  ALLOCATE ( casdiags % dseds_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dseds_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!638
IF (casdiags % l_dseds_am8) THEN
  ALLOCATE ( casdiags % dseds_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dseds_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!639
IF (casdiags % l_dseds_an11) THEN
  ALLOCATE ( casdiags % dseds_an11(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dseds_an11(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!640
IF (casdiags % l_dseds_an12) THEN
  ALLOCATE ( casdiags % dseds_an12(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dseds_an12(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!641
IF (casdiags % l_dsedg_am7) THEN
  ALLOCATE ( casdiags % dsedg_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsedg_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!642
IF (casdiags % l_dsedg_am8) THEN
  ALLOCATE ( casdiags % dsedg_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsedg_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!643
IF (casdiags % l_dsedg_an11) THEN
  ALLOCATE ( casdiags % dsedg_an11(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsedg_an11(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!644
IF (casdiags % l_dsedg_an12) THEN
  ALLOCATE ( casdiags % dsedg_an12(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsedg_an12(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!645
IF (casdiags % l_dssub_am2) THEN
  ALLOCATE ( casdiags % dssub_am2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dssub_am2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!646
IF (casdiags % l_dssub_an2) THEN
  ALLOCATE ( casdiags % dssub_an2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dssub_an2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!647
IF (casdiags % l_dssub_am6) THEN
  ALLOCATE ( casdiags % dssub_am6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dssub_am6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!648
IF (casdiags % l_dssub_an6) THEN
  ALLOCATE ( casdiags % dssub_an6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dssub_an6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!649
IF (casdiags % l_dgsub_am2) THEN
  ALLOCATE ( casdiags % dgsub_am2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgsub_am2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!650
IF (casdiags % l_dgsub_an2) THEN
  ALLOCATE ( casdiags % dgsub_an2(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgsub_an2(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!651
IF (casdiags % l_dgsub_am6) THEN
  ALLOCATE ( casdiags % dgsub_am6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgsub_am6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!652
IF (casdiags % l_dgsub_an6) THEN
  ALLOCATE ( casdiags % dgsub_an6(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgsub_an6(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!653
IF (casdiags % l_dhomc_am8) THEN
  ALLOCATE ( casdiags % dhomc_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dhomc_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!654
IF (casdiags % l_dhomc_am7) THEN
  ALLOCATE ( casdiags % dhomc_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dhomc_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!655
IF (casdiags % l_dhomr_am8) THEN
  ALLOCATE ( casdiags % dhomr_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dhomr_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!656
IF (casdiags % l_dhomr_am7) THEN
  ALLOCATE ( casdiags % dhomr_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dhomr_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!657
IF (casdiags % l_dimlt_am4) THEN
  ALLOCATE ( casdiags % dimlt_am4(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dimlt_am4(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!658
IF (casdiags % l_dimlt_am9) THEN
  ALLOCATE ( casdiags % dimlt_am9(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dimlt_am9(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!659
IF (casdiags % l_dsmlt_am4) THEN
  ALLOCATE ( casdiags % dsmlt_am4(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsmlt_am4(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!660
IF (casdiags % l_dsmlt_am9) THEN
  ALLOCATE ( casdiags % dsmlt_am9(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsmlt_am9(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!661
IF (casdiags % l_dgmlt_am4) THEN
  ALLOCATE ( casdiags % dgmlt_am4(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgmlt_am4(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!662
IF (casdiags % l_dgmlt_am9) THEN
  ALLOCATE ( casdiags % dgmlt_am9(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgmlt_am9(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!663
IF (casdiags % l_diacw_am8) THEN
  ALLOCATE ( casdiags % diacw_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % diacw_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!664
IF (casdiags % l_diacw_am7) THEN
  ALLOCATE ( casdiags % diacw_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % diacw_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!665
IF (casdiags % l_dsacw_am8) THEN
  ALLOCATE ( casdiags % dsacw_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsacw_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!666
IF (casdiags % l_dsacw_am7) THEN
  ALLOCATE ( casdiags % dsacw_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsacw_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!667
IF (casdiags % l_dgacw_am8) THEN
  ALLOCATE ( casdiags % dgacw_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgacw_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!668
IF (casdiags % l_dgacw_am7) THEN
  ALLOCATE ( casdiags % dgacw_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgacw_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!669
IF (casdiags % l_dsacr_am8) THEN
  ALLOCATE ( casdiags % dsacr_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsacr_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!670
IF (casdiags % l_dsacr_am7) THEN
  ALLOCATE ( casdiags % dsacr_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dsacr_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!671
IF (casdiags % l_dgacr_am8) THEN
  ALLOCATE ( casdiags % dgacr_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgacr_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!672
IF (casdiags % l_dgacr_am7) THEN
  ALLOCATE ( casdiags % dgacr_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % dgacr_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!673
IF (casdiags % l_draci_am8) THEN
  ALLOCATE ( casdiags % draci_am8(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % draci_am8(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!674
IF (casdiags % l_draci_am7) THEN
  ALLOCATE ( casdiags % draci_am7(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % draci_am7(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!675
IF (casdiags % l_asedl_am9) THEN
  ALLOCATE ( casdiags % asedl_am9(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % asedl_am9(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

!676
IF (casdiags % l_asedr_am9) THEN
  ALLOCATE ( casdiags % asedr_am9(i_start:i_end, j_start:j_end, k_start:k_end) )
  casdiags % asedr_am9(:,:,:) = zero_real_wp
  casdiags % l_process_rates = .TRUE.
END IF

! --------------------------------------

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



!-----------------------------------------
! aerosol diagnostics
!-----------------------------------------

IF ( ALLOCATED ( casdiags % asedr_am9 ) ) THEN
  DEALLOCATE ( casdiags % asedr_am9)
END IF !676

IF ( ALLOCATED ( casdiags % asedl_am9 ) ) THEN
  DEALLOCATE ( casdiags % asedl_am9)
END IF !675

IF ( ALLOCATED ( casdiags % draci_am7 ) ) THEN
  DEALLOCATE ( casdiags % draci_am7 )
END IF !674

IF ( ALLOCATED ( casdiags % draci_am8 ) ) THEN
  DEALLOCATE ( casdiags % draci_am8 )
END IF !673

IF ( ALLOCATED ( casdiags % dgacr_am7 ) ) THEN
  DEALLOCATE ( casdiags % dgacr_am7 )
END IF !672

IF ( ALLOCATED ( casdiags % dgacr_am8 ) ) THEN
  DEALLOCATE ( casdiags % dgacr_am8 )
END IF !671

IF ( ALLOCATED ( casdiags % dsacr_am7 ) ) THEN
  DEALLOCATE ( casdiags % dsacr_am7 )
END IF !670

IF ( ALLOCATED ( casdiags % dsacr_am8 ) ) THEN
  DEALLOCATE ( casdiags % dsacr_am8 )
END IF !669

IF ( ALLOCATED ( casdiags % dgacw_am7 ) ) THEN
  DEALLOCATE ( casdiags % dgacw_am7 )
END IF !668

IF ( ALLOCATED ( casdiags % dgacw_am8 ) ) THEN
  DEALLOCATE ( casdiags % dgacw_am8 )
END IF !667

IF ( ALLOCATED ( casdiags % dsacw_am7 ) ) THEN
  DEALLOCATE ( casdiags % dsacw_am7 )
END IF !666

IF ( ALLOCATED ( casdiags % dsacw_am8 ) ) THEN
  DEALLOCATE ( casdiags % dsacw_am8 )
END IF !665

IF ( ALLOCATED ( casdiags % diacw_am7 ) ) THEN
  DEALLOCATE ( casdiags % diacw_am7 )
END IF !664

IF ( ALLOCATED ( casdiags % diacw_am8 ) ) THEN
  DEALLOCATE ( casdiags % diacw_am8 )
END IF !663

IF ( ALLOCATED ( casdiags % dgmlt_am9 ) ) THEN
  DEALLOCATE ( casdiags % dgmlt_am9 )
END IF !662

IF ( ALLOCATED ( casdiags % dgmlt_am4 ) ) THEN
  DEALLOCATE ( casdiags % dgmlt_am4 )
END IF !661

IF ( ALLOCATED ( casdiags % dsmlt_am9 ) ) THEN
  DEALLOCATE ( casdiags % dsmlt_am9 )
END IF !660

IF ( ALLOCATED ( casdiags % dsmlt_am4 ) ) THEN
  DEALLOCATE ( casdiags % dsmlt_am4 )
END IF !659

IF ( ALLOCATED ( casdiags % dimlt_am9 ) ) THEN
  DEALLOCATE ( casdiags % dimlt_am9 )
END IF !658

IF ( ALLOCATED ( casdiags % dimlt_am4 ) ) THEN
  DEALLOCATE ( casdiags % dimlt_am4 )
END IF !657

IF ( ALLOCATED ( casdiags % dhomr_am7 ) ) THEN
  DEALLOCATE ( casdiags % dhomr_am7 )
END IF !656

IF ( ALLOCATED ( casdiags % dhomr_am8 ) ) THEN
  DEALLOCATE ( casdiags % dhomr_am8 )
END IF !655

IF ( ALLOCATED ( casdiags % dhomc_am7 ) ) THEN
  DEALLOCATE ( casdiags % dhomc_am7 )
END IF !654

IF ( ALLOCATED ( casdiags % dhomc_am8 ) ) THEN
  DEALLOCATE ( casdiags % dhomc_am8 )
END IF !653

IF ( ALLOCATED ( casdiags % dgsub_an6 ) ) THEN
  DEALLOCATE ( casdiags % dgsub_an6 )
END IF !652

IF ( ALLOCATED ( casdiags % dgsub_am6 ) ) THEN
  DEALLOCATE ( casdiags % dgsub_am6 )
END IF !651

IF ( ALLOCATED ( casdiags % dgsub_an2 ) ) THEN
  DEALLOCATE ( casdiags % dgsub_an2 )
END IF !650

IF ( ALLOCATED ( casdiags % dgsub_am2 ) ) THEN
  DEALLOCATE ( casdiags % dgsub_am2 )
END IF !649

IF ( ALLOCATED ( casdiags % dssub_an6 ) ) THEN
  DEALLOCATE ( casdiags % dssub_an6 )
END IF !648

IF ( ALLOCATED ( casdiags % dssub_am6 ) ) THEN
  DEALLOCATE ( casdiags % dssub_am6 )
END IF !647

IF ( ALLOCATED ( casdiags % dssub_an2 ) ) THEN
  DEALLOCATE ( casdiags % dssub_an2 )
END IF !646

IF ( ALLOCATED ( casdiags % dssub_am2 ) ) THEN
  DEALLOCATE ( casdiags % dssub_am2 )
END IF !645

IF ( ALLOCATED ( casdiags % dsedg_an12 ) ) THEN
  DEALLOCATE ( casdiags % dsedg_an12 )
END IF !644

IF ( ALLOCATED ( casdiags % dsedg_an11 ) ) THEN
  DEALLOCATE ( casdiags % dsedg_an11 )
END IF !643

IF ( ALLOCATED ( casdiags % dsedg_am8 ) ) THEN
  DEALLOCATE ( casdiags % dsedg_am8 )
END IF !642

IF ( ALLOCATED ( casdiags % dsedg_am7 ) ) THEN
  DEALLOCATE ( casdiags % dsedg_am7 )
END IF !641

IF ( ALLOCATED ( casdiags % dseds_an12 ) ) THEN
  DEALLOCATE ( casdiags % dseds_an12 )
END IF !640

IF ( ALLOCATED ( casdiags % dseds_an11 ) ) THEN
  DEALLOCATE ( casdiags % dseds_an11 )
END IF !639

IF ( ALLOCATED ( casdiags % dseds_am8 ) ) THEN
  DEALLOCATE ( casdiags % dseds_am8 )
END IF !638

IF ( ALLOCATED ( casdiags % dseds_am7 ) ) THEN
  DEALLOCATE ( casdiags % dseds_am7 )
END IF !637

IF ( ALLOCATED ( casdiags % dsedi_an12 ) ) THEN
  DEALLOCATE ( casdiags % dsedi_an12 )
END IF !636

IF ( ALLOCATED ( casdiags % dsedi_an11 ) ) THEN
  DEALLOCATE ( casdiags % dsedi_an11 )
END IF !635

IF ( ALLOCATED ( casdiags % dsedi_am8 ) ) THEN
  DEALLOCATE ( casdiags % dsedi_am8 )
END IF !634

IF ( ALLOCATED ( casdiags % dsedi_am7 ) ) THEN
  DEALLOCATE ( casdiags % dsedi_am7 )
END IF !633

IF ( ALLOCATED ( casdiags % dsub_an6 ) ) THEN
  DEALLOCATE ( casdiags % dsub_an6 )
END IF !632

IF ( ALLOCATED ( casdiags % dsub_am6 ) ) THEN
  DEALLOCATE ( casdiags % dsub_am6 )
END IF !631

IF ( ALLOCATED ( casdiags % dsub_an2 ) ) THEN
  DEALLOCATE ( casdiags % dsub_an2 )
END IF !630

IF ( ALLOCATED ( casdiags % dsub_am2 ) ) THEN
  DEALLOCATE ( casdiags % dsub_am2 )
END IF !629

IF ( ALLOCATED ( casdiags % dnuc_an6 ) ) THEN
  DEALLOCATE ( casdiags % dnuc_an6 )
END IF !628

IF ( ALLOCATED ( casdiags % dnuc_am9 ) ) THEN
  DEALLOCATE ( casdiags % dnuc_am9 )
END IF !627

IF ( ALLOCATED ( casdiags % dnuc_am6 ) ) THEN
  DEALLOCATE ( casdiags % dnuc_am6 )
END IF !626

IF ( ALLOCATED ( casdiags % dnuc_am8 ) ) THEN
  DEALLOCATE ( casdiags % dnuc_am8 )
END IF !625

IF ( ALLOCATED ( casdiags % asedl_an12 ) ) THEN
  DEALLOCATE ( casdiags % asedl_an12 )
END IF !624

IF ( ALLOCATED ( casdiags % asedl_an11 ) ) THEN
  DEALLOCATE ( casdiags % asedl_an11 )
END IF !623

IF ( ALLOCATED ( casdiags % asedl_am4 ) ) THEN
  DEALLOCATE ( casdiags % asedl_am4 )
END IF !622

IF ( ALLOCATED ( casdiags % arevp_an6 ) ) THEN
  DEALLOCATE ( casdiags % arevp_an6 )
END IF !621

IF ( ALLOCATED ( casdiags % arevp_am6 ) ) THEN
  DEALLOCATE ( casdiags % arevp_am6 )
END IF !620

IF ( ALLOCATED ( casdiags % arevp_am5 ) ) THEN
  DEALLOCATE ( casdiags % arevp_am5 )
END IF !619

IF ( ALLOCATED ( casdiags % arevp_am4 ) ) THEN
  DEALLOCATE ( casdiags % arevp_am4 )
END IF !618

IF ( ALLOCATED ( casdiags % arevp_an3 ) ) THEN
  DEALLOCATE ( casdiags % arevp_an3 )
END IF !617

IF ( ALLOCATED ( casdiags % arevp_am3 ) ) THEN
  DEALLOCATE ( casdiags % arevp_am3 )
END IF !616

IF ( ALLOCATED ( casdiags % arevp_an2 ) ) THEN
  DEALLOCATE ( casdiags % arevp_an2 )
END IF !615

IF ( ALLOCATED ( casdiags % arevp_am2 ) ) THEN
  DEALLOCATE ( casdiags % arevp_am2 )
END IF !614

IF ( ALLOCATED ( casdiags % asedr_an12 ) ) THEN
  DEALLOCATE ( casdiags % asedr_an12 )
END IF !613

IF ( ALLOCATED ( casdiags % asedr_an11 ) ) THEN
  DEALLOCATE ( casdiags % asedr_an11 )
END IF !612

IF ( ALLOCATED ( casdiags % asedr_am ) ) THEN
  DEALLOCATE ( casdiags % asedr_am )
END IF !611

IF ( ALLOCATED ( casdiags % aacw ) ) THEN
  DEALLOCATE ( casdiags % aacw )
END IF !610

IF ( ALLOCATED ( casdiags % aaut ) ) THEN
  DEALLOCATE ( casdiags % aaut )
END IF !609

IF ( ALLOCATED ( casdiags % aact_an6 ) ) THEN
  DEALLOCATE ( casdiags % aact_an6 )
END IF !608

IF ( ALLOCATED ( casdiags % aact_am9 ) ) THEN
  DEALLOCATE ( casdiags % aact_am9 )
END IF !607

IF ( ALLOCATED ( casdiags % aact_an3 ) ) THEN
  DEALLOCATE ( casdiags % aact_an3 )
END IF !606

IF ( ALLOCATED ( casdiags % aact_am3 ) ) THEN
  DEALLOCATE ( casdiags % aact_am3 )
END IF !605

IF ( ALLOCATED ( casdiags % aact_an2 ) ) THEN
  DEALLOCATE ( casdiags % aact_an2 )
END IF !604

IF ( ALLOCATED ( casdiags % aact_am2 ) ) THEN
  DEALLOCATE ( casdiags % aact_am2 )
END IF !603

IF ( ALLOCATED ( casdiags % aact_an1 ) ) THEN
  DEALLOCATE ( casdiags % aact_an1 )
END IF !602

IF ( ALLOCATED ( casdiags % aact_am1 ) ) THEN
  DEALLOCATE ( casdiags % aact_am1 )
END IF !601

IF ( ALLOCATED ( casdiags % niics_i ) ) THEN
  DEALLOCATE ( casdiags % niics_i )
END IF ! 531

IF ( ALLOCATED ( casdiags % niics_s ) ) THEN
  DEALLOCATE ( casdiags % niics_s )
END IF ! 530

IF ( ALLOCATED ( casdiags % ngaci ) ) THEN
  DEALLOCATE ( casdiags % ngaci )
END IF ! 529

IF ( ALLOCATED ( casdiags % pgaci ) ) THEN
  DEALLOCATE ( casdiags % pgaci )
END IF ! 528

IF ( ALLOCATED ( casdiags % nidps ) ) THEN
  DEALLOCATE ( casdiags % nidps )
END IF ! 527

IF ( ALLOCATED ( casdiags % pidps ) ) THEN
  DEALLOCATE ( casdiags % pidps )
END IF ! 526

IF ( ALLOCATED ( casdiags % nraci_i ) ) THEN
  DEALLOCATE ( casdiags % nraci_i )
END IF ! 525

IF ( ALLOCATED ( casdiags % nraci_r ) ) THEN
  DEALLOCATE ( casdiags % nraci_r )
END IF ! 524

IF ( ALLOCATED ( casdiags % nraci_g ) ) THEN
  DEALLOCATE ( casdiags % nraci_g )
END IF ! 523

IF ( ALLOCATED ( casdiags % praci_i ) ) THEN
  DEALLOCATE ( casdiags % praci_i )
END IF ! 522

IF ( ALLOCATED ( casdiags % praci_r ) ) THEN
  DEALLOCATE ( casdiags % praci_r )
END IF ! 521

IF ( ALLOCATED ( casdiags % praci_g ) ) THEN
  DEALLOCATE ( casdiags % praci_g )
END IF ! 520

IF ( ALLOCATED ( casdiags % pihal ) ) THEN
  DEALLOCATE ( casdiags % pihal )
END IF ! 519

IF ( ALLOCATED ( casdiags % ngmlt ) ) THEN
  DEALLOCATE ( casdiags % ngmlt )
END IF ! 518

IF ( ALLOCATED ( casdiags % ngacw ) ) THEN
  DEALLOCATE ( casdiags % ngacw )
END IF ! 517

IF ( ALLOCATED ( casdiags % nsmlt ) ) THEN
  DEALLOCATE ( casdiags % nsmlt )
END IF ! 516

IF ( ALLOCATED ( casdiags % nimlt ) ) THEN
  DEALLOCATE ( casdiags % nimlt )
END IF ! 515

IF ( ALLOCATED ( casdiags % nsacr ) ) THEN
  DEALLOCATE ( casdiags % nsacr )
END IF ! 514

IF ( ALLOCATED ( casdiags % nsacw ) ) THEN
  DEALLOCATE ( casdiags % nsacw )
END IF ! 513

IF ( ALLOCATED ( casdiags % niacw ) ) THEN
  DEALLOCATE ( casdiags % niacw )
END IF ! 512

IF ( ALLOCATED ( casdiags % ngsub ) ) THEN
  DEALLOCATE ( casdiags % ngsub )
END IF ! 511

IF ( ALLOCATED ( casdiags % ngacs ) ) THEN
  DEALLOCATE ( casdiags % ngacs )
END IF ! 510

IF ( ALLOCATED ( casdiags % nsaci ) ) THEN
  DEALLOCATE ( casdiags % nsaci )
END IF ! 509

IF ( ALLOCATED ( casdiags % nsaut ) ) THEN
  DEALLOCATE ( casdiags % nsaut )
END IF ! 508

IF ( ALLOCATED ( casdiags % nssub ) ) THEN
  DEALLOCATE ( casdiags % nssub )
END IF ! 507

IF ( ALLOCATED ( casdiags % nisub ) ) THEN
  DEALLOCATE ( casdiags % nisub )
END IF ! 506

IF ( ALLOCATED ( casdiags % nrevp ) ) THEN
  DEALLOCATE ( casdiags % nrevp )
END IF ! 505

IF ( ALLOCATED ( casdiags % nsedr ) ) THEN
  DEALLOCATE ( casdiags % nsedr )
END IF ! 504

IF ( ALLOCATED ( casdiags % nracr ) ) THEN
  DEALLOCATE ( casdiags % nracr )
END IF ! 503

IF ( ALLOCATED ( casdiags % nracw ) ) THEN
  DEALLOCATE ( casdiags % nracw )
END IF ! 502

IF ( ALLOCATED ( casdiags % nsedl ) ) THEN
  DEALLOCATE ( casdiags % nsedl )
END IF ! 501

IF ( ALLOCATED ( casdiags % nraut ) ) THEN
  DEALLOCATE ( casdiags % nraut )
END IF ! 500

!-----------------------------------------

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
