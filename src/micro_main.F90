module micro_main
  use variable_precision, only: wp
  use mphys_parameters, only: nz, nq, rain_params, cloud_params, ice_params, &
       snow_params, graupel_params, nspecies, ZERO_REAL_WP, a_s, b_s
  use process_routines, only: process_rate, zero_procs, allocate_procs, deallocate_procs, i_cond, i_praut, &
       i_pracw, i_pracr, i_prevp, i_psedr, i_psedl, i_aact, i_aaut, i_aacw, i_aevp, i_asedr, i_asedl, i_arevp, &
       i_tidy2, i_atidy2, i_inuc, i_idep, i_dnuc, i_dsub, i_saut, i_iacw, i_sacw, i_pseds, &
       i_sdep, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs, i_gdep, i_psedg, i_sagg, &
       i_gshd, i_ihal, i_smlt, i_gmlt, i_psedi, i_homr, i_homc, i_imlt, i_isub, i_ssub, i_gsub, i_sbrk, i_dssub, &
       i_dgsub, i_dsedi, i_dseds, i_dsedg, i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr, &
       i_dgacr, i_draci
  use sum_process, only: sum_procs, sum_aprocs, tend_temp, aerosol_tend_temp
  use aerosol_routines, only: examine_aerosol, aerosol_phys, aerosol_chem, aerosol_active, allocate_aerosol, &
       deallocate_aerosol
  use mphys_switches, only: hydro_complexity, aero_complexity, i_qv, i_ql, i_nl, i_qr, i_nr, i_m3r, i_th, i_qi, &
       i_qs, i_qg, i_ni, i_ns, i_ng, i_m3s, i_m3g, i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am4, i_am5, i_am6, &
       i_an6, i_am7, i_am8 , i_am9, i_am10, i_an10, i_an11, i_an12, i_ak1, i_ak2, i_ak3, &
       aerosol_option, l_warm, &
       l_sed, l_idep, aero_index, nq_l, nq_r, nq_i, nq_s, nq_g, &
       l_sg, l_g, l_process, max_sed_length, max_step_length, l_harrington, l_passive, ntotala, ntotalq, &
       l_onlycollect, pswitch, l_isub, l_pos1, l_pos2, l_pos3, l_pos4, &
       l_pos5, l_pos6, i_hstart, l_tidy_negonly,  &
       iopt_act, iopt_shipway_act, l_prf_cfrac, l_kfsm, l_gamma_online, l_subseds_maxv, &
       i_cfl, i_cfr, i_cfi, i_cfs, i_cfg 
! use mphys_switches, only: l_rain
  use passive_fields, only: rexner, min_dz
  use mphys_constants, only: cp, Lv
  use distributions, only: query_distributions, initialise_distributions, dist_lambda, dist_mu, dist_n0, dist_lams
  use passive_fields, only: initialise_passive_fields, set_passive_fields, TdegK, rhcrit_1d
  use autoconversion, only: raut
  use evaporation, only: revp
  use condensation, only: condevp_initialise, condevp_finalise, condevp
  use accretion, only: racw
  use aggregation, only: racr, ice_aggregation
  use sedimentation, only: sedr, sedr_1M_2M, terminal_velocity_CFL
  use ice_nucleation, only: inuc
  use ice_deposition, only: idep
  use ice_accretion, only: iacc
  use breakup, only: ice_breakup
  use snow_autoconversion, only: saut
  use ice_multiplication, only: hallet_mossop
  use graupel_wetgrowth, only: wetgrowth
  use graupel_embryo, only: graupel_embryos
  use ice_melting, only: melting
  use homogeneous, only: ihom_rain, ihom_droplets
  use adjust_deposition, only: adjust_dep
  use mphys_constants, only: fixed_aerosol_sigma, fixed_aerosol_density

  use mphys_tidy, only: initialise_mphystidy, finalise_mphystidy, qtidy, ensure_positive, &
       ensure_saturated, tidy_qin, tidy_ain, ensure_positive_aerosol
  use preconditioning, only: precondition, preconditioner

  use casim_reflec_mod, only: casim_reflec, setup_reflec_constants
  ! For initialization of Shipway (2015) activation scheme
  use shipway_lookup, only: generate_tables

  use generic_diagnostic_variables, only: casdiags

  use casim_runtime, only: casim_time
  use casim_parent_mod, only: casim_parent, parent_monc

! #if DEF_MODEL==MODEL_KiD
!   ! Kid modules
!   use diagnostics, only: save_dg, i_dgtime, n_sub, n_subsed
!   use runtime, only: time
!   use parameters, only: nx
!   Use namelists, only : no_precip_time, l_sediment
! #endif


  implicit none

  private

  character(len=*), parameter, private :: ModuleName='MICRO_MAIN'

  logical :: l_tendency_loc
  logical :: l_warm_loc

!$OMP THREADPRIVATE(l_tendency_loc, l_warm_loc)

  integer :: is, ie ! upper and lower i levels which are to be used
  integer :: js, je ! upper and lower j levels
  integer :: ks, ke ! upper and lower k levels

!$OMP THREADPRIVATE(is,ie,js,je,ks,ke)
!  integer :: nxny
  real(wp), allocatable, save :: precip(:,:) ! diagnostic for surface precip rate
  real(wp), allocatable :: dqfields(:,:), qfields(:,:), tend(:,:)
  real(wp), allocatable :: daerofields(:,:), aerofields(:,:), aerosol_tend(:,:)
  real(wp), allocatable :: cffields(:,:) !cloudfraction fields

!$OMP THREADPRIVATE(precip, dqfields, qfields, cffields, tend,                   &
!$OMP               daerofields, aerofields, aerosol_tend)

  type(process_rate), allocatable :: procs(:,:)
  type(process_rate), allocatable :: aerosol_procs(:,:)

!$OMP THREADPRIVATE(procs, aerosol_procs)

  type(aerosol_active), allocatable :: aeroact(:)
  type(aerosol_phys), allocatable   :: aerophys(:)
  type(aerosol_chem), allocatable   :: aerochem(:)

!$OMP THREADPRIVATE(aeroact, aerophys, aerochem)

  type(aerosol_active), allocatable :: dustact(:)
  type(aerosol_phys), allocatable   :: dustphys(:)
  type(aerosol_chem), allocatable   :: dustchem(:)

!$OMP THREADPRIVATE(dustact, dustphys, dustchem)

  type(aerosol_active), allocatable :: aeroice(:)  ! Soluble aerosol in ice
  type(aerosol_active), allocatable :: dustliq(:)! Insoluble aerosol in liquid

!$OMP THREADPRIVATE(aeroice, dustliq)

  real(wp), allocatable :: qfields_in(:,:)
  real(wp), allocatable :: qfields_mod(:,:)
  real(wp), allocatable :: aerofields_in(:,:)
  real(wp), allocatable :: aerofields_mod(:,:)

!$OMP THREADPRIVATE(qfields_in, qfields_mod, aerofields_in, aerofields_mod)

  logical, allocatable :: l_Tcold(:) ! temperature below freezing, i.e. .not. l_Twarm
  logical, allocatable :: l_sigevap(:) ! logical to determine significant evaporation

!$OMP THREADPRIVATE(l_Tcold, l_sigevap)

  real :: DTPUD ! Time step for puddle diagnostic

  public initialise_micromain, finalise_micromain, shipway_microphysics, DTPUD

!PRF
  public qfields_in, qfields_mod, aerofields_in, aerofields_mod, dqfields, qfields, tend, cffields
  public daerofields, aerofields, aerosol_tend, procs, aerosol_procs
  public aerophys, aeroact, aerochem, dustact, dustphys, dustchem, aeroice, dustliq
!PRF
contains

  subroutine initialise_micromain(il, iu, jl, ju, kl, ku,                 &
       is_in, ie_in, js_in, je_in, ks_in, ke_in, l_tendency)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='INITIALISE_MICROMAIN'

    integer, intent(in) :: il, iu ! upper and lower i levels
    integer, intent(in) :: jl, ju ! upper and lower j levels
    integer, intent(in) :: kl, ku ! upper and lower k levels

    integer, intent(in) :: is_in, ie_in ! upper and lower i levels which are to be used
    integer, intent(in) :: js_in, je_in ! upper and lower j levels
    integer, intent(in) :: ks_in, ke_in ! upper and lower k levels

    ! New optional l_tendency logical added...
    ! if true then a tendency is returned (i.e. units/s)
    ! if false then an increment is returned (i.e. units/timestep)
    logical, intent(in) :: l_tendency

    ! Local variables

    integer :: k

    integer :: nprocs     ! number of process rates stored
    integer :: naeroprocs ! number of process rates stored
    integer :: naero      ! number of aerosol fields

    real(wp) :: beta_init

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_warm_loc=l_warm ! Original setting

    is=is_in
    ie=ie_in
    js=js_in
    je=je_in
    ks=ks_in
    ke=ke_in

    allocate(rhcrit_1d(kl:ku))
    ! Set RHCrit to 1.0 as default; parent model can then overwrite
    ! this if needed
    rhcrit_1d(:) = 1.0

    l_tendency_loc = l_tendency

    nq=sum(hydro_complexity%nmoments)+2 ! also includes vapour and theta
    nz=ke-ks+1
    nprocs = hydro_complexity%nprocesses
    !nxny=(ie-is+1)*(je-js+1)

    allocate(precondition(nz))
    precondition=.true. ! Assume all points need to be considered
    allocate(l_Tcold(nz)) 
    l_Tcold =.false. ! Assumes no cold points, this is set at the beginning of microphysics_common 
    allocate(l_sigevap(nz))
    l_sigevap = .false.
    allocate(qfields(nz, nq))
    allocate(dqfields(nz, nq))
    qfields=ZERO_REAL_WP
    dqfields=ZERO_REAL_WP
    !allocate(procs(nz, nprocs))
    allocate(procs(ntotalq, nprocs))
    allocate(tend(nz, nq))
    allocate(tend_temp(nz,nq))
    allocate(cffields(nz,5)) !5 'cloud' fractions
    cffields=ZERO_REAL_WP

    ! Allocate aerosol storage
    if (aerosol_option > 0) then
      naero=ntotala
      naeroprocs=aero_complexity%nprocesses
      allocate(aerofields(nz, naero))
      allocate(daerofields(nz, naero))
      aerofields=ZERO_REAL_WP
      daerofields=ZERO_REAL_WP
      ! allocate(aerosol_procs(nz, naeroprocs))
      allocate(aerosol_procs(ntotala, naeroprocs))
      allocate(aerosol_tend(nz, naero))
      allocate(aerosol_tend_temp(nz, naero))
    else
      ! Dummy arrays required
      allocate(aerofields(1,1))
      allocate(daerofields(1,1))
      allocate(aerosol_procs(1,1))
      allocate(aerosol_tend(1,1))
      allocate(aerosol_tend_temp(1,1))
    end if

    allocate(aerophys(nz))
    allocate(aerochem(nz))
    allocate(aeroact(nz))

    call allocate_aerosol(aerophys, aerochem, aero_index%nccn)
    allocate(dustphys(nz))
    allocate(dustchem(nz))
    allocate(dustact(nz))
    call allocate_aerosol(dustphys, dustchem, aero_index%nin)

    allocate(aeroice(nz))
    allocate(dustliq(nz))

    ! Preserve initial values for non-Shipway activation
    if ( iopt_act == iopt_shipway_act ) then
      beta_init = 0.5
    else
      beta_init = 1.0
    end if

    ! Temporary initialization of chem and sigma
    do k =1,size(aerophys)
      aerophys(k)%sigma(:)=fixed_aerosol_sigma
      aerophys(k)%rpart(:)=0.0
      aerochem(k)%vantHoff(:)=3.0
      aerochem(k)%massMole(:)=132.0e-3
      aerochem(k)%density(:)=fixed_aerosol_density
      aerochem(k)%epsv(:)=1.0
      aerochem(k)%beta(:)=beta_init
    end do
    do k =1,size(dustphys)
      dustphys(k)%sigma(:)=fixed_aerosol_sigma
      dustphys(k)%rpart(:)=0.0
      dustchem(k)%vantHoff(:)=3.0
      dustchem(k)%massMole(:)=132.0e-3
      dustchem(k)%density(:)=fixed_aerosol_density
      dustchem(k)%epsv(:)=1.0
      dustchem(k)%beta(:)=beta_init
    end do

    !allocate space for the process rates
    call allocate_procs(procs, nz, nprocs, ntotalq)
    if (l_process) call allocate_procs(aerosol_procs, nz, naeroprocs, ntotala)

    ! allocate diagnostics
    allocate(precip(il:iu,jl:ju))

    call initialise_passive_fields(ks, ke)

    allocate(qfields_in(nz, nq))
    allocate(qfields_mod(nz, nq))
    if (aerosol_option > 0) then
      allocate(aerofields_in(nz, naero))
      allocate(aerofields_mod(nz, naero))
    end if
    call initialise_distributions(nz, nspecies)
    call initialise_mphystidy()
    call condevp_initialise()
! Here we initialise some things for the Shipway(2015) aerosol activation.

!$OMP SINGLE
    if ( iopt_act == iopt_shipway_act ) then
      call generate_tables()
    end if

    call setup_reflec_constants()
!$OMP END SINGLE

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine initialise_micromain

  subroutine finalise_micromain()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='FINALISE_MICROMAIN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! deallocate diagnostics
    deallocate(precip)

    ! deallocate process rates
    call deallocate_procs(procs)

    deallocate(procs)
    deallocate(qfields)
    deallocate(tend)
    deallocate(tend_temp)
    deallocate(precondition)
    deallocate(cffields)
    
    ! aerosol fields
    if (l_process) call deallocate_procs(aerosol_procs)

    deallocate(dustliq)
    deallocate(aeroice)

    call deallocate_aerosol(aerophys, aerochem)
    deallocate(aerophys)
    deallocate(aerochem)
    deallocate(aeroact)
    call deallocate_aerosol(dustphys, dustchem)
    deallocate(dustphys)
    deallocate(dustchem)
    deallocate(dustact)
    deallocate(aerosol_procs)
    deallocate(aerosol_tend)
    deallocate(aerosol_tend_temp)
    deallocate(aerofields)
    deallocate(daerofields)
    deallocate(rhcrit_1d)
    deallocate(qfields_in)
    deallocate(qfields_mod)
    if (allocated(aerofields_in)) deallocate(aerofields_in)
    if (allocated (aerofields_mod)) deallocate(aerofields_mod)
    deallocate(dist_lambda)
    deallocate(dist_mu)
    deallocate(dist_n0)
    deallocate(dist_lams)
    deallocate(a_s)
    deallocate(b_s)
    call finalise_mphystidy()
    call condevp_finalise()

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine finalise_micromain

  subroutine shipway_microphysics(il, iu, jl, ju, kl, ku, dt,               &
       qv, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13,          &
       theta, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13,       &
       a14, a15, a16, a17, a18, a19, a20,                                   &
       exner, pressure, rho, w, tke, z_half, z_centre, dz,                  &
       cfliq, cfice, cfsnow, cfrain, cfgr,    &
       dqv, dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11, dq12,  &
       dq13, dth, da1, da2, da3, da4, da5, da6, da7, da8, da9, da10, da11,  &
       da12, da13, da14, da15, da16, da17,                                  &
       is_in, ie_in, js_in, je_in, ks_in, ke_in,                            &
       l_tendency)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SHIPWAY_MICROPHYSICS'

    integer, intent(in) :: il, iu ! upper and lower i levels
    integer, intent(in) :: jl, ju ! upper and lower j levels
    integer, intent(in) :: kl, ku ! upper and lower k levels

    real(wp), intent(in) :: dt    ! parent model timestep (s)

    ! hydro fields in... 1-5 should be warm rain, 6+ are ice
    ! see mphys_casim for details of what is passed in
    real(wp), intent(in) :: q1( kl:ku, il:iu, jl:ju ), q2( kl:ku, il:iu, jl:ju )   &
         , q3( kl:ku, il:iu, jl:ju ), q4( kl:ku, il:iu, jl:ju ), q5( kl:ku, il:iu, jl:ju ) &
         , q6( kl:ku, il:iu, jl:ju ), q7( kl:ku, il:iu, jl:ju ), q8( kl:ku, il:iu, jl:ju ) &
         , q9( kl:ku, il:iu, jl:ju ), q10( kl:ku, il:iu, jl:ju ), q11( kl:ku, il:iu, jl:ju ) &
         , q12( kl:ku, il:iu, jl:ju ), q13( kl:ku, il:iu, jl:ju )

    real(wp) :: cfliq(kl:ku, il:iu, jl:ju ), cfrain(kl:ku, il:iu, jl:ju ), cfice(kl:ku, il:iu, jl:ju ), &
                cfsnow(kl:ku, il:iu, jl:ju ), cfgr(kl:ku, il:iu, jl:ju )



    real(wp), intent(in) :: qv( kl:ku, il:iu, jl:ju )
    real(wp), intent(in) :: theta( kl:ku, il:iu, jl:ju )
    real(wp), intent(in) :: exner( kl:ku, il:iu, jl:ju )
    real(wp), intent(in) :: pressure( kl:ku, il:iu, jl:ju )
    real(wp), intent(in) :: rho( kl:ku, il:iu, jl:ju )
    real(wp), intent(in) :: w( kl:ku, il:iu, jl:ju )
    real(wp), intent(in) :: tke( kl:ku, il:iu, jl:ju )
    real(wp), intent(in) :: dz( kl:ku, il:iu, jl:ju )
    real(wp), intent(in) :: z_half( kl-1:ku, il:iu, jl:ju )
    real(wp), intent(in) :: z_centre( kl:ku, il:iu, jl:ju )


    ! Aerosol fields in
    real(wp), intent(in) :: a1( kl:ku, il:iu, jl:ju ), a2( kl:ku, il:iu, jl:ju )   &
         , a3( kl:ku, il:iu, jl:ju ), a4( kl:ku, il:iu, jl:ju ), a5( kl:ku, il:iu, jl:ju ) &
         , a6( kl:ku, il:iu, jl:ju ), a7( kl:ku, il:iu, jl:ju ), a8( kl:ku, il:iu, jl:ju ) &
         , a9( kl:ku, il:iu, jl:ju ), a10( kl:ku, il:iu, jl:ju ), a11(kl:ku, il:iu, jl:ju ) &
         , a12( kl:ku, il:iu, jl:ju ), a13( kl:ku, il:iu, jl:ju ), a14( kl:ku, il:iu, jl:ju ) &
         , a15( kl:ku, il:iu, jl:ju ), a16( kl:ku, il:iu, jl:ju ), a17(kl:ku,il:iu, jl:ju )  &
         , a18( kl:ku, il:iu, jl:ju ), a19( kl:ku, il:iu, jl:ju ), a20( kl:ku,il:iu, jl:ju )

    ! hydro tendencies in:  from parent model forcing i.e. advection
    ! hydro tendencies out: from microphysics only...
    real(wp), intent(inout) :: dq1( kl:ku, il:iu, jl:ju ), dq2( kl:ku, il:iu, jl:ju ) &
         , dq3( kl:ku, il:iu, jl:ju ), dq4( kl:ku, il:iu, jl:ju ), dq5( kl:ku, il:iu, jl:ju ) &
         , dq6( kl:ku, il:iu, jl:ju ), dq7( kl:ku, il:iu, jl:ju ), dq8( kl:ku, il:iu, jl:ju ) &
         , dq9( kl:ku, il:iu, jl:ju ), dq10( kl:ku, il:iu, jl:ju ), dq11( kl:ku, il:iu, jl:ju ) &
         , dq12( kl:ku, il:iu, jl:ju ), dq13( kl:ku, il:iu, jl:ju )

    ! qv/theta tendencies in:  from parent model forcing i.e. advection
    ! qv/theta tendencies out: from microphysics only
    real(wp), intent(inout) :: dqv( kl:ku, il:iu, jl:ju ), dth( kl:ku, il:iu, jl:ju )

    ! aerosol tendencies in:  from parent model forcing i.e. advection
    ! aerosol tendencies out: from microphysics only
    real(wp), intent(inout) :: da1( kl:ku, il:iu, jl:ju ), da2( kl:ku, il:iu, jl:ju ) &
         , da3( kl:ku, il:iu, jl:ju ), da4( kl:ku, il:iu, jl:ju ), da5( kl:ku, il:iu, jl:ju ) &
         , da6( kl:ku, il:iu, jl:ju ), da7( kl:ku, il:iu, jl:ju ), da8( kl:ku, il:iu, jl:ju ) &
         , da9( kl:ku, il:iu, jl:ju ), da10( kl:ku, il:iu, jl:ju ), da11( kl:ku, il:iu, jl:ju ) &
         , da12( kl:ku, il:iu, jl:ju ), da13( kl:ku, il:iu, jl:ju ), da14( kl:ku, il:iu, jl:ju ) &
         , da15( kl:ku, il:iu, jl:ju ), da16( kl:ku, il:iu, jl:ju ), da17( kl:ku, il:iu, jl:ju )

    integer, intent(in), optional :: is_in, ie_in ! upper and lower i levels which are to be used
    integer, intent(in), optional :: js_in, je_in ! upper and lower j levels
    integer, intent(in), optional :: ks_in, ke_in ! upper and lower k levels
    !! ks_in and ke_in are not used in this routine so generate compile warnings

    ! New optional l_tendency logical added...
    ! if true then a tendency is returned (i.e. units/s)
    ! if false then an increment is returned (i.e. units/timestep)
    logical, intent(in), optional :: l_tendency
    !! l_tendency is not used in this routine so generates a compile warning

    ! Local variables

    integer :: k, i, j

    !diagnostics

    real(wp) :: precip_l
    real(wp) :: precip_r
    real(wp) :: precip_i
    real(wp) :: precip_s
    real(wp) :: precip_g

    real(wp) :: precip_l1d(nz)
    real(wp) :: precip_r1d(nz)
    real(wp) :: precip_i1d(nz)
    real(wp) :: precip_s1d(nz)
    real(wp) :: precip_so1d(nz)
    real(wp) :: precip_g1d(nz)

    real(wp) :: waterpath

    real(wp) :: dbz_tot_c(nz), dbz_g_c(nz), dbz_i_c(nz), &
                dbz_s_c(nz),   dbz_l_c(nz), dbz_r_c(nz)

    INTEGER :: kc ! Casim Z-level

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (casim_parent == parent_monc) casim_time = casim_time + dt
    ! (In the KiD model, casim_time is set to variable 'time')

! #if DEF_MODEL==MODEL_KiD
!     !AH - following code limits stops precip processes and sedimentation for no_precip_time. This 
!     !     is needed for the KiD-A 2d Sc case
!     if (time <= no_precip_time) then 
!        pswitch%l_praut=.false.
!        pswitch%l_pracw=.false.
!        pswitch%l_pracr=.false.
!        pswitch%l_prevp=.false.
!        l_sed=.false.
!     endif
!     ! if (time > no_precip_time .and. l_rain .and. l_sediment) then 
!     !    pswitch%l_praut=.true.
!     !    pswitch%l_pracw=.true.
!     !    pswitch%l_pracr=.true.
!     !    pswitch%l_prevp=.true.
!     !    l_sed=.true.
!     ! endif
! #endif

    do j=js_in,je_in

      do i=is_in,ie_in

        precip_l = 0.0
        precip_r = 0.0
        precip_i = 0.0
        precip_s = 0.0
        precip_g = 0.0

        precip_l1d(:)  = 0.0
        precip_r1d(:)  = 0.0
        precip_i1d(:)  = 0.0 
        precip_s1d(:)  = 0.0
        precip_g1d(:)  = 0.0
        precip_so1d(:) = 0.0

        tend=ZERO_REAL_WP
        call zero_procs(procs)
        aerosol_tend=ZERO_REAL_WP
        if (l_process) then
          call zero_procs(aerosol_procs)
        end if

        !set cloud fraction fields
        cffields(:,i_cfl)=cfliq(ks:ke,i,j)
        cffields(:,i_cfr)=cfrain(ks:ke,i,j)
        cffields(:,i_cfi)=cfice(ks:ke,i,j)
        cffields(:,i_cfs)=cfsnow(ks:ke,i,j)
        cffields(:,i_cfg)=cfgr(ks:ke,i,j)

        ! Set the qfields
        qfields(:, i_qv)=qv(ks:ke,i,j)
        qfields(:, i_th)=theta(ks:ke,i,j)

        if (nq_l > 0) qfields(:,i_ql)=q1(ks:ke,i,j)
        if (nq_r > 0) qfields(:,i_qr)=q2(ks:ke,i,j)
        if (nq_l > 1) qfields(:,i_nl)=q3(ks:ke,i,j)
        if (nq_r > 1) qfields(:,i_nr)=q4(ks:ke,i,j)
        if (nq_r > 2) qfields(:,i_m3r)=q5(ks:ke,i,j)
        if (nq_i > 0) qfields(:,i_qi)=q6(ks:ke,i,j)
        if (nq_s > 0) qfields(:,i_qs)=q7(ks:ke,i,j)
        if (nq_g > 0) qfields(:,i_qg)=q8(ks:ke,i,j)
        if (nq_i > 1) qfields(:,i_ni)=q9(ks:ke,i,j)
        if (nq_s > 1) qfields(:,i_ns)=q10(ks:ke,i,j)
        if (nq_g > 1) qfields(:,i_ng)=q11(ks:ke,i,j)
        if (nq_s > 2) qfields(:,i_m3s)=q12(ks:ke,i,j)
        if (nq_g > 2) qfields(:,i_m3g)=q13(ks:ke,i,j)
        dqfields(:, i_qv)=dqv(ks:ke,i,j)
        dqfields(:, i_th)=dth(ks:ke,i,j)

        if (nq_l > 0) dqfields(:,i_ql)=dq1(ks:ke,i,j)
        if (nq_r > 0) dqfields(:,i_qr)=dq2(ks:ke,i,j)
        if (nq_l > 1) dqfields(:,i_nl)=dq3(ks:ke,i,j)
        if (nq_r > 1) dqfields(:,i_nr)=dq4(ks:ke,i,j)
        if (nq_r > 2) dqfields(:,i_m3r)=dq5(ks:ke,i,j)
        if (nq_i > 0) dqfields(:,i_qi)=dq6(ks:ke,i,j)
        if (nq_s > 0) dqfields(:,i_qs)=dq7(ks:ke,i,j)
        if (nq_g > 0) dqfields(:,i_qg)=dq8(ks:ke,i,j)
        if (nq_i > 1) dqfields(:,i_ni)=dq9(ks:ke,i,j)
        if (nq_s > 1) dqfields(:,i_ns)=dq10(ks:ke,i,j)
        if (nq_g > 1) dqfields(:,i_ng)=dq11(ks:ke,i,j)
        if (nq_s > 2) dqfields(:,i_m3s)=dq12(ks:ke,i,j)
        if (nq_g > 2) dqfields(:,i_m3g)=dq13(ks:ke,i,j)

        if (aerosol_option > 0) then
          if (i_am1 >0) aerofields(:, i_am1)=a1(ks:ke,i,j)
          if (i_an1 >0) aerofields(:, i_an1)=a2(ks:ke,i,j)
          if (i_am2 >0) aerofields(:, i_am2)=a3(ks:ke,i,j)
          if (i_an2 >0) aerofields(:, i_an2)=a4(ks:ke,i,j)
          if (i_am3 >0) aerofields(:, i_am3)=a5(ks:ke,i,j)
          if (i_an3 >0) aerofields(:, i_an3)=a6(ks:ke,i,j)
          if (i_am4 >0) aerofields(:, i_am4)=a7(ks:ke,i,j)
          if (i_am5 >0) aerofields(:, i_am5)=a8(ks:ke,i,j)
          if (i_am6 >0) aerofields(:, i_am6)=a9(ks:ke,i,j)
          if (i_an6 >0) aerofields(:, i_an6)=a10(ks:ke,i,j)
          if (i_am7 >0) aerofields(:, i_am7)=a11(ks:ke,i,j)
          if (i_am8 >0) aerofields(:, i_am8)=a12(ks:ke,i,j)
          if (i_am9 >0) aerofields(:, i_am9)=a13(ks:ke,i,j)
          if (i_am10 >0) aerofields(:, i_am10)=a14(ks:ke,i,j)
          if (i_an10 >0) aerofields(:, i_an10)=a15(ks:ke,i,j)
          if (i_an11 >0) aerofields(:, i_an11)=a16(ks:ke,i,j)
          if (i_an12 >0) aerofields(:, i_an12)=a17(ks:ke,i,j)
          if (i_ak1 >0) aerofields(:, i_ak1)=a18(ks:ke,i,j)
          if (i_ak2 >0) aerofields(:, i_ak2)=a19(ks:ke,i,j)
          if (i_ak3 >0) aerofields(:, i_ak3)=a20(ks:ke,i,j)

          if (i_am1 >0) daerofields(:, i_am1)=da1(ks:ke,i,j)
          if (i_an1 >0) daerofields(:, i_an1)=da2(ks:ke,i,j)
          if (i_am2 >0) daerofields(:, i_am2)=da3(ks:ke,i,j)
          if (i_an2 >0) daerofields(:, i_an2)=da4(ks:ke,i,j)
          if (i_am3 >0) daerofields(:, i_am3)=da5(ks:ke,i,j)
          if (i_an3 >0) daerofields(:, i_an3)=da6(ks:ke,i,j)
          if (i_am4 >0) daerofields(:, i_am4)=da7(ks:ke,i,j)
          if (i_am5 >0) daerofields(:, i_am5)=da8(ks:ke,i,j)
          if (i_am6 >0) daerofields(:, i_am6)=da9(ks:ke,i,j)
          if (i_an6 >0) daerofields(:, i_an6)=da10(ks:ke,i,j)
          if (i_am7 >0) daerofields(:, i_am7)=da11(ks:ke,i,j)
          if (i_am8 >0) daerofields(:, i_am8)=da12(ks:ke,i,j)
          if (i_am9 >0) daerofields(:, i_am9)=da13(ks:ke,i,j)
          if (i_am10 >0) daerofields(:, i_am10)=da14(ks:ke,i,j)
          if (i_an10 >0) daerofields(:, i_an10)=da15(ks:ke,i,j)
          if (i_an11 >0) daerofields(:, i_an11)=da16(ks:ke,i,j)
          if (i_an12 >0) daerofields(:, i_an12)=da17(ks:ke,i,j)
        end if

        !--------------------------------------------------
        ! set fields which will not be modified
        !--------------------------------------------------

        call set_passive_fields(dt, rho(ks:ke,i,j),    &
             pressure(ks:ke,i,j), exner(ks:ke,i,j),            &
             z_half(ks-1:ke,i,j), z_centre(ks:ke,i,j), dz(ks:ke,i,j),     &
             w(ks:ke,i,j), tke(ks:ke,i,j), qfields )
        
        !--------------------------------------------------
        ! Do the business...
        !--------------------------------------------------
        call microphysics_common(dt, i , j, qfields, cffields, dqfields, tend, procs &
             , precip(i,j), precip_l, precip_r, precip_i, precip_s, precip_g       &
             , precip_r1d, precip_s1d, precip_so1d, precip_g1d                     &
             , aerophys, aerochem, aeroact                                         &
             , dustphys, dustchem, dustact                                         &
             , aeroice, dustliq                                                    &
             , aerofields, daerofields, aerosol_tend, aerosol_procs                &
             , rhcrit_1d)

        !--------------------------------------------------
        ! Relate back tendencies
        ! Check indices in mphys_switches that the appropriate
        ! fields are being passed back to mphys_casim
        !--------------------------------------------------
        dqv(ks:ke,i,j)=tend(:,i_qv)
        dth(ks:ke,i,j)=tend(:,i_th)
        dq1(ks:ke,i,j)=tend(:,i_ql)
        dq2(ks:ke,i,j)=tend(:,i_qr)

        if (cloud_params%l_2m) dq3(ks:ke,i,j)=tend(:,i_nl)
        if (rain_params%l_2m) dq4(ks:ke,i,j)=tend(:,i_nr)
        if (rain_params%l_3m) dq5(ks:ke,i,j)=tend(:,i_m3r)

        if (.not. l_warm) then
          if (ice_params%l_1m) dq6(ks:ke,i,j)=tend(:,i_qi)
          if (snow_params%l_1m) dq7(ks:ke,i,j)=tend(:,i_qs)
          if (graupel_params%l_1m) dq8(ks:ke,i,j)=tend(:,i_qg)
          if (ice_params%l_2m) dq9(ks:ke,i,j)=tend(:,i_ni)
          if (snow_params%l_2m) dq10(ks:ke,i,j)=tend(:,i_ns)
          if (graupel_params%l_2m) dq11(ks:ke,i,j)=tend(:,i_ng)
          if (snow_params%l_3m) dq12(ks:ke,i,j)=tend(:,i_m3s)
          if (graupel_params%l_3m) dq13(ks:ke,i,j)=tend(:,i_m3g)
        end if


        if (l_process) then
          if (i_am1 >0) da1(ks:ke,i,j)=aerosol_tend(:,i_am1)
          if (i_an1 >0) da2(ks:ke,i,j)=aerosol_tend(:,i_an1)
          if (i_am2 >0) da3(ks:ke,i,j)=aerosol_tend(:,i_am2)
          if (i_an2 >0) da4(ks:ke,i,j)=aerosol_tend(:,i_an2)
          if (i_am3 >0) da5(ks:ke,i,j)=aerosol_tend(:,i_am3)
          if (i_an3 >0) da6(ks:ke,i,j)=aerosol_tend(:,i_an3)
          if (i_am4 >0) da7(ks:ke,i,j)=aerosol_tend(:,i_am4)
          if (i_am5 >0) da8(ks:ke,i,j)=aerosol_tend(:,i_am5)
          if (i_am6 >0) da9(ks:ke,i,j)=aerosol_tend(:,i_am6)
          if (i_an6 >0) da10(ks:ke,i,j)=aerosol_tend(:,i_an6)
          if (i_am7 >0) da11(ks:ke,i,j)=aerosol_tend(:,i_am7)
          if (i_am8 >0) da12(ks:ke,i,j)=aerosol_tend(:,i_am8)
          if (i_am9 >0) da13(ks:ke,i,j)=aerosol_tend(:,i_am9)
          if (i_am10 >0) da14(ks:ke,i,j)=aerosol_tend(:,i_am10)
          if (i_an10 >0) da15(ks:ke,i,j)=aerosol_tend(:,i_an10)
          if (i_an11 >0) da16(ks:ke,i,j)=aerosol_tend(:,i_an11)
          if (i_an12 >0) da17(ks:ke,i,j)=aerosol_tend(:,i_an12)
        else
          da1(ks:ke,i,j)=0.0
          da2(ks:ke,i,j)=0.0
          da3(ks:ke,i,j)=0.0
          da4(ks:ke,i,j)=0.0
          da5(ks:ke,i,j)=0.0
          da6(ks:ke,i,j)=0.0
          da7(ks:ke,i,j)=0.0
          da9(ks:ke,i,j)=0.0
          da10(ks:ke,i,j)=0.0
          da11(ks:ke,i,j)=0.0
          da12(ks:ke,i,j)=0.0
          da13(ks:ke,i,j)=0.0
          da14(ks:ke,i,j)=0.0
          da15(ks:ke,i,j)=0.0
          da16(ks:ke,i,j)=0.0
          da17(ks:ke,i,j)=0.0
        end if


        if ( l_warm ) then

          if ( casdiags % l_surface_cloud ) casdiags % SurfaceCloudR(i,j) = precip_l
          if ( casdiags % l_surface_rain ) casdiags % SurfaceRainR(i,j)  = precip_r
          if ( casdiags % l_surface_snow ) casdiags % SurfaceSnowR(i,j)  = 0.0
          if ( casdiags % l_surface_graup) casdiags % SurfaceGraupR(i,j) = 0.0

          if ( casdiags % l_rainfall_3d ) casdiags % rainfall_3d(i,j,ks:ke)  = precip_r1d(:)
          if ( casdiags % l_snowfall_3d ) casdiags % snowfall_3d(i,j,ks:ke)  = 0.0
          if ( casdiags % l_snowonly_3d ) casdiags % snowonly_3d(i,j,ks:ke)  = 0.0
          if ( casdiags % l_graupfall_3d) casdiags % graupfall_3d(i,j,ks:ke) = 0.0


        else ! l_warm

          if ( casdiags % l_surface_rain ) casdiags % SurfaceRainR(i,j)  = precip_r
          if ( casdiags % l_surface_snow ) casdiags % SurfaceSnowR(i,j)  = precip_s
          if ( casdiags % l_surface_graup) casdiags % SurfaceGraupR(i,j) = precip_g

          if ( casdiags % l_rainfall_3d ) casdiags % rainfall_3d(i,j,ks:ke)  = precip_r1d(:)
          if ( casdiags % l_snowfall_3d ) casdiags % snowfall_3d(i,j,ks:ke)  = precip_s1d(:)
          if ( casdiags % l_snowonly_3d ) casdiags % snowonly_3d(i,j,ks:ke)  = precip_so1d(:)
          if ( casdiags % l_graupfall_3d) casdiags % graupfall_3d(i,j,ks:ke) = precip_g1d(:)

        end if ! l_warm

        if ( casdiags % l_radar ) then

          call tidy_qin(qfields)  !check this is conserving. If i do this here do we need a tidy_ain?
          call casim_reflec( nz, nq, rho(ks:ke,i,j), qfields, cffields,           &
                             dbz_tot_c, dbz_g_c, dbz_i_c,                      &
                             dbz_s_c,   dbz_l_c, dbz_r_c  )

          casdiags % dbz_tot(i,j, ks:ke) = dbz_tot_c(:)
          casdiags % dbz_g(i,j,   ks:ke) = dbz_g_c(:)
          casdiags % dbz_s(i,j,   ks:ke) = dbz_s_c(:)
          casdiags % dbz_i(i,j,   ks:ke) = dbz_i_c(:)
          casdiags % dbz_l(i,j,   ks:ke) = dbz_l_c(:)
          casdiags % dbz_r(i,j,   ks:ke) = dbz_r_c(:)

        end if ! casdiags % l_radar

        if ( casdiags % l_tendency_dg ) then
           DO k = ks, ke
              kc = k - ks + 1
              casdiags % dth_cond_evap(i,j,k) = procs(cloud_params%i_1m,i_cond%id)%column_data(kc) * &
                   Lv/cp * rexner(kc)
              casdiags % dqv_cond_evap(i,j,k) = -(procs(cloud_params%i_1m,i_cond%id)%column_data(kc))
              casdiags % dth_total(i,j,k) = tend(kc,i_th)
              casdiags % dqv_total(i,j,k) = tend(kc,i_qv)
              casdiags % dqc(i,j,k) = tend(kc,i_ql)
              casdiags % dqr(i,j,k) = tend(kc,i_qr)

              if (.not. l_warm) then
                 casdiags % dqi(i,j,k) = tend(kc,i_qi)
                 casdiags % dqs(i,j,k) = tend(kc,i_qs)
                 casdiags % dqg(i,j,k) = tend(kc,i_qg)
              endif
           enddo
        endif

        if ( casdiags % l_lwp ) then
           waterpath=0.0
           DO k = ks, ke
              waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_ql) )
           enddo
           casdiags % lwp(i,j)=waterpath
        endif
        if ( casdiags % l_rwp ) then
           waterpath=0.0
           DO k = ks, ke
              waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qr) )
           enddo
           casdiags % rwp(i,j)=waterpath
        endif
        if ( casdiags % l_iwp ) then
           waterpath=0.0
           DO k = ks, ke
              waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qi) )
           enddo
           casdiags % iwp(i,j)=waterpath
        endif
        if ( casdiags % l_swp ) then
           waterpath=0.0
           DO k = ks, ke
              waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qs) )
           enddo
           casdiags % swp(i,j)=waterpath
        endif
        if ( casdiags % l_gwp ) then
           waterpath=0.0
           DO k = ks, ke
              waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qg) )
           enddo
           casdiags % gwp(i,j)=waterpath
        endif

     end do ! i
  end do   ! j

! #if DEF_MODEL==MODEL_KiD
!     call save_dg(sum(casdiags % SurfaceRainR(:, :))/nxny, 'precip', i_dgtime)
!     call save_dg(sum(casdiags % SurfaceRainR(:, :))/nxny*3600.0, 'surface_precip_mmhr', i_dgtime)
! #endif

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine shipway_microphysics

  subroutine microphysics_common(dt, ix, jy, qfields, cffields, dqfields, tend & 
       , procs, precip, precip_l, precip_r, precip_i, precip_s, precip_g      & 
       , precip_r1d, precip_s1d, precip_so1d, precip_g1d                      &
       , aerophys, aerochem, aeroact                                          &
       , dustphys, dustchem, dustact                                          &
       , aeroice, dustliq                                                     &
       , aerofields, daerofields, aerosol_tend, aerosol_procs                 &
       , rhcrit_1d)


    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    use casim_parent_mod, only: casim_parent, parent_um

    implicit none

    real(wp), intent(in) :: dt  ! timestep from parent model
    integer, intent(in) :: ix, jy
    real(wp), intent(in) :: rhcrit_1d(:)
    real(wp), intent(inout) :: qfields(:,:), dqfields(:,:), tend(:,:)
    real(wp), intent(in) :: cffields(:,:)
    
    type(process_rate), intent(inout) :: procs(:,:)
    real(wp), intent(out) :: precip
    real(wp), intent(INOUT) :: precip_l
    real(wp), INTENT(INOUT) :: precip_r
    real(wp), intent(INOUT) :: precip_i
    real(wp), INTENT(INOUT) :: precip_s
    real(wp), intent(inout) :: precip_g

    real(wp), INTENT(INOUT) :: precip_r1d(:)
    real(wp), INTENT(INOUT) :: precip_s1d(:)
    real(wp), INTENT(INOUT) :: precip_so1d(:)
    real(wp), INTENT(INOUT) :: precip_g1d(:)

    real(wp) :: mindz

    ! Aerosol fields
    type(aerosol_phys), intent(inout)   :: aerophys(:)
    type(aerosol_chem), intent(in)      :: aerochem(:)
    type(aerosol_active), intent(inout) :: aeroact(:)
    type(aerosol_phys), intent(inout)   :: dustphys(:)
    type(aerosol_chem), intent(in)      :: dustchem(:)
    type(aerosol_active), intent(inout) :: dustact(:)

    type(aerosol_active), intent(inout) :: aeroice(:)
    type(aerosol_active), intent(inout) :: dustliq(:)

    real(wp), intent(inout) :: aerofields(:,:), daerofields(:,:), aerosol_tend(:,:)
    type(process_rate), intent(inout), optional :: aerosol_procs(:,:)

    real(wp) :: step_length

    real(wp) :: sed_length, sed_length_cloud, sed_length_rain, sed_length_ice, sed_length_snow &
      , sed_length_graupel

    integer :: n, k, nsed, iq

    logical :: l_Twarm   ! temperature above freezing

    integer, parameter :: level1 = 1

    ! Local working precipitation rates
    real(wp) :: precip_l_w ! Liquid cloud precip
    real(wp) :: precip_r_w ! Rain precip
    real(wp) :: precip_i_w ! Ice precip
    real(wp) :: precip_g_w ! Graupel precip
    real(wp) :: precip_s_w ! Snow precip

    !AH - note that nz is derived in mphys_init and accounts for the lowest level
    !     not equal to 1
    real(wp) :: precip1d(nz) ! local working precip rate
    real(wp) :: precip_l_w1d(nz) ! Liquid cloud precip 1D
    real(wp) :: precip_r_w1d(nz) ! Rain precip 1D
    real(wp) :: precip_i_w1d(nz) ! Ice precip 1D
    real(wp) :: precip_g_w1d(nz) ! Graupel precip 1D
    real(wp) :: precip_s_w1d(nz) ! Snow precip

    integer :: nsubsteps, nsubseds

    real :: inv_nsubsteps, inv_nsubseds, inv_allsubs
    ! inverse number of substeps for each hydrometor
    real :: inv_nsubseds_cloud, inv_nsubseds_ice, inv_nsubseds_rain,           &
      inv_nsubseds_snow, inv_nsubseds_graupel
    real :: inv_allsubs_cloud, inv_allsubs_ice, inv_allsubs_rain,              &
      inv_allsubs_snow, inv_allsubs_graupel

    ! number of substeps for each hydrometor
    integer :: nsubseds_cloud, nsubseds_ice, nsubseds_rain,                    &
      nsubseds_snow, nsubseds_graupel

    character(len=*), parameter :: RoutineName='MICROPHYSICS_COMMON'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    qfields_in=qfields ! Initial values of q
    qfields_mod=qfields ! Modified initial values of q (may be modified if bad values sent in)

    !! AH - Derive the number of microphysical substeps. The default 
    !!      max_step_length = 10.0 s in MONC (default) and 120 s in the UM. 
    !!      This is set in mphys_switches 

    nsubsteps=max(1, ceiling(dt/max_step_length))
    step_length=dt/nsubsteps
    inv_nsubsteps = 1.0 / real(nsubsteps)

    
    !! AH - Derive the maximum number of sedimentation substeps, which
    !!      are performed every microphysics step. max_sed_length = 2.0 s for MONC 
    !!      and 120 s for the UM and is set in mphys_switches. 
    !!      If step_length (the microphysical timestep) is longer than max_sed_length 
    !!      then the sedimentation will substep

    nsubseds=max(1, ceiling(step_length/max_sed_length))
    sed_length=step_length/nsubseds
    inv_nsubseds = 1.0 / real(nsubseds)

    inv_allsubs = 1.0 / real( nsubseds * nsubsteps)
    
    if (l_subseds_maxv) then
       if (casim_parent == parent_um) then 
          mindz = 20.0
       else
          mindz = min_dz ! derived in passive_fields
       endif
       ! cloud
       !print *, 'terminal vt called'
       call terminal_velocity_CFL(step_length, cloud_params%maxv, nsubseds_cloud, &
            sed_length_cloud, nsubseds, sed_length, mindz)
       inv_nsubseds_cloud = 1.0 / real(nsubseds_cloud)
       inv_allsubs_cloud = 1.0 / real(nsubseds_cloud * nsubsteps)
       ! rain
       call terminal_velocity_CFL(step_length, rain_params%maxv, nsubseds_rain, &
            sed_length_rain, nsubseds, sed_length, mindz)
       inv_nsubseds_rain = 1.0 / real(nsubseds_rain) 
       inv_allsubs_rain = 1.0 / real(nsubseds_rain * nsubsteps)
       if (.not. l_warm_loc) then
          ! ice
          call terminal_velocity_CFL(step_length, ice_params%maxv, nsubseds_ice, &
               sed_length_ice, nsubseds, sed_length, mindz)
          inv_nsubseds_ice = 1.0 / real(nsubseds_ice)
          inv_allsubs_ice = 1.0 / real(nsubseds_ice * nsubsteps)
          ! snow
          call terminal_velocity_CFL(step_length, snow_params%maxv, nsubseds_snow, &
               sed_length_snow, nsubseds, sed_length, mindz)
          inv_nsubseds_snow = 1.0 / real(nsubseds_snow)
          inv_allsubs_snow = 1.0 / real(nsubseds_snow * nsubsteps)
          ! graupel
          call terminal_velocity_CFL(step_length, graupel_params%maxv, nsubseds_graupel, &
               sed_length_graupel, nsubseds, sed_length, mindz)
          inv_nsubseds_graupel = 1.0 / real(nsubseds_graupel)
          inv_allsubs_graupel = 1.0 / real(nsubseds_graupel * nsubsteps)
          !print *, nsubseds_cloud,nsubseds_rain,nsubseds_ice,nsubseds_snow,nsubseds_graupel
          !print *, sed_length_cloud,sed_length_rain,sed_length_ice,sed_length_snow,sed_length_graupel
       endif
    else
       nsubseds_cloud = nsubseds
       sed_length_cloud = sed_length 
       inv_nsubseds_cloud = inv_nsubseds
       inv_allsubs_cloud = inv_allsubs
       ! rain
       nsubseds_rain = nsubseds
       sed_length_rain = sed_length 
       inv_nsubseds_rain = inv_nsubseds 
       inv_allsubs_rain = inv_allsubs
       if (.not. l_warm_loc) then
          ! ice
          nsubseds_ice = nsubseds
          sed_length_ice = sed_length 
          inv_nsubseds_ice = inv_nsubseds
          inv_allsubs_ice = inv_allsubs
          ! snow
          nsubseds_snow = nsubseds
          sed_length_snow = sed_length 
          inv_nsubseds_snow = inv_nsubseds
          inv_allsubs_snow = inv_allsubs
          ! graupel
          nsubseds_graupel = nsubseds
          sed_length_graupel = sed_length 
          inv_nsubseds_graupel = inv_nsubseds
          inv_allsubs_graupel = inv_allsubs
       endif
    endif
     
    if (l_tendency_loc) then! Parent model uses tendencies
      qfields_mod=qfields_in+dt*dqfields
    else! Parent model uses increments
      qfields_mod=qfields_in+dqfields
    end if

    if (.not. l_passive) then
      call tidy_qin(qfields_mod)
    end if

    !---------------------------------------------------------------
    ! Determine (and possibly limit) size distribution
    !---------------------------------------------------------------
    call query_distributions(cloud_params, qfields_mod, cffields, icall=1)
    call query_distributions(rain_params, qfields_mod, cffields, icall=1)
    if (.not. l_warm_loc) then
      call query_distributions(ice_params, qfields_mod, cffields, icall=1)
      call query_distributions(snow_params, qfields_mod, cffields, icall=1)
      call query_distributions(graupel_params, qfields_mod, cffields, icall=1)
    end if

    qfields=qfields_mod

    if (aerosol_option > 0) then
      aerofields_in=aerofields ! Initial values of aerosol
      aerofields_mod=aerofields ! Modified initial values  (may be modified if bad values sent in)

      if (l_tendency_loc) then! Parent model uses tendencies
        aerofields_mod=aerofields_in+dt*daerofields
      else! Parent model uses increments
        aerofields_mod=aerofields_in+daerofields
      end if

      if (l_process) call tidy_ain(qfields_mod, aerofields_mod)

      aerofields=aerofields_mod
    end if

    do n=1,nsubsteps

      call preconditioner(qfields)

      if ( casdiags % l_mphys_pts ) then
        ! Set microphysics points flag based on precondition
        do k = 1, nz
          casdiags % mphys_pts(ix, jy, k) = precondition(k)
        end do
      end if ! casdiags % l_mphys_pts

      !------------------------------------------------------
      ! Early exit if we will have nothing to do.
      ! (i.e. no hydrometeors and subsaturated)
      !------------------------------------------------------
      if (.not. any(precondition(:))) exit

      !-------------------------------
      ! Derive aerosol distribution
      ! parameters
      !-------------------------------
      if (aerosol_option > 0)  call examine_aerosol(aerofields, qfields, aerophys, aerochem, &
           aeroact, dustphys, dustchem, dustact, aeroice, dustliq, icall=1)

      do k=1,nz

          l_Twarm=TdegK(k) > 273.15
          l_Tcold(k)=.not. l_Twarm
       
       enddo
       !
       !=================================
       !
       ! WARM MICROPHYSICAL PROCESSES....
       !
       !=================================
       !
       !-------------------------------
       ! Do the autoconversion to rain
       !-------------------------------
       if (pswitch%l_praut) then
          call raut(step_length, qfields, cffields, aerofields, procs, aerosol_procs)
       end if

       !-------------------------------
       ! Do the rain accreting cloud
       !-------------------------------
       if (pswitch%l_pracw) then
          call racw(step_length, qfields, cffields, aerofields, procs, rain_params, aerosol_procs)
       end if

       !-------------------------------
       ! Do the rain self-collection
       !-------------------------------
       if (pswitch%l_pracr) then
          call racr(step_length, qfields, procs)
       end if

       !-------------------------------
       ! Do the evaporation of rain
       !-------------------------------
       if (pswitch%l_prevp) then
          !! initialise l_sigevap to false for all levels so that previous 
          !! "trues" are not included when rain_precondition is false
          l_sigevap(:) = .false.
          call revp(step_length, nz, qfields, aerophys, &
               aerochem, aeroact, dustliq, procs, aerosol_procs, l_sigevap)
       endif

       !=================================
       !
       ! ICE MICROPHYSICAL PROCESSES....
       !
       !=================================
       ! Start of all ice processes that occur at T < 0C
       !
       if (.not. l_warm_loc) then
          
          !------------------------------------------------------
          ! Condensation/immersion/contact nucleation of cloud ice
          !------------------------------------------------------
          if (pswitch%l_pinuc) then 
             call inuc(step_length, nz, l_Tcold, qfields, cffields, procs,     &
                  dustphys, aeroact, dustliq, aerosol_procs)
          end if

          !------------------------------------------------------
          ! Autoconverion to snow
          !------------------------------------------------------
          if (pswitch%l_psaut .and. .not. l_kfsm) then
             call saut(step_length, nz, l_Tcold, qfields, procs)
          end if
          !------------------------------------------------------
          ! Accretion processes
          !------------------------------------------------------
          ! Ice -> Cloud -> Ice
          if (pswitch%l_piacw) then 
             call iacc(step_length, nz, l_Tcold, ice_params, cloud_params, ice_params, qfields, &
                  cffields, procs, l_sigevap, aeroact, dustliq, aerosol_procs)
          end if
          ! Snow -> Cloud -> Snow
          if (l_sg) then
             if (pswitch%l_psacw .and. .not. l_kfsm) then
                call iacc(step_length, nz,  l_Tcold,snow_params, cloud_params, snow_params, &
                     qfields, cffields, procs, l_sigevap, aeroact, dustliq, aerosol_procs)
             end if
             !
             ! Snow -> Ice -> Snow
             if (pswitch%l_psaci .and. .not. l_kfsm) then
                call iacc(step_length, nz, l_Tcold, snow_params, ice_params, snow_params, qfields, &
                     cffields, procs, l_sigevap, aeroact, dustliq, aerosol_procs)
             end if
             !
             if (pswitch%l_praci) then
                ! Rain -> Ice -> Graupel AND Rain -> Ice -> snow, decision made in iacc
                call iacc(step_length, nz,  l_Tcold, rain_params, ice_params, graupel_params, &
                     qfields, cffields, procs, l_sigevap, aeroact, dustliq, aerosol_procs, snow_params)
                ! only one call needed and the decision of graupel or snow is made within iacc
                ! NOTE: this will break kfsm!!
                
             end if
             if (pswitch%l_psacr .and. .not. l_kfsm) then 
                ! Snow -> Rain -> Graupel AND Snow -> Rain -> Snow, decision made in iacc
                call iacc(step_length, nz,  l_Tcold, snow_params, rain_params, graupel_params, &
                     qfields, cffields, procs, l_sigevap, aeroact, dustliq, aerosol_procs, &
                     snow_params)
             end if
             if (l_g) then
                ! Graupel -> Cloud -> Graupel
                if (pswitch%l_pgacw) then
                   call iacc(step_length, nz,  l_Tcold, graupel_params, cloud_params, &
                        graupel_params, qfields, cffields, procs, l_sigevap, aeroact, dustliq, &
                        aerosol_procs)
                end if
                  ! Graupel -> Rain -> Graupel
                if (pswitch%l_pgacr) then
                   call iacc(step_length, nz,  l_Tcold, graupel_params, rain_params, &
                        graupel_params, qfields, cffields,  &
                        procs, l_sigevap, aeroact, dustliq, aerosol_procs)
                end if
                ! Graupel -> Ice -> Graupel
                !                   if(pswitch%l_gsaci)call iacc(step_length, k, graupel_params, ice_params, graupel_params, qfields, &
                !                       procs, aeroact, dustliq, aerosol_procs)
                ! Graupel -> Snow -> Graupel
                !                   if(pswitch%l_gsacs)call iacc(step_length, k, graupel_params, snow_params, graupel_params, qfields, &
                !                       procs, aeroact, dustliq, aerosol_procs)
             end if
          end if
          
          !------------------------------------------------------
          ! Small snow accreting cloud should be sent to graupel
          ! (Ikawa & Saito 1991)
          !------------------------------------------------------
          if (.not. l_kfsm) then
             ! Only do this process when Kalli's single moment code is
             ! not in use; otherwise we ignore it.
             if (l_g .and. .not. l_onlycollect) then
                call graupel_embryos(step_length, nz, l_Tcold, qfields, cffields, & 
                     procs)
             end if ! l_g
          end if ! not l_kfsm and precondition
          
          !------------------------------------------------------
          ! Wet deposition/shedding (resulting from graupel
          ! accretion processes)
          ! NB This alters some of the accretion processes, so
          ! must come after their calculation and before they
          ! are used/rescaled elsewhere
          !------------------------------------------------------
          if (l_g .and. .not. l_onlycollect) then  
             call wetgrowth(nz, l_Tcold, qfields, cffields, &
                  procs, l_sigevap)
          end if
          
          !------------------------------------------------------
          ! Aggregation (self-collection)
          !------------------------------------------------------
          if (pswitch%l_psagg .and. .not. l_kfsm) then 
             call ice_aggregation(nz, l_Tcold, snow_params, qfields, procs)
          end if

          !------------------------------------------------------
          ! Break up (snow only)
          !------------------------------------------------------
          if (pswitch%l_psbrk .and. .not. l_kfsm) then 
             call ice_breakup(nz, l_Tcold, snow_params, qfields, procs)
          end if

          !------------------------------------------------------
          ! Ice multiplication (Hallet-mossop)
          !------------------------------------------------------
          if (pswitch%l_pihal .and. .not. l_kfsm) then 
             call hallet_mossop(step_length, nz, cffields, &
                  procs)
          end if
          
          !------------------------------------------------------
          ! Homogeneous freezing (rain and cloud)
          !------------------------------------------------------
          if (pswitch%l_phomr) then 
             call ihom_rain(step_length, nz, l_Tcold, qfields, &
                  l_sigevap, aeroact, dustliq, procs, aerosol_procs)
          end if

          if (pswitch%l_phomc) then
             call ihom_droplets(step_length, nz, l_Tcold, qfields, aeroact, & 
                  dustliq, procs, aerosol_procs)
          endif
          !------------------------------------------------------
          ! Deposition/sublimation of ice/snow/graupel
          !------------------------------------------------------
          if (pswitch%l_pidep) then
             call idep(step_length, nz,  l_Tcold, ice_params, qfields,  &
                  procs, dustact, aeroice, aerosol_procs)
          end if

          if (pswitch%l_psdep .and. .not. l_kfsm ) then
             call idep(step_length, nz, l_Tcold,  snow_params, qfields, &
                  procs, dustact, aeroice, aerosol_procs)
          end if

          if (pswitch%l_pgdep) then 
             call idep(step_length, nz,  l_Tcold, graupel_params, qfields, &
                  procs, dustact, aeroice, aerosol_procs)
          end if

          if (l_harrington .and. .not. l_onlycollect) then  
             call adjust_dep(nz,  l_Tcold, procs)
          end if

          !-----------------------------------------------------------
          ! Make sure we don't remove more than saturation allows
          !-----------------------------------------------------------
          if (l_idep) then 
             call ensure_saturated(nz, l_Tcold, step_length, qfields, procs, &
               (/i_idep, i_sdep, i_gdep/))
          end if
          !-----------------------------------------------------------
          ! Make sure we don't put back more than saturation allows
          !-----------------------------------------------------------
          if (l_isub) then 
             call ensure_saturated(nz, l_Tcold,  step_length, qfields, procs, &
                  (/i_isub, i_ssub, i_gsub/))
          end if
          ! END all processes at T < 0C
          !
          ! start all ice processes that occur T > 0C, i.e. melting
              !------------------------------------------------------
              ! Melting of ice/snow/graupel
              !------------------------------------------------------
          if (pswitch%l_psmlt .and. .not. l_kfsm) then
             call melting(step_length, nz, snow_params, qfields, &
                  procs, l_sigevap, aeroice, dustact, aerosol_procs)
          end if
          if (pswitch%l_pgmlt) then
             call melting(step_length, nz, graupel_params, qfields, &
                  procs, l_sigevap, aeroice, dustact, aerosol_procs)
          end if
          if (pswitch%l_pimlt) then  
             call melting(step_length, nz, ice_params, qfields, &
                  procs, l_sigevap, aeroice, dustact, aerosol_procs)
          end if
      end if ! end if .not. warm_loc
      !-----------------------------------------------------------
      ! Make sure we don't remove more than we have to start with
      !-----------------------------------------------------------
        if (.not. l_warm_loc) then
         if (l_pos1) call ensure_positive(nz, step_length, qfields, procs, cloud_params, &
                (/i_praut, i_pracw, i_iacw, i_sacw, i_gacw, i_homc, i_inuc/))
              
         if (l_pos2) call ensure_positive(nz, step_length, qfields, procs, ice_params, &
                (/i_raci, i_saci, i_gaci, i_saut, i_isub, i_imlt/), &
                (/i_ihal, i_gshd, i_inuc, i_homc, i_iacw, i_idep/))
           
         if (l_pos3) call ensure_positive(nz, step_length, qfields, procs, rain_params, &
                (/i_prevp, i_sacr, i_gacr, i_homr/), &
                (/i_praut, i_pracw, i_raci, i_gshd, i_smlt, i_gmlt/))
              
         if (l_pos4) call ensure_positive(nz, step_length, qfields, procs, snow_params, &
                (/i_gacs, i_smlt, i_sacr, i_ssub /), &
                (/i_sdep, i_sacw, i_saut, i_saci, i_raci, i_gshd, i_ihal/))              
        else
           if (pswitch%l_praut .and. pswitch%l_pracw) then
            if (l_pos5) call ensure_positive(nz, step_length, qfields, procs, cloud_params, (/i_praut, i_pracw/))  
                   
           end if

           if (pswitch%l_prevp) then
            if (l_pos6) call ensure_positive(nz, step_length, qfields, procs, rain_params, (/i_prevp/), (/i_praut, i_pracw/))
           end if

        end if
      
      !-------------------------------
      ! Collect terms we have so far
      !-------------------------------

      call sum_procs(step_length, nz, procs, tend, (/i_praut, i_pracw, i_pracr, i_prevp/)     &
           , l_thermalexchange=.true., qfields=qfields , l_passive=l_passive)


      if (.not. l_warm_loc) then

        call sum_procs(step_length, nz, procs, tend,      &
             (/i_idep, i_sdep, i_gdep, i_iacw, i_sacr, i_sacw, i_saci, i_raci,&
             i_gacw, i_gacr, i_gaci, i_gacs, i_ihal, i_gshd, i_sbrk,&
             i_saut, i_sagg, i_isub, i_ssub, i_gsub/),        &
             l_thermalexchange=.true., qfields=qfields,&
             l_passive=l_passive, i_thirdmoment=2)
        call sum_procs(step_length, nz, procs, tend,      &
             (/i_inuc, i_imlt, i_smlt, i_gmlt, i_homr, i_homc/),     &
             l_thermalexchange=.true., qfields=qfields , l_passive=l_passive)
      end if

      call update_q(qfields_mod, qfields, tend,l_fixneg=.true.)

      if (l_process) then
        do k=1,nz
            call ensure_positive_aerosol(k, step_length, aerofields, aerosol_procs,&
                 (/i_aaut, i_aacw, i_aevp, i_arevp, i_dnuc, i_dsub, i_dssub, i_dgsub, &
                 i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr,    &
                 i_dgacr, i_draci  /) )
        end do

        call sum_aprocs(step_length, aerosol_procs, aerosol_tend,      &
             (/i_aaut, i_aacw, i_aevp, i_arevp/), aerofields)

        if (.not. l_warm) then
          call sum_aprocs(step_length, aerosol_procs, aerosol_tend,      &
               (/i_dnuc, i_dsub, i_dssub, i_dgsub, i_dimlt, i_dsmlt, i_dgmlt,     &
               i_diacw, i_dsacw, i_dgacw, i_dsacr, i_dgacr, i_draci /), aerofields)
        end if

        call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.,l_fixneg=.true.)

        !-------------------------------
        ! Re-Derive aerosol distribution
        ! parameters
        !-------------------------------
        call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact, &
             dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
      end if

      !-------------------------------
      ! Do the condensation/evaporation
      ! of cloud and activation of new
      ! drops
      !-------------------------------

       if (casim_parent == parent_um .and. l_prf_cfrac) then
        ! In um and using Paul's cloud fraction so just do update of fields


         !call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)
         call query_distributions(cloud_params, qfields, cffields, icall=2)
         call query_distributions(rain_params, qfields, cffields, icall=2)
         if (.not. l_warm_loc) then
           call query_distributions(ice_params, qfields, cffields, icall=2)
           call query_distributions(snow_params, qfields, cffields, icall=2)
           call query_distributions(graupel_params, qfields, cffields, icall=2)
         end if

       else ! Not in UM or not using Paul's cloud fraction in UM
        if (pswitch%l_pcond)then
            call condevp(step_length, nz, qfields,                 &
               procs, aerophys, aerochem, aeroact, dustphys, dustchem, dustliq, &
               aerosol_procs, rhcrit_1d)
        !-------------------------------
        ! Collect terms we have so far
        !-------------------------------
        call sum_procs(step_length, nz, procs, tend,      &
             (/i_cond/)     &
             , l_thermalexchange=.true., qfields=qfields     &
             , l_passive=l_passive)

          call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)
        end if ! pswitch%l_pcond

      !---------------------------------------------------------------
      ! Re-Determine (and possibly limit) size distribution
      !---------------------------------------------------------------
        call query_distributions(cloud_params, qfields, cffields, icall=2)
        call query_distributions(rain_params, qfields, cffields, icall=2)
        if (.not. l_warm_loc) then
          call query_distributions(ice_params, qfields, cffields, icall=2)
          call query_distributions(snow_params, qfields, cffields, icall=2)
          call query_distributions(graupel_params, qfields, cffields, icall=2)
        end if

        if (l_process) then
          call sum_aprocs(step_length, aerosol_procs, aerosol_tend,      &
             (/i_aact /), aerofields)

          call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)

        !-------------------------------
        ! Re-Derive aerosol distribution
        ! parameters
        !-------------------------------
          call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact     &
             , dustphys, dustchem, dustact, aeroice, dustliq, icall=3)
        end if

      endif ! (casim_parent == parent_um .and. l_prf_cfrac)

!!PRF

      !-------------------------------
      ! Do the sedimentation
      !-------------------------------

      precip_l_w = 0.0
      precip_r_w = 0.0
      precip_i_w = 0.0
      precip_g_w = 0.0
      precip_s_w = 0.0

      do k = 1, nz
        precip_l_w1d(k) = 0.0
        precip_r_w1d(k) = 0.0
        precip_i_w1d(k) = 0.0
        precip_g_w1d(k) = 0.0
        precip_s_w1d(k) = 0.0
      end do

      if ( casdiags % l_process_rates ) then
         call gather_process_diagnostics(ix, jy, ks, ke,ncall=0)
      end if

      if (l_sed) then
         if (.not. l_subseds_maxv) then ! need to add check for 3M code 

      !! AH - Following block of code performs sedimentation using the standard (original) 
      !!      method, where number of substeps for all hydrometeors are derived using 
      !!      max_sed_length and this substep is applied to all hydrometeors
      !!
            do nsed=1,nsubseds

               if (nsed > 1) then
                  !-------------------------------
                  ! Reset process rates if they
                  ! are to be re-used
                  !-------------------------------
                  !call zero_procs_exp(procs)
                  call zero_procs(procs)
                  if (l_process) call zero_procs(aerosol_procs)
                  !---------------------------------------------------------------
                  ! Re-Determine (and possibly limit) size distribution
                  !---------------------------------------------------------------
                  call query_distributions(cloud_params, qfields, cffields, icall=nsed+2)
                  call query_distributions(rain_params, qfields, cffields, icall=nsed+2)
                  if (.not. l_warm_loc) then
                     call query_distributions(ice_params, qfields, cffields, icall=nsed+2)
                     call query_distributions(snow_params, qfields, cffields, icall=nsed+2)
                     call query_distributions(graupel_params, qfields, cffields, icall=nsed+2)
                  end if

                  !-------------------------------
                  ! Re-Derive aerosol distribution
                  ! parameters
                  !-------------------------------
                  if (l_process) call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact &
                       , dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
               end if
               
               ! NOTE: if l_gamma_online is true then the original CASIM sedimentation will be used, 
               !       which will calculate the gamma function every timestep. This has to be done 
               !       when 3M or diagnostic shape is used, as shape will change. For single and 
               !       and double moment simulations, it is recommended that l_gamma_online is false. 
               !       This is computationally much more efficient (on CPU and GPU).
               if (pswitch%l_psedl) then
                  if (l_gamma_online) then
                     call sedr(qfields, aeroact, dustliq, cloud_params, &
                          procs, aerosol_procs, precip1d, l_process)
                  else
                     call sedr_1M_2M(sed_length, qfields, aeroact, dustliq, cloud_params, &
                          procs, aerosol_procs, precip1d, l_process)
                  endif
                  
                  precip_l_w = precip_l_w + precip1d(level1)
                  
                  do k = 1, nz
                     precip_l_w1d(k) = precip_l_w1d(k) + precip1d(k)
                  end do
                  
                  call sum_procs(sed_length, nz, procs, tend, (/i_psedl/), qfields=qfields)
                  
               end if

               if (pswitch%l_psedr) then
                  if (l_gamma_online) then
                     call sedr(qfields, aeroact, dustliq, rain_params, &
                          procs, aerosol_procs, precip1d, l_process)
                  else
                     call sedr_1M_2M(sed_length, qfields, aeroact, dustliq, rain_params, &
                          procs, aerosol_procs, precip1d, l_process)
                  endif
                  
                  precip_r_w = precip_r_w + precip1d(level1)
                  
                  do k = 1, nz
                     precip_r_w1d(k) = precip_r_w1d(k) + precip1d(k)
                  end do

                  call sum_procs(sed_length, nz, procs, tend, (/i_psedr/), qfields=qfields)

               end if

               if (.not. l_warm_loc) then
                  
                  if (pswitch%l_psedi) then
                     if (l_gamma_online) then 
                        call sedr(qfields, aeroice, dustact, ice_params, &
                             procs, aerosol_procs, precip1d, l_process)
                     else
                        call sedr_1M_2M(sed_length, qfields, aeroice, dustact, ice_params, &
                             procs, aerosol_procs, precip1d, l_process)
                     endif
                     precip_i_w = precip_i_w + precip1d(level1)
                     
                     do k = 1, nz
                        precip_i_w1d(k) = precip_i_w1d(k) + precip1d(k)
                     end do
                     
                     call sum_procs(sed_length, nz, procs, tend, (/i_psedi/), qfields=qfields)

                  end if
                  if (pswitch%l_pseds .and. .not. l_kfsm) then
                     if (l_gamma_online) then 
                        call sedr(qfields, aeroice, dustact, snow_params, &
                             procs, aerosol_procs, precip1d, l_process)
                     else
                        call sedr_1M_2M(sed_length, qfields, aeroice, dustact, snow_params, &
                             procs, aerosol_procs, precip1d, l_process) 
                     endif
                     precip_s_w = precip_s_w + precip1d(level1)

                     do k = 1, nz
                        precip_s_w1d(k) = precip_s_w1d(k) + precip1d(k)
                     end do
                     
                     call sum_procs(sed_length, nz, procs, tend, (/ i_pseds /), qfields=qfields)

                  end if
                  if (pswitch%l_psedg) then
                     if (l_gamma_online) then
                        call sedr(qfields,  aeroice, dustact, graupel_params, &
                             procs, aerosol_procs, precip1d, l_process)
                     else
                        call sedr_1M_2M(sed_length, qfields, aeroice, dustact, graupel_params, &
                             procs, aerosol_procs, precip1d, l_process)
                     endif
                     precip_g_w = precip_g_w + precip1d(level1)
                     
                     do k = 1, nz
                        precip_g_w1d(k) = precip_g_w1d(k) + precip1d(k)
                     end do
                     
                     call sum_procs(sed_length, nz, procs, tend, (/i_psedg/), qfields=qfields)

                  end if
                  
                  !!call sum_procs(sed_length, nz, procs, tend, (/i_psedi, i_pseds, i_psedg/), qfields=qfields)
               end if

               call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)

               if (l_process) then
                  call sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_asedr, i_asedl/),&
                       aerofields)
                  
                  call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
                  
                  if (l_process .and. .not. l_warm) then
                     call sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_dsedi, i_dseds, i_dsedg/),&
                          aerofields)
                     
                     call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
                  end if
               end if
               
               if ( casdiags % l_process_rates ) then
                  call gather_process_diagnostics(ix, jy, ks, ke,ncall=1)
               end if

            end do

         else ! l_subseds_maxv = true
      !! AH - Following block of code performs sedimentation using a CFL based on the 
      !!      prescribed max terminal velocity (mphys_params) for each hydrometeor This 
      !!      creates a number of substeps for each  hydrometeor and hence a loop for 
      !!      each hydrometeor
      !! 
            if (pswitch%l_psedl) then
               do nsed=1,nsubseds_cloud 
                  if (nsed > 1) then
                     !---------------------------------------------------------------
                     ! Re-Determine (and possibly limit) size distribution
                     !---------------------------------------------------------------
                     call query_distributions(cloud_params, qfields, cffields, icall=nsed+2)
                     !-------------------------------
                     ! Re-Derive aerosol distribution
                     ! parameters
                     !-------------------------------
                     if (l_process) call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact &
                          , dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
                     
                  endif
                  if (l_gamma_online) then
                     call sedr(qfields, aeroact, dustliq, cloud_params, &
                          procs, aerosol_procs, precip1d, l_process)
                  else
                     call sedr_1M_2M(sed_length_cloud, qfields, aeroact, dustliq, cloud_params, &
                       procs, aerosol_procs, precip1d, l_process)
                  endif
                  
                  precip_l_w = precip_l_w + precip1d(level1)
                  
                  do k = 1, nz
                     precip_l_w1d(k) = precip_l_w1d(k) + precip1d(k)
                  end do
                  
                  if ( casdiags % l_process_rates ) then
                     call gather_process_diagnostics(ix, jy, ks, ke,ncall=1)
                  end if
                  
                  call sum_procs(sed_length_cloud, nz, procs, tend, (/i_psedl/), qfields=qfields)

                  call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)

                  if (l_process) then
                     call sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_asedl/),&
                          aerofields)                  
                     call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
                  endif
                  ! !-------------------------------
                  ! ! Reset process rates if they
                  ! ! are to be re-used
                  ! !-------------------------------
                  call zero_procs(procs, (/i_psedl/)) 
                  if (l_process) call zero_procs(aerosol_procs, (/i_asedl/))
               enddo
            end if
               
            if (pswitch%l_psedr) then
               do nsed=1,nsubseds_rain                 
                  if (nsed > 1) then
                     !---------------------------------------------------------------
                     ! Re-Determine (and possibly limit) size distribution
                     !---------------------------------------------------------------
                     call query_distributions(rain_params, qfields, cffields, icall=nsed+2)
                     !-------------------------------
                     ! Re-Derive aerosol distribution
                     ! parameters
                     !-------------------------------
                     if (l_process) call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact &
                          , dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
                  endif
                  if (l_gamma_online) then
                     call sedr(qfields, aeroact, dustliq, rain_params, &
                          procs, aerosol_procs, precip1d, l_process)
                  else
                     call sedr_1M_2M(sed_length_rain, qfields, aeroact, dustliq, rain_params, &
                          procs, aerosol_procs, precip1d, l_process)
                  endif

                  precip_r_w = precip_r_w + precip1d(level1)
                  
                  do k = 1, nz
                     precip_r_w1d(k) = precip_r_w1d(k) + precip1d(k)
                  end do
                  
                  if ( casdiags % l_process_rates ) then
                     call gather_process_diagnostics(ix, jy, ks, ke,ncall=1)
                  end if
                  
                  call sum_procs(sed_length_rain, nz, procs, tend, (/i_psedr/), qfields=qfields)
                  call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)
                  if (l_process) then
                     call sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_asedr, i_asedl/),&
                          aerofields)
                     call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
                  endif
                  !-------------------------------
                  ! Reset process rates if they
                  ! are to be re-used
                  !-------------------------------
                  call zero_procs(procs, (/i_psedr/))
                  if (l_process) call zero_procs(aerosol_procs, (/i_asedr/))
                  
               enddo
            endif
                        
            if (.not. l_warm_loc) then
                  
               if (pswitch%l_psedi) then
                  do nsed=1,nsubseds_ice                 
                     if (nsed > 1) then
                        !---------------------------------------------------------------
                        ! Re-Determine (and possibly limit) size distribution
                        !---------------------------------------------------------------
                        call query_distributions(ice_params, qfields, cffields, icall=nsed+2)
                        !-------------------------------
                        ! Re-Derive aerosol distribution
                        ! parameters
                        !-------------------------------
                        if (l_process) call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact &
                              , dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
                     endif
                     if (l_gamma_online) then
                        call sedr(qfields, aeroice, dustact, ice_params, &
                             procs, aerosol_procs, precip1d, l_process)
                     else
                        call sedr_1M_2M(sed_length_ice, qfields, aeroice, dustact, ice_params, &
                             procs, aerosol_procs, precip1d, l_process)
                     endif
                     
                     precip_i_w = precip_i_w + precip1d(level1)
                     
                     do k = 1, nz
                        precip_i_w1d(k) = precip_i_w1d(k) + precip1d(k)
                     end do
                     
                     if ( casdiags % l_process_rates ) then
                        call gather_process_diagnostics(ix, jy, ks, ke,ncall=1)
                     end if
                     
                     call sum_procs(sed_length_ice, nz, procs, tend, (/i_psedi/), qfields=qfields)
                     call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)
                     if (l_process) then
                        call sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_dsedi/),&
                             aerofields)
                        call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
                     endif
                     !-------------------------------
                     ! Reset process rates if they
                     ! are to be re-used
                     !-------------------------------
                     call zero_procs(procs, (/i_psedi/))
                     if (l_process) call zero_procs(aerosol_procs, (/i_dsedi/))
                  enddo
               end if
               if (pswitch%l_pseds) then
                  do nsed=1,nsubseds_snow                 
                     if (nsed > 1) then
                        !---------------------------------------------------------------
                        ! Re-Determine (and possibly limit) size distribution
                        !---------------------------------------------------------------
                        call query_distributions(snow_params, qfields, cffields, icall=nsed+2)
                        ! !-------------------------------
                        ! ! Re-Derive aerosol distribution
                        ! ! parameters
                        ! !-------------------------------
                        if (l_process) call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact &
                              , dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
                     endif
                     
                     if (l_gamma_online) then
                        call sedr(qfields, aeroice, dustact, snow_params, &
                             procs, aerosol_procs, precip1d, l_process) 
                     else
                        call sedr_1M_2M(sed_length_snow, qfields, aeroice, dustact, snow_params, &
                             procs, aerosol_procs, precip1d, l_process) 
                     endif

                     precip_s_w = precip_s_w + precip1d(level1)
                     
                     do k = 1, nz
                        precip_s_w1d(k) = precip_s_w1d(k) + precip1d(k)
                     end do
                     
                     if ( casdiags % l_process_rates ) then
                        call gather_process_diagnostics(ix, jy, ks, ke,ncall=1)
                     end if
                     
                     call sum_procs(sed_length_snow, nz, procs, tend, (/i_pseds/), qfields=qfields)
                     call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)
                     if (l_process) then
                        call sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_dseds/),&
                             aerofields)
                        call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
                     endif
                     !-------------------------------
                     ! Reset process rates if they
                     ! are to be re-used
                     !-------------------------------
                     call zero_procs(procs, (/i_pseds/))
                     if (l_process) call zero_procs(aerosol_procs, (/i_dseds/))
                  enddo
               end if
               
               if (pswitch%l_psedg) then
                  do nsed=1,nsubseds_graupel                 
                     if (nsed > 1) then
                        !---------------------------------------------------------------
                        ! Re-Determine (and possibly limit) size distribution
                        !---------------------------------------------------------------
                        call query_distributions(graupel_params, qfields, cffields, icall=nsed+2)
                        !-------------------------------
                        ! Re-Derive aerosol distribution
                        ! parameters
                        !-------------------------------
                        if (l_process) call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact &
                             , dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
                     endif
                     if (l_gamma_online) then
                        call sedr(qfields, aeroice, dustact, &
                             graupel_params, procs, aerosol_procs, precip1d, l_process)
                     else
                        call sedr_1M_2M(sed_length_graupel, qfields, aeroice, dustact, &
                             graupel_params, procs, aerosol_procs, precip1d, l_process)
                     endif
                     
                     precip_g_w = precip_g_w + precip1d(level1)
                     
                     do k = 1, nz
                        precip_g_w1d(k) = precip_g_w1d(k) + precip1d(k)
                     end do
                     
                     if ( casdiags % l_process_rates ) then
                        call gather_process_diagnostics(ix, jy, ks, ke,ncall=1)
                     end if

                     call sum_procs(sed_length_graupel, nz, procs, tend, (/i_psedg/), &
                       qfields=qfields)
                     call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)
                     if (l_process) then
                        call sum_aprocs(sed_length_graupel, aerosol_procs, aerosol_tend, (/i_dsedg/),&
                             aerofields)
                     
                        call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
                     endif
                     !-------------------------------
                     ! Reset process rates if they
                     ! are to be re-used
                     !-------------------------------
                     call zero_procs(procs, (/i_psedg/))
                     if (l_process) call zero_procs(aerosol_procs, (/i_dsedg/))
                  enddo
               end if
            endif
         endif
      endif
 
      precip_l = precip_l + precip_l_w
      ! For diagnostic purposes, set precip_r, precip_s and precip to pass out
      ! For the UM, rainfall rate is assumed as sum of all liquid components
      ! (so includes sedimentation of rain and liquid cloud) 
      precip_r = precip_r + precip_r_w

      ! For the UM, snowfall rate is assumed to be a sum of all solid components
      ! (so includes ice, snow and graupel) 
      precip_s = precip_s + precip_s_w 
      precip_i = precip_i + precip_i_w 
      ! For the UM, graupel rate is just itself
      precip_g = precip_g + precip_g_w

      do k = 1, nz
        precip_r1d(k)  = precip_r1d(k)  + precip_l_w1d(k) + precip_r_w1d(k)
        precip_g1d(k)  = precip_g1d(k)  + precip_g_w1d(k)
        precip_s1d(k)  = precip_s1d(k)  + precip_s_w1d(k) + precip_i_w1d(k) + precip_g_w1d(k)
        precip_so1d(k) = precip_so1d(k) + precip_s_w1d(k) + precip_i_w1d(k)
      end do

      if (nsubsteps>1)then
        !-------------------------------
        ! Reset process rates if they
        ! are to be re-used
        !-------------------------------
        !call zero_procs_exp(procs) 
        call zero_procs(procs)
        if (l_process) call zero_procs(aerosol_procs)
        !---------------------------------------------------------------
        ! Re-Determine (and possibly limit) size distribution
        !---------------------------------------------------------------
        call query_distributions(cloud_params, qfields, cffields, icall=nsubseds+3)
        call query_distributions(rain_params, qfields, cffields, icall=nsubseds+3)
        if (.not. l_warm_loc) then
          call query_distributions(ice_params, qfields, cffields, icall=nsubseds+3)
          call query_distributions(snow_params, qfields, cffields, icall=nsubseds+3)
          call query_distributions(graupel_params, qfields, cffields, icall=nsubseds+3)
        end if
     end if
  end do

    ! We want the mean precipitation over the parent timestep - so divide
    ! by the total number of substeps - multiply by inv_allsubs is quicker
    precip_l = precip_l * inv_allsubs_cloud
    precip_r = precip_r * inv_allsubs_rain
    precip_i = precip_i * inv_allsubs_ice
    precip_s = precip_s * inv_allsubs_snow
    precip_g = precip_g * inv_allsubs_graupel

    ! UM precip rates are 
    precip_r = precip_l + precip_r
    precip_s = precip_i + precip_s + precip_g

    do k = 1, nz
      precip_r1d(k)  = precip_r1d(k)  * inv_allsubs
      precip_s1d(k)  = precip_s1d(k)  * inv_allsubs
      precip_so1d(k) = precip_so1d(k) * inv_allsubs
      precip_g1d(k)  = precip_g1d(k)  * inv_allsubs
    end do ! k

    ! Precip is a sum of everything, so just add rain and snow together which
    ! has all components added.
    ! Do not add precip_g, otherwise graupel contributions will be double-counted
    precip = precip_r   + precip_s

    !--------------------------------------------------
    ! Tidy up any small/negative numbers
    ! we may have generated.
    !--------------------------------------------------
    if (pswitch%l_tidy2) then

       call qtidy(step_length, nz, qfields, procs, aerofields, aeroact, dustact, &
            aeroice, dustliq , aerosol_procs, i_tidy2, i_atidy2, l_negonly=l_tidy_negonly)

      call sum_procs(step_length, nz, procs, tend, (/i_tidy2/), qfields=qfields, l_passive=l_passive)

      call update_q(qfields_mod, qfields, tend)

      if (l_process) then
        call sum_aprocs(step_length, aerosol_procs, aerosol_tend, (/i_atidy2/), aerofields )
        call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
      end if
    end if

    !
    ! Add on initial adjustments that may have been made
    !

    if (l_tendency_loc) then! Convert back from cumulative value to tendency
      tend=tend+qfields_mod-qfields_in-dqfields*dt
      tend=tend/dt
    else
      tend=tend+qfields_mod-qfields_in-dqfields
      ! prevent negative values
      do iq=i_hstart,ntotalq
        do k=1,nz
          tend(k,iq)=max(tend(k,iq), -(qfields_in(k,iq)-dqfields(k,iq)))
        end do
      end do
    end if

    if (aerosol_option > 0) then
      ! processing
      if (l_process) then
        if (l_tendency_loc) then! Convert back from cumulative value to tendency
          aerosol_tend=aerosol_tend+aerofields_mod-aerofields_in-daerofields*dt
          aerosol_tend=aerosol_tend/dt
        else
          aerosol_tend=aerosol_tend+aerofields_mod-aerofields_in-daerofields
          ! prevent negative values
          do iq=1,ntotala
            do k=1,nz
              aerosol_tend(k,iq)=max(aerosol_tend(k,iq), -(aerofields_in(k,iq)-daerofields(k,iq)))
            end do
          end do
        end if
      end if
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine microphysics_common

  subroutine update_q(qfields_in, qfields, tend, l_aerosol, l_fixneg)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(in) :: qfields_in(:,:)
    real(wp), intent(inout) :: qfields(:,:)
    real(wp), intent(in) :: tend(:,:)
    logical, intent(in), optional :: l_aerosol ! flag to indicate updating of aerosol
    logical, intent(in), optional :: l_fixneg  ! Flag to use cludge to bypass negative/zero numbers
    integer :: k, iqx
    logical :: l_fix

    character(len=*), parameter :: RoutineName='UPDATE_Q'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_fix=.false.
    if (present(l_fixneg)) l_fix=l_fixneg

    do iqx=1, ubound(tend,2)
       do k=lbound(tend,1), ubound(tend,1)
          qfields(k,iqx)=qfields_in(k,iqx)+tend(k,iqx)
       enddo
    enddo
    
     if (.not. present(l_aerosol) .and. l_fix) then
      !quick lem fixes  - this code should never be used ?
        do iqx=1, ubound(tend,2)
          do k=lbound(tend,1), ubound(tend,1)
             if (iqx==i_ni .and. qfields(k,iqx)<=0) then
                qfields(k,iqx)=0
                qfields(k,i_qi)=0
             end if
             if (iqx==i_nr .and. qfields(k,iqx)<=0) then
                qfields(k,iqx)=0
                qfields(k,i_qr)=0
                if (i_m3r/=0)qfields(k,i_m3r)=0
             end if
             if (iqx==i_nl .and. qfields(k,iqx)<=0) then
                qfields(k,iqx)=0
                qfields(k,i_ql)=0
             end if
             if (iqx==i_ns .and. qfields(k,iqx)<=0) then
                qfields(k,iqx)=0
                qfields(k,i_qs)=0
                if (i_m3s/=0)qfields(k,i_m3s)=0
             end if
             if (iqx==i_ng .and. qfields(k,iqx)<=0) then
                qfields(k,iqx)=0
                qfields(k,i_qg)=0
                if (i_m3g/=0)qfields(k,i_m3g)=0
             end if
          end do
       end do
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine update_q

  subroutine gather_process_diagnostics(i, j, ks, ke,ncall)

    ! Gathers all process rate diagnostics if in use and outputs them to the
    ! CASIM generic diagnostic fields, ready for use in any model.

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Indices of this particular grid square
    integer, intent(in) :: i, j,ncall
    integer, intent(in) :: ks, ke ! Start/end points of grid

    ! Local variables

    character(len=*), parameter :: RoutineName='GATHER_PROCESS_DIAGNOSTICS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    INTEGER :: kc ! Casim Z-level
    INTEGER :: k  ! Loop counter in z-direction

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Based on code from RGS: fill in process rates:

    if (ncall==0) THEN
    IF (casdiags % l_phomc) THEN
      IF (pswitch%l_phomc) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % phomc(i,j,k) = procs(ice_params%i_1m,i_homc%id)%column_data(kc)
        END DO
      ELSE
        casdiags % phomc(i,j,:) = ZERO_REAL_WP
      END IF
    END iF ! casdiags % l_phomc

    IF (casdiags % l_nhomc) THEN
      IF ((pswitch%l_phomc) .and. (ice_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % nhomc(i,j,k) = procs(ice_params%i_2m,i_homc%id)%column_data(kc)
        END DO
      ELSE
        casdiags % nhomc(i,j,:) = ZERO_REAL_WP
      END IF
    END iF ! casdiags % l_phomc

    IF (casdiags % l_pinuc) THEN
      IF (pswitch%l_pinuc) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pinuc(i,j,k) = procs(ice_params%i_1m,i_inuc%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pinuc(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! casdiags % l_pinuc

    IF (casdiags % l_ninuc) THEN
      IF ((pswitch%l_pinuc) .and. (ice_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % ninuc(i,j,k) = procs(ice_params%i_2m,i_inuc%id)%column_data(kc)
        END DO
      ELSE
        casdiags % ninuc (i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! casdiags % l_pinuc

    IF (casdiags % l_pidep) THEN
      IF (pswitch%l_pidep) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pidep(i,j,k) = procs(ice_params%i_1m,i_idep%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pidep(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psdep) THEN
      IF (pswitch%l_psdep) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psdep(i,j,k) = procs(snow_params%i_1m,i_sdep%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psdep(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_piacw) THEN
      IF (pswitch%l_piacw) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % piacw(i,j,k) = procs(ice_params%i_1m,i_iacw%id)%column_data(kc)
        END DO
      ELSE
        casdiags % piacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psacw) THEN
      IF (pswitch%l_psacw) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psacw(i,j,k) = procs(snow_params%i_1m,i_sacw%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psacr) THEN
      IF (pswitch%l_psacr) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psacr(i,j,k) = procs(snow_params%i_1m,i_sacr%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psacr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pisub) THEN
      IF (pswitch%l_pisub) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pisub(i,j,k) = -1.0 * procs(ice_params%i_1m,i_isub%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pisub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pssub) THEN
      IF (pswitch%l_pssub) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pssub(i,j,k) = -1.0 * procs(snow_params%i_1m,i_ssub%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pssub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF


    IF (casdiags % l_pimlt) THEN
      IF (pswitch%l_pimlt) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pimlt(i,j,k) = procs(rain_params%i_1m,i_imlt%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pimlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psmlt) THEN
      IF (pswitch%l_psmlt) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psmlt(i,j,k) = procs(rain_params%i_1m,i_smlt%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psaut) THEN
      IF (pswitch%l_psaut) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psaut(i,j,k) = procs(snow_params%i_1m,i_saut%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psaut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psaci) THEN
      IF (pswitch%l_psaci) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psaci(i,j,k) = procs(snow_params%i_1m,i_saci%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psaci(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_praut) THEN
      IF (pswitch%l_praut) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % praut(i,j,k) = procs(rain_params%i_1m,i_praut%id)%column_data(kc)
        END DO
      ELSE
        casdiags % praut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pracw) THEN
      IF (pswitch%l_pracw) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pracw(i,j,k) = procs(rain_params%i_1m,i_pracw%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pracw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF
    IF (casdiags % l_prevp) THEN
      IF (pswitch%l_prevp) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % prevp(i,j,k) = -1.0 * procs(rain_params%i_1m,i_prevp%id)%column_data(kc)
        END DO
      ELSE
        casdiags % prevp(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

! #if DEF_MODEL==MODEL_KiD
!      IF (casdiags % l_praut) THEN
!        IF (pswitch%l_praut) THEN
!           DO k = ks, ke
!              kc = k - ks + 1
!              call save_dg(k, casdiags % praut(i,j,k) , 'praut', i_dgtime)
!           END DO
!        ENDIF
!     ENDIF
    
!     IF (casdiags % l_pracw) THEN
!         IF (pswitch%l_pracw) THEN
!            DO k = ks, ke
!               kc = k - ks + 1
!               call save_dg(k, casdiags % pracw(i,j,k) , 'pracw', i_dgtime) 
!            END DO
!         END IF
!      END IF

!      IF (casdiags % l_prevp) THEN
!         IF (pswitch%l_prevp) THEN
!            DO k = ks, ke
!               kc = k - ks + 1
!               call save_dg(k, casdiags % prevp(i,j,k) , 'prevp', i_dgtime) 
!             END DO
!         END IF
!      END IF

!      IF (casdiags % l_psedr) THEN
!         IF (pswitch%l_psedr) THEN
!            DO k = ks, ke
!               kc = k - ks + 1
!               call save_dg(k, procs(rain_params%i_1m,i_psedr%id)%column_data(kc) , 'psedr', i_dgtime) 
!            END DO
!         END IF
!      END IF

!      IF (casdiags % l_psedl) THEN
!         IF (pswitch%l_psedl) THEN
!            DO k = ks, ke
!               kc = k - ks + 1
!               call save_dg(k, procs(cloud_params%i_1m,i_psedl%id)%column_data(kc) , 'psedl', i_dgtime) 
!            END DO
!         END IF
!      END IF
     
!      IF (casdiags % l_pracr) THEN
!         IF (pswitch%l_pracr) THEN
!            DO k = ks, ke
!               kc = k - ks + 1
!               call save_dg(k, procs(rain_params%i_1m,i_pracr%id)%column_data(kc), 'pracr', i_dgtime) 
!            END DO
!         END IF
!      END IF
! #endif
    IF (casdiags % l_pgacw) THEN
      IF (pswitch%l_pgacw) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pgacw(i,j,k) = procs(graupel_params%i_1m, i_gacw%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pgacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgacs) THEN
      IF (pswitch%l_pgacs) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pgacs(i,j,k) = procs(graupel_params%i_1m, i_gacs%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pgacs(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgmlt) THEN
      IF (pswitch%l_pgmlt) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pgmlt(i,j,k) = procs(rain_params%i_1m, i_gmlt%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pgmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgsub) THEN
      IF (pswitch%l_pgsub) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pgsub(i,j,k) = -1.0 * procs(graupel_params%i_1m,i_gsub%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pgsub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedi) THEN
      IF (pswitch%l_psedi) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psedi(i,j,k) = procs(ice_params%i_1m, i_psedi%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psedi(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedi) THEN
      IF ((pswitch%l_psedi) .and. (snow_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % nsedi(i,j,k) = procs(ice_params%i_2m,i_psedi%id)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedi(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pseds) THEN
      IF (pswitch%l_pseds) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pseds(i,j,k) = procs(snow_params%i_1m, i_pseds%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pseds(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nseds) THEN
      IF ((pswitch%l_pseds) .and. (snow_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % nseds(i,j,k) = procs(snow_params%i_2m, i_pseds%id)%column_data(kc)
        END DO
      ELSE
        casdiags % nseds(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedr) THEN
      IF (pswitch%l_psedr) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psedr(i,j,k) = procs(rain_params%i_1m,i_psedr%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psedr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedg) THEN
      IF (pswitch%l_psedg) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psedg(i,j,k) = procs(graupel_params%i_1m,i_psedg%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psedg(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedg) THEN
      IF ((pswitch%l_psedg) .and. (graupel_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % nsedg(i,j,k) = procs(graupel_params%i_2m,i_psedg%id)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedg(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedl) THEN
      IF (pswitch%l_psedl) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psedl(i,j,k) = procs(cloud_params%i_1m,i_psedl%id)%column_data(kc)
        END DO
      ELSE
        casdiags % psedl(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pcond) THEN
      IF (pswitch%l_pcond) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pcond(i,j,k) = procs(cloud_params%i_1m,i_cond%id)%column_data(kc)
        END DO
      ELSE
        casdiags % pcond(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_phomr) THEN
      IF (pswitch%l_phomr) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % phomr(i,j,k) = procs(graupel_params%i_1m,i_homr%id)%column_data(kc)
        END DO
      ELSE
        casdiags % phomr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF


    IF (casdiags % l_nihal) THEN
      IF ((pswitch%l_pihal) .and. (ice_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % nihal(i,j,k) = procs(ice_params%i_2m,i_ihal%id)%column_data(kc)
        END DO
      ELSE
        casdiags % nihal(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nhomr) THEN
      IF ((pswitch%l_phomr) .and. (ice_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % nhomr(i,j,k) = procs(graupel_params%i_2m,i_homr%id)%column_data(kc)
        END DO
      ELSE
        casdiags % nhomr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    else !ncall > 0

    IF (casdiags % l_psedi) THEN
      IF (pswitch%l_psedi) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psedi(i,j,k) = casdiags % psedi(i,j,k)+procs(ice_params%i_1m,i_psedi%id)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nsedi) THEN
      IF ((pswitch%l_psedi) .and. (snow_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % nsedi(i,j,k) = casdiags % nsedi(i,j,k)+procs(ice_params%i_2m,i_psedi%id)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_pseds) THEN
      IF (pswitch%l_pseds) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % pseds(i,j,k) = casdiags % pseds(i,j,k)+procs(snow_params%i_1m,i_pseds%id)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nseds) THEN
      IF ((pswitch%l_pseds) .and. (snow_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % nseds(i,j,k) = casdiags % nseds(i,j,k)+procs(snow_params%i_2m,i_pseds%id)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_psedr) THEN
      IF (pswitch%l_psedr) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psedr(i,j,k) = casdiags % psedr(i,j,k)+procs(rain_params%i_1m,i_psedr%id)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_psedg) THEN
      IF (pswitch%l_psedg) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psedg(i,j,k) = casdiags % psedg(i,j,k)+procs(graupel_params%i_1m,i_psedg%id)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nsedg) THEN
      IF ((pswitch%l_psedg) .and. (graupel_params%l_2m)) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % nsedg(i,j,k) = casdiags % nsedg(i,j,k)+procs(graupel_params%i_2m,i_psedg%id)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_psedl) THEN
      IF (pswitch%l_psedl) THEN
        DO k = ks, ke
          kc = k - ks + 1
          casdiags % psedl(i,j,k) = casdiags % psedl(i,j,k)+procs(cloud_params%i_1m,i_psedl%id)%column_data(kc)
        END DO
      END IF
    END IF

    end if ! ncall

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine gather_process_diagnostics
end module micro_main
