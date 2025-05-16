module micro_main
  use variable_precision, only: wp
  use mphys_parameters, only: nz, nq, rain_params, cloud_params, ice_params, &
       snow_params, graupel_params, nspecies, ZERO_REAL_WP, a_s, b_s, &
       nxy_inner
  use process_routines, only: process_rate, zero_procs, allocate_procs, deallocate_procs, i_cond, i_praut, &
       i_pracw, i_pracr, i_prevp, i_psedr, i_psedl, i_aact, i_aaut, i_aacw, i_aevp, i_asedr, i_asedl, i_arevp, &
       i_tidy2, i_atidy2, i_inuc, i_idep, i_dnuc, i_dsub, i_saut, i_iacw, i_sacw, i_pseds, &
       i_sdep, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs, i_gdep, i_psedg, i_sagg, &
       i_gshd, i_ihal, i_smlt, i_gmlt, i_psedi, i_homr, i_homc, i_imlt, i_isub, i_ssub, i_gsub, i_sbrk, i_dssub, &
       i_dgsub, i_dsedi, i_dseds, i_dsedg, i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr, &
       i_dgacr, i_draci, i_dhomr, i_dhomc, i_idps, i_iics
  use sum_process, only: sum_procs, sum_aprocs, tend_temp, aerosol_tend_temp
  use aerosol_routines, only: examine_aerosol, aerosol_phys, aerosol_chem, aerosol_active, allocate_aerosol, &
       deallocate_aerosol
  use mphys_switches, only: hydro_complexity, aero_complexity, i_qv, i_ql, i_nl, i_qr, i_nr, i_m3r, i_th, i_qi, &
       i_qs, i_qg, i_ni, i_ns, i_ng, i_m3s, i_m3g, i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am4, i_am5, i_am6, &
       i_an6, i_am7, i_am8 , i_am9, i_am10, i_an10, i_an11, i_an12, i_ak1, i_ak2, i_ak3, &
       aerosol_option, l_warm, &
       l_sed, l_idep, aero_index, nq_l, nq_r, nq_i, nq_s, nq_g, &
       l_sg, l_g, l_process, max_sed_length, max_step_length, l_harrington, l_passive, ntotala, ntotalq, &
       l_onlycollect, pswitch, l_isub, l_pos1, l_pos2, l_pos3, l_pos4, l_no_pgacs_in_sumprocs, &
       l_pos5, l_pos6, i_hstart, l_tidy_negonly,  &
       iopt_act, iopt_shipway_act, l_prf_cfrac, l_kfsm, l_gamma_online, l_subseds_maxv, &
       i_cfl, i_cfr, i_cfi, i_cfs, i_cfg, l_reisner_graupel_embryo
! use mphys_switches, only: l_rain,
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
  use ice_multiplication, only: hallet_mossop, droplet_shattering, ice_collision
  use graupel_wetgrowth, only: wetgrowth
  use graupel_embryo, only: graupel_embryos
  use ice_melting, only: melting
  use homogeneous, only: ihom_rain, ihom_droplets
  use adjust_deposition, only: adjust_dep
  use mphys_constants, only: fixed_aerosol_sigma, fixed_aerosol_density
!AJM  removing line below causes model to fail
   use lookup, only: get_slope_generic

  use mphys_tidy, only: initialise_mphystidy, finalise_mphystidy, qtidy, ensure_positive, &
       ensure_saturated, tidy_qin, tidy_ain, ensure_positive_aerosol
  use preconditioning, only: precondition, preconditioner

  use casim_reflec_mod, only: casim_reflec, setup_reflec_constants
  ! For initialization of Shipway (2015) activation scheme
  use shipway_lookup, only: generate_tables

  use generic_diagnostic_variables, only: casdiags

  use casim_runtime, only: casim_time
  use casim_parent_mod, only: casim_parent, parent_monc

  use casim_stph, only: l_rp2_casim, snow_a_x_rp, ice_a_x_rp

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

  integer :: i_start, i_end ! upper and lower i levels which are to be used
  integer :: j_start, j_end ! upper and lower j levels
  integer :: k_start, k_end ! upper and lower k levels

!$OMP THREADPRIVATE(i_start,i_end,j_start,j_end,k_start,k_end)
!  integer :: nxny
  real(wp), allocatable, save :: precip(:,:) ! diagnostic for surface precip rate

!--Add one more dimention for arrays to be used in ixy_inner loop--!
  real(wp), allocatable :: dqfields(:,:,:), qfields(:,:,:), tend(:,:,:)
  real(wp), allocatable :: daerofields(:,:,:), aerofields(:,:,:), aerosol_tend(:,:,:)
  real(wp), allocatable :: cffields(:,:,:) !cloudfraction fields

!$OMP THREADPRIVATE(precip, dqfields, qfields, cffields, tend,                   &
!$OMP               daerofields, aerofields, aerosol_tend)

  type(process_rate), allocatable :: procs(:,:,:)
  type(process_rate), allocatable :: aerosol_procs(:,:,:)

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

  real(wp), allocatable :: qfields_in(:,:,:)
  real(wp), allocatable :: qfields_mod(:,:,:)
  real(wp), allocatable :: aerofields_in(:,:,:)
  real(wp), allocatable :: aerofields_mod(:,:,:)

!$OMP THREADPRIVATE(qfields_in, qfields_mod, aerofields_in, aerofields_mod)

  logical, allocatable :: l_Tcold(:,:) ! temperature below freezing, i.e. .not. l_Twarm
  logical, allocatable :: l_sigevap(:,:) ! logical to determine significant evaporation

!$OMP THREADPRIVATE(l_Tcold, l_sigevap)

!Arrays needed by inner loop

  integer, allocatable :: columns(:)

!$OMP THREADPRIVATE(columns)

  real :: DTPUD ! Time step for puddle diagnostic

  public initialise_micromain, finalise_micromain, shipway_microphysics, DTPUD

!PRF
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

    i_start=is_in
    i_end=ie_in
    j_start=js_in
    j_end=je_in
    k_start=ks_in
    k_end=ke_in

    !nxy_inner = (i_end-i_start+1)*(j_end-j_start+1) ! test

    allocate(rhcrit_1d(kl:ku))
    ! Set RHCrit to 1.0 as default; parent model can then overwrite
    ! this if needed
    rhcrit_1d(:) = 1.0

    l_tendency_loc = l_tendency

    nq=sum(hydro_complexity%nmoments)+2 ! also includes vapour and theta
    nz=k_end-k_start+1
    nprocs = hydro_complexity%nprocesses


    allocate(precondition(nz, nxy_inner))
    precondition=.true. ! Assume all points need to be considered
    allocate(l_Tcold(nz, nxy_inner))
    l_Tcold =.false. ! Assumes no cold points, this is set at the beginning of microphysics_common
    allocate(l_sigevap(nz, nxy_inner))
    l_sigevap = .false.
    allocate(qfields(nz, nq, nxy_inner))
    allocate(dqfields(nz, nq, nxy_inner))
    qfields=ZERO_REAL_WP
    dqfields=ZERO_REAL_WP
    !allocate(procs(nz, nprocs))
    allocate(procs(ntotalq, nprocs, nxy_inner))
    allocate(tend(nz, nq, nxy_inner))
    allocate(tend_temp(nz,nq))
    allocate(cffields(nz,5, nxy_inner)) !5 'cloud' fractions
    cffields=ZERO_REAL_WP

    ! Allocate aerosol storage
    if (aerosol_option > 0) then
      naero=ntotala
      naeroprocs=aero_complexity%nprocesses
      allocate(aerofields(nz, naero, nxy_inner))
      allocate(daerofields(nz, naero, nxy_inner))
      aerofields=ZERO_REAL_WP
      daerofields=ZERO_REAL_WP
      ! allocate(aerosol_procs(nz, naeroprocs))
      allocate(aerosol_procs(ntotala, naeroprocs, nxy_inner))
      allocate(aerosol_tend(nz, naero, nxy_inner))
      allocate(aerosol_tend_temp(nz, naero))
    else
      ! Dummy arrays required
      allocate(aerofields(1,1,nxy_inner))
      allocate(daerofields(1,1,nxy_inner))
      allocate(aerosol_procs(1,1,nxy_inner))
      allocate(aerosol_tend(1,1,nxy_inner))
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
!    do ixy = 1, nxy_inner ! push ixy loop into allocate_procs
    call allocate_procs(nxy_inner, procs, nz, nprocs, ntotalq)
    if (l_process) call allocate_procs(nxy_inner, aerosol_procs, nz, naeroprocs, ntotala)
!    end do


    ! allocate diagnostics
    allocate(precip(il:iu,jl:ju))

    call initialise_passive_fields(k_start, k_end)

    allocate(qfields_in(nz, nq, nxy_inner))
    allocate(qfields_mod(nz, nq, nxy_inner))
    if (aerosol_option > 0) then
      allocate(aerofields_in(nz, naero, nxy_inner))
      allocate(aerofields_mod(nz, naero, nxy_inner))
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

! Array needed for inner loop
    allocate(columns(nxy_inner))


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

!    nxy_inner = 2

    ! deallocate diagnostics
    deallocate(precip)

    ! deallocate process rates
    call deallocate_procs(nxy_inner, procs)

    deallocate(procs)
    deallocate(qfields)
    deallocate(tend)
    deallocate(tend_temp)
    deallocate(precondition)
    deallocate(cffields)

    ! aerosol fields
    if (l_process) call deallocate_procs(nxy_inner, aerosol_procs)

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

! Array needed for inner loop
    deallocate(columns)


    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine finalise_micromain

  subroutine shipway_microphysics(il, iu, jl, ju, kl, ku, dt,               &
       qv, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13,          &
       theta, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13,       &
       a14, a15, a16, a17, a18, a19, a20,                                   &
       exner, pressure, rho, w, tke, dz,                                    &
       cfliq, cfice, cfsnow, cfrain, cfgr,    &
       dqv, dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11, dq12,  &
       dq13, dth, da1, da2, da3, da4, da5, da6, da7, da8, da9, da10, da11,  &
       da12, da13, da14, da15, da16, da17,                                  &
       is_in, ie_in, js_in, je_in)

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

    ! Local variables

    integer :: k, i, j, ixy
    integer :: nxy_all, ixy_outer, ixy_inner, nxy_inner_loop

    real(wp) :: precip_l(nxy_inner)
    real(wp) :: precip_r(nxy_inner)
    real(wp) :: precip_i(nxy_inner)
    real(wp) :: precip_s(nxy_inner)
    real(wp) :: precip_g(nxy_inner)

    real(wp) :: precip_l1d(nz, nxy_inner)
    real(wp) :: precip_r1d(nz, nxy_inner)
    real(wp) :: precip_i1d(nz, nxy_inner)
    real(wp) :: precip_s1d(nz, nxy_inner)
    real(wp) :: precip_so1d(nz, nxy_inner)
    real(wp) :: precip_g1d(nz, nxy_inner)

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



    nxy_all = (je_in-js_in+1)*(ie_in-is_in+1)
    do ixy_outer=1, int(ceiling(dble(nxy_all)/nxy_inner))
       nxy_inner_loop = min(nxy_inner, &
                        nxy_all-(ixy_outer-1)*nxy_inner)
       do ixy_inner=1, nxy_inner_loop
          ixy = (ixy_outer-1)*nxy_inner + ixy_inner
          j = modulo(ixy-1,(je_in-js_in+1))+js_in
          i = (ixy-1)/(je_in-js_in+1)+is_in

          precip_l(ixy_inner) = 0.0
          precip_r(ixy_inner) = 0.0
          precip_i(ixy_inner) = 0.0
          precip_s(ixy_inner) = 0.0
          precip_g(ixy_inner) = 0.0

          precip_l1d(:, ixy_inner)  = 0.0
          precip_r1d(:, ixy_inner)  = 0.0
          precip_i1d(:, ixy_inner)  = 0.0
          precip_s1d(:, ixy_inner)  = 0.0
          precip_g1d(:, ixy_inner)  = 0.0
          precip_so1d(:, ixy_inner) = 0.0
          tend(:,:,ixy_inner)=ZERO_REAL_WP
          call zero_procs(procs(:,:,ixy_inner))
          aerosol_tend(:,:,ixy_inner)=ZERO_REAL_WP
          if (l_process) then
             call zero_procs(aerosol_procs(:,:,ixy_inner))
          end if
          !set cloud fraction fields
          cffields(:,i_cfl,ixy_inner)=cfliq(k_start:k_end,i,j)
          cffields(:,i_cfr,ixy_inner)=cfrain(k_start:k_end,i,j)
          cffields(:,i_cfi,ixy_inner)=cfice(k_start:k_end,i,j)
          cffields(:,i_cfs,ixy_inner)=cfsnow(k_start:k_end,i,j)
          cffields(:,i_cfg,ixy_inner)=cfgr(k_start:k_end,i,j)
          ! Set the qfields
          qfields(:, i_qv, ixy_inner)=qv(k_start:k_end,i,j)
          qfields(:, i_th, ixy_inner)=theta(k_start:k_end,i,j)
          if (nq_l > 0) qfields(:,i_ql,ixy_inner)=q1(k_start:k_end,i,j)
          if (nq_r > 0) qfields(:,i_qr,ixy_inner)=q2(k_start:k_end,i,j)
          if (nq_l > 1) qfields(:,i_nl,ixy_inner)=q3(k_start:k_end,i,j)
          if (nq_r > 1) qfields(:,i_nr,ixy_inner)=q4(k_start:k_end,i,j)
          if (nq_r > 2) qfields(:,i_m3r,ixy_inner)=q5(k_start:k_end,i,j)
          if (nq_i > 0) qfields(:,i_qi,ixy_inner)=q6(k_start:k_end,i,j)
          if (nq_s > 0) qfields(:,i_qs,ixy_inner)=q7(k_start:k_end,i,j)
          if (nq_g > 0) qfields(:,i_qg,ixy_inner)=q8(k_start:k_end,i,j)
          if (nq_i > 1) qfields(:,i_ni,ixy_inner)=q9(k_start:k_end,i,j)
          if (nq_s > 1) qfields(:,i_ns,ixy_inner)=q10(k_start:k_end,i,j)
          if (nq_g > 1) qfields(:,i_ng,ixy_inner)=q11(k_start:k_end,i,j)
          if (nq_s > 2) qfields(:,i_m3s,ixy_inner)=q12(k_start:k_end,i,j)
          if (nq_g > 2) qfields(:,i_m3g,ixy_inner)=q13(k_start:k_end,i,j)
          dqfields(:, i_qv, ixy_inner)=dqv(k_start:k_end,i,j)
          dqfields(:, i_th, ixy_inner)=dth(k_start:k_end,i,j)
          if (nq_l > 0) dqfields(:,i_ql,ixy_inner)=dq1(k_start:k_end,i,j)
          if (nq_r > 0) dqfields(:,i_qr,ixy_inner)=dq2(k_start:k_end,i,j)
          if (nq_l > 1) dqfields(:,i_nl,ixy_inner)=dq3(k_start:k_end,i,j)
          if (nq_r > 1) dqfields(:,i_nr,ixy_inner)=dq4(k_start:k_end,i,j)
          if (nq_r > 2) dqfields(:,i_m3r,ixy_inner)=dq5(k_start:k_end,i,j)
          if (nq_i > 0) dqfields(:,i_qi,ixy_inner)=dq6(k_start:k_end,i,j)
          if (nq_s > 0) dqfields(:,i_qs,ixy_inner)=dq7(k_start:k_end,i,j)
          if (nq_g > 0) dqfields(:,i_qg,ixy_inner)=dq8(k_start:k_end,i,j)
          if (nq_i > 1) dqfields(:,i_ni,ixy_inner)=dq9(k_start:k_end,i,j)
          if (nq_s > 1) dqfields(:,i_ns,ixy_inner)=dq10(k_start:k_end,i,j)
          if (nq_g > 1) dqfields(:,i_ng,ixy_inner)=dq11(k_start:k_end,i,j)
          if (nq_s > 2) dqfields(:,i_m3s,ixy_inner)=dq12(k_start:k_end,i,j)
          if (nq_g > 2) dqfields(:,i_m3g,ixy_inner)=dq13(k_start:k_end,i,j)
          if (aerosol_option > 0) then
             if (i_am1 >0) aerofields(:, i_am1, ixy_inner)=a1(k_start:k_end,i,j)
             if (i_an1 >0) aerofields(:, i_an1, ixy_inner)=a2(k_start:k_end,i,j)
             if (i_am2 >0) aerofields(:, i_am2, ixy_inner)=a3(k_start:k_end,i,j)
             if (i_an2 >0) aerofields(:, i_an2, ixy_inner)=a4(k_start:k_end,i,j)
             if (i_am3 >0) aerofields(:, i_am3, ixy_inner)=a5(k_start:k_end,i,j)
             if (i_an3 >0) aerofields(:, i_an3, ixy_inner)=a6(k_start:k_end,i,j)
             if (i_am4 >0) aerofields(:, i_am4, ixy_inner)=a7(k_start:k_end,i,j)
             if (i_am5 >0) aerofields(:, i_am5, ixy_inner)=a8(k_start:k_end,i,j)
             if (i_am6 >0) aerofields(:, i_am6, ixy_inner)=a9(k_start:k_end,i,j)
             if (i_an6 >0) aerofields(:, i_an6, ixy_inner)=a10(k_start:k_end,i,j)
             if (i_am7 >0) aerofields(:, i_am7, ixy_inner)=a11(k_start:k_end,i,j)
             if (i_am8 >0) aerofields(:, i_am8, ixy_inner)=a12(k_start:k_end,i,j)
             if (i_am9 >0) aerofields(:, i_am9, ixy_inner)=a13(k_start:k_end,i,j)
             if (i_am10 >0) aerofields(:, i_am10, ixy_inner)=a14(k_start:k_end,i,j)
             if (i_an10 >0) aerofields(:, i_an10, ixy_inner)=a15(k_start:k_end,i,j)
             if (i_an11 >0) aerofields(:, i_an11, ixy_inner)=a16(k_start:k_end,i,j)
             if (i_an12 >0) aerofields(:, i_an12, ixy_inner)=a17(k_start:k_end,i,j)
             if (i_ak1 >0) aerofields(:, i_ak1, ixy_inner)=a18(k_start:k_end,i,j)
             if (i_ak2 >0) aerofields(:, i_ak2, ixy_inner)=a19(k_start:k_end,i,j)
             if (i_ak3 >0) aerofields(:, i_ak3, ixy_inner)=a20(k_start:k_end,i,j)
             if (i_am1 >0) daerofields(:, i_am1, ixy_inner)=da1(k_start:k_end,i,j)
             if (i_an1 >0) daerofields(:, i_an1, ixy_inner)=da2(k_start:k_end,i,j)
             if (i_am2 >0) daerofields(:, i_am2, ixy_inner)=da3(k_start:k_end,i,j)
             if (i_an2 >0) daerofields(:, i_an2, ixy_inner)=da4(k_start:k_end,i,j)
             if (i_am3 >0) daerofields(:, i_am3, ixy_inner)=da5(k_start:k_end,i,j)
             if (i_an3 >0) daerofields(:, i_an3, ixy_inner)=da6(k_start:k_end,i,j)
             if (i_am4 >0) daerofields(:, i_am4, ixy_inner)=da7(k_start:k_end,i,j)
             if (i_am5 >0) daerofields(:, i_am5, ixy_inner)=da8(k_start:k_end,i,j)
             if (i_am6 >0) daerofields(:, i_am6, ixy_inner)=da9(k_start:k_end,i,j)
             if (i_an6 >0) daerofields(:, i_an6, ixy_inner)=da10(k_start:k_end,i,j)
             if (i_am7 >0) daerofields(:, i_am7, ixy_inner)=da11(k_start:k_end,i,j)
             if (i_am8 >0) daerofields(:, i_am8, ixy_inner)=da12(k_start:k_end,i,j)
             if (i_am9 >0) daerofields(:, i_am9, ixy_inner)=da13(k_start:k_end,i,j)
             if (i_am10 >0) daerofields(:, i_am10, ixy_inner)=da14(k_start:k_end,i,j)
             if (i_an10 >0) daerofields(:, i_an10, ixy_inner)=da15(k_start:k_end,i,j)
             if (i_an11 >0) daerofields(:, i_an11, ixy_inner)=da16(k_start:k_end,i,j)
             if (i_an12 >0) daerofields(:, i_an12, ixy_inner)=da17(k_start:k_end,i,j)
          end if
       end do ! ixy_inner

       !inner loop pushed into set_passive_fields
       !--------------------------------------------------
       ! set fields which will not be modified
       !--------------------------------------------------
       call set_passive_fields(nxy_inner_loop, ixy_outer, is_in, js_in, je_in, &
            dt, rho,    &
            pressure, exner,    &
            dz,                            &
            w, tke, qfields)

       !--------------------------------------------------
       ! Do the business...
       !--------------------------------------------------
       call microphysics_common(&
            nxy_inner_loop, &
            ! Inner loop size
            ixy_outer, is_in, js_in, je_in, &
            ! To calculate i, j
            dt, &
            !i , j,
            ! To be calculated so no need anymore
            qfields, cffields, dqfields, tend, procs &
            !, precip(i,j)
            , precip &
            , precip_l, precip_r, precip_i, precip_s, precip_g       &
            , precip_r1d, precip_s1d, precip_so1d, precip_g1d                     &
            , aerophys, aerochem, aeroact                                         &
            , dustphys, dustchem, dustact                                         &
            , aeroice, dustliq                                                    &
            , aerofields, daerofields, aerosol_tend, aerosol_procs                &
            , rhcrit_1d)

       !end do ! ixy_inner
       !--------------------------------------------------
       ! Relate back tendencies
       ! Check indices in mphys_switches that the appropriate
       ! fields are being passed back to mphys_casim
       !--------------------------------------------------

       do ixy_inner=1, nxy_inner_loop
           ixy = (ixy_outer-1)*nxy_inner + ixy_inner
           j = modulo(ixy-1,(je_in-js_in+1))+js_in
           i = (ixy-1)/(je_in-js_in+1)+is_in

           dqv(k_start:k_end,i,j)=tend(:,i_qv,ixy_inner)
           dth(k_start:k_end,i,j)=tend(:,i_th,ixy_inner)
           dq1(k_start:k_end,i,j)=tend(:,i_ql,ixy_inner)
           dq2(k_start:k_end,i,j)=tend(:,i_qr,ixy_inner)
           if (cloud_params%l_2m) dq3(k_start:k_end,i,j)=tend(:,i_nl,ixy_inner)
           if (rain_params%l_2m) dq4(k_start:k_end,i,j)=tend(:,i_nr,ixy_inner)
           if (rain_params%l_3m) dq5(k_start:k_end,i,j)=tend(:,i_m3r,ixy_inner)

           if (.not. l_warm) then
              if (ice_params%l_1m) dq6(k_start:k_end,i,j)=tend(:,i_qi,ixy_inner)
              if (snow_params%l_1m) dq7(k_start:k_end,i,j)=tend(:,i_qs,ixy_inner)
              if (graupel_params%l_1m) dq8(k_start:k_end,i,j)=tend(:,i_qg,ixy_inner)
              if (ice_params%l_2m) dq9(k_start:k_end,i,j)=tend(:,i_ni,ixy_inner)
              if (snow_params%l_2m) dq10(k_start:k_end,i,j)=tend(:,i_ns,ixy_inner)
              if (graupel_params%l_2m) dq11(k_start:k_end,i,j)=tend(:,i_ng,ixy_inner)
              if (snow_params%l_3m) dq12(k_start:k_end,i,j)=tend(:,i_m3s,ixy_inner)
              if (graupel_params%l_3m) dq13(k_start:k_end,i,j)=tend(:,i_m3g,ixy_inner)
           end if

           if (l_process) then
              if (i_am1 >0) da1(k_start:k_end,i,j)=aerosol_tend(:,i_am1,ixy_inner)
              if (i_an1 >0) da2(k_start:k_end,i,j)=aerosol_tend(:,i_an1,ixy_inner)
              if (i_am2 >0) da3(k_start:k_end,i,j)=aerosol_tend(:,i_am2,ixy_inner)
              if (i_an2 >0) da4(k_start:k_end,i,j)=aerosol_tend(:,i_an2,ixy_inner)
              if (i_am3 >0) da5(k_start:k_end,i,j)=aerosol_tend(:,i_am3,ixy_inner)
              if (i_an3 >0) da6(k_start:k_end,i,j)=aerosol_tend(:,i_an3,ixy_inner)
              if (i_am4 >0) da7(k_start:k_end,i,j)=aerosol_tend(:,i_am4,ixy_inner)
              if (i_am5 >0) da8(k_start:k_end,i,j)=aerosol_tend(:,i_am5,ixy_inner)
              if (i_am6 >0) da9(k_start:k_end,i,j)=aerosol_tend(:,i_am6,ixy_inner)
              if (i_an6 >0) da10(k_start:k_end,i,j)=aerosol_tend(:,i_an6,ixy_inner)
              if (i_am7 >0) da11(k_start:k_end,i,j)=aerosol_tend(:,i_am7,ixy_inner)
              if (i_am8 >0) da12(k_start:k_end,i,j)=aerosol_tend(:,i_am8,ixy_inner)
              if (i_am9 >0) da13(k_start:k_end,i,j)=aerosol_tend(:,i_am9,ixy_inner)
              if (i_am10 >0) da14(k_start:k_end,i,j)=aerosol_tend(:,i_am10,ixy_inner)
              if (i_an10 >0) da15(k_start:k_end,i,j)=aerosol_tend(:,i_an10,ixy_inner)
              if (i_an11 >0) da16(k_start:k_end,i,j)=aerosol_tend(:,i_an11,ixy_inner)
              if (i_an12 >0) da17(k_start:k_end,i,j)=aerosol_tend(:,i_an12,ixy_inner)
           else
              da1(k_start:k_end,i,j)=0.0
              da2(k_start:k_end,i,j)=0.0
              da3(k_start:k_end,i,j)=0.0
              da4(k_start:k_end,i,j)=0.0
              da5(k_start:k_end,i,j)=0.0
              da6(k_start:k_end,i,j)=0.0
              da7(k_start:k_end,i,j)=0.0
              da9(k_start:k_end,i,j)=0.0
              da10(k_start:k_end,i,j)=0.0
              da11(k_start:k_end,i,j)=0.0
              da12(k_start:k_end,i,j)=0.0
              da13(k_start:k_end,i,j)=0.0
              da14(k_start:k_end,i,j)=0.0
              da15(k_start:k_end,i,j)=0.0
              da16(k_start:k_end,i,j)=0.0
              da17(k_start:k_end,i,j)=0.0
           end if

           if ( l_warm ) then
              if ( casdiags % l_surface_cloud ) casdiags % SurfaceCloudR(i,j) = precip_l(ixy_inner)
              if ( casdiags % l_surface_rain ) casdiags % SurfaceRainR(i,j)  = precip_r(ixy_inner)
              if ( casdiags % l_surface_snow ) casdiags % SurfaceSnowR(i,j)  = 0.0
              if ( casdiags % l_surface_graup) casdiags % SurfaceGraupR(i,j) = 0.0
              if ( casdiags % l_rainfall_3d ) casdiags % rainfall_3d(i,j,k_start:k_end)  = precip_r1d(:,ixy_inner)
              if ( casdiags % l_snowfall_3d ) casdiags % snowfall_3d(i,j,k_start:k_end)  = 0.0
              if ( casdiags % l_snowonly_3d ) casdiags % snowonly_3d(i,j,k_start:k_end)  = 0.0
              if ( casdiags % l_graupfall_3d) casdiags % graupfall_3d(i,j,k_start:k_end) = 0.0
           else ! l_warm

              if ( casdiags % l_surface_rain ) casdiags % SurfaceRainR(i,j)  = precip_r(ixy_inner)
              if ( casdiags % l_surface_snow ) casdiags % SurfaceSnowR(i,j)  = precip_s(ixy_inner)
              if ( casdiags % l_surface_graup) casdiags % SurfaceGraupR(i,j) = precip_g(ixy_inner)
              if ( casdiags % l_rainfall_3d ) casdiags % rainfall_3d(i,j,k_start:k_end)  = precip_r1d(:,ixy_inner)
              if ( casdiags % l_snowfall_3d ) casdiags % snowfall_3d(i,j,k_start:k_end)  = precip_s1d(:,ixy_inner)
              if ( casdiags % l_snowonly_3d ) casdiags % snowonly_3d(i,j,k_start:k_end)  = precip_so1d(:,ixy_inner)
              if ( casdiags % l_graupfall_3d) casdiags % graupfall_3d(i,j,k_start:k_end) = precip_g1d(:,ixy_inner)
           end if ! l_warm

           if ( casdiags % l_radar ) then

              call tidy_qin(ixy_inner, qfields(:,:,ixy_inner))  !check this is conserving. If i do this here do we need a tidy_ain?
              call casim_reflec(ixy_inner, nz, nq, rho(k_start:k_end,i,j), qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),  &
                               dbz_tot_c, dbz_g_c, dbz_i_c,                      &
                               dbz_s_c,   dbz_l_c, dbz_r_c  )

              casdiags % dbz_tot(i,j, k_start:k_end) = dbz_tot_c(:)
              casdiags % dbz_g(i,j,   k_start:k_end) = dbz_g_c(:)
              casdiags % dbz_s(i,j,   k_start:k_end) = dbz_s_c(:)
              casdiags % dbz_i(i,j,   k_start:k_end) = dbz_i_c(:)
              casdiags % dbz_l(i,j,   k_start:k_end) = dbz_l_c(:)
              casdiags % dbz_r(i,j,   k_start:k_end) = dbz_r_c(:)

           end if ! casdiags % l_radar

           if ( casdiags % l_tendency_dg ) then
              DO k = k_start, k_end
                 kc = k - k_start + 1
                 casdiags % dth_cond_evap(i,j,k) = procs(cloud_params%i_1m,i_cond%id,ixy_inner)%column_data(kc) * &
                                                   Lv/cp * rexner(kc,ixy_inner)
                 casdiags % dqv_cond_evap(i,j,k) = -(procs(cloud_params%i_1m,i_cond%id,ixy_inner)%column_data(kc))
                 casdiags % dth_total(i,j,k) = tend(kc,i_th,ixy_inner)
                 casdiags % dqv_total(i,j,k) = tend(kc,i_qv,ixy_inner)
                 casdiags % dqc(i,j,k) = tend(kc,i_ql,ixy_inner)
                 casdiags % dqr(i,j,k) = tend(kc,i_qr,ixy_inner)

                 if (.not. l_warm) then
                    casdiags % dqi(i,j,k) = tend(kc,i_qi,ixy_inner)
                    casdiags % dqs(i,j,k) = tend(kc,i_qs,ixy_inner)
                    casdiags % dqg(i,j,k) = tend(kc,i_qg,ixy_inner)
                 endif
              enddo
           endif

           if ( casdiags % l_lwp ) then
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_ql,ixy_inner) )
              enddo
              casdiags % lwp(i,j)=waterpath
           endif
           if ( casdiags % l_rwp ) then
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qr,ixy_inner) )
              enddo
              casdiags % rwp(i,j)=waterpath
           endif
           if ( casdiags % l_iwp ) then
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qi,ixy_inner) )
              enddo
              casdiags % iwp(i,j)=waterpath
           endif
           if ( casdiags % l_swp ) then
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qs,ixy_inner) )
              enddo
              casdiags % swp(i,j)=waterpath
           endif
           if ( casdiags % l_gwp ) then
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qg,ixy_inner) )
              enddo
              casdiags % gwp(i,j)=waterpath
           endif
       !end do ! i
    !end do   ! j
        end do ! ixy_inner
     end do ! ixy_outer


! #if DEF_MODEL==MODEL_KiD
!     call save_dg(sum(casdiags % SurfaceRainR(:, :))/nxny, 'precip', i_dgtime)
!     call save_dg(sum(casdiags % SurfaceRainR(:, :))/nxny*3600.0, 'surface_precip_mmhr', i_dgtime)
! #endif

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine shipway_microphysics

  subroutine microphysics_common(&
         nxy_inner_loop &
       , ixy_outer, is_in, js_in, je_in &
       , dt &
       !ix, jy,
       , qfields, cffields, dqfields, tend &
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

    integer, intent(in) :: nxy_inner_loop
    integer, intent(in) :: ixy_outer
    integer, intent(in) :: is_in, js_in, je_in

    real(wp), intent(in) :: dt  ! timestep from parent model
    ! integer, intent(in) :: ix, jy

    real(wp), intent(in) :: rhcrit_1d(:)
    real(wp), intent(inout) :: qfields(:,:,:), dqfields(:,:,:), tend(:,:,:)
    real(wp), intent(in) :: cffields(:,:,:)

    type(process_rate), intent(inout) :: procs(:,:,:)
    ! real(wp), intent(out) :: precip
    real(wp), intent(out) :: precip(:,:)
    real(wp), intent(INOUT) :: precip_l(:)
    real(wp), INTENT(INOUT) :: precip_r(:)
    real(wp), intent(INOUT) :: precip_i(:)
    real(wp), INTENT(INOUT) :: precip_s(:)
    real(wp), intent(inout) :: precip_g(:)

    real(wp), INTENT(INOUT) :: precip_r1d(:,:)
    real(wp), INTENT(INOUT) :: precip_s1d(:,:)
    real(wp), INTENT(INOUT) :: precip_so1d(:,:)
    real(wp), INTENT(INOUT) :: precip_g1d(:,:)

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

    real(wp), intent(inout) :: aerofields(:,:,:), daerofields(:,:,:), aerosol_tend(:,:,:)
    type(process_rate), intent(inout), optional :: aerosol_procs(:,:,:)

    real(wp) :: step_length

    real(wp) :: sed_length, sed_length_cloud, sed_length_rain, sed_length_ice, sed_length_snow &
      , sed_length_graupel

    !--not input anymore--
    integer :: ix, jy, ixy_inner, ixy, i_column

    integer :: n, k, nsed, iq

    logical :: l_Twarm   ! temperature above freezing

    integer, parameter :: level1 = 1

    ! Local working precipitation rates
    real(wp) :: precip_l_w(nxy_inner) ! Liquid cloud precip
    real(wp) :: precip_r_w(nxy_inner) ! Rain precip
    real(wp) :: precip_i_w(nxy_inner) ! Ice precip
    real(wp) :: precip_g_w(nxy_inner) ! Graupel precip
    real(wp) :: precip_s_w(nxy_inner) ! Snow precip

    !AH - note that nz is derived in mphys_init and accounts for the lowest level
    !     not equal to 1
    real(wp) :: precip1d(nz,nxy_inner) ! local working precip rate
    real(wp) :: precip_l_w1d(nz,nxy_inner) ! Liquid cloud precip 1D
    real(wp) :: precip_r_w1d(nz,nxy_inner) ! Rain precip 1D
    real(wp) :: precip_i_w1d(nz,nxy_inner) ! Ice precip 1D
    real(wp) :: precip_g_w1d(nz,nxy_inner) ! Graupel precip 1D
    real(wp) :: precip_s_w1d(nz,nxy_inner) ! Snow precip

    integer :: nsubsteps, nsubseds, n_inner

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

    ! Apply RP scheme
    if ( l_rp2_casim ) then
      snow_params%a_x = snow_a_x_rp
      ice_params%a_x = ice_a_x_rp
    end if

! Calculations for substeps do not need to stay inside inner loop
! So moving this part after those calculations
!------
!    do ixy_inner=1, nxy_inner_loop

!       ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
!       jy = modulo(ixy-1,(je_in-js_in+1))+js_in
!       ix = (ixy-1)/(je_in-js_in+1)+is_in

!       qfields_in(:,:,ixy_inner)=qfields(:,:,ixy_inner) ! Initial values of q
!       qfields_mod(:,:,ixy_inner)=qfields(:,:,ixy_inner) ! Modified initial values of q (may be modified if bad values sent in)
!-----
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


   ! starting inner loop here --
    do ixy_inner=1, nxy_inner_loop

       ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
       jy = modulo(ixy-1,(je_in-js_in+1))+js_in
       ix = (ixy-1)/(je_in-js_in+1)+is_in

       if (l_subseds_maxv) then

          if (casim_parent == parent_um) then
             mindz = 20.0
          else
             !mindz = min_dz ! derived in passive_fields
             mindz = min_dz(ixy_inner)
          endif
          ! cloud
          !print *, 'terminal vt called'                    ! be careful of mindz, need to have inner dimension too?
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
       end if !! V - Do we need to calculate the subseds in every timestep?

       qfields_in(:,:,ixy_inner)=qfields(:,:,ixy_inner) ! Initial values of q
       !V -- Don't need to give qields_mod value two times (?)
       !V qfields_mod(:,:,ixy_inner)=qfields(:,:,ixy_inner) ! Modified initial values of q (may be modified if bad values sent in)

       if (l_tendency_loc) then! Parent model uses tendencies
          qfields_mod(:,:,ixy_inner)=qfields_in(:,:,ixy_inner)+dt*dqfields(:,:,ixy_inner)
       else! Parent model uses increments
          qfields_mod(:,:,ixy_inner)=qfields_in(:,:,ixy_inner)+dqfields(:,:,ixy_inner)
       end if

       if (.not. l_passive) then
          call tidy_qin(ixy_inner, qfields_mod(:,:,ixy_inner))
       end if

       !---------------------------------------------------------------
       ! Determine (and possibly limit) size distribution
       !---------------------------------------------------------------
       call query_distributions(ixy_inner, cloud_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
       call query_distributions(ixy_inner, rain_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
       if (.not. l_warm_loc) then
          call query_distributions(ixy_inner, ice_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
          call query_distributions(ixy_inner, snow_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
          call query_distributions(ixy_inner, graupel_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
       end if

       qfields(:,:,ixy_inner)=qfields_mod(:,:,ixy_inner)

       if (aerosol_option > 0) then
          aerofields_in(:,:,ixy_inner)=aerofields(:,:,ixy_inner) ! Initial values of aerosol
          aerofields_mod(:,:,ixy_inner)=aerofields(:,:,ixy_inner) ! Modified initial values  (may be modified if bad values sent in)

          if (l_tendency_loc) then! Parent model uses tendencies
             aerofields_mod(:,:,ixy_inner)=aerofields_in(:,:,ixy_inner)+dt*daerofields(:,:,ixy_inner)
          else! Parent model uses increments
             aerofields_mod(:,:,ixy_inner)=aerofields_in(:,:,ixy_inner)+daerofields(:,:,ixy_inner)
          end if

          if (l_process) call tidy_ain(qfields_mod(:,:,ixy_inner), aerofields_mod(:,:,ixy_inner))

          aerofields(:,:,ixy_inner)=aerofields_mod(:,:,ixy_inner)
       end if

    end do ! ixy_inner


    ! switch n and ixy_inner loop using an "if" condition
    do n = 1, nsubsteps

       n_inner = 0 !(how many columns has precondition=true under n-loop)

       do ixy_inner=1, nxy_inner_loop

          ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
          jy = modulo(ixy-1,(je_in-js_in+1))+js_in
          ix = (ixy-1)/(je_in-js_in+1)+is_in

          !do n=1,nsubsteps

          call preconditioner(ixy_inner, qfields(:,:,ixy_inner))

          if ( casdiags % l_mphys_pts ) then
             ! Set microphysics points flag based on precondition
             do k = 1, nz
                casdiags % mphys_pts(ix, jy, k) = precondition(k,ixy_inner)
             end do
          end if ! casdiags % l_mphys_pts

         !!------------------------------------------------------
         !! Early exit if we will have nothing to do.
         !! (i.e. no hydrometeors and subsaturated)
         !!------------------------------------------------------
         !if (.not. any(precondition(:,ixy_inner))) exit
         if (any(precondition(:,ixy_inner))) then

            !-------------------------------
            ! Derive aerosol distribution
            ! parameters
            !-------------------------------
            if (aerosol_option > 0)                                           &
               call examine_aerosol(aerofields(:,:,ixy_inner),                &
                    qfields(:,:,ixy_inner), aerophys, aerochem, aeroact,      &
                    dustphys, dustchem, dustact, aeroice, dustliq, icall=1)

            ! In order to get rid of "if precondition", here collect the columns that
            ! have "if precondition" to be "true"
            !--------------------------------------------------------------------------------
            n_inner = n_inner + 1
            columns(n_inner) = ixy_inner

         end if
       end do ! ixy_inner


       ! The microphysics will only do in columns with precondtion = true
       if (n_inner .eq. 0) exit ! nothing to do for all columns

       do i_column = 1, n_inner

          ixy_inner = columns(i_column)

          ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
          jy = modulo(ixy-1,(je_in-js_in+1))+js_in
          ix = (ixy-1)/(je_in-js_in+1)+is_in

          ! Later on the microphysics will only do in columns with precondition=true

          do k=1,nz

             l_Twarm=TdegK(k,ixy_inner) > 273.15
             l_Tcold(k,ixy_inner)=.not. l_Twarm

          end do
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
             call raut(ixy_inner, step_length, qfields(:,:,ixy_inner),         &
               cffields(:,:,ixy_inner), aerofields(:,:,ixy_inner),             &
               procs(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner))
          end if

          !-------------------------------
          ! Do the rain accreting cloud
          !-------------------------------
          if (pswitch%l_pracw) then
             call racw(ixy_inner, step_length, qfields(:,:,ixy_inner),         &
               cffields(:,:,ixy_inner), aerofields(:,:,ixy_inner),             &
               procs(:,:,ixy_inner), rain_params, aerosol_procs(:,:,ixy_inner))
          end if

          !-------------------------------
          ! Do the rain self-collection
          !-------------------------------
          if (pswitch%l_pracr) then
              call racr(ixy_inner, step_length, qfields(:,:,ixy_inner),        &
                procs(:,:,ixy_inner))
          end if

          !-------------------------------
          ! Do the evaporation of rain
          !-------------------------------
          if (pswitch%l_prevp) then
             !! initialise l_sigevap to false for all levels so that previous
             !! "trues" are not included when rain_precondition is false
             l_sigevap(:,ixy_inner) = .false.
             call revp(ixy_inner, step_length, nz, qfields(:,:,ixy_inner),     &
               cffields(:,:,ixy_inner), aerophys, aerochem, aeroact, dustliq,  &
               procs(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner),             &
               l_sigevap(:,ixy_inner))
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
                call inuc(ixy_inner,step_length, nz, l_Tcold(:,ixy_inner),     &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), dustphys, aeroact, dustliq,            &
                  aerosol_procs(:,:,ixy_inner))
             end if

             !------------------------------------------------------
             ! Autoconverion to snow
             !------------------------------------------------------
             if (pswitch%l_psaut .and. .not. l_kfsm) then
                call saut(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner),    &
                  qfields(:,:,ixy_inner), procs(:,:,ixy_inner))
             end if
             !------------------------------------------------------
             ! Accretion processes
             !------------------------------------------------------
             ! Ice -> Cloud -> Ice
             if (pswitch%l_piacw) then
                call iacc(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner),    &
                  ice_params, cloud_params, ice_params, qfields(:,:,ixy_inner),&
                  cffields(:,:,ixy_inner), procs(:,:,ixy_inner),               &
                  l_sigevap(:,ixy_inner), aeroact, dustliq,                    &
                  aerosol_procs(:,:,ixy_inner))
             end if
             ! Snow -> Cloud -> Snow
             if (l_sg) then
                if (pswitch%l_psacw .and. .not. l_kfsm) then
                   call iacc(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner), &
                     snow_params, cloud_params, snow_params,                   &
                     qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),          &
                     procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact,    &
                     dustliq, aerosol_procs(:,:,ixy_inner))
                end if
                !
                ! Snow -> Ice -> Snow
                if (pswitch%l_psaci .and. .not. l_kfsm) then
                   call iacc(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner), &
                     snow_params, ice_params, snow_params,                     &
                     qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),          &
                     procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact,    &
                     dustliq, aerosol_procs(:,:,ixy_inner))
                end if
                !
                if (pswitch%l_praci) then
                   ! Rain -> Ice -> Graupel AND Rain -> Ice -> snow, decision made in iacc
                   call iacc(ixy_inner, step_length, nz,  l_Tcold(:,ixy_inner),&
                     rain_params, ice_params, graupel_params,                  &
                     qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),          &
                     procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact,    &
                     dustliq, aerosol_procs(:,:,ixy_inner), snow_params)
                   ! only one call needed and the decision of graupel or snow is made within iacc
                   ! NOTE: this will break kfsm!!

                end if
                if (pswitch%l_psacr .and. .not. l_kfsm) then
                   ! Snow -> Rain -> Graupel AND Snow -> Rain -> Snow, decision made in iacc
                   call iacc(ixy_inner, step_length, nz,  l_Tcold(:,ixy_inner),&
                     snow_params, rain_params, graupel_params,                 &
                     qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),          &
                     procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact,    &
                     dustliq, aerosol_procs(:,:,ixy_inner), snow_params)
                end if
                if (l_g) then
                   ! Graupel -> Cloud -> Graupel
                   if (pswitch%l_pgacw) then
                      call iacc(ixy_inner, step_length, nz,                    &
                        l_Tcold(:,ixy_inner), graupel_params, cloud_params,    &
                        graupel_params, qfields(:,:,ixy_inner),                &
                        cffields(:,:,ixy_inner), procs(:,:,ixy_inner),         &
                        l_sigevap(:,ixy_inner), aeroact, dustliq,              &
                        aerosol_procs(:,:,ixy_inner))
                   end if
                   ! Graupel -> Rain -> Graupel
                   if (pswitch%l_pgacr) then
                      call iacc(ixy_inner, step_length, nz,                    &
                        l_Tcold(:,ixy_inner), graupel_params, rain_params,     &
                        graupel_params, qfields(:,:,ixy_inner),                &
                        cffields(:,:,ixy_inner), procs(:,:,ixy_inner),         &
                        l_sigevap(:,ixy_inner), aeroact, dustliq,              &
                        aerosol_procs(:,:,ixy_inner))
                   end if
                   ! Graupel -> Ice -> Graupel
                   !                   if(pswitch%l_gsaci)call iacc(step_length, k, graupel_params, ice_params, graupel_params, qfields, &
                   !                       procs, aeroact, dustliq, aerosol_procs)
                   ! Graupel -> Snow -> Graupel
                   !                   if(pswitch%l_gsacs)call iacc(step_length, k, graupel_params, snow_params, graupel_params, qfields, &
                   !                       procs, aeroact, dustliq, aerosol_procs)

                ! Graupel -> Ice -> Graupel
                if (pswitch%l_pgaci) then
                   call iacc(ixy_inner, step_length, nz,  l_Tcold(:,ixy_inner), graupel_params, ice_params, &
                        graupel_params, qfields(:,:,ixy_inner), cffields(:,:,ixy_inner), & 
                        procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact, dustliq, aerosol_procs(:,:,ixy_inner))
                end if

                ! Graupel -> Snow -> Graupel
                if (pswitch%l_pgacs) then
                   call iacc(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner),  graupel_params, snow_params, &
                        graupel_params, qfields(:,:,ixy_inner), cffields(:,:,ixy_inner), &
                        procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact, dustliq, aerosol_procs(:,:,ixy_inner))
                end if
                end if
             end if

             !------------------------------------------------------
             ! Small snow accreting cloud should be sent to graupel
             ! (Ikawa & Saito 1991)
             !------------------------------------------------------
             if (.not. l_kfsm .and. l_reisner_graupel_embryo) then
                ! Only do this process when Kalli's single moment code is
                ! not in use and l_reisner_graupel_embryo is true; otherwise we ignore it.
                if (l_g .and. .not. l_onlycollect) then
                call graupel_embryos(ixy_inner, step_length, nz,               &
                  l_Tcold(:,ixy_inner), qfields(:,:,ixy_inner),                &
                  cffields(:,:,ixy_inner), procs(:,:,ixy_inner))
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
                call wetgrowth(ixy_inner, nz, l_Tcold(:,ixy_inner),            &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), l_sigevap(:,ixy_inner))
             end if

             !------------------------------------------------------
             ! Aggregation (self-collection)
             !------------------------------------------------------
             if (pswitch%l_psagg .and. .not. l_kfsm) then
                call ice_aggregation(ixy_inner, step_length, nz,               &
                  l_Tcold(:,ixy_inner), snow_params, qfields(:,:,ixy_inner),   &
                  procs(:,:,ixy_inner))
             end if

             !------------------------------------------------------
             ! Break up (snow only)
             !------------------------------------------------------
             if (pswitch%l_psbrk .and. .not. l_kfsm) then
                call ice_breakup(nz, l_Tcold(:,ixy_inner), snow_params,        &
                  qfields(:,:,ixy_inner), procs(:,:,ixy_inner))
             end if

             !------------------------------------------------------
             ! Ice multiplication (Hallet-mossop)
             !------------------------------------------------------
             if (pswitch%l_pihal .and. .not. l_kfsm) then
                call hallet_mossop(ixy_inner, step_length, nz,                 &
                  cffields(:,:,ixy_inner), procs(:,:,ixy_inner))
             end if

             !------------------------------------------------------
             ! Homogeneous freezing (rain and cloud)
             !------------------------------------------------------
             if (pswitch%l_phomr) then
                call ihom_rain(ixy_inner, step_length, nz,                     &
                  l_Tcold(:,ixy_inner), qfields(:,:,ixy_inner),                &
                  l_sigevap(:,ixy_inner), aeroact, dustliq,                    &
                  procs(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner))
             end if

             if (pswitch%l_phomc) then
                call ihom_droplets(ixy_inner, step_length, nz,                 &
                  l_Tcold(:,ixy_inner), qfields(:,:,ixy_inner), aeroact,       &
                  dustliq, procs(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner))
             endif
          
          !------------------------------------------------------
          ! Droplet shattering
          !------------------------------------------------------
          if (pswitch%l_pidps .and. .not. l_kfsm) then
             call droplet_shattering(ixy_inner, step_length, nz, cffields(:,:,ixy_inner), &
                  qfields(:,:,ixy_inner), procs(:,:,ixy_inner))
          end if

          !------------------------------------------------------
          ! Ice-ice collision (breakup)
          if (pswitch%l_piics .and. .not. l_kfsm) then
             call ice_collision(ixy_inner, step_length, nz, cffields(:,:,ixy_inner), &
                  procs(:,:,ixy_inner))
          end if

             !------------------------------------------------------
             ! Deposition/sublimation of ice/snow/graupel
             !------------------------------------------------------
             if (pswitch%l_pidep) then
                call idep(ixy_inner, step_length, nz,  l_Tcold(:,ixy_inner),   &
                  ice_params, qfields(:,:,ixy_inner), cffields(:,:,ixy_inner), &
                  procs(:,:,ixy_inner), dustact, aeroice,                      &
                  aerosol_procs(:,:,ixy_inner))
             end if

             if (pswitch%l_psdep .and. .not. l_kfsm ) then
                call idep(ixy_inner,step_length, nz, l_Tcold(:,ixy_inner),     &
                  snow_params, qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),&
                  procs(:,:,ixy_inner), dustact, aeroice,                      &
                  aerosol_procs(:,:,ixy_inner))
             end if

             if (pswitch%l_pgdep) then
                call idep(ixy_inner,step_length, nz,  l_Tcold(:,ixy_inner),    &
                  graupel_params, qfields(:,:,ixy_inner),                      &
                  cffields(:,:,ixy_inner), procs(:,:,ixy_inner), dustact,      &
                  aeroice, aerosol_procs(:,:,ixy_inner))
             end if

             if (l_harrington .and. .not. l_onlycollect) then
                call adjust_dep(nz, l_Tcold(:,ixy_inner), procs(:,:,ixy_inner))
             end if

             !-----------------------------------------------------------
             ! Make sure we don't remove more than saturation allows
             !-----------------------------------------------------------
             if (l_idep) then
                call ensure_saturated(ixy_inner, nz, l_Tcold(:,ixy_inner),     &
                  step_length, qfields(:,:,ixy_inner), procs(:,:,ixy_inner),   &
                  (/i_idep, i_sdep, i_gdep/))
             end if
             !-----------------------------------------------------------
             ! Make sure we don't put back more than saturation allows
             !-----------------------------------------------------------
             if (l_isub) then
                call ensure_saturated(ixy_inner, nz, l_Tcold(:,ixy_inner),     &
                  step_length, qfields(:,:,ixy_inner), procs(:,:,ixy_inner),   &
                  (/i_isub, i_ssub, i_gsub/))
             end if
             ! END all processes at T < 0C
             !
             ! start all ice processes that occur T > 0C, i.e. melting
             !------------------------------------------------------
             ! Melting of ice/snow/graupel
             !------------------------------------------------------
             if (pswitch%l_psmlt .and. .not. l_kfsm) then
                call melting(ixy_inner, step_length, nz, snow_params,          &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroice,       &
                  dustact, aerosol_procs(:,:,ixy_inner))
             end if
             if (pswitch%l_pgmlt) then
                call melting(ixy_inner, step_length, nz, graupel_params,       &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroice,       &
                  dustact, aerosol_procs(:,:,ixy_inner))
             end if
             if (pswitch%l_pimlt) then
                call melting(ixy_inner, step_length, nz, ice_params,           &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroice,       &
                  dustact, aerosol_procs(:,:,ixy_inner))
             end if
         end if ! end if .not. warm_loc
         !-----------------------------------------------------------
         ! Make sure we don't remove more than we have to start with
         !-----------------------------------------------------------
         if (.not. l_warm_loc) then
            if (l_pos1) call ensure_positive(nz, step_length,                  &
              qfields(:,:,ixy_inner), procs(:,:,ixy_inner), cloud_params,      &
              (/i_praut, i_pracw, i_iacw, i_sacw, i_gacw, i_homc, i_inuc/),    &
              aeroprocs=aerosol_procs(:,:,ixy_inner),                          &
              iprocs_dependent=(/i_aaut, i_aacw/))

            if (l_pos2) call ensure_positive(nz, step_length,                  &
              qfields(:,:,ixy_inner), procs(:,:,ixy_inner), ice_params,        &
              (/i_raci, i_saci, i_gaci, i_saut, i_isub, i_imlt/),              &
              (/i_ihal, i_idps, i_iics, i_gshd, i_inuc, i_homc, i_iacw, i_idep/))

            if (l_pos3) call ensure_positive(nz, step_length,                  &
              qfields(:,:,ixy_inner), procs(:,:,ixy_inner), rain_params,       &
              (/i_prevp, i_sacr, i_gacr, i_homr/),                             &
              (/i_praut, i_pracw, i_raci, i_gshd, i_smlt, i_gmlt/),            &
              aeroprocs=aerosol_procs(:,:,ixy_inner),                          &
              iprocs_dependent=(/i_arevp/))

            if (l_pos4) call ensure_positive(nz, step_length,                  &
              qfields(:,:,ixy_inner), procs(:,:,ixy_inner), snow_params,       &
              (/i_gacs, i_smlt, i_sacr, i_ssub /),                             &
              (/i_sdep, i_sacw, i_saut, i_saci, i_raci, i_gshd, i_ihal, i_iics/)) 
         else
            if (pswitch%l_praut .and. pswitch%l_pracw) then
                if (l_pos5) call ensure_positive(nz, step_length,              &
                  qfields(:,:,ixy_inner), procs(:,:,ixy_inner), cloud_params,  &
                  (/i_praut, i_pracw/),                                        &
                  aeroprocs=aerosol_procs(:,:,ixy_inner),                      &
                  iprocs_dependent=(/i_aaut, i_aacw/))
            end if

            if (pswitch%l_prevp) then
               if (l_pos6) call ensure_positive(nz, step_length,               &
                 qfields(:,:,ixy_inner), procs(:,:,ixy_inner),                 &
                 rain_params, (/i_prevp/), (/i_praut, i_pracw/),               &
                 aerosol_procs(:,:,ixy_inner), (/i_arevp/))
            end if

         end if

         !-------------------------------
         ! Collect terms we have so far
         !-------------------------------

         call sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner),      &
           tend(:,:,ixy_inner), (/i_praut, i_pracw, i_pracr, i_prevp/),        &
           l_thermalexchange=.true., qfields=qfields(:,:,ixy_inner),           &
           l_passive=l_passive)


         if (.not. l_warm_loc) then
           if (.not. l_no_pgacs_in_sumprocs) then
             call sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                (/i_idep, i_sdep, i_gdep, i_iacw, i_sacr, i_sacw, i_saci, i_raci,&
                i_gacw, i_gacr, i_gaci, i_gacs, i_ihal, i_iics, i_idps,  i_gshd, i_sbrk,&
                i_saut, i_sagg, i_isub, i_ssub, i_gsub/),        &
                l_thermalexchange=.true., qfields=qfields(:,:,ixy_inner),&
                l_passive=l_passive, i_thirdmoment=2)
           else
             call sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                (/i_idep, i_sdep, i_gdep, i_iacw, i_sacr, i_sacw, i_saci, i_raci,&
                i_gacw, i_gacr, i_gaci, i_ihal, i_iics, i_idps, i_gshd, i_sbrk,&
                i_saut, i_sagg, i_isub, i_ssub, i_gsub/),        &
                l_thermalexchange=.true., qfields=qfields(:,:,ixy_inner),&
                l_passive=l_passive, i_thirdmoment=2)
           end if

            call sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner),   &
                 tend(:,:,ixy_inner),                                          &
                 (/i_inuc, i_imlt, i_smlt, i_gmlt, i_homr, i_homc/),           &
                 l_thermalexchange=.true., qfields=qfields(:,:,ixy_inner),     &
                 l_passive=l_passive)
         end if

         call update_q(qfields_mod(:,:,ixy_inner), qfields(:,:,ixy_inner),     &
                                          tend(:,:,ixy_inner), l_fixneg=.true.)

         if (l_process) then
             call ensure_positive_aerosol(nz, step_length,                     &
                  aerofields(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner),     &
                  (/i_aaut, i_aacw, i_aevp, i_arevp, i_dnuc, i_dsub, i_dssub,  &
                  i_dgsub, i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw,        &
                  i_dgacw, i_dsacr, i_dgacr, i_draci, i_dhomr, i_dhomc /) )

            call sum_aprocs(step_length, nz, aerosol_procs(:,:,ixy_inner),     &
                 aerosol_tend(:,:,ixy_inner),                                  &
                 (/i_aaut, i_aacw, i_aevp, i_arevp/) )

            if (.not. l_warm) then
               call sum_aprocs(step_length, nz, aerosol_procs(:,:,ixy_inner),  &
                    aerosol_tend(:,:,ixy_inner),                               &
                    (/i_dnuc, i_dsub, i_dssub, i_dgsub, i_dimlt, i_dsmlt,      &
                    i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr, i_dgacr,      &
                    i_draci, i_dhomr, i_dhomc /) )
            end if

            call update_q(aerofields_mod(:,:,ixy_inner),                       &
                 aerofields(:,:,ixy_inner), aerosol_tend(:,:,ixy_inner),       &
                 l_aerosol=.true.,l_fixneg=.true.)

            !-------------------------------
            ! Re-Derive aerosol distribution
            ! parameters
            !-------------------------------
            call examine_aerosol(aerofields(:,:,ixy_inner),                    &
                 qfields(:,:,ixy_inner), aerophys, aerochem, aeroact,          &
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
            call query_distributions(ixy_inner, cloud_params,                  &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            call query_distributions(ixy_inner, rain_params,                   &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            if (.not. l_warm_loc) then
               call query_distributions(ixy_inner, ice_params,                 &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               call query_distributions(ixy_inner, snow_params,                &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               call query_distributions(ixy_inner, graupel_params,             &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            end if

         else ! Not in UM or not using Paul's cloud fraction in UM
            if (pswitch%l_pcond)then
               call condevp(ixy_inner, step_length, nz, qfields(:,:,ixy_inner),&
                    procs(:,:,ixy_inner), aerophys, aerochem, aeroact,         &
                    dustphys, dustchem, dustliq,                               &
                    aerosol_procs(:,:,ixy_inner), rhcrit_1d)
               !-------------------------------
               ! Collect terms we have so far
               !-------------------------------
               call sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner),&
                    tend(:,:,ixy_inner), (/i_cond/),                           &
                    l_thermalexchange=.true., qfields=qfields(:,:,ixy_inner),  &
                    l_passive=l_passive)

               call update_q(qfields_mod(:,:,ixy_inner),                       &
                    qfields(:,:,ixy_inner), tend(:,:,ixy_inner),               &
                    l_fixneg=.true.)
            end if ! pswitch%l_pcond

            !---------------------------------------------------------------
            ! Re-Determine (and possibly limit) size distribution
            !---------------------------------------------------------------
            call query_distributions(ixy_inner, cloud_params,                  &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            call query_distributions(ixy_inner, rain_params,                   &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            if (.not. l_warm_loc) then
               call query_distributions(ixy_inner, ice_params,                 &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               call query_distributions(ixy_inner, snow_params,                &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               call query_distributions(ixy_inner, graupel_params,             &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            end if

            if (l_process) then
               call sum_aprocs(step_length, nz, aerosol_procs(:,:,ixy_inner),  &
                    aerosol_tend(:,:,ixy_inner), (/i_aact /) )

               call update_q(aerofields_mod(:,:,ixy_inner),                    &
                    aerofields(:,:,ixy_inner), aerosol_tend(:,:,ixy_inner),    &
                    l_aerosol=.true.)

               !-------------------------------
               ! Re-Derive aerosol distribution
               ! parameters
               !-------------------------------
               call examine_aerosol(aerofields(:,:,ixy_inner),                 &
                    qfields(:,:,ixy_inner), aerophys, aerochem, aeroact,       &
                    dustphys, dustchem, dustact, aeroice, dustliq, icall=3)
            end if

         endif ! (casim_parent == parent_um .and. l_prf_cfrac)

!!PRF

         !-------------------------------
         ! Do the sedimentation
         !-------------------------------

         precip_l_w(ixy_inner) = 0.0
         precip_r_w(ixy_inner) = 0.0
         precip_i_w(ixy_inner) = 0.0
         precip_g_w(ixy_inner) = 0.0
         precip_s_w(ixy_inner) = 0.0

         do k = 1, nz
            precip_l_w1d(k, ixy_inner) = 0.0
            precip_r_w1d(k, ixy_inner) = 0.0
            precip_i_w1d(k, ixy_inner) = 0.0
            precip_g_w1d(k, ixy_inner) = 0.0
            precip_s_w1d(k, ixy_inner) = 0.0
         end do

         if ( casdiags % l_process_rates ) then
            call gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=0)
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
                     call zero_procs(procs(:,:,ixy_inner))
                     if (l_process) call zero_procs(aerosol_procs(:,:,ixy_inner))
                     !---------------------------------------------------------------
                     ! Re-Determine (and possibly limit) size distribution
                     !---------------------------------------------------------------
                     call query_distributions(ixy_inner, cloud_params,         &
                          qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                     call query_distributions(ixy_inner, rain_params,          &
                          qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                     if (.not. l_warm_loc) then
                        call query_distributions(ixy_inner, ice_params,        &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                        call query_distributions(ixy_inner, snow_params,       &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                        call query_distributions(ixy_inner, graupel_params,    &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                     end if

                     !-------------------------------
                     ! Re-Derive aerosol distribution
                     ! parameters
                     !-------------------------------
                     if (l_process) call examine_aerosol(                      &
                             aerofields(:,:,ixy_inner), qfields(:,:,ixy_inner),&
                             aerophys, aerochem, aeroact, dustphys, dustchem,  &
                             dustact, aeroice, dustliq, icall=2)
                  end if ! nsed > 1

                  ! NOTE: if l_gamma_online is true then the original CASIM sedimentation will be used,
                  !       which will calculate the gamma function every timestep. This has to be done
                  !       when 3M or diagnostic shape is used, as shape will change. For single and
                  !       and double moment simulations, it is recommended that l_gamma_online is false.
                  !       This is computationally much more efficient (on CPU and GPU).
                  if (pswitch%l_psedl) then
                     if (l_gamma_online) then
                        call sedr(ixy_inner, qfields(:,:,ixy_inner), aeroact,  &
                             dustliq, cloud_params, procs(:,:,ixy_inner),      &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     else
                        call sedr_1M_2M(ixy_inner, sed_length,                 &
                              qfields(:,:,ixy_inner), aeroact, dustliq,        &
                              cloud_params, procs(:,:,ixy_inner),              &
                              aerosol_procs(:,:,ixy_inner),                    &
                              precip1d(:,ixy_inner), l_process)
                     endif

                     precip_l_w(ixy_inner) = precip_l_w(ixy_inner) +           &
                                                     precip1d(level1,ixy_inner)

                     do k = 1, nz
                        precip_l_w1d(k,ixy_inner) = precip_l_w1d(k,ixy_inner)  &
                                                        + precip1d(k,ixy_inner)
                     end do

                     call sum_procs(ixy_inner, sed_length, nz,                 &
                          procs(:,:,ixy_inner), tend(:,:,ixy_inner),           &
                          (/i_psedl/), qfields=qfields(:,:,ixy_inner))

                  end if

                  if (pswitch%l_psedr) then
                     if (l_gamma_online) then
                        call sedr(ixy_inner, qfields(:,:,ixy_inner), aeroact,  &
                             dustliq, rain_params, procs(:,:,ixy_inner),       &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     else
                        call sedr_1M_2M(ixy_inner, sed_length,                 &
                             qfields(:,:,ixy_inner), aeroact, dustliq,         &
                             rain_params, procs(:,:,ixy_inner),                &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     endif

                     precip_r_w(ixy_inner) = precip_r_w(ixy_inner) +           &
                                                     precip1d(level1,ixy_inner)

                     do k = 1, nz
                        precip_r_w1d(k,ixy_inner) = precip_r_w1d(k,ixy_inner) +&
                                                          precip1d(k,ixy_inner)
                     end do

                     call sum_procs(ixy_inner, sed_length, nz,                 &
                          procs(:,:,ixy_inner), tend(:,:,ixy_inner),           &
                          (/i_psedr/), qfields=qfields(:,:,ixy_inner))

                  end if

                  if (.not. l_warm_loc) then

                     if (pswitch%l_psedi) then
                        if (l_gamma_online) then
                           call sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, ice_params,                  &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        else
                           call sedr_1M_2M(ixy_inner, sed_length,              &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                ice_params, procs(:,:,ixy_inner),              &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        endif
                        precip_i_w(ixy_inner) = precip_i_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        do k = 1, nz
                           precip_i_w1d(k,ixy_inner) =                         &
                              precip_i_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        end do

                        call sum_procs(ixy_inner, sed_length, nz,              &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_psedi/), qfields=qfields(:,:,ixy_inner))

                     end if

                     if (pswitch%l_pseds .and. .not. l_kfsm) then
                        if (l_gamma_online) then
                           call sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, snow_params,                 &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        else
                           call sedr_1M_2M(ixy_inner, sed_length,              &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                snow_params, procs(:,:,ixy_inner),             &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        endif
                        precip_s_w(ixy_inner) = precip_s_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        do k = 1, nz
                           precip_s_w1d(k,ixy_inner) =                         &
                              precip_s_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        end do

                        call sum_procs(ixy_inner, sed_length, nz,              &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/ i_pseds /), qfields=qfields(:,:,ixy_inner))

                     end if

                     if (pswitch%l_psedg) then
                        if (l_gamma_online) then
                           call sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, graupel_params,              &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        else
                           call sedr_1M_2M(ixy_inner, sed_length,              &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                graupel_params, procs(:,:,ixy_inner),          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        endif
                        precip_g_w(ixy_inner) = precip_g_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        do k = 1, nz
                           precip_g_w1d(k,ixy_inner) =                         &
                              precip_g_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        end do

                        call sum_procs(ixy_inner, sed_length, nz,              &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_psedg/), qfields=qfields(:,:,ixy_inner))

                     end if

                     !!call sum_procs(sed_length, nz, procs, tend, (/i_psedi, i_pseds, i_psedg/), qfields=qfields)
                  end if !.not. l_warm_loc

                  call update_q(qfields_mod(:,:,ixy_inner),                    &
                       qfields(:,:,ixy_inner), tend(:,:,ixy_inner),            &
                       l_fixneg=.true.)

                  if (l_process) then
                     if (l_warm) then
                        call ensure_positive_aerosol(nz, step_length,          &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             (/i_asedr, i_asedl/) )
                        call sum_aprocs(sed_length, nz,                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             aerosol_tend(:,:,ixy_inner), (/i_asedr, i_asedl/))
                        call update_q(aerofields_mod(:,:,ixy_inner),           &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_tend(:,:,ixy_inner), l_aerosol=.true.)
                    else ! not l_warm - includes ice procs
                        call ensure_positive_aerosol(nz, step_length,          &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             (/i_asedr, i_asedl,i_dsedi, i_dseds, i_dsedg/) )
                        call sum_aprocs(sed_length, nz,                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             aerosol_tend(:,:,ixy_inner),                      &
                             (/i_asedr, i_asedl,i_dsedi, i_dseds, i_dsedg/) )
                        call update_q(aerofields_mod(:,:,ixy_inner),           &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_tend(:,:,ixy_inner), l_aerosol=.true.)
                    end if ! l_warm
                  end if ! l_process

                  if ( casdiags % l_process_rates ) then
                     call gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end, ncall=1)
                  end if

               end do ! nsed

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
                        call query_distributions(ixy_inner, cloud_params,      &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                        !-------------------------------
                        ! Re-Derive aerosol distribution
                        ! parameters
                        !-------------------------------
                        if (l_process) call examine_aerosol(                   &
                                            aerofields(:,:,ixy_inner),         &
                                            qfields(:,:,ixy_inner), aerophys,  &
                                            aerochem, aeroact , dustphys,      &
                                            dustchem, dustact, aeroice,        &
                                            dustliq, icall=2)

                     endif

                     if (l_gamma_online) then
                        call sedr(ixy_inner, qfields(:,:,ixy_inner), aeroact,  &
                             dustliq, cloud_params, procs(:,:,ixy_inner),      &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     else
                        call sedr_1M_2M(ixy_inner, sed_length_cloud,           &
                             qfields(:,:,ixy_inner), aeroact, dustliq,         &
                             cloud_params, procs(:,:,ixy_inner),               &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     endif

                     precip_l_w(ixy_inner) = precip_l_w(ixy_inner) +           &
                                                     precip1d(level1,ixy_inner)

                     do k = 1, nz
                        precip_l_w1d(k,ixy_inner) = precip_l_w1d(k,ixy_inner) +&
                                                          precip1d(k,ixy_inner)
                     end do

                     if ( casdiags % l_process_rates ) then
                        call gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                     end if

                     call sum_procs(ixy_inner, sed_length_cloud, nz,           &
                          procs(:,:,ixy_inner), tend(:,:,ixy_inner),           &
                          (/i_psedl/), qfields=qfields(:,:,ixy_inner))

                     call update_q(qfields_mod(:,:,ixy_inner),                 &
                          qfields(:,:,ixy_inner), tend(:,:,ixy_inner),         &
                          l_fixneg=.true.)

                     if (l_process) then
                        call sum_aprocs(sed_length, nz,                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             aerosol_tend(:,:,ixy_inner), (/i_asedl/) )
                        call update_q(aerofields_mod(:,:,ixy_inner),           &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_tend(:,:,ixy_inner), l_aerosol=.true.)
                     endif
                     ! !-------------------------------
                     ! ! Reset process rates if they
                     ! ! are to be re-used
                     ! !-------------------------------
                     call zero_procs(procs(:,:,ixy_inner), (/i_psedl/))
                     if (l_process) call zero_procs(aerosol_procs(:,:,ixy_inner), (/i_asedl/))
                  end do ! k
               end if ! pswitch%l_psedl

               if (pswitch%l_psedr) then
                  do nsed=1,nsubseds_rain
                     if (nsed > 1) then
                        !---------------------------------------------------------------
                        ! Re-Determine (and possibly limit) size distribution
                        !---------------------------------------------------------------
                        call query_distributions(ixy_inner, rain_params,       &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                        !-------------------------------
                        ! Re-Derive aerosol distribution
                        ! parameters
                        !-------------------------------
                        if (l_process) call examine_aerosol(                   &
                             aerofields(:,:,ixy_inner), qfields(:,:,ixy_inner),&
                             aerophys, aerochem, aeroact, dustphys, dustchem,  &
                             dustact, aeroice, dustliq, icall=2)
                     endif
                     if (l_gamma_online) then
                        call sedr(ixy_inner, qfields(:,:,ixy_inner), aeroact,  &
                             dustliq, rain_params, procs(:,:,ixy_inner),       &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     else
                        call sedr_1M_2M(ixy_inner, sed_length_rain,            &
                             qfields(:,:,ixy_inner), aeroact, dustliq,         &
                             rain_params, procs(:,:,ixy_inner),                &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     endif

                     precip_r_w(ixy_inner) = precip_r_w(ixy_inner) +           &
                                                     precip1d(level1,ixy_inner)

                     do k = 1, nz
                        precip_r_w1d(k,ixy_inner) = precip_r_w1d(k,ixy_inner) +&
                                                          precip1d(k,ixy_inner)
                     end do

                     if ( casdiags % l_process_rates ) then
                        call gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                     end if

                     call sum_procs(ixy_inner, sed_length_rain, nz,            &
                          procs(:,:,ixy_inner), tend(:,:,ixy_inner),           &
                          (/i_psedr/), qfields=qfields(:,:,ixy_inner))
                     call update_q(qfields_mod(:,:,ixy_inner),                 &
                          qfields(:,:,ixy_inner), tend(:,:,ixy_inner),         &
                          l_fixneg=.true.)
                     if (l_process) then
                        call sum_aprocs(sed_length, nz,                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             aerosol_tend(:,:,ixy_inner), (/i_asedr, i_asedl/))
                        call update_q(aerofields_mod(:,:,ixy_inner),           &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_tend(:,:,ixy_inner), l_aerosol=.true.)
                     end if
                     !-------------------------------
                     ! Reset process rates if they
                     ! are to be re-used
                     !-------------------------------
                     call zero_procs(procs(:,:,ixy_inner), (/i_psedr/))
                     if (l_process) call zero_procs(aerosol_procs(:,:,ixy_inner), (/i_asedr/))

                  end do ! nsed
               end if ! pswitch%l_psedr

               if (.not. l_warm_loc) then

                  if (pswitch%l_psedi) then
                     do nsed=1,nsubseds_ice
                        if (nsed > 1) then
                           !---------------------------------------------------------------
                           ! Re-Determine (and possibly limit) size distribution
                           !---------------------------------------------------------------
                           call query_distributions(ixy_inner, ice_params,     &
                               qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                           !-------------------------------
                           ! Re-Derive aerosol distribution
                           ! parameters
                           !-------------------------------
                           if (l_process) call examine_aerosol(                &
                                        aerofields(:,:,ixy_inner),             &
                                        qfields(:,:,ixy_inner), aerophys,      &
                                        aerochem, aeroact , dustphys, dustchem,&
                                        dustact, aeroice, dustliq, icall=2)
                        end if
                        if (l_gamma_online) then
                           call sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, ice_params,                  &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        else
                           call sedr_1M_2M(ixy_inner, sed_length_ice,          &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                ice_params, procs(:,:,ixy_inner),              &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        end if

                        precip_i_w(ixy_inner) = precip_i_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        do k = 1, nz
                           precip_i_w1d(k,ixy_inner) =                         &
                              precip_i_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        end do

                        if ( casdiags % l_process_rates ) then
                           call gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                        end if

                        call sum_procs(ixy_inner, sed_length_ice, nz,          &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_psedi/), qfields=qfields(:,:,ixy_inner))
                        call update_q(qfields_mod(:,:,ixy_inner),              &
                             qfields(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                             l_fixneg=.true.)
                        if (l_process) then
                           call sum_aprocs(sed_length, nz,                     &
                                aerosol_procs(:,:,ixy_inner),                  &
                                aerosol_tend(:,:,ixy_inner), (/i_dsedi/))
                           call update_q(aerofields_mod(:,:,ixy_inner),        &
                                aerofields(:,:,ixy_inner),                     &
                                aerosol_tend(:,:,ixy_inner), l_aerosol=.true.)
                        end if
                        !-------------------------------
                        ! Reset process rates if they
                        ! are to be re-used
                        !-------------------------------
                        call zero_procs(procs(:,:,ixy_inner), (/i_psedi/))
                        if (l_process) call zero_procs(aerosol_procs(:,:,ixy_inner), (/i_dsedi/))
                     end do ! nsed
                  end if !pswitch%l_psedi

                  if (pswitch%l_pseds) then
                     do nsed=1,nsubseds_snow
                        if (nsed > 1) then
                           !---------------------------------------------------------------
                           ! Re-Determine (and possibly limit) size distribution
                           !---------------------------------------------------------------
                           call query_distributions(ixy_inner, snow_params,    &
                               qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                           ! !-------------------------------
                           ! ! Re-Derive aerosol distribution
                           ! ! parameters
                           ! !-------------------------------
                           if (l_process) call examine_aerosol(                &
                                    aerofields(:,:,ixy_inner),                 &
                                    qfields(:,:,ixy_inner), aerophys, aerochem,&
                                    aeroact, dustphys, dustchem, dustact,      &
                                    aeroice, dustliq, icall=2)
                        end if

                        if (l_gamma_online) then
                           call sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, snow_params,                 &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        else
                           call sedr_1M_2M(ixy_inner, sed_length_snow,         &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                snow_params, procs(:,:,ixy_inner),             &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        end if

                        precip_s_w(ixy_inner) = precip_s_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        do k = 1, nz
                           precip_s_w1d(k,ixy_inner) =                         &
                              precip_s_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        end do

                        if ( casdiags % l_process_rates ) then
                           call gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                        end if

                        call sum_procs(ixy_inner, sed_length_snow, nz,         &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_pseds/), qfields=qfields(:,:,ixy_inner))
                        call update_q(qfields_mod(:,:,ixy_inner),              &
                             qfields(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                             l_fixneg=.true.)
                        if (l_process) then
                           call sum_aprocs(sed_length, nz,                     &
                                aerosol_procs(:,:,ixy_inner),                  &
                                aerosol_tend(:,:,ixy_inner), (/i_dseds/) )
                           call update_q(aerofields_mod(:,:,ixy_inner),        &
                                aerofields(:,:,ixy_inner),                     &
                                aerosol_tend(:,:,ixy_inner), l_aerosol=.true.)
                        end if
                        !-------------------------------
                        ! Reset process rates if they
                        ! are to be re-used
                        !-------------------------------
                        call zero_procs(procs(:,:,ixy_inner), (/i_pseds/))
                        if (l_process) call zero_procs(aerosol_procs(:,:,ixy_inner), (/i_dseds/))
                     end do ! nsed
                  end if ! pswitch%l_pseds

                  if (pswitch%l_psedg) then
                     do nsed=1,nsubseds_graupel
                        if (nsed > 1) then
                           !---------------------------------------------------------------
                           ! Re-Determine (and possibly limit) size distribution
                           !---------------------------------------------------------------
                           call query_distributions(ixy_inner, graupel_params, &
                               qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                           !-------------------------------
                           ! Re-Derive aerosol distribution
                           ! parameters
                           !-------------------------------
                           if (l_process) call examine_aerosol(                &
                                   aerofields(:,:,ixy_inner),                  &
                                   qfields(:,:,ixy_inner), aerophys, aerochem, &
                                   aeroact, dustphys, dustchem, dustact,       &
                                   aeroice, dustliq, icall=2)
                        end if
                        if (l_gamma_online) then
                           call sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, graupel_params,              &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        else
                           call sedr_1M_2M(ixy_inner, sed_length_graupel,      &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                graupel_params, procs(:,:,ixy_inner),          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        end if

                        precip_g_w(ixy_inner) = precip_g_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        do k = 1, nz
                           precip_g_w1d(k,ixy_inner) =                         &
                              precip_g_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        end do

                        if ( casdiags % l_process_rates ) then
                           call gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                        end if

                        call sum_procs(ixy_inner, sed_length_graupel, nz,      &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_psedg/), qfields=qfields(:,:,ixy_inner))
                        call update_q(qfields_mod(:,:,ixy_inner),              &
                             qfields(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                             l_fixneg=.true.)
                        if (l_process) then
                           call sum_aprocs(sed_length_graupel, nz,             &
                                aerosol_procs(:,:,ixy_inner),                  &
                                aerosol_tend(:,:,ixy_inner), (/i_dsedg/))

                           call update_q(aerofields_mod(:,:,ixy_inner),        &
                                aerofields(:,:,ixy_inner),                     &
                                aerosol_tend(:,:,ixy_inner), l_aerosol=.true.)
                        end if
                        !-------------------------------
                        ! Reset process rates if they
                        ! are to be re-used
                        !-------------------------------
                        call zero_procs(procs(:,:,ixy_inner), (/i_psedg/))
                        if (l_process) call zero_procs(aerosol_procs(:,:,ixy_inner), (/i_dsedg/))
                     end do ! nsed
                  end if ! pswitch%l_psedg
               end if  ! l_warm_loc
              end if !.not. l_subseds_maxv
         end if ! l_sed

         precip_l(ixy_inner) = precip_l(ixy_inner) + precip_l_w(ixy_inner)
         ! For diagnostic purposes, set precip_r, precip_s and precip to pass out
         ! For the UM, rainfall rate is assumed as sum of all liquid components
         ! (so includes sedimentation of rain and liquid cloud)
         precip_r(ixy_inner) = precip_r(ixy_inner) + precip_r_w(ixy_inner)

         ! For the UM, snowfall rate is assumed to be a sum of all solid components
         ! (so includes ice, snow and graupel)
         precip_s(ixy_inner) = precip_s(ixy_inner) + precip_s_w(ixy_inner)
         precip_i(ixy_inner) = precip_i(ixy_inner) + precip_i_w(ixy_inner)
         ! For the UM, graupel rate is just itself
         precip_g(ixy_inner) = precip_g(ixy_inner) + precip_g_w(ixy_inner)

         do k = 1, nz
            precip_r1d(k,ixy_inner)  = precip_r1d(k,ixy_inner)  +              &
                          precip_l_w1d(k,ixy_inner) + precip_r_w1d(k,ixy_inner)
            precip_g1d(k,ixy_inner)  = precip_g1d(k,ixy_inner)  +              &
                                                      precip_g_w1d(k,ixy_inner)
            precip_s1d(k,ixy_inner)  = precip_s1d(k,ixy_inner)  +              &
                       precip_s_w1d(k,ixy_inner) + precip_i_w1d(k,ixy_inner) + &
                       precip_g_w1d(k,ixy_inner)
            precip_so1d(k,ixy_inner) = precip_so1d(k,ixy_inner) +              &
                          precip_s_w1d(k,ixy_inner) + precip_i_w1d(k,ixy_inner)
         end do

         if (nsubsteps>1)then
            !-------------------------------
            ! Reset process rates if they
            ! are to be re-used
            !-------------------------------
            !call zero_procs_exp(procs)
            call zero_procs(procs(:,:,ixy_inner))
            if (l_process) call zero_procs(aerosol_procs(:,:,ixy_inner))
            !---------------------------------------------------------------
            ! Re-Determine (and possibly limit) size distribution
            !---------------------------------------------------------------
            call query_distributions(ixy_inner, cloud_params,                  &
                 qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            call query_distributions(ixy_inner, rain_params,                   &
                 qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            if (.not. l_warm_loc) then
               call query_distributions(ixy_inner, ice_params,                 &
                    qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               call query_distributions(ixy_inner, snow_params,                &
                    qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               call query_distributions(ixy_inner, graupel_params,             &
                    qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            end if
         end if
       !end if ! precondition
    end do ! i_column
    end do ! nsubsteps

    do ixy_inner=1, nxy_inner_loop

       ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
       jy = modulo(ixy-1,(je_in-js_in+1))+js_in
       ix = (ixy-1)/(je_in-js_in+1)+is_in

       ! We want the mean precipitation over the parent timestep - so divide
       ! by the total number of substeps - multiply by inv_allsubs is quicker
       precip_l(ixy_inner) = precip_l(ixy_inner) * inv_allsubs_cloud
       precip_r(ixy_inner) = precip_r(ixy_inner) * inv_allsubs_rain
       precip_i(ixy_inner) = precip_i(ixy_inner) * inv_allsubs_ice
       precip_s(ixy_inner) = precip_s(ixy_inner) * inv_allsubs_snow
       precip_g(ixy_inner) = precip_g(ixy_inner) * inv_allsubs_graupel

       ! UM precip rates are
       precip_r(ixy_inner) = precip_l(ixy_inner) + precip_r(ixy_inner)
       precip_s(ixy_inner) = precip_i(ixy_inner) + precip_s(ixy_inner) + precip_g(ixy_inner)

       do k = 1, nz
          precip_r1d(k,ixy_inner)  = precip_r1d(k,ixy_inner)  * inv_allsubs
          precip_s1d(k,ixy_inner)  = precip_s1d(k,ixy_inner)  * inv_allsubs
          precip_so1d(k,ixy_inner) = precip_so1d(k,ixy_inner) * inv_allsubs
          precip_g1d(k,ixy_inner)  = precip_g1d(k,ixy_inner)  * inv_allsubs
       end do ! k

       ! Precip is a sum of everything, so just add rain and snow together which
       ! has all components added.
       ! Do not add precip_g, otherwise graupel contributions will be double-counted
       precip(ix,jy) = precip_r(ixy_inner)   + precip_s(ixy_inner)

       !--------------------------------------------------
       ! Tidy up any small/negative numbers
       ! we may have generated.
       !--------------------------------------------------
       if (pswitch%l_tidy2) then

          call qtidy(ixy_inner, step_length, nz, qfields(:,:,ixy_inner),       &
               procs(:,:,ixy_inner), aerofields(:,:,ixy_inner), aeroact,       &
               dustact, aeroice, dustliq , aerosol_procs(:,:,ixy_inner),       &
               i_tidy2, i_atidy2, l_negonly=l_tidy_negonly)

          call sum_procs(ixy_inner, step_length, nz,                           &
               procs(:,:,ixy_inner), tend(:,:,ixy_inner), (/i_tidy2/),         &
               qfields=qfields(:,:,ixy_inner), l_passive=l_passive)

          call update_q(qfields_mod(:,:,ixy_inner), qfields(:,:,ixy_inner), tend(:,:,ixy_inner))

          if (l_process) then
             call sum_aprocs(step_length, nz, aerosol_procs(:,:,ixy_inner),    &
                  aerosol_tend(:,:,ixy_inner), (/i_atidy2/) )
             call update_q(aerofields_mod(:,:,ixy_inner),                      &
                  aerofields(:,:,ixy_inner), aerosol_tend(:,:,ixy_inner),      &
                  l_aerosol=.true.)
          end if
       end if

       !
       ! Add on initial adjustments that may have been made
       !

       if (l_tendency_loc) then! Convert back from cumulative value to tendency
          tend(:,:,ixy_inner)=tend(:,:,ixy_inner)+qfields_mod(:,:,ixy_inner)-  &
                           qfields_in(:,:,ixy_inner)-dqfields(:,:,ixy_inner)*dt
          tend(:,:,ixy_inner)=tend(:,:,ixy_inner)/dt
       else
          tend(:,:,ixy_inner)=tend(:,:,ixy_inner)+qfields_mod(:,:,ixy_inner)-  &
                              qfields_in(:,:,ixy_inner)-dqfields(:,:,ixy_inner)
          ! prevent negative values
          do iq=i_hstart,ntotalq
             do k=1,nz
                tend(k,iq,ixy_inner)=max(tend(k,iq,ixy_inner),                 &
                        -(qfields_in(k,iq,ixy_inner)-dqfields(k,iq,ixy_inner)))
            end do
          end do
       end if

       if (aerosol_option > 0) then
          ! processing
          if (l_process) then
             if (l_tendency_loc) then! Convert back from cumulative value to tendency
                aerosol_tend(:,:,ixy_inner)=aerosol_tend(:,:,ixy_inner)+       &
                    aerofields_mod(:,:,ixy_inner)-aerofields_in(:,:,ixy_inner)-&
                    daerofields(:,:,ixy_inner)*dt
                aerosol_tend(:,:,ixy_inner)=aerosol_tend(:,:,ixy_inner)/dt
             else
                aerosol_tend(:,:,ixy_inner)=aerosol_tend(:,:,ixy_inner)+       &
                    aerofields_mod(:,:,ixy_inner)-aerofields_in(:,:,ixy_inner)-&
                    daerofields(:,:,ixy_inner)
                ! prevent negative values
                do iq=1,ntotala
                   do k=1,nz
                      aerosol_tend(k,iq,ixy_inner)=                            &
                           max(aerosol_tend(k,iq,ixy_inner),                   &
                           -(aerofields_in(k,iq,ixy_inner)-daerofields(k,iq,ixy_inner)))
                   end do
                end do
             end if
          end if
       end if

    end do ! ixy_inner

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
             if (iqx==i_ni .and. qfields(k,iqx)<=0.0) then
                qfields(k,iqx)=0.0
                qfields(k,i_qi)=0.0
             end if
             if (iqx==i_nr .and. qfields(k,iqx)<=0.0) then
                qfields(k,iqx)=0.0
                qfields(k,i_qr)=0.0
                if (i_m3r/=0) qfields(k,i_m3r)=0.0
             end if
             if (iqx==i_nl .and. qfields(k,iqx)<=0.0) then
                qfields(k,iqx)=0.0
                qfields(k,i_ql)=0.0
             end if
             if (iqx==i_ns .and. qfields(k,iqx)<=0.0) then
                qfields(k,iqx)=0.0
                qfields(k,i_qs)=0.0
                if (i_m3s/=0) qfields(k,i_m3s)=0.0
             end if
             if (iqx==i_ng .and. qfields(k,iqx)<=0.0) then
                qfields(k,iqx)=0.0
                qfields(k,i_qg)=0.0
                if (i_m3g/=0) qfields(k,i_m3g)=0.0
             end if
          end do
       end do
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine update_q

  subroutine gather_process_diagnostics(ixy_inner, i, j, k_start, k_end,ncall)

    ! Gathers all process rate diagnostics if in use and outputs them to the
    ! CASIM generic diagnostic fields, ready for use in any model.

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Indices of this particular grid square
    integer, intent(in) :: ixy_inner, i, j,ncall
    integer, intent(in) :: k_start, k_end ! Start/end points of grid

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
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % phomc(i,j,k) = procs(ice_params%i_1m,i_homc%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % phomc(i,j,:) = ZERO_REAL_WP
      END IF
    END iF ! casdiags % l_phomc

    IF (casdiags % l_nhomc) THEN
      IF ((pswitch%l_phomc) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nhomc(i,j,k) = procs(ice_params%i_2m,i_homc%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nhomc(i,j,:) = ZERO_REAL_WP
      END IF
    END iF ! casdiags % l_phomc

    IF (casdiags % l_pinuc) THEN
      IF (pswitch%l_pinuc) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pinuc(i,j,k) = procs(ice_params%i_1m,i_inuc%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pinuc(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! casdiags % l_pinuc

    IF (casdiags % l_ninuc) THEN
      IF ((pswitch%l_pinuc) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ninuc(i,j,k) = procs(ice_params%i_2m,i_inuc%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ninuc (i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! casdiags % l_pinuc

    IF (casdiags % l_pidep) THEN
      IF (pswitch%l_pidep) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pidep(i,j,k) = procs(ice_params%i_1m,i_idep%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pidep(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psdep) THEN
      IF (pswitch%l_psdep) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psdep(i,j,k) = procs(snow_params%i_1m,i_sdep%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psdep(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_piacw) THEN
      IF (pswitch%l_piacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % piacw(i,j,k) = procs(ice_params%i_1m,i_iacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % piacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_niacw) THEN
      IF ((pswitch%l_piacw) .and. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % niacw(i,j,k) = procs(cloud_params%i_2m,i_iacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % niacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF  

    IF (casdiags % l_psacw) THEN
      IF (pswitch%l_psacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psacw(i,j,k) = procs(snow_params%i_1m,i_sacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsacw) THEN
      IF ((pswitch%l_psacw) .and. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsacw(i,j,k) = procs(cloud_params%i_2m,i_sacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psacr) THEN
      IF (pswitch%l_psacr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psacr(i,j,k) = procs(snow_params%i_1m,i_sacr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psacr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsacr) THEN
      IF ((pswitch%l_psacr) .and. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsacr(i,j,k) = procs(rain_params%i_2m,i_sacr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsacr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pisub) THEN
      IF (pswitch%l_pisub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pisub(i,j,k) = -1.0 * procs(ice_params%i_1m,i_isub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pisub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nisub) THEN
      IF ((pswitch%l_pisub) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nisub(i,j,k) = -1.0 * procs(ice_params%i_2m,i_isub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nisub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pssub) THEN
      IF (pswitch%l_pssub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pssub(i,j,k) = -1.0 * procs(snow_params%i_1m,i_ssub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pssub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nssub) THEN
      IF ((pswitch%l_pssub) .and. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nssub(i,j,k) = -1.0 * procs(snow_params%i_2m,i_ssub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nssub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pimlt) THEN
      IF (pswitch%l_pimlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pimlt(i,j,k) = procs(rain_params%i_1m,i_imlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pimlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nimlt) THEN
      IF ((pswitch%l_pimlt) .and. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nimlt(i,j,k) = procs(rain_params%i_2m,i_imlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nimlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psmlt) THEN
      IF (pswitch%l_psmlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psmlt(i,j,k) = procs(rain_params%i_1m,i_smlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsmlt) THEN
      IF ((pswitch%l_psmlt) .and. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsmlt(i,j,k) = procs(rain_params%i_2m,i_smlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psaut) THEN
      IF (pswitch%l_psaut) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psaut(i,j,k) = procs(snow_params%i_1m,i_saut%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psaut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsaut) THEN
      IF ((pswitch%l_psaut) .and. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsaut(i,j,k) = procs(snow_params%i_2m,i_saut%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsaut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psaci) THEN
      IF (pswitch%l_psaci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psaci(i,j,k) = procs(snow_params%i_1m,i_saci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psaci(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsaci) THEN
      IF ((pswitch%l_psaci) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsaci(i,j,k) = procs(ice_params%i_2m,i_saci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsaci(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_praut) THEN
      IF (pswitch%l_praut) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % praut(i,j,k) = procs(rain_params%i_1m,i_praut%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % praut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nraut) THEN
      IF ((pswitch%l_praut) .and. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nraut(i,j,k) = procs(rain_params%i_2m,i_praut%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nraut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pracw) THEN
      IF (pswitch%l_pracw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pracw(i,j,k) = procs(rain_params%i_1m,i_pracw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pracw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nracw) THEN
      IF ((pswitch%l_pracw) .and. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nracw(i,j,k) = procs(cloud_params%i_2m,i_pracw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nracw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nracr) THEN
      IF ((pswitch%l_pracr) .and. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nracr(i,j,k) = procs(rain_params%i_2m,i_pracr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nracr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_prevp) THEN
      IF (pswitch%l_prevp) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % prevp(i,j,k) = -1.0 * procs(rain_params%i_1m,i_prevp%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % prevp(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nrevp) THEN
      IF ((pswitch%l_prevp) .and. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nrevp(i,j,k) = -1.0 * procs(rain_params%i_2m,i_prevp%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nrevp(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

! #if DEF_MODEL==MODEL_KiD
!      IF (casdiags % l_praut) THEN
!        IF (pswitch%l_praut) THEN
!           DO k = k_start, k_end
!              kc = k - k_start + 1
!              call save_dg(k, casdiags % praut(i,j,k) , 'praut', i_dgtime)
!           END DO
!        ENDIF
!     ENDIF

!     IF (casdiags % l_pracw) THEN
!         IF (pswitch%l_pracw) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, casdiags % pracw(i,j,k) , 'pracw', i_dgtime)
!            END DO
!         END IF
!      END IF

!      IF (casdiags % l_prevp) THEN
!         IF (pswitch%l_prevp) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, casdiags % prevp(i,j,k) , 'prevp', i_dgtime)
!             END DO
!         END IF
!      END IF

!      IF (casdiags % l_psedr) THEN
!         IF (pswitch%l_psedr) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, procs(rain_params%i_1m,i_psedr%id,ixy_inner)%column_data(kc) , 'psedr', i_dgtime)
!            END DO
!         END IF
!      END IF

!      IF (casdiags % l_psedl) THEN
!         IF (pswitch%l_psedl) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, procs(cloud_params%i_1m,i_psedl%id,ixy_inner)%column_data(kc) , 'psedl', i_dgtime)
!            END DO
!         END IF
!      END IF

!      IF (casdiags % l_pracr) THEN
!         IF (pswitch%l_pracr) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, procs(rain_params%i_1m,i_pracr%id,ixy_inner)%column_data(kc), 'pracr', i_dgtime)
!            END DO
!         END IF
!      END IF
! #endif
    IF (casdiags % l_pgacw) THEN
      IF (pswitch%l_pgacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgacw(i,j,k) = procs(graupel_params%i_1m, i_gacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngacw) THEN
      IF ((pswitch%l_pgacw) .and. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngacw(i,j,k) = procs(cloud_params%i_2m, i_gacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgacs) THEN
      IF (pswitch%l_pgacs) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgacs(i,j,k) = procs(graupel_params%i_1m, i_gacs%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgacs(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngacs) THEN
      IF ((pswitch%l_pgacs) .and. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngacs(i,j,k) = procs(snow_params%i_2m, i_gacs%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngacs(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgmlt) THEN
      IF (pswitch%l_pgmlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgmlt(i,j,k) = procs(rain_params%i_1m, i_gmlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngmlt) THEN
      IF ((pswitch%l_pgmlt) .and. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngmlt(i,j,k) = procs(rain_params%i_2m, i_gmlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgsub) THEN
      IF (pswitch%l_pgsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgsub(i,j,k) = -1.0 * procs(graupel_params%i_1m,i_gsub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgsub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngsub) THEN
      IF ((pswitch%l_pgsub) .and. (graupel_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngsub(i,j,k) = -1.0 * procs(graupel_params%i_2m, i_gsub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngsub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedi) THEN
      IF (pswitch%l_psedi) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedi(i,j,k) = procs(ice_params%i_1m, i_psedi%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psedi(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedi) THEN
      IF ((pswitch%l_psedi) .and. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedi(i,j,k) = procs(ice_params%i_2m,i_psedi%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedi(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pseds) THEN
      IF (pswitch%l_pseds) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pseds(i,j,k) = procs(snow_params%i_1m, i_pseds%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pseds(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nseds) THEN
      IF ((pswitch%l_pseds) .and. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nseds(i,j,k) = procs(snow_params%i_2m, i_pseds%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nseds(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedr) THEN
      IF (pswitch%l_psedr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedr(i,j,k) = procs(rain_params%i_1m,i_psedr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psedr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedr) THEN
      IF ((pswitch%l_psedr) .and. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedr(i,j,k) = procs(rain_params%i_2m,i_psedr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedg) THEN
      IF (pswitch%l_psedg) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedg(i,j,k) = procs(graupel_params%i_1m,i_psedg%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psedg(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedg) THEN
      IF ((pswitch%l_psedg) .and. (graupel_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedg(i,j,k) = procs(graupel_params%i_2m,i_psedg%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedg(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedl) THEN
      IF (pswitch%l_psedl) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedl(i,j,k) = procs(cloud_params%i_1m,i_psedl%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psedl(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedl) THEN
      IF ((pswitch%l_psedl) .and. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedl(i,j,k) = procs(cloud_params%i_2m,i_psedl%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedl(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pcond) THEN
      IF (pswitch%l_pcond) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pcond(i,j,k) = procs(cloud_params%i_1m,i_cond%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pcond(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_phomr) THEN
      IF (pswitch%l_phomr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % phomr(i,j,k) = procs(graupel_params%i_1m,i_homr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % phomr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

   IF (casdiags % l_pihal) THEN
      IF (pswitch%l_pihal) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pihal(i,j,k) = procs(ice_params%i_1m,i_ihal%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pihal(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nihal) THEN
      IF ((pswitch%l_pihal) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nihal(i,j,k) = procs(ice_params%i_2m,i_ihal%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nihal(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nhomr) THEN
      IF ((pswitch%l_phomr) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nhomr(i,j,k) = procs(graupel_params%i_2m,i_homr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nhomr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_praci_g) THEN
      IF (pswitch%l_praci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % praci_g(i,j,k) = procs(graupel_params%i_1m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % praci_g(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_praci_r) THEN
      IF (pswitch%l_praci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % praci_r(i,j,k) = procs(rain_params%i_1m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % praci_r(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_praci_i) THEN
      IF (pswitch%l_praci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % praci_i(i,j,k) = procs(ice_params%i_1m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % praci_i(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nraci_g) THEN
      IF ((pswitch%l_praci) .and. (graupel_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nraci_g(i,j,k) = procs(graupel_params%i_2m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nraci_g(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nraci_r) THEN
      IF ((pswitch%l_praci) .and. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nraci_r(i,j,k) = procs(rain_params%i_2m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nraci_r(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nraci_i) THEN
      IF ((pswitch%l_praci) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nraci_i(i,j,k) = procs(ice_params%i_2m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nraci_i(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pidps) THEN
      IF ((pswitch%l_pidps) .and. (ice_params%l_1m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pidps(i,j,k) = procs(ice_params%i_1m,i_idps%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pidps(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nidps) THEN
      IF ((pswitch%l_pidps) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nidps(i,j,k) = procs(ice_params%i_2m,i_idps%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nidps(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgaci) THEN
      IF ((pswitch%l_pgaci) .and. (ice_params%l_1m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgaci(i,j,k) = procs(ice_params%i_1m,i_gaci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgaci(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngaci) THEN
      IF ((pswitch%l_pgaci) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngaci(i,j,k) = procs(ice_params%i_2m,i_gaci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngaci(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_niics_s) THEN
      IF ((pswitch%l_piics) .and. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % niics_s(i,j,k) = procs(snow_params%i_2m,i_iics%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % niics_s(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_niics_i) THEN
      IF ((pswitch%l_piics) .and. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % niics_i(i,j,k) = procs(ice_params%i_2m,i_iics%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % niics_i(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    else !ncall > 0

    IF (casdiags % l_psedi) THEN
      IF (pswitch%l_psedi) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedi(i,j,k) = casdiags % psedi(i,j,k)+                   &
                    procs(ice_params%i_1m,i_psedi%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nsedi) THEN
      IF ((pswitch%l_psedi) .and. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedi(i,j,k) = casdiags % nsedi(i,j,k)+                  &
                   procs(ice_params%i_2m,i_psedi%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_pseds) THEN
      IF (pswitch%l_pseds) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pseds(i,j,k) = casdiags % pseds(i,j,k)+                  &
                  procs(snow_params%i_1m,i_pseds%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nseds) THEN
      IF ((pswitch%l_pseds) .and. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nseds(i,j,k) = casdiags % nseds(i,j,k)+                   &
                   procs(snow_params%i_2m,i_pseds%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_psedr) THEN
      IF (pswitch%l_psedr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedr(i,j,k) = casdiags % psedr(i,j,k)+                   &
                   procs(rain_params%i_1m,i_psedr%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_psedg) THEN
      IF (pswitch%l_psedg) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedg(i,j,k) = casdiags % psedg(i,j,k)+                   &
                procs(graupel_params%i_1m,i_psedg%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nsedg) THEN
      IF ((pswitch%l_psedg) .and. (graupel_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedg(i,j,k) = casdiags % nsedg(i,j,k)+                   &
                procs(graupel_params%i_2m,i_psedg%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_psedl) THEN
      IF (pswitch%l_psedl) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedl(i,j,k) = casdiags % psedl(i,j,k)+                   &
                  procs(cloud_params%i_1m,i_psedl%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nsedl) THEN
      IF ((pswitch%l_psedl) .and. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedl(i,j,k) = procs(cloud_params%i_2m,i_psedl%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedl(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    end if ! ncall

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine gather_process_diagnostics
end module micro_main
