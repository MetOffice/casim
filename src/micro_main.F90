module micro_main
  use mphys_die, only: throw_mphys_error
  use variable_precision, only: wp
  use mphys_parameters, only: nz, nq, naero, nprocs, naeroprocs, rain_params, cloud_params, ice_params, &
       snow_params, graupel_params, nspecies, parent_dt, ZERO_REAL_WP
  use process_routines, only: process_rate, zero_procs, allocate_procs, deallocate_procs, i_cond, i_praut, &
       i_pracw, i_pracr, i_prevp, i_psedr, i_psedl, i_aact, i_aaut, i_aacw, i_aevp, i_asedr, i_asedl, i_arevp, &
       i_tidy, i_tidy2, i_atidy, i_atidy2, i_inuc, i_idep, i_dnuc, i_dsub, i_saut, i_iacw, i_sacw, i_pseds, &
       i_sdep, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs, i_gdep, i_psedg, i_iagg, i_sagg, i_gagg, &
       i_gshd, i_ihal, i_smlt, i_gmlt, i_psedi, i_homr, i_homc, i_imlt, i_isub, i_ssub, i_gsub, i_sbrk, i_dssub, &
       i_dgsub, i_dsedi, i_dseds, i_dsedg, i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr, &
       i_dgacr, i_draci, process_name
  use sum_process, only: sum_procs, sum_aprocs
  use aerosol_routines, only: examine_aerosol, aerosol_phys, aerosol_chem, aerosol_active, allocate_aerosol, &
       deallocate_aerosol, MNtoRm
  use mphys_switches, only: hydro_complexity, aero_complexity, i_qv, i_ql, i_nl, i_qr, i_nr, i_m3r, i_th, i_qi, &
       i_qs, i_qg, i_ni, i_ns, i_ng, i_m3s, i_m3g, i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am4, i_am5, i_am6, &
       i_an6, i_am7, i_am8 , i_am9, i_am10, i_an10, i_an11, i_an12, hydro_names, aero_names,  aerosol_option, l_warm, &
       l_inuc, l_sed, l_condensation, l_iaut, l_imelt, l_iacw, l_idep, aero_index, l_sedl, nq_l, nq_r, nq_i, nq_s, nq_g, &
       l_sg, l_g, l_process, l_halletmossop, max_sed_length, max_step_length, l_harrington, l_passive, ntotala, ntotalq, &
       active_number, isol, iinsol, l_raci_g, l_onlycollect, l_pracr, pswitch, l_isub, l_pos1, l_pos2, l_pos3, l_pos4, &
       l_pos5, l_pos6, i_hstart, nsubsteps, nsubseds, l_tidy_negonly

  use distributions, only: query_distributions, dist_lambda, dist_mu, dist_n0
  use passive_fields, only: initialise_passive_fields, set_passive_fields, TdegK
  use autoconversion, only: raut
  use evaporation, only: revp
  use condensation, only: condevp
  use accretion, only: racw
  use aggregation, only: racr, ice_aggregation
  use sedimentation, only: sedr
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
  use lookup, only: get_slope_generic, moment
  use thresholds, only: qr_small, nr_small, m3r_small, thresh_sig, thresh_tidy
  use mphys_constants, only:  fixed_aerosol_number, fixed_aerosol_rm, &
       fixed_aerosol_sigma, fixed_aerosol_density

  use mphys_tidy, only: qtidy, ensure_positive, ensure_saturated, tidy_qin, tidy_ain, ensure_positive_aerosol
  use preconditioning, only: precondition, preconditioner

#if DEF_MODEL==MODEL_LEM
  use diaghelp_lem, only: i_here, j_here, k_here, koff, n_sub, n_subsed
#endif

#if DEF_MODEL==MODEL_KiD
  ! Kid modules
  use diagnostics, only: save_dg, i_dgtime, i_here, j_here, k_here, n_sub, n_subsed
  use runtime, only: time
  use parameters, only: nx
#elif DEF_MODEL==MODEL_LEM_DIAG
  use initialize, only: DTPUD

  use diags_lem, only: lem_procdgs
  use com_prametr, only: ilatentdgp, nprc_prametr => nprc, iforscalp
  use com_dgstore, only: procrate, fallq
  use com_dgchars, only: procchar
  use com_rain, only: puddle
  use com_classy, only: precip_lem => precip
#elif DEF_MODEL==MODEL_UM
  use diaghelp_um, only: i_here, j_here, k_here, l_debug_um, debug_i, debug_j, debug_k, debug_pe, debug_i2, debug_j2, debug_k2, &
       n_sub, n_subsed
  use timestep_mod, only: timestep_number
  use mphys_casim_diagnostics, only: SurfacePrecip, ProcessRates
#elif  DEF_MODEL==MODEL_MONC
  use diaghelp_monc, only: i_here, j_here, k_here, n_sub, n_subsed, mype
#endif

  implicit none

  private

  logical :: l_tendency_loc
  logical :: l_warm_loc

  integer :: is, ie ! upper and lower i levels which are to be used
  integer :: js, je ! upper and lower j levels
  integer :: ks, ke ! upper and lower k levels    

  real(wp), allocatable :: rhcrit(:)
  real(wp), allocatable, save :: precip(:,:) ! diagnostic for surface precip rate
  real(wp), allocatable :: dqfields(:,:), qfields(:,:), tend(:,:)
  real(wp), allocatable :: daerofields(:,:), aerofields(:,:), aerosol_tend(:,:)

  type(process_rate), allocatable :: procs(:,:)
  type(process_rate), allocatable :: aerosol_procs(:,:)

  type(aerosol_active), allocatable :: aeroact(:)
  type(aerosol_phys), allocatable   :: aerophys(:)
  type(aerosol_chem), allocatable   :: aerochem(:)

  type(aerosol_active), allocatable :: dustact(:)
  type(aerosol_phys), allocatable   :: dustphys(:)
  type(aerosol_chem), allocatable   :: dustchem(:)

  type(aerosol_active), allocatable :: aeroice(:)  ! Soluble aerosol in ice
  type(aerosol_active), allocatable :: dustliq(:)! Insoluble aerosol in liquid

  public initialise_casim, finalise_casim, shipway_microphysics
contains

  subroutine initialise_casim(il, iu, jl, ju, kl, ku,                 &       
       is_in, ie_in, js_in, je_in, ks_in, ke_in, l_tendency, rhcrit_in)

    integer, intent(in) :: il, iu ! upper and lower i levels
    integer, intent(in) :: jl, ju ! upper and lower j levels
    integer, intent(in) :: kl, ku ! upper and lower k levels

    integer, intent(in), optional :: is_in, ie_in ! upper and lower i levels which are to be used
    integer, intent(in), optional :: js_in, je_in ! upper and lower j levels
    integer, intent(in), optional :: ks_in, ke_in ! upper and lower k levels

    ! New optional l_tendency logical added...
    ! if true then a tendency is returned (i.e. units/s)
    ! if false then an increment is returned (i.e. units/timestep)
    logical, intent(in), optional :: l_tendency
    real(wp), intent(in), optional :: rhcrit_in(kl:ku)

    ! Local variables
    real(wp) :: p1, p2, p3

    integer :: n, k, i, j, nxny, imode

    integer :: NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

    integer, parameter :: CHUNKSIZE=10

    real(wp) :: rain_number
    real(wp) :: rain_mass
    real(wp) :: rain_m3
    real(wp) :: n0,lam,mu
    real(wp) :: m1,m2,m3

    !diagnostics
    real(wp) :: field, res
    character(100) :: units, name
    integer :: im
    character(1) :: char

    !    if (timestep_number*dt < 3600)then
    !      l_warm_loc=.true. ! Force warm in spinup
    !    else
    l_warm_loc=l_warm ! Original setting
    !    end if

    ! Set grid extents to operate on
    if (present(is_in)) is=is_in
    if (present(ie_in)) ie=ie_in
    if (present(js_in)) js=js_in
    if (present(je_in)) je=je_in
    if (present(ks_in)) ks=ks_in
    if (present(ke_in)) ke=ke_in

    ! if not passed in, then default to full grid
    if (.not. present(is_in)) is=il
    if (.not. present(ie_in)) ie=iu
    if (.not. present(js_in)) js=jl
    if (.not. present(je_in)) je=ju
    if (.not. present(ks_in)) ks=kl
    if (.not. present(ke_in)) ke=ku

    allocate(rhcrit(kl:ku))

    l_tendency_loc=.true.
    if (present(l_tendency)) l_tendency_loc=l_tendency
    if (present(rhcrit_in)) then
      rhcrit(:)=rhcrit_in(:)
    else
      rhcrit(:)=1.0
    end if

    nq=sum(hydro_complexity%nmoments)+2 ! also includes vapour and theta
    nz=ke-ks+1
    nprocs = hydro_complexity%nprocesses
    nxny=(ie-is+1)*(je-js+1)

    allocate(precondition(nz))
    precondition=.true. ! Assume all points need to be considered
    allocate(qfields(nz, nq))
    allocate(dqfields(nz, nq))
    qfields=ZERO_REAL_WP
    dqfields=ZERO_REAL_WP
    allocate(procs(nz, nprocs))
    allocate(tend(nz, nq))

    ! Allocate aerosol storage
    if (aerosol_option > 0) then
      naero=ntotala
      naeroprocs=aero_complexity%nprocesses
      allocate(aerofields(nz, naero))
      allocate(daerofields(nz, naero))
      aerofields=ZERO_REAL_WP
      daerofields=ZERO_REAL_WP
      allocate(aerosol_procs(nz, naeroprocs))
      allocate(aerosol_tend(nz, naero))
    else
      ! Dummy arrays required
      allocate(aerofields(1,1))
      allocate(daerofields(1,1))
      allocate(aerosol_procs(1,1))
      allocate(aerosol_tend(1,1))
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

    ! Temporary initialization of chem and sigma
    do k =1,size(aerophys)
      aerophys(k)%sigma(:)=fixed_aerosol_sigma
      aerophys(k)%rpart(:)=0.0
      aerochem(k)%vantHoff(:)=3.0
      aerochem(k)%massMole(:)=132.0e-3
      aerochem(k)%density(:)=fixed_aerosol_density
      aerochem(k)%epsv(:)=1.0
      aerochem(k)%beta(:)=1.0
    end do
    do k =1,size(dustphys)
      dustphys(k)%sigma(:)=fixed_aerosol_sigma
      dustphys(k)%rpart(:)=0.0
      dustchem(k)%vantHoff(:)=3.0
      dustchem(k)%massMole(:)=132.0e-3
      dustchem(k)%density(:)=fixed_aerosol_density
      dustchem(k)%epsv(:)=1.0
      dustchem(k)%beta(:)=1.0
    end do

    !allocate space for the process rates
    call allocate_procs(procs, nz, nprocs, ntotalq)
    if (l_process) call allocate_procs(aerosol_procs, nz, naeroprocs, ntotala)

    ! allocate diagnostics
    allocate(precip(il:iu,jl:ju))   

    call initialise_passive_fields(ks, ke)
  end subroutine initialise_casim

  subroutine finalise_casim()
    ! deallocate diagnostics
    deallocate(precip)

    ! deallocate process rates
    call deallocate_procs(procs)

    deallocate(procs)
    deallocate(qfields)
    deallocate(tend)
    deallocate(precondition)

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
    deallocate(aerofields)
    deallocate(daerofields)
    deallocate(rhcrit)
  end subroutine finalise_casim  

  subroutine shipway_microphysics(il, iu, jl, ju, kl, ku, dt,               &
       qv, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13,          &
       theta, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13,       &
       a14, a15, a16, a17,                                                  &
       exner, pressure, rho, w, tke, z_half, z_centre, dz,                  &
       dqv, dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11, dq12,  &
       dq13, dth, da1, da2, da3, da4, da5, da6, da7, da8, da9, da10, da11,  &
       da12, da13, da14, da15, da16, da17,                                  &
       is_in, ie_in, js_in, je_in, ks_in, ke_in,                            &
       l_tendency)

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
         , a15( kl:ku, il:iu, jl:ju ), a16( kl:ku, il:iu, jl:ju ), a17( kl:ku, il:iu, jl:ju )

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

    ! New optional l_tendency logical added...
    ! if true then a tendency is returned (i.e. units/s)
    ! if false then an increment is returned (i.e. units/timestep)
    logical, intent(in), optional :: l_tendency

    ! Local variables

    real(wp) :: p1, p2, p3

    integer :: n, k, i, j, nxny, imode

    integer :: NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

    integer, parameter :: CHUNKSIZE=10

    real(wp) :: rain_number
    real(wp) :: rain_mass
    real(wp) :: rain_m3
    real(wp) :: n0,lam,mu
    real(wp) :: m1,m2,m3

    !diagnostics
    real(wp) :: field, res
    character(100) :: units, name
    integer :: im
    character(1) :: char   

    ! Save parent model timestep for later use (e.g. in diagnostics)
    parent_dt=dt

    do i=is,ie
#if DEF_MODEL!=MODEL_MONC
      i_here=i
#endif

      do j=js,je
#if DEF_MODEL!=MODEL_MONC
        j_here=j
#endif

        tend=ZERO_REAL_WP
        call zero_procs(procs)
        aerosol_tend=ZERO_REAL_WP
        if (l_process) then
          call zero_procs(aerosol_procs)
        end if

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
             w(ks:ke,i,j), tke(ks:ke,i,j), qfields)
        
        !--------------------------------------------------
        ! Do the business...
        !--------------------------------------------------
        call microphysics_common(dt, ks, ke, qfields, dqfields, tend, procs, precip(i,j) &
             , aerophys, aerochem, aeroact                                         &
             , dustphys, dustchem, dustact                                         &
             , aeroice, dustliq                                                    &
             , aerofields, daerofields, aerosol_tend, aerosol_procs                &
             , rhcrit)

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

#if DEF_MODEL==MODEL_LEM_DIAG
        PUDDLE(J_here,1,I_here) = PUDDLE(J_here,1,I_here) +precip(i,j)*DTPUD
        precip_lem(j_here,i_here) = precip_lem(j_here,i_here) + precip(i,j)*3600.0
#endif
      end do
    end do

#if DEF_MODEL==MODEL_KiD
    call save_dg(sum(precip(:, :))/nxny, 'precip', i_dgtime)
    call save_dg(sum(precip(:, :))/nxny*3600.0, 'surface_precip_mmhr', i_dgtime)
#endif
  end subroutine shipway_microphysics

  subroutine microphysics_common(dt, kl, ku, qfields, dqfields, tend, procs, precip &
       , aerophys, aerochem, aeroact                                          &
       , dustphys, dustchem, dustact                                          &
       , aeroice, dustliq                                                     &
       , aerofields, daerofields, aerosol_tend, aerosol_procs                 &
       , rhcrit)

    real(wp), intent(in) :: dt  ! timestep from parent model
    integer, intent(in) :: kl, ku
    real(wp), intent(in) :: rhcrit(:)
    real(wp), intent(inout) :: qfields(:,:), dqfields(:,:), tend(:,:)
    type(process_rate), intent(inout) :: procs(:,:)
    real(wp), intent(out) :: precip

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


    real(wp), allocatable :: qfields_in(:,:)
    real(wp), allocatable :: qfields_mod(:,:)
    real(wp), allocatable :: aerofields_in(:,:)
    real(wp), allocatable :: aerofields_mod(:,:)
    real(wp) :: step_length
    real(wp) :: sed_length
    integer :: n, k, nsed, iq

    logical :: l_Twarm   ! temperature above freezing
    logical :: l_Tcold   ! temperature below freezing, i.e. .not. l_Twarm

    logical :: l_sigevap ! Is there significant evaporation of rain

    allocate(qfields_in(nz, nq))
    allocate(qfields_mod(nz, nq))
    qfields_in=qfields ! Initial values of q
    qfields_mod=qfields ! Modified initial values of q (may be modified if bad values sent in)

    allocate(dist_lambda(nz,nspecies))
    allocate(dist_mu(nz,nspecies))
    allocate(dist_n0(nz,nspecies))

    nsubsteps=max(1, ceiling(dt/max_step_length))
    step_length=dt/nsubsteps

    nsubseds=max(1, ceiling(step_length/max_sed_length))
    sed_length=step_length/nsubseds

    n_sub=1
    n_subsed=1

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
    call query_distributions(cloud_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
    call query_distributions(rain_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
    if (.not. l_warm_loc) then
      call query_distributions(ice_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
      call query_distributions(snow_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
      call query_distributions(graupel_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
    end if

    qfields=qfields_mod

    if (aerosol_option > 0) then
      allocate(aerofields_in(nz, naero))
      allocate(aerofields_mod(nz, naero))
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
      n_sub=n

      call preconditioner(qfields)

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

      do k=nz,1,-1
        k_here=k

        if (precondition(k)) then
          l_Twarm=TdegK(k) > 273.15
          l_Tcold=.not. l_Twarm
          !=================================
          !
          ! WARM MICROPHYSICAL PROCESSES....
          !
          !=================================
          !-------------------------------
          ! Do the autoconversion to rain
          !-------------------------------
          if (pswitch%l_praut) call raut(step_length, k, qfields, aerofields, procs, aerosol_procs)

          !-------------------------------
          ! Do the rain accreting cloud
          !-------------------------------
          if (pswitch%l_pracw) call racw(step_length, k, qfields, aerofields, procs, rain_params, aerosol_procs)

          !-------------------------------
          ! Do the rain self-collection
          !-------------------------------
          if (pswitch%l_pracr) call racr(step_length, k, qfields, procs)

          !-------------------------------
          ! Do the evaporation of rain
          !-------------------------------
          if (pswitch%l_prevp) call revp(step_length, k, qfields, aerofields, aerophys, aerochem, &
               aeroact, dustliq, procs, aerosol_procs, l_sigevap)
          !=================================
          !
          ! ICE MICROPHYSICAL PROCESSES....
          !
          !=================================
          if (.not. l_warm_loc) then
            if (l_Tcold) then
              !------------------------------------------------------
              ! Autoconverion to snow
              !------------------------------------------------------
              if (pswitch%l_psaut) call saut(step_length, k, qfields, aerofields, procs, aerosol_procs)

              !------------------------------------------------------
              ! Accretion processes
              !------------------------------------------------------
              ! Ice -> Cloud -> Ice
              if (pswitch%l_piacw) call iacc(step_length, k, ice_params, cloud_params, ice_params, qfields, &
                   procs, aeroact, dustliq, aerosol_procs)
              ! Snow -> Cloud -> Snow
              if (l_sg) then
                if (pswitch%l_psacw) call iacc(step_length, k, snow_params, cloud_params, snow_params, &
                     qfields, procs, aeroact, dustliq, aerosol_procs)
                ! Snow -> Ice -> Snow
                if (pswitch%l_psaci) call iacc(step_length, k, snow_params, ice_params, snow_params, qfields, &
                     procs, aeroact, dustliq, aerosol_procs)
                if (pswitch%l_praci .and. (.not. l_sigevap)) then
                  ! Rain mass dependence to decide if we produce snow or graupel
                  if (qfields(k,i_qr) > thresh_sig(i_qr) .and. l_g .and. l_raci_g) then
                    ! Rain -> Ice -> Graupel
                    call iacc(step_length, k, rain_params, ice_params, graupel_params, &
                         qfields, procs, aeroact, dustliq, aerosol_procs)
                  else
                    ! Rain -> Ice -> Snow
                    call iacc(step_length, k, rain_params, ice_params, snow_params, qfields, &
                         procs, aeroact, dustliq, aerosol_procs)
                  end if
                end if
                if (pswitch%l_psacr .and. (.not. l_sigevap)) then
                  ! Temperature dependence to decide if we produce snow or graupel
                  if (TdegK(k) < 268.15 .and. l_g) then
                    ! Snow -> Rain -> Graupel
                    call iacc(step_length, k, snow_params, rain_params, graupel_params, &
                         qfields, procs, aeroact, dustliq, aerosol_procs)
                  else
                    ! Snow -> Rain -> Snow
                    call iacc(step_length, k, snow_params, rain_params, snow_params, &
                         qfields, procs, aeroact, dustliq, aerosol_procs)
                  end if
                end if
                if (l_g) then
                  ! Graupel -> Cloud -> Graupel
                  if (pswitch%l_pgacw)call iacc(step_length, k, graupel_params, cloud_params, &
                       graupel_params, qfields, procs, aeroact, dustliq, aerosol_procs)
                  ! Graupel -> Rain -> Graupel
                  if (pswitch%l_pgacr .and. (.not. l_sigevap))call iacc(step_length, k,       &
                       graupel_params, rain_params, graupel_params, qfields,                    &
                       procs, aeroact, dustliq, aerosol_procs)
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
              if (l_g .and. .not. l_onlycollect) call graupel_embryos(step_length, k, qfields, &
                   procs, aerophys, aerochem, aeroact, aerosol_procs)

              !------------------------------------------------------
              ! Wet deposition/shedding (resulting from graupel
              ! accretion processes)
              ! NB This alters some of the accretion processes, so
              ! must come after their calculation and before they
              ! are used/rescaled elsewhere
              !------------------------------------------------------
              if (l_g .and. .not. l_onlycollect .and. (.not. l_sigevap)) call wetgrowth(step_length, k, qfields, &
                   procs, aerophys, aerochem, aeroact, aerosol_procs)

              !------------------------------------------------------
              ! Aggregation (self-collection)
              !------------------------------------------------------
              if (pswitch%l_psagg) call ice_aggregation(step_length, k, snow_params, qfields, procs)

              !------------------------------------------------------
              ! Break up (snow only)
              !------------------------------------------------------
              if (pswitch%l_psbrk) call ice_breakup(step_length, k, snow_params, qfields, procs)

              !------------------------------------------------------
              ! Ice multiplication (Hallet-mossop)
              !------------------------------------------------------
              if (pswitch%l_pihal) call hallet_mossop(step_length, k, qfields,     &
                   procs, aerophys, aerochem, aeroact, aerosol_procs)

              !------------------------------------------------------
              ! Homogeneous freezing (rain and cloud)
              !------------------------------------------------------
              if (pswitch%l_phomr .and. (.not. l_sigevap)) call ihom_rain(step_length, k, qfields, &
                   aeroact, dustliq, procs, aerosol_procs)
              if (pswitch%l_phomc) call ihom_droplets(step_length, k, qfields, aeroact, dustliq, procs, aerosol_procs)

              !------------------------------------------------------
              ! Condensation/immersion/contact nucleation of cloud ice
              !------------------------------------------------------
              if (pswitch%l_pinuc) call inuc(step_length, k, qfields, procs,     &
                   dustphys, dustchem, aeroact, dustliq, aerosol_procs)

              !------------------------------------------------------
              ! Deposition/sublimation of ice/snow/graupel
              !------------------------------------------------------
              if (pswitch%l_pidep) call idep(step_length, k, ice_params, qfields,  &
                   procs, dustact, aeroice, aerosol_procs)

              if (pswitch%l_psdep) call idep(step_length, k, snow_params, qfields, &
                   procs, dustact, aeroice, aerosol_procs)

              if (pswitch%l_pgdep) call idep(step_length, k, graupel_params, qfields, &
                   procs, dustact, aeroice, aerosol_procs)

              if (l_harrington .and. .not. l_onlycollect) call adjust_dep(dt, k, procs, qfields)

              !-----------------------------------------------------------
              ! Make sure we don't remove more than saturation allows
              !-----------------------------------------------------------
              if (l_idep) call ensure_saturated(k, step_length, qfields, procs, (/i_idep, i_sdep, i_gdep/))
              !-----------------------------------------------------------
              ! Make sure we don't put back more than saturation allows
              !-----------------------------------------------------------
              if (l_isub) call ensure_saturated(k, step_length, qfields, procs, (/i_isub, i_ssub, i_gsub/))

            else ! T is  above freezing
              !------------------------------------------------------
              ! Melting of ice/snow/graupel
              !------------------------------------------------------
              if (pswitch%l_psmlt .and. (.not. l_sigevap)) call melting(step_length, k, snow_params, qfields, &
                   procs, aeroice, dustact, aerosol_procs)
              if (pswitch%l_pgmlt .and. (.not. l_sigevap)) call melting(step_length, k, graupel_params, qfields, &
                   procs, aeroice, dustact, aerosol_procs)
              if (pswitch%l_pimlt .and. (.not. l_sigevap)) call melting(step_length, k, ice_params, qfields, &
                   procs, aeroice, dustact, aerosol_procs)

            end if

            !-----------------------------------------------------------
            ! Make sure we don't remove more than we have to start with
            !-----------------------------------------------------------
            if (l_pos1) call ensure_positive(k, step_length, qfields, procs, cloud_params, &
                 (/i_praut, i_pracw, i_iacw, i_sacw, i_gacw, i_homc/),  &
                 aeroprocs=aerosol_procs, iprocs_dependent=(/i_aaut, i_aacw/))

            if (.not. l_warm_loc) then
              if (l_pos2) call ensure_positive(k, step_length, qfields, procs, ice_params, &
                   (/i_raci, i_saci, i_gaci, i_saut, i_isub, i_imlt/), &
                   (/i_ihal, i_gshd, i_inuc, i_homc, i_iacw, i_idep/))
            end if

            if (l_pos3) call ensure_positive(k, step_length, qfields, procs, rain_params, &
                 (/i_prevp, i_sacr, i_gacr, i_homr/), &
                 (/i_praut, i_pracw, i_raci, i_gshd, i_smlt, i_gmlt/), &
                 aeroprocs=aerosol_procs, iprocs_dependent=(/i_arevp/))

            if (.not. l_warm_loc) then
              if (l_pos4) call ensure_positive(k, step_length, qfields, procs, snow_params, &
                   (/i_gacs, i_smlt, i_sacr, i_ssub /), &
                   (/i_sdep, i_sacw, i_saut, i_saci, i_raci, i_gshd, i_ihal/))
            end if

          else
            if (pswitch%l_praut .and. pswitch%l_pracw) then
              if (l_pos5) call ensure_positive(k, step_length, qfields, procs, cloud_params, (/i_praut, i_pracw/), &
                   aeroprocs=aerosol_procs, iprocs_dependent=(/i_aaut, i_aacw/))
            end if

            if (pswitch%l_prevp) then
              if (l_pos6) call ensure_positive(k, step_length, qfields, procs, rain_params, (/i_prevp/), (/i_praut, i_pracw/) &
                   , aerosol_procs, (/i_arevp/))
            end if

          end if

        end if
      end do
      !-------------------------------
      ! Collect terms we have so far
      !-------------------------------
      call sum_procs(step_length, procs, tend, (/i_praut, i_pracw, i_pracr, i_prevp/)     &
           , hydro_names, l_thermalexchange=.true., qfields=qfields , l_passive=l_passive)

      if (.not. l_warm_loc) then
        call sum_procs(step_length, procs, tend,      &
             (/i_idep, i_sdep, i_gdep, i_iacw, i_sacr, i_sacw, i_saci, i_raci,&
             i_gacw, i_gacr, i_gaci, i_gacs, i_ihal, i_gshd, i_sbrk,&
             i_saut, i_iagg, i_sagg, i_gagg, i_isub, i_ssub, i_gsub/),        &
             hydro_names, l_thermalexchange=.true., qfields=qfields,&
             l_passive=l_passive, i_thirdmoment=2)
        call sum_procs(step_length, procs, tend,      &
             (/i_inuc, i_imlt, i_smlt, i_gmlt, i_homr, i_homc/), hydro_names,     &
             l_thermalexchange=.true., qfields=qfields , l_passive=l_passive)
      end if

      call update_q(qfields_mod, qfields, tend)

      if (l_process) then
        do k=1,nz
          if (precondition(k)) then
            call ensure_positive_aerosol(k, step_length, aerofields, aerosol_procs,&
                 (/i_aaut, i_aacw, i_aevp, i_arevp, i_dnuc, i_dsub, i_dssub, i_dgsub, &
                 i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr,    &
                 i_dgacr, i_draci  /) )
          end if
        end do

        call sum_aprocs(step_length, aerosol_procs, aerosol_tend,      &
             (/i_aaut, i_aacw, i_aevp, i_arevp/), aero_names, aerofields)

        if (.not. l_warm) then
          call sum_aprocs(step_length, aerosol_procs, aerosol_tend,      &
               (/i_dnuc, i_dsub, i_dssub, i_dgsub, i_dimlt, i_dsmlt, i_dgmlt,     &
               i_diacw, i_dsacw, i_dgacw, i_dsacr, i_dgacr, i_draci /), aero_names, aerofields)
        end if

        call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)

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
      if (pswitch%l_pcond)then
        do k=nz,1,-1
#if DEF_MODEL==MODEL_KiD
          k_here=k
#endif
#if DEF_MODEL==MODEL_LEM_DIAG
          k_here=k
#endif
          call condevp(step_length, k, qfields, aerofields,     &
               procs, aerophys, aerochem, aeroact, dustphys, dustchem, dustliq, &
               aerosol_procs, rhcrit(k))
        end do

        !-------------------------------
        ! Collect terms we have so far
        !-------------------------------
        call sum_procs(step_length, procs, tend,      &
             (/i_cond/)     &
             , hydro_names, l_thermalexchange=.true., qfields=qfields     &
             , l_passive=l_passive)

        call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)
      end if

      !---------------------------------------------------------------
      ! Re-Determine (and possibly limit) size distribution
      !---------------------------------------------------------------
      call query_distributions(cloud_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
      call query_distributions(rain_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
      if (.not. l_warm_loc) then
        call query_distributions(ice_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
        call query_distributions(snow_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
        call query_distributions(graupel_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
      end if

      if (l_process) then
        call sum_aprocs(step_length, aerosol_procs, aerosol_tend,      &
             (/i_aact /), aero_names, aerofields)

        call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)

        !-------------------------------
        ! Re-Derive aerosol distribution
        ! parameters
        !-------------------------------
        call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact     &
             , dustphys, dustchem, dustact, aeroice, dustliq, icall=3)
      end if

#if DEF_MODEL==MODEL_LEM_DIAG
      ! Need to recode these with new process rates
      ! Save the process rates
      !call lem_procdgs(n, kl, ku, procs)
#endif
      !-------------------------------
      ! Do the sedimentation
      !-------------------------------

      if (l_sed) then
        do nsed=1,nsubseds
          n_subsed=nsed

          if (nsed > 1) then
            !-------------------------------
            ! Reset process rates if they
            ! are to be re-used
            !-------------------------------
            call zero_procs(procs)
            if (l_process) call zero_procs(aerosol_procs)
            !---------------------------------------------------------------
            ! Re-Determine (and possibly limit) size distribution
            !---------------------------------------------------------------
            call query_distributions(cloud_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
            call query_distributions(rain_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
            if (.not. l_warm_loc) then
              call query_distributions(ice_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
              call query_distributions(snow_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
              call query_distributions(graupel_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
            end if

            !-------------------------------
            ! Re-Derive aerosol distribution
            ! parameters
            !-------------------------------
            if (l_process) call examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact &
                 , dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
          end if

          if (pswitch%l_psedl) then
            call sedr(sed_length, qfields, aerofields, aeroact, dustliq, tend, cloud_params, &
                 procs, aerosol_procs, precip, l_process)

            call sum_procs(sed_length, procs, tend, (/i_psedl/), hydro_names, qfields=qfields)
          end if

          if (pswitch%l_psedr) then
            call sedr(sed_length, qfields, aerofields, aeroact, dustliq, tend, rain_params, &
                 procs, aerosol_procs, precip, l_process)

            call sum_procs(sed_length, procs, tend, (/i_psedr/), hydro_names, qfields=qfields)
          end if

          if (.not. l_warm_loc) then

            if (pswitch%l_psedi) then
              call sedr(sed_length, qfields, aerofields, aeroice, dustact, tend, ice_params, &
                   procs, aerosol_procs, precip, l_process)
            end if
            if (pswitch%l_pseds) then
              call sedr(sed_length, qfields, aerofields, aeroice, dustact, tend, snow_params, &
                   procs, aerosol_procs, precip, l_process)
            end if
            if (pswitch%l_psedg) then
              call sedr(sed_length, qfields, aerofields, aeroice, dustact, tend, graupel_params, &
                   procs, aerosol_procs, precip, l_process)
            end if

            call sum_procs(sed_length, procs, tend, (/i_psedi, i_pseds, i_psedg/), hydro_names, qfields=qfields)
          end if

          call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)

          if (l_process) then
            call sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_asedr, i_asedl/),&
                 aero_names, aerofields)

            call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)

            if (l_process .and. .not. l_warm) then
              call sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_dsedi, i_dseds, i_dsedg/),&
                   aero_names, aerofields)

              call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
            end if
          end if

        end do

      end if

      if (nsubsteps>1)then
        !-------------------------------
        ! Reset process rates if they
        ! are to be re-used
        !-------------------------------
        call zero_procs(procs)
        if (l_process) call zero_procs(aerosol_procs)
        !---------------------------------------------------------------
        ! Re-Determine (and possibly limit) size distribution
        !---------------------------------------------------------------
        call query_distributions(cloud_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsubseds+3)
        call query_distributions(rain_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsubseds+3)
        if (.not. l_warm_loc) then
          call query_distributions(ice_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsubseds+3)
          call query_distributions(snow_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsubseds+3)
          call query_distributions(graupel_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsubseds+3)
        end if
      end if
    end do

    !--------------------------------------------------
    ! Tidy up any small/negative numbers
    ! we may have generated.
    !--------------------------------------------------
    if (pswitch%l_tidy2) then
      do k=nz,1,-1

        k_here=k

        if (precondition(k)) then
          call qtidy(step_length, k, qfields, procs, aerofields, aeroact, dustact, &
               aeroice, dustliq , aerosol_procs, i_tidy2, i_atidy2, l_negonly=l_tidy_negonly)
        end if

      end do

      call sum_procs(step_length, procs, tend, (/i_tidy2/), hydro_names, qfields=qfields, l_passive=l_passive)

      call update_q(qfields_mod, qfields, tend)

      if (l_process) then
        call sum_aprocs(step_length, aerosol_procs, aerosol_tend, (/i_atidy2/), aero_names, aerofields )
        call update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.true.)
      end if
    end if

    !
    ! Add on initial adjustments that may have been made
    !

    if (l_tendency_loc) then! Convert back from cummulative value to tendency
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
        if (l_tendency_loc) then! Convert back from cummulative value to tendency
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
      deallocate(aerofields_in)
      deallocate(aerofields_mod)
    end if

    deallocate(qfields_in)
    deallocate(qfields_mod)
    deallocate(dist_lambda)
    deallocate(dist_mu)
    deallocate(dist_n0)
  end subroutine microphysics_common

  subroutine update_q(qfields_in, qfields, tend, l_aerosol, l_fixneg)
    real(wp), intent(in) :: qfields_in(:,:)
    real(wp), intent(inout) :: qfields(:,:)
    real(wp), intent(in) :: tend(:,:)
    logical, intent(in), optional :: l_aerosol ! flag to indicate updating of aerosol
    logical, intent(in), optional :: l_fixneg  ! Flag to use cludge to bypass negative/zero numbers
    integer :: k, iqx
    logical :: l_fix

    l_fix=.false.
    if (present(l_fixneg)) l_fix=l_fixneg

    do iqx=1, ubound(tend,2)
      do k=lbound(tend,1), ubound(tend,1)
        qfields(k,iqx)=qfields_in(k,iqx)+tend(k,iqx)
        !quick lem fixes  - this code should never be used ?
        if (.not. present(l_aerosol) .and. l_fix) then
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
        end if
      end do
    end do
  end subroutine update_q

  subroutine update_tend(dt, qfields, dqfields, tend)
    real(wp) :: dt
    real(wp), intent(in) :: qfields(:,:)
    real(wp), intent(in) :: dqfields(:,:)
    real(wp), intent(inout) :: tend(:,:)

    real(wp), allocatable :: tmp(:,:)

    integer :: k, iqx

    ! Using tmp for premultiplication seems to optimize
    ! the loop much better for some compilers
    allocate(tmp(ubound(dqfields,1), lbound(dqfields,2):ubound(dqfields,2)))
    tmp=dt*dqfields

    do iqx=1, ubound(tend,2)
      do k=lbound(tend,1), ubound(tend,1)
        tend(k,iqx) = tend(k,iqx) + tmp(k,iqx)
      end do
    end do
    deallocate(tmp)
  end subroutine update_tend
end module micro_main
