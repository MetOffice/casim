MODULE micro_main

USE mphys_die, ONLY: throw_mphys_error
USE variable_precision, ONLY: wp
USE mphys_parameters, ONLY: nz, nq, naero, nprocs, naeroprocs               &
     , rain_params, cloud_params, ice_params, snow_params, graupel_params   &
     , nspecies, parent_dt
USE process_routines, ONLY: process_rate, zero_procs,                          &
     allocate_procs, deallocate_procs                                       &
     , i_cond, i_praut, i_pracw, i_pracr, i_prevp, i_psedr, i_psedl         &
     , i_aact, i_aaut, i_aacw, i_aevp, i_asedr, i_asedl, i_arevp            &
     , i_tidy, i_tidy2, i_atidy, i_atidy2                                   &
     , i_inuc, i_idep, i_dnuc, i_dsub, i_saut, i_iacw, i_sacw, i_pseds      &
     , i_sdep, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs       &
     , i_gdep, i_psedg, i_iagg, i_sagg, i_gagg, i_gshd, i_ihal              &
     , i_smlt, i_gmlt, i_psedi, i_homr, i_homc, i_imlt                      &
     , i_isub, i_ssub, i_gsub, i_sbrk                                       &
     , i_dssub, i_dgsub, i_dsedi, i_dseds, i_dsedg                          &
     , i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw, i_dgacw                 &
     , i_dsacr, i_dgacr, i_draci, process_name

USE sum_process, ONLY: sum_procs, sum_aprocs
USE aerosol_routines, ONLY:                                                    &
       examine_aerosol, aerosol_phys, aerosol_chem, aerosol_active          &
     , allocate_aerosol, deallocate_aerosol, MNtoRm

USE mphys_switches, ONLY: hydro_complexity, aero_complexity,                &
     i_qv, i_ql, i_nl, i_qr, i_nr, i_m3r, i_th                              &
     , i_qi, i_qs, i_qg, i_ni, i_ns, i_ng, i_m3s, i_m3g                     &
     , i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am4, i_am5               &
     , i_am6, i_an6, i_am7, i_am8 , i_am9, i_am10, i_an10, i_an11, i_an12   &
     , hydro_names, aero_names,  aerosol_option                             &
     , l_warm, l_inuc, l_sed, l_condensation                                &
     , l_iaut, l_imelt, l_iacw, l_idep, aero_index, l_sedl                  &
     , nq_l, nq_r, nq_i, nq_s, nq_g                                         &
     , l_sg, l_g, l_process                                                 &
     , l_halletmossop, max_sed_length, max_step_length, l_harrington        &
     , l_passive, ntotala, ntotalq                                          &
     , active_number, isol, iinsol, l_raci_g, l_onlycollect, l_pracr        &
     , pswitch, l_isub, l_pos1, l_pos2, l_pos3, l_pos4, l_pos5, l_pos6      &
     , i_hstart, nsubsteps, nsubseds

USE distributions, ONLY: query_distributions, dist_lambda, dist_mu, dist_n0
USE passive_fields, ONLY: set_passive_fields, TdegK
USE autoconversion, ONLY: raut
USE evaporation, ONLY: revp
USE condensation, ONLY: condevp
USE accretion, ONLY: racw
USE aggregation, ONLY: racr, ice_aggregation
USE sedimentation, ONLY: sedr
USE ice_nucleation, ONLY: inuc
USE ice_deposition, ONLY: idep
USE ice_accretion, ONLY: iacc
USE breakup, ONLY: ice_breakup
USE snow_autoconversion, ONLY: saut
USE ice_multiplication, ONLY: hallet_mossop
USE graupel_wetgrowth, ONLY: wetgrowth
USE graupel_embryo, ONLY: graupel_embryos
USE ice_melting, ONLY: melting
USE homogeneous, ONLY: ihom_rain, ihom_droplets
USE adjust_deposition, ONLY: adjust_dep

USE lookup, ONLY: get_slope_generic, moment
USE thresholds, ONLY: qr_small, nr_small, m3r_small, thresh_sig, thresh_tidy
USE mphys_constants, ONLY:  fixed_aerosol_number, fixed_aerosol_rm          &
     , fixed_aerosol_sigma, fixed_aerosol_density

USE mphys_tidy, ONLY: qtidy, ensure_positive, ensure_saturated, tidy_qin, tidy_ain &
       , ensure_positive_aerosol
USE preconditioning, ONLY: precondition, preconditioner


#if DEF_MODEL==MODEL_LEM
USE diaghelp_lem, ONLY: i_here, j_here, k_here, koff,   &
       n_sub, n_subsed
#endif

#if DEF_MODEL==MODEL_KiD
  ! Kid modules
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, j_here, k_here,   &
       n_sub, n_subsed
USE runtime, ONLY: time
USE parameters, ONLY: nx
#elif DEF_MODEL==MODEL_LEM_DIAG
USE initialize, ONLY: DTPUD

USE diags_lem, ONLY: lem_procdgs
USE com_prametr, ONLY: ilatentdgp, nprc_prametr => nprc, iforscalp
USE com_dgstore, ONLY: procrate, fallq
USE com_dgchars, ONLY: procchar
USE com_rain, ONLY: puddle
USE com_classy, ONLY: precip_lem => precip
#elif DEF_MODEL==MODEL_UM
USE diaghelp_um, ONLY: i_here, j_here, k_here, l_debug_um, debug_i, debug_j, debug_k, debug_pe, debug_i2, debug_j2, debug_k2, &
       n_sub, n_subsed
USE timestep_mod, ONLY: timestep_number
USE mphys_casim_diagnostics, ONLY: SurfacePrecip, ProcessRates
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here, k_here, n_sub, n_subsed, mype
#endif

IMPLICIT NONE

  ! Debug switches
LOGICAL :: l_debug_main=.FALSE., l_debug_common=.FALSE.

LOGICAL :: l_tendency_loc
LOGICAL :: l_warm_loc

CONTAINS

SUBROUTINE shipway_microphysics(il, iu, jl, ju, kl, ku, dt,               &
     qv, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13,          &
     theta, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13,       &
     a14, a15, a16, a17,                                                  &
     exner, pressure, rho, w, tke, z_half, z_centre, dz,                  &
     dqv, dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11, dq12,  &
     dq13, dth, da1, da2, da3, da4, da5, da6, da7, da8, da9, da10, da11,  &
     da12, da13, da14, da15, da16, da17,                                  &
     is_in, ie_in, js_in, je_in, ks_in, ke_in,                            &
     l_tendency, rhcrit_in                                                &
     )

INTEGER, INTENT(IN) :: il, iu ! upper and lower i levels
INTEGER, INTENT(IN) :: jl, ju ! upper and lower j levels
INTEGER, INTENT(IN) :: kl, ku ! upper and lower k levels

REAL(wp), INTENT(IN) :: dt    ! parent model timestep (s)

    ! hydro fields in... 1-5 should be warm rain, 6+ are ice
    ! see mphys_casim for details of what is passed in
REAL(wp), INTENT(IN) :: q1( kl:ku, il:iu, jl:ju ), q2( kl:ku, il:iu, jl:ju )   &
       , q3( kl:ku, il:iu, jl:ju ), q4( kl:ku, il:iu, jl:ju ), q5( kl:ku, il:iu, jl:ju ) &
       , q6( kl:ku, il:iu, jl:ju ), q7( kl:ku, il:iu, jl:ju ), q8( kl:ku, il:iu, jl:ju ) &
       , q9( kl:ku, il:iu, jl:ju ), q10( kl:ku, il:iu, jl:ju ), q11( kl:ku, il:iu, jl:ju ) &
       , q12( kl:ku, il:iu, jl:ju ), q13( kl:ku, il:iu, jl:ju )


REAL(wp), INTENT(IN) :: qv( kl:ku, il:iu, jl:ju )
REAL(wp), INTENT(IN) :: theta( kl:ku, il:iu, jl:ju )
REAL(wp), INTENT(IN) :: exner( kl:ku, il:iu, jl:ju )
REAL(wp), INTENT(IN) :: pressure( kl:ku, il:iu, jl:ju )
REAL(wp), INTENT(IN) :: rho( kl:ku, il:iu, jl:ju )
REAL(wp), INTENT(IN) :: w( kl:ku, il:iu, jl:ju )
REAL(wp), INTENT(IN) :: tke( kl:ku, il:iu, jl:ju )
REAL(wp), INTENT(IN) :: dz( kl:ku, il:iu, jl:ju )
REAL(wp), INTENT(IN) :: z_half( kl-1:ku, il:iu, jl:ju )
REAL(wp), INTENT(IN) :: z_centre( kl:ku, il:iu, jl:ju )

    ! Aerosol fields in
REAL(wp), INTENT(IN) :: a1( kl:ku, il:iu, jl:ju ), a2( kl:ku, il:iu, jl:ju )   &
       , a3( kl:ku, il:iu, jl:ju ), a4( kl:ku, il:iu, jl:ju ), a5( kl:ku, il:iu, jl:ju ) &
       , a6( kl:ku, il:iu, jl:ju ), a7( kl:ku, il:iu, jl:ju ), a8( kl:ku, il:iu, jl:ju ) &
       , a9( kl:ku, il:iu, jl:ju ), a10( kl:ku, il:iu, jl:ju ), a11(kl:ku, il:iu, jl:ju ) &
       , a12( kl:ku, il:iu, jl:ju ), a13( kl:ku, il:iu, jl:ju ), a14( kl:ku, il:iu, jl:ju ) &
       , a15( kl:ku, il:iu, jl:ju ), a16( kl:ku, il:iu, jl:ju ), a17( kl:ku, il:iu, jl:ju )

    ! hydro tendencies in:  from parent model forcing i.e. advection
    ! hydro tendencies out: from microphysics only...
REAL(wp), INTENT(INOUT) :: dq1( kl:ku, il:iu, jl:ju ), dq2( kl:ku, il:iu, jl:ju ) &
       , dq3( kl:ku, il:iu, jl:ju ), dq4( kl:ku, il:iu, jl:ju ), dq5( kl:ku, il:iu, jl:ju ) &
       , dq6( kl:ku, il:iu, jl:ju ), dq7( kl:ku, il:iu, jl:ju ), dq8( kl:ku, il:iu, jl:ju ) &
       , dq9( kl:ku, il:iu, jl:ju ), dq10( kl:ku, il:iu, jl:ju ), dq11( kl:ku, il:iu, jl:ju ) &
       , dq12( kl:ku, il:iu, jl:ju ), dq13( kl:ku, il:iu, jl:ju )

    ! qv/theta tendencies in:  from parent model forcing i.e. advection
    ! qv/theta tendencies out: from microphysics only
REAL(wp), INTENT(INOUT) :: dqv( kl:ku, il:iu, jl:ju ), dth( kl:ku, il:iu, jl:ju )

    ! aerosol tendencies in:  from parent model forcing i.e. advection
    ! aerosol tendencies out: from microphysics only
REAL(wp), INTENT(INOUT) :: da1( kl:ku, il:iu, jl:ju ), da2( kl:ku, il:iu, jl:ju ) &
       , da3( kl:ku, il:iu, jl:ju ), da4( kl:ku, il:iu, jl:ju ), da5( kl:ku, il:iu, jl:ju ) &
       , da6( kl:ku, il:iu, jl:ju ), da7( kl:ku, il:iu, jl:ju ), da8( kl:ku, il:iu, jl:ju ) &
       , da9( kl:ku, il:iu, jl:ju ), da10( kl:ku, il:iu, jl:ju ), da11( kl:ku, il:iu, jl:ju ) &
       , da12( kl:ku, il:iu, jl:ju ), da13( kl:ku, il:iu, jl:ju ), da14( kl:ku, il:iu, jl:ju ) &
       , da15( kl:ku, il:iu, jl:ju ), da16( kl:ku, il:iu, jl:ju ), da17( kl:ku, il:iu, jl:ju )

INTEGER, INTENT(IN), OPTIONAL :: is_in, ie_in ! upper and lower i levels which are to be used
INTEGER, INTENT(IN), OPTIONAL :: js_in, je_in ! upper and lower j levels
INTEGER, INTENT(IN), OPTIONAL :: ks_in, ke_in ! upper and lower k levels

    ! New optional l_tendency logical added...
    ! if true then a tendency is returned (i.e. units/s)
    ! if false then an increment is returned (i.e. units/timestep)
LOGICAL, INTENT(IN), OPTIONAL :: l_tendency
REAL(wp), INTENT(IN), OPTIONAL :: rhcrit_in(kl:ku)

    ! Local variables
INTEGER :: is, ie ! upper and lower i levels which are to be used
INTEGER :: js, je ! upper and lower j levels
INTEGER :: ks, ke ! upper and lower k levels

REAL(wp) :: rhcrit(kl:ku)
REAL(wp), ALLOCATABLE :: dqfields(:,:), qfields(:,:), tend(:,:)
REAL(wp), ALLOCATABLE :: daerofields(:,:), aerofields(:,:), aerosol_tend(:,:)

TYPE(process_rate), ALLOCATABLE :: procs(:,:)
TYPE(process_rate), ALLOCATABLE :: aerosol_procs(:,:)

TYPE(aerosol_active), ALLOCATABLE :: aeroact(:)
TYPE(aerosol_phys), ALLOCATABLE   :: aerophys(:)
TYPE(aerosol_chem), ALLOCATABLE   :: aerochem(:)

TYPE(aerosol_active), ALLOCATABLE :: dustact(:)
TYPE(aerosol_phys), ALLOCATABLE   :: dustphys(:)
TYPE(aerosol_chem), ALLOCATABLE   :: dustchem(:)

TYPE(aerosol_active), ALLOCATABLE :: aeroice(:)  ! Soluble aerosol in ice
TYPE(aerosol_active), ALLOCATABLE :: dustliq(:)! Insoluble aerosol in liquid

REAL(wp), ALLOCATABLE, SAVE :: precip(:,:) ! diagnostic for surface precip rate

REAL(wp) :: p1, p2, p3

INTEGER :: n, k, i, j, nxny, imode

INTEGER :: NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

INTEGER, PARAMETER :: CHUNKSIZE=10

REAL(wp) :: rain_number
REAL(wp) :: rain_mass
REAL(wp) :: rain_m3
REAL(wp) :: n0,lam,mu
REAL(wp) :: m1,m2,m3

    !diagnostics
REAL(wp) :: field, res
CHARACTER(100) :: units, name
INTEGER :: im
CHARACTER(1) :: char

IF (l_debug_main)PRINT*,'DEBUG-M1'

    ! Save parent model timestep for later use (e.g. in diagnostics)
parent_dt=dt

!    if (timestep_number*dt < 3600)then
!      l_warm_loc=.true. ! Force warm in spinup
!    else
l_warm_loc=l_warm ! Original setting
!    end if

      ! Set grid extents to operate on
IF (PRESENT(is_in))is=is_in
IF (PRESENT(ie_in))ie=ie_in
IF (PRESENT(js_in))js=js_in
IF (PRESENT(je_in))je=je_in
IF (PRESENT(ks_in))ks=ks_in
IF (PRESENT(ke_in))ke=ke_in

      ! if not passed in, then default to full grid
IF (.NOT. PRESENT(is_in))is=il
IF (.NOT. PRESENT(ie_in))ie=iu
IF (.NOT. PRESENT(js_in))js=jl
IF (.NOT. PRESENT(je_in))je=ju
IF (.NOT. PRESENT(ks_in))ks=kl
IF (.NOT. PRESENT(ke_in))ke=ku

l_tendency_loc=.TRUE.
IF (PRESENT(l_tendency))l_tendency_loc=l_tendency
IF (PRESENT(rhcrit_in)) THEN
  rhcrit(:)=rhcrit_in(:)
ELSE
  rhcrit(:)=1.0
END IF

nq = SUM(hydro_complexity%nmoments) + 2 ! also includes vapour and theta
nz = ke - ks + 1
nprocs = hydro_complexity%nprocesses
nxny=(ie-is+1)*(je-js+1)

ALLOCATE(precondition(nz))
precondition = .TRUE. ! Assume all points need to be considered
ALLOCATE(qfields(nz, nq))
ALLOCATE(dqfields(nz, nq))
qfields  = 0.0
dqfields = 0.0
ALLOCATE(procs(nz, nprocs))
ALLOCATE(tend(nz, nq))

    ! Allocate aerosol storage
IF (aerosol_option > 0) THEN
  naero = ntotala
  naeroprocs=aero_complexity%nprocesses
  ALLOCATE(aerofields(nz, naero))
  ALLOCATE(daerofields(nz, naero))
  aerofields  = 0.0
  daerofields = 0.0
  ALLOCATE(aerosol_procs(nz, naeroprocs))
  ALLOCATE(aerosol_tend(nz, naero))
ELSE
      ! Dummy arrays required
  ALLOCATE(aerofields(1,1))
  ALLOCATE(daerofields(1,1))
  ALLOCATE(aerosol_procs(1,1))
  ALLOCATE(aerosol_tend(1,1))
END IF

ALLOCATE(aerophys(nz))
ALLOCATE(aerochem(nz))
ALLOCATE(aeroact(nz))
CALL allocate_aerosol(aerophys, aerochem, aero_index%nccn)
ALLOCATE(dustphys(nz))
ALLOCATE(dustchem(nz))
ALLOCATE(dustact(nz))
CALL allocate_aerosol(dustphys, dustchem, aero_index%nin)

ALLOCATE(aeroice(nz))
ALLOCATE(dustliq(nz))

IF (l_debug_main)PRINT*,'DEBUG-M2'

    ! Temporary initialization of chem and sigma
DO k =1,SIZE(aerophys)
  aerophys(k)%sigma(:)=fixed_aerosol_sigma
  aerophys(k)%rpart(:)=0.0
  aerochem(k)%vantHoff(:)=3.0
  aerochem(k)%massMole(:)=132.0e-3
  aerochem(k)%density(:)=fixed_aerosol_density
  aerochem(k)%epsv(:)=1.0
  aerochem(k)%beta(:)=1.0
END DO
DO k =1,SIZE(dustphys)
  dustphys(k)%sigma(:)=fixed_aerosol_sigma
  dustphys(k)%rpart(:)=0.0
  dustchem(k)%vantHoff(:)=3.0
  dustchem(k)%massMole(:)=132.0e-3
  dustchem(k)%density(:)=fixed_aerosol_density
  dustchem(k)%epsv(:)=1.0
  dustchem(k)%beta(:)=1.0
END DO

    !allocate space for the process rates
CALL allocate_procs(procs, nz, nprocs, ntotalq)
IF (l_process)CALL allocate_procs(aerosol_procs, nz, naeroprocs, ntotala)

    ! allocate diagnostics
ALLOCATE(precip(il:iu,jl:ju))

IF (l_debug_main)PRINT*,'DEBUG-M3'

DO i=is,ie
#if DEF_MODEL!=MODEL_MONC
  i_here=i
#endif

  DO j=js,je
#if DEF_MODEL!=MODEL_MONC
  j_here=j
#endif



    tend = 0.0
    CALL zero_procs(procs)
    aerosol_tend = 0.0
    IF (l_process) THEN
      CALL zero_procs(aerosol_procs)
    END IF

        ! Set the qfields
    qfields(:, i_qv) = qv(ks:ke,i,j)
    qfields(:, i_th) = theta(ks:ke,i,j)

    IF (nq_l > 0)qfields(:,i_ql) = q1(ks:ke,i,j)
    IF (nq_r > 0)qfields(:,i_qr) = q2(ks:ke,i,j)
    IF (nq_l > 1)qfields(:,i_nl) = q3(ks:ke,i,j)
    IF (nq_r > 1)qfields(:,i_nr) = q4(ks:ke,i,j)
    IF (nq_r > 2)qfields(:,i_m3r) = q5(ks:ke,i,j)
    IF (nq_i > 0)qfields(:,i_qi) = q6(ks:ke,i,j)
    IF (nq_s > 0)qfields(:,i_qs) = q7(ks:ke,i,j)
    IF (nq_g > 0)qfields(:,i_qg) = q8(ks:ke,i,j)
    IF (nq_i > 1)qfields(:,i_ni) = q9(ks:ke,i,j)
    IF (nq_s > 1)qfields(:,i_ns) = q10(ks:ke,i,j)
    IF (nq_g > 1)qfields(:,i_ng) = q11(ks:ke,i,j)
    IF (nq_s > 2)qfields(:,i_m3s) = q12(ks:ke,i,j)
    IF (nq_g > 2)qfields(:,i_m3g) = q13(ks:ke,i,j)
    dqfields(:, i_qv) = dqv(ks:ke,i,j)
    dqfields(:, i_th) = dth(ks:ke,i,j)

    IF (nq_l > 0)dqfields(:,i_ql) = dq1(ks:ke,i,j)
    IF (nq_r > 0)dqfields(:,i_qr) = dq2(ks:ke,i,j)
    IF (nq_l > 1)dqfields(:,i_nl) = dq3(ks:ke,i,j)
    IF (nq_r > 1)dqfields(:,i_nr) = dq4(ks:ke,i,j)
    IF (nq_r > 2)dqfields(:,i_m3r) = dq5(ks:ke,i,j)
    IF (nq_i > 0)dqfields(:,i_qi) = dq6(ks:ke,i,j)
    IF (nq_s > 0)dqfields(:,i_qs) = dq7(ks:ke,i,j)
    IF (nq_g > 0)dqfields(:,i_qg) = dq8(ks:ke,i,j)
    IF (nq_i > 1)dqfields(:,i_ni) = dq9(ks:ke,i,j)
    IF (nq_s > 1)dqfields(:,i_ns) = dq10(ks:ke,i,j)
    IF (nq_g > 1)dqfields(:,i_ng) = dq11(ks:ke,i,j)
    IF (nq_s > 2)dqfields(:,i_m3s) = dq12(ks:ke,i,j)
    IF (nq_g > 2)dqfields(:,i_m3g) = dq13(ks:ke,i,j)

    IF (aerosol_option > 0) THEN
      IF (i_am1 >0)aerofields(:, i_am1) = a1(ks:ke,i,j)
      IF (i_an1 >0)aerofields(:, i_an1) = a2(ks:ke,i,j)
      IF (i_am2 >0)aerofields(:, i_am2) = a3(ks:ke,i,j)
      IF (i_an2 >0)aerofields(:, i_an2) = a4(ks:ke,i,j)
      IF (i_am3 >0)aerofields(:, i_am3) = a5(ks:ke,i,j)
      IF (i_an3 >0)aerofields(:, i_an3) = a6(ks:ke,i,j)
      IF (i_am4 >0)aerofields(:, i_am4) = a7(ks:ke,i,j)
      IF (i_am5 >0)aerofields(:, i_am5) = a8(ks:ke,i,j)
      IF (i_am6 >0)aerofields(:, i_am6) = a9(ks:ke,i,j)
      IF (i_an6 >0)aerofields(:, i_an6) = a10(ks:ke,i,j)
      IF (i_am7 >0)aerofields(:, i_am7) = a11(ks:ke,i,j)
      IF (i_am8 >0)aerofields(:, i_am8) = a12(ks:ke,i,j)
      IF (i_am9 >0)aerofields(:, i_am9) = a13(ks:ke,i,j)
      IF (i_am10 >0)aerofields(:, i_am10) = a14(ks:ke,i,j)
      IF (i_an10 >0)aerofields(:, i_an10) = a15(ks:ke,i,j)
      IF (i_an11 >0)aerofields(:, i_an11) = a16(ks:ke,i,j)
      IF (i_an12 >0)aerofields(:, i_an12) = a17(ks:ke,i,j)

      IF (i_am1 >0)daerofields(:, i_am1) = da1(ks:ke,i,j)
      IF (i_an1 >0)daerofields(:, i_an1) = da2(ks:ke,i,j)
      IF (i_am2 >0)daerofields(:, i_am2) = da3(ks:ke,i,j)
      IF (i_an2 >0)daerofields(:, i_an2) = da4(ks:ke,i,j)
      IF (i_am3 >0)daerofields(:, i_am3) = da5(ks:ke,i,j)
      IF (i_an3 >0)daerofields(:, i_an3) = da6(ks:ke,i,j)
      IF (i_am4 >0)daerofields(:, i_am4) = da7(ks:ke,i,j)
      IF (i_am5 >0)daerofields(:, i_am5) = da8(ks:ke,i,j)
      IF (i_am6 >0)daerofields(:, i_am6) = da9(ks:ke,i,j)
      IF (i_an6 >0)daerofields(:, i_an6) = da10(ks:ke,i,j)
      IF (i_am7 >0)daerofields(:, i_am7) = da11(ks:ke,i,j)
      IF (i_am8 >0)daerofields(:, i_am8) = da12(ks:ke,i,j)
      IF (i_am9 >0)daerofields(:, i_am9) = da13(ks:ke,i,j)
      IF (i_am10 >0)daerofields(:, i_am10) = da14(ks:ke,i,j)
      IF (i_an10 >0)daerofields(:, i_an10) = da15(ks:ke,i,j)
      IF (i_an11 >0)daerofields(:, i_an11) = da16(ks:ke,i,j)
      IF (i_an12 >0)daerofields(:, i_an12) = da17(ks:ke,i,j)
    END IF

        !--------------------------------------------------
        ! set fields which will not be modified
        !--------------------------------------------------
    CALL set_passive_fields(ks, ke, dt, rho(ks:ke,i,j),    &
         pressure(ks:ke,i,j), exner(ks:ke,i,j),            &
         z_half(ks-1:ke,i,j), z_centre(ks:ke,i,j), dz(ks:ke,i,j),     &
         w(ks:ke,i,j), tke(ks:ke,i,j), qfields)

        !--------------------------------------------------
        ! Do the business...
        !--------------------------------------------------
    CALL microphysics_common(dt, ks, ke, qfields, dqfields, tend, procs, precip(i,j) &
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
    dqv(ks:ke,i,j) = tend(:,i_qv)
    dth(ks:ke,i,j) = tend(:,i_th)
    dq1(ks:ke,i,j) = tend(:,i_ql)
    dq2(ks:ke,i,j) = tend(:,i_qr)

    IF (cloud_params%l_2m)dq3(ks:ke,i,j) = tend(:,i_nl)
    IF (rain_params%l_2m)dq4(ks:ke,i,j) = tend(:,i_nr)
    IF (rain_params%l_3m)dq5(ks:ke,i,j) = tend(:,i_m3r)

    IF (.NOT. l_warm) THEN
      IF (ice_params%l_1m)    dq6(ks:ke,i,j)  = tend(:,i_qi)
      IF (snow_params%l_1m)   dq7(ks:ke,i,j)  = tend(:,i_qs)
      IF (graupel_params%l_1m)dq8(ks:ke,i,j)  = tend(:,i_qg)
      IF (ice_params%l_2m)    dq9(ks:ke,i,j)  = tend(:,i_ni)
      IF (snow_params%l_2m)   dq10(ks:ke,i,j) = tend(:,i_ns)
      IF (graupel_params%l_2m)dq11(ks:ke,i,j) = tend(:,i_ng)
      IF (snow_params%l_3m)   dq12(ks:ke,i,j) = tend(:,i_m3s)
      IF (graupel_params%l_3m)dq13(ks:ke,i,j) = tend(:,i_m3g)
    END IF


    IF (l_process) THEN
      IF (i_am1 >0)da1(ks:ke,i,j) = aerosol_tend(:,i_am1)
      IF (i_an1 >0)da2(ks:ke,i,j) = aerosol_tend(:,i_an1)
      IF (i_am2 >0)da3(ks:ke,i,j) = aerosol_tend(:,i_am2)
      IF (i_an2 >0)da4(ks:ke,i,j) = aerosol_tend(:,i_an2)
      IF (i_am3 >0)da5(ks:ke,i,j) = aerosol_tend(:,i_am3)
      IF (i_an3 >0)da6(ks:ke,i,j) = aerosol_tend(:,i_an3)
      IF (i_am4 >0)da7(ks:ke,i,j) = aerosol_tend(:,i_am4)
      IF (i_am5 >0)da8(ks:ke,i,j) = aerosol_tend(:,i_am5)
      IF (i_am6 >0)da9(ks:ke,i,j) = aerosol_tend(:,i_am6)
      IF (i_an6 >0)da10(ks:ke,i,j) = aerosol_tend(:,i_an6)
      IF (i_am7 >0)da11(ks:ke,i,j) = aerosol_tend(:,i_am7)
      IF (i_am8 >0)da12(ks:ke,i,j) = aerosol_tend(:,i_am8)
      IF (i_am9 >0)da13(ks:ke,i,j) = aerosol_tend(:,i_am9)
      IF (i_am10 >0)da14(ks:ke,i,j) = aerosol_tend(:,i_am10)
      IF (i_an10 >0)da15(ks:ke,i,j) = aerosol_tend(:,i_an10)
      IF (i_an11 >0)da16(ks:ke,i,j) = aerosol_tend(:,i_an11)
      IF (i_an12 >0)da17(ks:ke,i,j) = aerosol_tend(:,i_an12)
    ELSE
      da1(ks:ke,i,j) = 0.0
      da2(ks:ke,i,j) = 0.0
      da3(ks:ke,i,j) = 0.0
      da4(ks:ke,i,j) = 0.0
      da5(ks:ke,i,j) = 0.0
      da6(ks:ke,i,j) = 0.0
      da7(ks:ke,i,j) = 0.0
      da9(ks:ke,i,j) = 0.0
      da10(ks:ke,i,j) = 0.0
      da11(ks:ke,i,j) = 0.0
      da12(ks:ke,i,j) = 0.0
      da13(ks:ke,i,j) = 0.0
      da14(ks:ke,i,j) = 0.0
      da15(ks:ke,i,j) = 0.0
      da16(ks:ke,i,j) = 0.0
      da17(ks:ke,i,j) = 0.0
    END IF

        ! !--------------------------------------------------
        ! ! Some diagnostic stuff
        ! !--------------------------------------------------

! #if DEF_MODEL==MODEL_KiD || DEF_MODEL==MODEL_LEM
!          do k=1,nz
!            rain_mass = qfields(k, i_qr)
!            if (rain_mass > qr_small)then
!              p1=rain_params%p1
!              p3=rain_params%p2
!              p3=rain_params%p3
!              if (rain_params%l_2m)then
!                m2 = qfields(k, i_nr)
!              end if
!              if (rain_params%l_3m)then
!                m3 = qfields(k, i_m3r)
!              end if

!              call get_slope_generic(k, rain_params, n0, lam, mu, rain_mass, m2, m3)

! !#if DEF_MODEL==MODEL_KiD
!              if (nx==1)then
! !               call save_dg(k, mu, 'mu', i_dgtime)
! !               call save_dg(k, lam, 'lam', i_dgtime)
! !               call save_dg(k, n0, 'n0', i_dgtime)
!               call save_dg(k, rho(k,i,j), 'rho', i_dgtime)
!                do im=0,7
!                  write(char,'(i1)') im
!                  call save_dg(k, moment(n0, lam, mu, real(im,wp)), 'M'//char,&
!                     i_dgtime)
!                end do
!              else
! !               call save_dg(k, i, mu, 'mu', i_dgtime)
! !               call save_dg(k, i, lam, 'lam', i_dgtime)
! !               call save_dg(k, i, n0, 'n0', i_dgtime)
!               call save_dg(k, i, rho(k,i,j), 'rho', i_dgtime)
!                do im=0,7
!                  write(char,'(i1)') im
!                  call save_dg(k, i, moment(n0, lam, mu, real(im,wp)), 'M'//char, &
!                     i_dgtime)
!                end do
!              endif
! ! #elif DEF_MODEL==MODEL_LEM
! !             do im=0,7
! !               dgfields(j_here, k+koff, im+1, i_here) = moment(n0, lam, mu, real(im,wp))
! !             end do
! ! #endif
!            end if
!          end do
! #endif

#if DEF_MODEL==MODEL_LEM_DIAG
    PUDDLE(J_here,1,I_here) = PUDDLE(J_here,1,I_here) +precip(i,j)*DTPUD
    precip_lem(j_here,i_here) = precip_lem(j_here,i_here) + precip(i,j)*3600.0
#endif

  END DO
END DO

IF (l_debug_main)PRINT*,'DEBUG-M4'

#if DEF_MODEL==MODEL_KiD
CALL save_dg(SUM(precip(:, :))/nxny, 'precip', i_dgtime)
CALL save_dg(SUM(precip(:, :))/nxny*3600.0, 'surface_precip_mmhr', i_dgtime)
#endif

    ! deallocate diagnostics
DEALLOCATE(precip)

    ! deallocate process rates
IF (l_debug_main)PRINT*,'DEBUG-M5'
CALL deallocate_procs(procs)
IF (l_debug_main)PRINT*,'DEBUG-M6'
DEALLOCATE(procs)
DEALLOCATE(qfields)
DEALLOCATE(tend)
DEALLOCATE(precondition)
IF (l_debug_main)PRINT*,'DEBUG-M6a'


    ! aerosol fields
IF (l_process)CALL deallocate_procs(aerosol_procs)

IF (l_debug_main)PRINT*,'DEBUG-M6d'

DEALLOCATE(dustliq)
DEALLOCATE(aeroice)

CALL deallocate_aerosol(aerophys, aerochem)
DEALLOCATE(aerophys)
DEALLOCATE(aerochem)
DEALLOCATE(aeroact)
CALL deallocate_aerosol(dustphys, dustchem)
DEALLOCATE(dustphys)
DEALLOCATE(dustchem)
DEALLOCATE(dustact)
DEALLOCATE(aerosol_procs)
DEALLOCATE(aerosol_tend)
DEALLOCATE(aerofields)
DEALLOCATE(daerofields)

IF (l_debug_main)PRINT*,'DEBUG-M7'

END SUBROUTINE shipway_microphysics

SUBROUTINE microphysics_common(dt, kl, ku, qfields, dqfields, tend, procs, precip &
     , aerophys, aerochem, aeroact                                          &
     , dustphys, dustchem, dustact                                          &
     , aeroice, dustliq                                                     &
     , aerofields, daerofields, aerosol_tend, aerosol_procs                 &
     , rhcrit )

REAL(wp), INTENT(IN) :: dt  ! timestep from parent model
INTEGER, INTENT(IN) :: kl, ku
REAL(wp), INTENT(IN) :: rhcrit(:)
REAL(wp), INTENT(INOUT) :: qfields(:,:), dqfields(:,:), tend(:,:)
TYPE(process_rate), INTENT(INOUT) :: procs(:,:)
REAL(wp), INTENT(OUT) :: precip

    ! Aerosol fields
TYPE(aerosol_phys), INTENT(INOUT)   :: aerophys(:)
TYPE(aerosol_chem), INTENT(IN)      :: aerochem(:)
TYPE(aerosol_active), INTENT(INOUT) :: aeroact(:)
TYPE(aerosol_phys), INTENT(INOUT)   :: dustphys(:)
TYPE(aerosol_chem), INTENT(IN)      :: dustchem(:)
TYPE(aerosol_active), INTENT(INOUT) :: dustact(:)

TYPE(aerosol_active), INTENT(INOUT) :: aeroice(:)
TYPE(aerosol_active), INTENT(INOUT) :: dustliq(:)

REAL(wp), INTENT(INOUT) :: aerofields(:,:), daerofields(:,:), aerosol_tend(:,:)
TYPE(process_rate), INTENT(INOUT), OPTIONAL :: aerosol_procs(:,:)


REAL(wp), ALLOCATABLE :: qfields_in(:,:)
REAL(wp), ALLOCATABLE :: qfields_mod(:,:)
REAL(wp), ALLOCATABLE :: aerofields_in(:,:)
REAL(wp), ALLOCATABLE :: aerofields_mod(:,:)
REAL(wp) :: step_length
REAL(wp) :: sed_length
INTEGER :: n, k, nsed, iq

LOGICAL :: l_Twarm   ! temperature above freezing
LOGICAL :: l_Tcold   ! temperature below freezing, i.e. .not. l_Twarm

LOGICAL :: l_sigevap ! Is there significant evaporation of rain

ALLOCATE(qfields_in(nz, nq))
ALLOCATE(qfields_mod(nz, nq))
qfields_in=qfields ! Initial values of q
qfields_mod=qfields ! Modified initial values of q (may be modified if bad values sent in)

ALLOCATE(dist_lambda(nz,nspecies))
ALLOCATE(dist_mu(nz,nspecies))
ALLOCATE(dist_n0(nz,nspecies))

nsubsteps = MAX(1, CEILING(dt/max_step_length))
step_length = dt/nsubsteps

nsubseds = MAX(1, CEILING(step_length/max_sed_length))
sed_length = step_length/nsubseds

n_sub=1
n_subsed=1

IF (l_tendency_loc) THEN! Parent model uses tendencies
  qfields_mod=qfields_in + dt*dqfields
ELSE! Parent model uses increments
  qfields_mod=qfields_in + dqfields
END IF

IF (.NOT. l_passive) THEN
  CALL tidy_qin(qfields_mod)
END IF

    !---------------------------------------------------------------
    ! Determine (and possibly limit) size distribution
    !---------------------------------------------------------------
CALL query_distributions(cloud_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
CALL query_distributions(rain_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
IF (.NOT. l_warm_loc) THEN
  CALL query_distributions(ice_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
  CALL query_distributions(snow_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
  CALL query_distributions(graupel_params, qfields_mod, dist_lambda, dist_mu, dist_n0, icall=1)
END IF

qfields=qfields_mod

IF (aerosol_option > 0) THEN
  ALLOCATE(aerofields_in(nz, naero))
  ALLOCATE(aerofields_mod(nz, naero))
  aerofields_in=aerofields ! Initial values of aerosol
  aerofields_mod=aerofields ! Modified initial values  (may be modified if bad values sent in)

  IF (l_tendency_loc) THEN! Parent model uses tendencies
    aerofields_mod=aerofields_in + dt*daerofields
  ELSE! Parent model uses increments
    aerofields_mod=aerofields_in + daerofields
  END IF

  IF (l_process)CALL tidy_ain(qfields_mod, aerofields_mod)

aerofields=aerofields_mod
END IF

IF (l_debug_common)PRINT*,'DEBUG1'

DO n=1,nsubsteps
  n_sub=n

  CALL preconditioner(qfields)

      !------------------------------------------------------
      ! Early exit if we will have nothing to do.
      ! (i.e. no hydrometeors and subsaturated)
      !------------------------------------------------------
  IF (.NOT. ANY(precondition(:)))EXIT

      !-------------------------------
      ! Derive aerosol distribution
      ! parameters
      !-------------------------------
  IF (aerosol_option > 0)  CALL examine_aerosol(aerofields, qfields, aerophys, aerochem, &
         aeroact, dustphys, dustchem, dustact, aeroice, dustliq, icall=1)


  IF (l_debug_common)PRINT*,'DEBUG2'
  DO k=nz,1,-1
    k_here=k

    IF (precondition(k)) THEN
      l_Twarm = TdegK(k) > 273.15
      l_Tcold = .NOT. l_Twarm
          !=================================
          !
          ! WARM MICROPHYSICAL PROCESSES....
          !
          !=================================
          !-------------------------------
          ! Do the autoconversion to rain
          !-------------------------------
      IF (pswitch%l_praut)CALL raut(step_length, k, qfields, aerofields, procs, aerosol_procs)

          !-------------------------------
          ! Do the rain accreting cloud
          !-------------------------------
      IF (pswitch%l_pracw)CALL racw(step_length, k, qfields, aerofields, procs, rain_params, aerosol_procs)

          !-------------------------------
          ! Do the rain self-collection
          !-------------------------------
      IF (pswitch%l_pracr)CALL racr(step_length, k, qfields, procs)

          !-------------------------------
          ! Do the evaporation of rain
          !-------------------------------
      IF (pswitch%l_prevp)CALL revp(step_length, k, qfields, aerofields, aerophys, aerochem, &
             aeroact, dustliq, procs, aerosol_procs, l_sigevap)
          !=================================
          !
          ! ICE MICROPHYSICAL PROCESSES....
          !
          !=================================
      IF (.NOT. l_warm_loc) THEN

        IF (l_Tcold) THEN
              !------------------------------------------------------
              ! Autoconverion to snow
              !------------------------------------------------------
          IF (pswitch%l_psaut) CALL saut(step_length, k, qfields, aerofields, procs, aerosol_procs)

              !------------------------------------------------------
              ! Accretion processes
              !------------------------------------------------------
              ! Ice -> Cloud -> Ice
          IF (pswitch%l_piacw)CALL iacc(step_length, k, ice_params, cloud_params, ice_params, qfields, &
                   procs, aeroact, dustliq, aerosol_procs)
              ! Snow -> Cloud -> Snow
          IF (l_sg) THEN
            IF (pswitch%l_psacw)CALL iacc(step_length, k, snow_params, cloud_params, snow_params, &
                      qfields, procs, aeroact, dustliq, aerosol_procs)
                 ! Snow -> Ice -> Snow
            IF (pswitch%l_psaci)CALL iacc(step_length, k, snow_params, ice_params, snow_params, qfields, &
                      procs, aeroact, dustliq, aerosol_procs)
            IF (pswitch%l_praci .AND. (.NOT. l_sigevap)) THEN
                    ! Rain mass dependence to decide if we produce snow or graupel
              IF (qfields(k,i_qr) > thresh_sig(i_qr) .AND. l_g .AND. l_raci_g) THEN
                       ! Rain -> Ice -> Graupel
                CALL iacc(step_length, k, rain_params, ice_params, graupel_params, &
                     qfields, procs, aeroact, dustliq, aerosol_procs)
              ELSE
                       ! Rain -> Ice -> Snow
                CALL iacc(step_length, k, rain_params, ice_params, snow_params, qfields, &
                     procs, aeroact, dustliq, aerosol_procs)
              END IF
            END IF
            IF (pswitch%l_psacr .AND. (.NOT. l_sigevap)) THEN
                    ! Temperature dependence to decide if we produce snow or graupel
              IF (TdegK(k) < 268.15 .AND. l_g) THEN
                       ! Snow -> Rain -> Graupel
                CALL iacc(step_length, k, snow_params, rain_params, graupel_params, &
                     qfields, procs, aeroact, dustliq, aerosol_procs)
              ELSE
                       ! Snow -> Rain -> Snow
                CALL iacc(step_length, k, snow_params, rain_params, snow_params, &
                     qfields, procs, aeroact, dustliq, aerosol_procs)
              END IF
            END IF
            IF (l_g) THEN
                    ! Graupel -> Cloud -> Graupel
              IF (pswitch%l_pgacw)CALL iacc(step_length, k, graupel_params, cloud_params, &
                 graupel_params, qfields, procs, aeroact, dustliq, aerosol_procs)
                    ! Graupel -> Rain -> Graupel
              IF (pswitch%l_pgacr .AND. (.NOT. l_sigevap))CALL iacc(step_length, k,       &
                 graupel_params, rain_params, graupel_params, qfields,                    &
                 procs, aeroact, dustliq, aerosol_procs)
                    ! Graupel -> Ice -> Graupel
                    !                   if(pswitch%l_gsaci)call iacc(step_length, k, graupel_params, ice_params, graupel_params, qfields, &
                    !                       procs, aeroact, dustliq, aerosol_procs)
                    ! Graupel -> Snow -> Graupel
                    !                   if(pswitch%l_gsacs)call iacc(step_length, k, graupel_params, snow_params, graupel_params, qfields, &
                    !                       procs, aeroact, dustliq, aerosol_procs)
            END IF
          END IF


              !------------------------------------------------------
              ! Small snow accreting cloud should be sent to graupel
              ! (Ikawa & Saito 1991)
              !------------------------------------------------------
          IF (l_g .AND. .NOT. l_onlycollect)CALL graupel_embryos(step_length, k, qfields, &
                 procs, aerophys, aerochem, aeroact, aerosol_procs)

              !------------------------------------------------------
              ! Wet deposition/shedding (resulting from graupel
              ! accretion processes)
              ! NB This alters some of the accretion processes, so
              ! must come after their calculation and before they
              ! are used/rescaled elsewhere
              !------------------------------------------------------
          IF (l_g .AND. .NOT. l_onlycollect .AND. (.NOT. l_sigevap))CALL wetgrowth(step_length, k, qfields, &
                 procs, aerophys, aerochem, aeroact, aerosol_procs)

              !------------------------------------------------------
              ! Aggregation (self-collection)
              !------------------------------------------------------
          IF (pswitch%l_psagg)CALL ice_aggregation(step_length, k, snow_params, qfields, procs)

              !------------------------------------------------------
              ! Break up (snow only)
              !------------------------------------------------------
          IF (pswitch%l_psbrk)CALL ice_breakup(step_length, k, snow_params, qfields, procs)

              !------------------------------------------------------
              ! Ice multiplication (Hallet-mossop)
              !------------------------------------------------------
          IF (pswitch%l_pihal)CALL hallet_mossop(step_length, k, qfields,     &
                 procs, aerophys, aerochem, aeroact, aerosol_procs)

              !------------------------------------------------------
              ! Homogeneous freezing (rain and cloud)
              !------------------------------------------------------
          IF (pswitch%l_phomr .AND. (.NOT. l_sigevap))CALL ihom_rain(step_length, k, qfields, &
             aeroact, dustliq, procs, aerosol_procs)
          IF (pswitch%l_phomc)CALL ihom_droplets(step_length, k, qfields, aeroact, dustliq, procs, aerosol_procs)

              !------------------------------------------------------
              ! Condensation/immersion/contact nucleation of cloud ice
              !------------------------------------------------------
          IF (pswitch%l_pinuc)CALL inuc(step_length, k, qfields, procs,     &
                 dustphys, dustchem, aeroact, dustliq, aerosol_procs)

              !------------------------------------------------------
              ! Deposition/sublimation of ice/snow/graupel
              !------------------------------------------------------
          IF (pswitch%l_pidep) CALL idep(step_length, k, ice_params, qfields,  &
                 procs, dustact, aeroice, aerosol_procs)

          IF (pswitch%l_psdep) CALL idep(step_length, k, snow_params, qfields, &
                 procs, dustact, aeroice, aerosol_procs)

          IF (pswitch%l_pgdep) CALL idep(step_length, k, graupel_params, qfields, &
                 procs, dustact, aeroice, aerosol_procs)

          IF (l_harrington .AND. .NOT. l_onlycollect) CALL adjust_dep(dt, k, procs, qfields)

              !-----------------------------------------------------------
              ! Make sure we don't remove more than saturation allows
              !-----------------------------------------------------------
          IF (l_idep)CALL ensure_saturated(k, step_length, qfields, procs, (/i_idep, i_sdep, i_gdep/))
              !-----------------------------------------------------------
              ! Make sure we don't put back more than saturation allows
              !-----------------------------------------------------------
          IF (l_isub)CALL ensure_saturated(k, step_length, qfields, procs, (/i_isub, i_ssub, i_gsub/))

        ELSE ! T is  above freezing
              !------------------------------------------------------
              ! Melting of ice/snow/graupel
              !------------------------------------------------------
          IF (pswitch%l_psmlt .AND. (.NOT. l_sigevap))CALL melting(step_length, k, snow_params, qfields, &
                 procs, aeroice, dustact, aerosol_procs)
          IF (pswitch%l_pgmlt .AND. (.NOT. l_sigevap))CALL melting(step_length, k, graupel_params, qfields, &
                   procs, aeroice, dustact, aerosol_procs)
          IF (pswitch%l_pimlt .AND. (.NOT. l_sigevap))CALL melting(step_length, k, ice_params, qfields, &
                 procs, aeroice, dustact, aerosol_procs)

        END IF

            !-----------------------------------------------------------
            ! Make sure we don't remove more than we have to start with
            !-----------------------------------------------------------
        IF (l_pos1)CALL ensure_positive(k, step_length, qfields, procs, cloud_params, &
               (/i_praut, i_pracw, i_iacw, i_sacw, i_gacw, i_homc/),  &
               aeroprocs=aerosol_procs, iprocs_dependent=(/i_aaut, i_aacw/))

        IF (.NOT. l_warm_loc) THEN
          IF (l_pos2)CALL ensure_positive(k, step_length, qfields, procs, ice_params, &
                 (/i_raci, i_saci, i_gaci, i_saut, i_isub, i_imlt/), &
                 (/i_ihal, i_gshd, i_inuc, i_homc, i_iacw, i_idep/))
        END IF

        IF (l_pos3) CALL ensure_positive(k, step_length, qfields, procs, rain_params, &
               (/i_prevp, i_sacr, i_gacr, i_homr/), &
               (/i_praut, i_pracw, i_raci, i_gshd, i_smlt, i_gmlt/), &
               aeroprocs=aerosol_procs, iprocs_dependent=(/i_arevp/))

        IF (.NOT. l_warm_loc) THEN
          IF (l_pos4)CALL ensure_positive(k, step_length, qfields, procs, snow_params, &
                 (/i_gacs, i_smlt, i_sacr, i_ssub /), &
                 (/i_sdep, i_sacw, i_saut, i_saci, i_raci, i_gshd, i_ihal/))
        END IF

      ELSE

        IF (pswitch%l_praut .AND. pswitch%l_pracw) THEN
          IF (l_pos5)CALL ensure_positive(k, step_length, qfields, procs, cloud_params, (/i_praut, i_pracw/), &
                     aeroprocs=aerosol_procs, iprocs_dependent=(/i_aaut, i_aacw/))
        END IF

        IF (pswitch%l_prevp) THEN
          IF (l_pos6)CALL ensure_positive(k, step_length, qfields, procs, rain_params, (/i_prevp/), (/i_praut, i_pracw/) &
                     , aerosol_procs, (/i_arevp/))
        END IF

      END IF

    END IF
  END DO
  IF (l_debug_common)PRINT*,'DEBUG3'
      !-------------------------------
      ! Collect terms we have so far
      !-------------------------------
  CALL sum_procs(step_length, procs, tend,      &
       (/i_praut, i_pracw, i_pracr, i_prevp/)     &
       , hydro_names, l_thermalexchange=.TRUE., qfields=qfields     &
       , l_passive=l_passive)

  IF (.NOT. l_warm_loc) THEN
    CALL sum_procs(step_length, procs, tend,      &
         (/i_idep, i_sdep, i_gdep, i_iacw, i_sacr, i_sacw, i_saci, i_raci     &
         , i_gacw, i_gacr, i_gaci, i_gacs, i_ihal, i_gshd, i_sbrk           &
         , i_saut, i_iagg, i_sagg, i_gagg, i_isub, i_ssub, i_gsub/),        &
         hydro_names, l_thermalexchange=.TRUE., qfields=qfields             &
         , l_passive=l_passive, i_thirdmoment=2)
    CALL sum_procs(step_length, procs, tend,      &
         (/i_inuc, i_imlt, i_smlt, i_gmlt, i_homr, i_homc/), hydro_names,     &
         l_thermalexchange=.TRUE., qfields=qfields     &
         , l_passive=l_passive)
  END IF

  CALL update_q(qfields_mod, qfields, tend)


  IF (l_process) THEN

    DO k=1,nz
      IF (precondition(k)) THEN
        CALL ensure_positive_aerosol(k, step_length, aerofields, aerosol_procs,&
             (/i_aaut, i_aacw, i_aevp, i_arevp, i_dnuc, i_dsub, i_dssub, i_dgsub, &
             i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr,    &
             i_dgacr, i_draci  /) )
      END IF
    END DO

    CALL sum_aprocs(step_length, aerosol_procs, aerosol_tend,      &
         (/i_aaut, i_aacw, i_aevp, i_arevp/), aero_names, aerofields)

    IF (.NOT. l_warm) THEN
      CALL sum_aprocs(step_length, aerosol_procs, aerosol_tend,      &
           (/i_dnuc, i_dsub, i_dssub, i_dgsub, i_dimlt, i_dsmlt, i_dgmlt,     &
           i_diacw, i_dsacw, i_dgacw, i_dsacr, i_dgacr, i_draci /), aero_names, aerofields)
    END IF

    CALL update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.TRUE.)

        !-------------------------------
        ! Re-Derive aerosol distribution
        ! parameters
        !-------------------------------
    CALL examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact     &
         , dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
  END IF

      !-------------------------------
      ! Do the condensation/evaporation
      ! of cloud and activation of new
      ! drops
      !-------------------------------
  IF (pswitch%l_pcond)THEN
    DO k=nz,1,-1
#if DEF_MODEL==MODEL_KiD
      k_here=k
#endif
#if DEF_MODEL==MODEL_LEM_DIAG
      k_here=k
#endif
      if (k==104.and.i_here==67.and.j_here==3) print*, qfields(k,i_qr), qfields(k,i_nr)
      CALL condevp(step_length, k, qfields, aerofields,     &
         procs, aerophys, aerochem, aeroact, dustphys, dustchem, dustliq, &
         aerosol_procs, rhcrit(k))
    END DO

    !-------------------------------
    ! Collect terms we have so far
    !-------------------------------
    CALL sum_procs(step_length, procs, tend,      &
       (/i_cond/)     &
       , hydro_names, l_thermalexchange=.TRUE., qfields=qfields     &
       , l_passive=l_passive)
    
    CALL update_q(qfields_mod, qfields, tend, l_fixneg=.TRUE.)

  END IF

      !---------------------------------------------------------------
      ! Re-Determine (and possibly limit) size distribution
      !---------------------------------------------------------------
  CALL query_distributions(cloud_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
  CALL query_distributions(rain_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
  IF (.NOT. l_warm_loc) THEN
    CALL query_distributions(ice_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
    CALL query_distributions(snow_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
    CALL query_distributions(graupel_params, qfields, dist_lambda, dist_mu, dist_n0, icall=2)
  END IF

  IF (l_process) THEN
    CALL sum_aprocs(step_length, aerosol_procs, aerosol_tend,      &
         (/i_aact /), aero_names, aerofields)

    CALL update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.TRUE.)

        !-------------------------------
        ! Re-Derive aerosol distribution
        ! parameters
        !-------------------------------
    CALL examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact     &
         , dustphys, dustchem, dustact, aeroice, dustliq, icall=3)
  END IF

  IF (l_debug_common)PRINT*,'DEBUG4'



#if DEF_MODEL==MODEL_LEM_DIAG
! Need to recode these with new process rates
      ! Save the process rates
      !call lem_procdgs(n, kl, ku, procs)
#endif
      !-------------------------------
      ! Do the sedimentation
      !-------------------------------

  IF (l_sed) THEN

    DO nsed=1,nsubseds
      n_subsed=nsed

      IF (nsed > 1) THEN
            !-------------------------------
            ! Reset process rates if they
            ! are to be re-used
            !-------------------------------
        CALL zero_procs(procs)
        IF (l_process)CALL zero_procs(aerosol_procs)
            !---------------------------------------------------------------
            ! Re-Determine (and possibly limit) size distribution
            !---------------------------------------------------------------
        CALL query_distributions(cloud_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
        CALL query_distributions(rain_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
        IF (.NOT. l_warm_loc) THEN
          CALL query_distributions(ice_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
          CALL query_distributions(snow_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
          CALL query_distributions(graupel_params, qfields, dist_lambda, dist_mu, dist_n0, icall=nsed+2)
        END IF

            !-------------------------------
            ! Re-Derive aerosol distribution
            ! parameters
            !-------------------------------
        IF (l_process)CALL examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact &
               , dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
      END IF

      IF (pswitch%l_psedl) THEN
        CALL sedr(sed_length, qfields, aerofields, aeroact, dustliq, tend, cloud_params, procs, aerosol_procs, precip, l_process)

        CALL sum_procs(sed_length, procs, tend, (/i_psedl/), hydro_names, qfields=qfields)
      END IF

      IF (pswitch%l_psedr) THEN
        CALL sedr(sed_length, qfields, aerofields, aeroact, dustliq, tend, rain_params, procs, aerosol_procs, precip, l_process)

        CALL sum_procs(sed_length, procs, tend, (/i_psedr/), hydro_names, qfields=qfields)
      END IF

      IF (.NOT. l_warm_loc) THEN

        IF (pswitch%l_psedi) THEN
          CALL sedr(sed_length, qfields, aerofields, aeroice, dustact, tend, ice_params, procs, aerosol_procs, precip, l_process)
        END IF
        IF (pswitch%l_pseds) THEN
          CALL sedr(sed_length, qfields, aerofields, aeroice, dustact, tend, snow_params, &
             procs, aerosol_procs, precip, l_process)
        END IF
        IF (pswitch%l_psedg) THEN
          CALL sedr(sed_length, qfields, aerofields, aeroice, dustact, tend, graupel_params, &
             procs, aerosol_procs, precip, l_process)
        END IF

        CALL sum_procs(sed_length, procs, tend, (/i_psedi, i_pseds, i_psedg/), hydro_names, qfields=qfields)
      END IF

      CALL update_q(qfields_mod, qfields, tend, l_fixneg=.TRUE.)

      IF (l_process) THEN
        CALL sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_asedr, i_asedl/) &
             , aero_names, aerofields)

        CALL update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.TRUE.)

        IF (l_process .AND. .NOT. l_warm) THEN
          CALL sum_aprocs(sed_length, aerosol_procs, aerosol_tend, (/i_dsedi, i_dseds, i_dsedg/) &
               , aero_names, aerofields)

          CALL update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.TRUE.)
        END IF
      END IF

    END DO

  END IF

      !-------------------------------
      ! Reset process rates if they
      ! are to be re-used
      !-------------------------------
  IF (nsubsteps>1)CALL zero_procs(procs)
  IF (nsubsteps>1 .AND. l_process)CALL zero_procs(aerosol_procs)


END DO
IF (l_debug_common)PRINT*,'DEBUG5'

    !--------------------------------------------------
    ! Tidy up any small/negative numbers
    ! we may have generated.
!--------------------------------------------------
IF (pswitch%l_tidy2) THEN
  DO k=nz,1,-1

    k_here=k

    IF (precondition(k)) THEN
      CALL qtidy(step_length, k, qfields, procs, aerofields     &
         , aeroact, dustact, aeroice, dustliq           &
         , aerosol_procs, i_tidy2, i_atidy2,     &
         l_negonly=.TRUE.)
    END IF

  END DO

  CALL sum_procs(step_length, procs, tend, (/i_tidy2/), hydro_names, qfields=qfields &
     , l_passive=l_passive)

  CALL update_q(qfields_mod, qfields, tend)

  IF (l_process) THEN
    CALL sum_aprocs(step_length, aerosol_procs, aerosol_tend, (/i_atidy2/), aero_names, aerofields )
    CALL update_q(aerofields_mod, aerofields, aerosol_tend, l_aerosol=.TRUE.)
  END IF

END IF

    !
    ! Add on initial adjustments that may have been made
    !

IF (l_tendency_loc) THEN! Convert back from cummulative value to tendency
  tend = tend + qfields_mod - qfields_in - dqfields*dt
  tend = tend/dt
ELSE
  tend = tend + qfields_mod - qfields_in - dqfields
      ! prevent negative values
  DO iq=i_hstart,ntotalq
    DO k=1,nz
      tend(k,iq) = MAX(tend(k,iq), -(qfields_in(k,iq) - dqfields(k,iq)))
    END DO
  END DO

END IF

    ! if (timestep_number*dt < 3600 .and. timestep_number > 0)then
    !   tend = (timestep_number*dt/3600.)*tend  ! Slow spin up period down after initial shock.
    ! end if


IF (aerosol_option > 0) THEN
      ! processing
  IF (l_process) THEN
    IF (l_tendency_loc) THEN! Convert back from cummulative value to tendency
      aerosol_tend = aerosol_tend + aerofields_mod - aerofields_in - daerofields*dt
      aerosol_tend = aerosol_tend/dt
    ELSE
      aerosol_tend = aerosol_tend + aerofields_mod - aerofields_in - daerofields
            ! prevent negative values
      DO iq=1,ntotala
        DO k=1,nz
          aerosol_tend(k,iq) = MAX(aerosol_tend(k,iq), -(aerofields_in(k,iq) - daerofields(k,iq)))
        END DO
      END DO
    END IF
  END IF
  DEALLOCATE(aerofields_in)
  DEALLOCATE(aerofields_mod)
END IF
IF (l_debug_common)PRINT*,'DEBUG6'

DEALLOCATE(qfields_in)
DEALLOCATE(qfields_mod)

DEALLOCATE(dist_lambda)
DEALLOCATE(dist_mu)
DEALLOCATE(dist_n0)

END SUBROUTINE microphysics_common

SUBROUTINE update_q(qfields_in, qfields, tend, l_aerosol, l_fixneg)

REAL(wp), INTENT(IN)    :: qfields_in(:,:)
REAL(wp), INTENT(INOUT) :: qfields(:,:)
REAL(wp), INTENT(IN)    :: tend(:,:)
LOGICAL, INTENT(IN), OPTIONAL :: l_aerosol ! flag to indicate updating of aerosol
                                               !  - may be used for diagnostics.
LOGICAL, INTENT(IN), OPTIONAL :: l_fixneg  ! Flag to use cludge to bypass negative/zero numbers
INTEGER :: k, iqx
LOGICAL :: l_fix

l_fix=.FALSE.
IF (PRESENT(l_fixneg))l_fix=l_fixneg

DO iqx=1,UBOUND(tend,2)
  DO k=LBOUND(tend,1), UBOUND(tend,1)
    qfields(k,iqx) = qfields_in(k,iqx) + tend(k,iqx)
!    IF (PRESENT(l_aerosol)) THEN
!      IF (l_aerosol .AND. qfields(k,iqx)< -1e-20) THEN
!        PRINT*, 'NEGA', k, i_here, iqx, qfields(k,iqx), SPACING(tend(k,iqx))
!      END IF
!    END IF
        !quick lem fixes  - this code should never be used ?
    IF (.NOT. PRESENT(l_aerosol) .AND. l_fix) THEN
      IF (iqx==i_ni .AND. qfields(k,iqx)<=0) THEN
        qfields(k,iqx)=0
        qfields(k,i_qi)=0
      END IF
      IF (iqx==i_nr .AND. qfields(k,iqx)<=0) THEN
        qfields(k,iqx)=0
        qfields(k,i_qr)=0
        IF (i_m3r/=0)qfields(k,i_m3r)=0
      END IF
      IF (iqx==i_nl .AND. qfields(k,iqx)<=0) THEN
        qfields(k,iqx)=0
        qfields(k,i_ql)=0
      END IF
      IF (iqx==i_ns .AND. qfields(k,iqx)<=0) THEN
        qfields(k,iqx)=0
        qfields(k,i_qs)=0
        IF (i_m3s/=0)qfields(k,i_m3s)=0
      END IF
      IF (iqx==i_ng .AND. qfields(k,iqx)<=0) THEN
        qfields(k,iqx)=0
        qfields(k,i_qg)=0
        IF (i_m3g/=0)qfields(k,i_m3g)=0
      END IF
    END IF
  END DO
END DO

END SUBROUTINE update_q

SUBROUTINE update_tend(dt, qfields, dqfields, tend)

REAL(wp) :: dt
REAL(wp), INTENT(IN) :: qfields(:,:)
REAL(wp), INTENT(IN) :: dqfields(:,:)
REAL(wp), INTENT(INOUT) :: tend(:,:)

REAL(wp), ALLOCATABLE :: tmp(:,:)

INTEGER :: k, iqx

    ! Using tmp for premultiplication seems to optimize
    ! the loop much better for some compilers
ALLOCATE(tmp(UBOUND(dqfields,1), LBOUND(dqfields,2):UBOUND(dqfields,2)))
tmp= dt*dqfields

DO iqx=1,UBOUND(tend,2)
  DO k=LBOUND(tend,1), UBOUND(tend,1)
    tend(k,iqx) = tend(k,iqx) + tmp(k,iqx)
  END DO
END DO

DEALLOCATE(tmp)

END SUBROUTINE update_tend

END MODULE micro_main
