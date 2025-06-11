module mphys_switches
  use variable_precision, only: wp
  use mphys_parameters, only: cloud_params, rain_params, ice_params,           &
                              rain_params_orig, rain_params_kf,                &
                              ice_params_orig, ice_params_kf,                  &
                              snow_params, graupel_params, p1,p2,p3
  use thresholds, only: thresh_small, th_small, qv_small, ql_small, qr_small, nl_small, nr_small, m3r_small, &
       qi_small, qs_small, ni_small, ns_small, m3s_small, qg_small, ng_small, m3g_small, thresh_tidy, &
       th_tidy, qv_tidy, ql_tidy, qr_tidy, nl_tidy, nr_tidy, m3r_tidy, qi_tidy, qs_tidy, ni_tidy, &
       ns_tidy, m3s_tidy, qg_tidy, ng_tidy, m3g_tidy, thresh_sig, qs_sig, ql_sig, qr_sig, qi_sig, qg_sig, &
       thresh_large, qs_large, ql_large, qr_large, qi_large, qg_large, ns_large, nl_large, nr_large, ni_large, &
       ng_large, thresh_atidy, aeromass_small, aeronumber_small
  use process_routines, only: i_cond, i_praut, &
       i_pracw, i_pracr, i_prevp, i_psedr, i_psedl, i_aact, i_aaut, i_aacw, i_aevp, i_asedr, i_asedl, i_arevp, &
       i_tidy, i_tidy2, i_atidy, i_atidy2, i_inuc, i_idep, i_dnuc, i_dsub, i_saut, i_iacw, i_sacw, i_pseds, &
       i_sdep, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs, i_gdep, i_psedg, i_iagg, i_sagg, i_gagg, &
       i_gshd, i_ihal, i_smlt, i_gmlt, i_psedi, i_homr, i_homc, i_imlt, i_isub, i_ssub, i_gsub, i_sbrk, i_dssub, &
       i_dgsub, i_dsedi, i_dseds, i_dsedg, i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr, &
       i_dgacr, i_draci, i_dhomc, i_dhomr, i_iics, i_idps, process_name

  implicit none

  character(len=*), parameter, private :: ModuleName='MPHYS_SWITCHES'

  logical :: mphys_is_set = .false.

  ! Microphysics option are now specified through a 5 digit integer, such that
  ! each digit represents the number of moments of cloud, rain, ice, snow, graupel
  ! respectively.  E.g. 23000 would use double moment cloud and triple moment rain
  integer :: option = 22222

  !-------------------------------------------------
  ! Switches to determine what aerosol we should use
  ! (replaces old aerosol_option)
  !-------------------------------------------------
  ! number of moments to use for soluble modes
  ! (/ aitken, accumulation, coarse /)
  integer :: soluble_modes(3) = (/ 0, 0, 0/)
  ! number of moments to use for insoluble modes
  ! (/ accumulation, coarse /)
  integer :: insoluble_modes(2) = (/ 0, 0/)
  ! do we want to carry activated aersols in cloud/rain?
  ! (/ for soluble aerosol, for insoluble aerosol /)
  logical :: active_cloud(2) = (/.false., .false./)
  ! do we want to carry activated aersols in ice/snow/graupel?
  ! (/ for soluble aerosol, for insoluble aerosol /)
  logical :: active_ice(2) = (/.false., .false./)
  ! do we want to carry a separate activated category in rain?
  ! (/ for soluble aerosol /)
  logical :: active_rain(1) = (/.false./)
  ! do we want to retain information amout activated number
  ! so that we can use this to replenish aerosol (rather than
  ! using the processed hydrometeor numbers)
  ! (/ for soluble aerosol, for insoluble aerosol /)
  logical :: active_number(2) = (/.false., .false./)
  logical :: l_warm = .false.  ! Only use warm rain microphysics

  logical :: l_cfrac_casim_diag_scheme = .false. ! Use diagnostic cloud scheme

  logical :: l_prf_cfrac = .true. ! Paul Field's new Smith Cloud Fraction
                                   ! scheme

  logical :: l_adjust_D0 = .true.  !adjust the psd in distributions
  
  logical :: l_kk00 = .true.  ! true=use KK2000 autoconv+accretion, false=use KK2013

! Flag to decide whether to transfer the evaporating aerosol based on whether 
! it is less or more than halfway between the sizes of accum and coarse modes
! This is True (Dan Grosvenor Bug) in the package branch
! Turn to false to turn off Dan's bug fix
  logical :: l_aeroproc_midway =.FALSE.

  logical :: l_aeroproc_no_coarse =.false.
  logical :: l_bypass_which_mode = .FALSE. ! DPG - bypass which_mode_to_use and
                                          !put all evaporated aeroosl in either
                                          !accum or coarse mode
!DPG - option of where to transfer aerosol from evap; 1=accum mode (default),
!2 = coarse, anything else = keep mass in activated mode
!(but transfer number to the mode selected by iact_mode below)  
  integer :: iopt_which_mode=1  

  logical :: l_ukca_casim = .false. ! CASIM is coupled to UKCA. 
  ! This is set to false, but needs setting to .true. if l_ukca
  ! is true in the UM. This is done in init_casim_run_mod.F90

  character(10), allocatable :: hydro_names(:)
  character(20), allocatable :: aero_names(:)

  ! index to determine if something is a mass or a number
  integer :: inumber=1
  integer :: imass=2
  integer :: ihyg=3
  ! index to determine if something is a soluble(ccn) or insoluble(IN) (NB currently can't be both)
  integer :: isol=1
  integer :: iinsol=2
  ! ONLY for ARG. 0: don't activate in-cloud, 1: activate in-cloud as per
  ! ARG default. 2: use smaller smax out of ARG smax and the smax calculated
  ! from existing cloud droplets
  integer :: activate_in_cloud=2

  ! standard hydrometeor indices
  integer :: i_qv  = 0 ! water vapour
  integer :: i_ql  = 0 ! cloud liquid mass
  integer :: i_nl  = 0 ! cloud liquid number
  integer :: i_qr  = 0 ! rain mass
  integer :: i_nr  = 0 ! rain number
  integer :: i_m3r = 0 ! rain 3rd moment
  integer :: i_qi  = 0 ! ice mass
  integer :: i_ni  = 0 ! ice number
  integer :: i_qs  = 0 ! snow mass
  integer :: i_ns  = 0 ! snow number
  integer :: i_m3s = 0 ! snow 3rd moment
  integer :: i_qg  = 0 ! graupel mass
  integer :: i_ng  = 0 ! graupel number
  integer :: i_m3g = 0 ! graupel 3rd moment
  integer :: i_qh  = 0 ! hail mass
  integer :: i_nh  = 0 ! hail number
  integer :: i_m3h = 0 ! graupel 3rd moment
  integer :: i_th  = 0 ! theta
  
  !cloud fraction indices (fixed)
  integer :: i_cfl=1 !liquid
  integer :: i_cfr=2 !rain
  integer :: i_cfi=3 !ice
  integer :: i_cfs=4 !snow
  integer :: i_cfg=5 !graupel
  
  

  ! information about location indices of different moments
  integer :: i_qstart = 1 ! First index in for q variables (always 1)
  integer :: i_nstart     ! First index in for n variables
  integer :: i_m3start    ! First index in for m3 variables
  integer :: i_hstart = 3 ! First index in for hydrometeors (always 3)

  integer :: ntotalq   ! total number of q variables (includes theta and qv)
  integer :: ntotala   ! total number of aerosol variables

  integer :: nactivea  ! total number of active aerosol(soluble) variables
  integer :: nactived  ! total number of active aerosol(insoluble) variables

  ! number of moments for each hydrometeor
  integer :: nq_l = 0  ! cloud liquid
  integer :: nq_r = 0  ! rain
  integer :: nq_i = 0  ! cloud ice
  integer :: nq_s = 0  ! snow
  integer :: nq_g = 0  ! graupel

  ! logicals for moments
  logical :: l_2mc = .false.
  logical :: l_2mr = .false.
  logical :: l_3mr = .false.

  logical :: l_2mi = .false.
  logical :: l_2ms = .false.
  logical :: l_3ms = .false.
  logical :: l_2mg = .false.
  logical :: l_3mg = .false.

  ! aerosol indices (soluble)
  integer :: i_am1 = 0    ! aitken aerosol mass
  integer :: i_an1 = 0    ! aitken aerosol number
  integer :: i_am2 = 0    ! accumulation aerosol mass
  integer :: i_an2 = 0    ! accumulation aerosol number
  integer :: i_am3 = 0    ! coarse aerosol mass
  integer :: i_an3 = 0    ! coarse aerosol number
  integer :: i_am4 = 0    ! activated aerosol mass
  integer :: i_am5 = 0    ! activated aerosol mass (in rain)

  ! aerosol indices (insoluble)
  integer :: i_am6 = 0    ! dust mass
  integer :: i_an6 = 0    ! dust number
  integer :: i_am7 = 0    ! activated dust mass

  ! aerosol indices (soluble/insoluble)
  integer :: i_am8 = 0    ! soluble mass in ice (incl. snow/graupel)
  integer :: i_am9 = 0    ! dust mass in water (cloud/rain)

  ! aerosol indices (accumulation mode insoluble)
  integer :: i_am10 = 0    ! dust mass (accum)
  integer :: i_an10 = 0    ! dust number (accum)

  ! aerosol indices (activated number)
  integer :: i_an11 = 0    ! soluble
  integer :: i_an12 = 0    ! insoluble

  ! aerosol indices (hygroscopicity of soluble modes)
  integer :: i_ak1 = 0 ! Aitken
  integer :: i_ak2 = 0 ! accumulation
  integer :: i_ak3 = 0 ! coarse

  type :: complexity
     integer :: nspecies
     integer, pointer :: nmoments(:)
     integer :: nprocesses = 0
  end type complexity

  integer, parameter :: maxmodes=3  !< maximum ccn/in modes

  type :: aerosol_index              !< contains indexing imformation for aerosol as ccn/in
     integer :: nccn            = 0  !< How many aerosol act as ccn
     integer :: nin             = 0  !< How many aerosol act as in
     integer :: ccn_m(maxmodes) = 0  !< list of ccn mass indices
     integer :: ccn_n(maxmodes) = 0  !< list of ccn number indices
     integer :: ccn_k(maxmodes) = 0  !< list of ccn hygroscopicity indices
     integer :: in_m(maxmodes)  = 0  !< list of in mass indices
     integer :: in_n(maxmodes)  = 0  !< list of in number indices
     ! Some useful indices for the different modes
     integer :: i_aitken = 0     ! atiken mode
     integer :: i_accum  = 0     ! accumulation mode
     integer :: i_coarse = 0     ! Coarse mode
     integer :: i_accum_dust  = 0     ! accumulation mode dust
     integer :: i_coarse_dust = 0     ! Coarse mode dust
  end type aerosol_index

  type(complexity), save :: hydro_complexity
  type(complexity), save :: aero_complexity

  type(aerosol_index), save :: aero_index

  ! substepping
  real(wp) :: max_step_length = 10.0
  real(wp) :: max_sed_length = 2.0
  
  ! Switch used to determine whether gamma functions  computed for each timestep (true) 
  ! or at initialisation. Default is false, since precomputing the gamma functions results in 
  ! significant compute speed up and reduces the need for timestep lookups, which is beneficial 
  ! to GPUs. In general this should only be true when 3-moment or diagnostic shape is used
  logical :: l_gamma_online = .false. 
  logical :: l_subseds_maxv = .false.

  ! Switch to use eulexp sedimentation, which is based on the sedimentation from the
  ! UM. It permits one step sedimentation on long timesteps. False will give the 
  ! standard CASIM sed with substep on timesteps greater than max_sed_length
  ! Note: This has only been tested with single moment and double moment, not sure how this 
  !       works with three moment stuff
  logical :: l_sed_eulexp = .true.

  real(wp) :: cfl_vt_max = 1.0
  
  !mixed-phase overlap factor (1=max overlap, 0=min overlap) for cloud fraction
  real(wp) :: mpof = 0.5 
  
  logical :: l_srg = .false. ! if true snow collecting rain makes graupel, otherwise makes snow
  logical :: l_reisner_graupel_embryo = .false. ! if true use Reisner et al. 1998 QJRMS approach making graupel embryos 

  ! process switches
  ! Some of these switches are obsolete or inactive - review these
  ! and then remove this message

  logical :: l_abelshipway=.true.
  logical :: l_sed_3mdiff=.false.
  ! logicals to set 1M sedimentation for 2M or 3M configurations
  logical :: l_sed_icecloud_as_1m=.false.
  logical :: l_sed_rain_1m = .false.
  logical :: l_sed_snow_1m = .false.
  logical :: l_sed_graupel_1m = .false.
  logical :: l_cons=.false.
  logical :: l_inuc=.true.

  logical :: l_sg             = .true.  ! run with snow and graupel
  logical :: l_g              = .true.  ! run with graupel
  logical :: l_halletmossop   = .true.
  logical :: l_sip_icebreakup = .false.
  logical :: l_sip_dropletshatter = .false.
  logical :: l_no_pgacs_in_sumprocs = .false. ! If running ice breakup then no pgacs added in sumprocs

  logical :: l_harrington     = .false.  ! Use Jerry's method for ice autoconvertion

  logical :: l_condensation   = .true.  ! condense water
  logical :: l_evaporation    = .true.  ! evaporate rain
  logical :: l_rain           = .true.  ! rain sources
  logical :: l_sed            = .true.  ! sedimentation
  logical :: l_boussinesq     = .false. ! set rho=1 everywhere
  integer :: diag_mu_option   = -999    ! select diagnostic mu           &

  logical :: l_iaut           = .true.  ! autoconversion of ice/snow
  logical :: l_imelt          = .true.  ! melting of ice/snow/graupel
  logical :: l_iacw           = .true.  ! ice/snow/graupel accumulation cloud water
  logical :: l_idep           = .true.  ! deposition of ice/snow/graupel
  logical :: l_isub           = .true.  ! sublimation of ice/snow/graupel

  logical :: l_pos1           = .true.   ! switch on positivity check
  logical :: l_pos2           = .true.   ! switch on positivity check
  logical :: l_pos3           = .true.   ! switch on positivity check
  logical :: l_pos4           = .true.   ! switch on positivity check
  logical :: l_pos5           = .true.   ! switch on positivity check
  logical :: l_pos6           = .true.   ! switch on positivity check

  logical, target :: l_pcond   = .true.  ! Condensation
  logical, target :: l_praut   = .true.  ! Autoconversion cloud -> rain
  logical, target :: l_pracw   = .true.  ! Accretion  cloud -> rain
  logical, target :: l_pracr   = .true.  ! aggregation of rain drops
  logical, target :: l_prevp   = .true.  ! evaporation of rain
  logical, target :: l_psedl   = .true.  ! sedimentation of cloud
  logical, target :: l_psedr   = .true.  ! sedimentation of rain
  logical, target :: l_ptidy   = .true.  ! tidying term 1
  logical, target :: l_ptidy2  = .true.  ! tidying term 2
  logical, target :: l_pinuc   = .true.  ! ice nucleation
  logical, target :: l_pidep   = .true.  ! ice deposition
  logical, target :: l_piacw   = .true.  ! ice accreting water
  logical, target :: l_psaut   = .true.  ! ice autoconversion ice -> snow
  logical, target :: l_psdep   = .true.  ! vapour deposition onto snow
  logical, target :: l_psacw   = .true.  ! snow accreting water
  logical, target :: l_pgdep   = .true.  ! vapour deposition onto graupel
  logical, target :: l_pseds   = .true.  ! snow sedimentation
  logical, target :: l_psedi   = .true.  ! ice sedimentation
  logical, target :: l_psedg   = .true.  ! graupel sedimentation
  logical, target :: l_psaci   = .true.  ! snow accreting ice
  logical, target :: l_praci   = .true.  ! rain accreting ice
  logical, target :: l_psacr   = .true.  ! snow accreting rain
  logical, target :: l_pgacr   = .false.  ! graupel accreting rain
  logical, target :: l_pgacw   = .true.  ! graupel accreting cloud water
  logical, target :: l_pgaci   = .false.  ! graupel accreting ice
  logical, target :: l_pgacs   = .false.  ! graupel accreting snow
  logical, target :: l_piagg   = .false. ! aggregation of ice particles
  logical, target :: l_psagg   = .true.  ! aggregation of snow particles
  logical, target :: l_pgagg   = .false. ! aggregation of graupel particles
  logical, target :: l_psbrk   = .true.  ! break up of snow flakes
  logical, target :: l_pgshd   = .true.  ! shedding of liquid from graupel
  logical, target :: l_pihal   = .true.  ! hallet mossop
  logical, target :: l_piics   = .false.  ! ice-ice collision
  logical, target :: l_pidps   = .false.  ! droplet shattering
  logical, target :: l_psmlt   = .true.  ! snow melting
  logical, target :: l_pgmlt   = .true.  ! graupel melting
  logical, target :: l_phomr   = .true.  ! homogeneous freezing of rain
  logical, target :: l_phomc   = .true.  ! homogeneous freezing of cloud droplets
  logical, target :: l_pssub   = .true.  ! sublimation of snow
  logical, target :: l_pgsub   = .true.  ! sublimation of graupel
  logical, target :: l_pisub   = .true.  ! sublimation of ice
  logical, target :: l_pimlt   = .true.  ! ice melting

  logical :: l_tidy_conserve_E = .true.  !fixed in qtidy (missing exner)
  logical :: l_tidy_conserve_q = .true.
  logical :: l_tidy_negonly    = .true.  ! Only tidy up negative (i.e. not small) values 

  logical :: l_preventsmall = .false.  ! Modify tendencies to prevent production of small quantities

  type :: process_switch
     logical, pointer :: l_pcond ! Condensation
     logical, pointer :: l_praut ! Autoconversion cloud -> rain
     logical, pointer :: l_pracw ! Accretion  cloud -> rain
     logical, pointer :: l_pracr ! aggregation of rain drops
     logical, pointer :: l_prevp ! evaporation of rain
     logical, pointer :: l_psedl ! sedimentation of cloud
     logical, pointer :: l_psedr ! sedimentation of rain
     logical, pointer :: l_pinuc ! ice nucleation
     logical, pointer :: l_pidep ! ice deposition
     logical, pointer :: l_piacw ! ice accreting water
     logical, pointer :: l_psaut ! ice autoconversion ice -> snow
     logical, pointer :: l_psdep ! vapour deposition onto snow
     logical, pointer :: l_psacw ! snow accreting water
     logical, pointer :: l_pgdep ! vapour deposition onto graupel
     logical, pointer :: l_pseds ! snow sedimentation
     logical, pointer :: l_psedi ! ice sedimentation
     logical, pointer :: l_psedg ! graupel sedimentation
     logical, pointer :: l_psaci ! snow accreting ice
     logical, pointer :: l_praci ! rain accreting ice
     logical, pointer :: l_psacr ! snow accreting rain
     logical, pointer :: l_pgacr ! graupel accreting rain
     logical, pointer :: l_pgacw ! graupel accreting cloud water
     logical, pointer :: l_pgaci ! graupel accreting ice
     logical, pointer :: l_pgacs ! graupel accreting snow
     logical, pointer :: l_piagg ! aggregation of ice particles
     logical, pointer :: l_psagg ! aggregation of snow particles
     logical, pointer :: l_pgagg ! aggregation of graupel particles
     logical, pointer :: l_psbrk ! break up of snow flakes
     logical, pointer :: l_pgshd ! shedding of liquid from graupel
     logical, pointer :: l_pihal ! hallet mossop
     logical, pointer :: l_piics ! ice-ice collision
     logical, pointer :: l_pidps ! droplet shattering
     logical, pointer :: l_psmlt ! snow melting
     logical, pointer :: l_pgmlt ! graupel melting
     logical, pointer :: l_phomr ! homogeneous freezing of rain
     logical, pointer :: l_phomc ! homogeneous freezing of cloud droplets
     logical, pointer :: l_pssub ! sublimation of snow
     logical, pointer :: l_pgsub ! sublimation of graupel
     logical, pointer :: l_pisub ! sublimation of ice
     logical, pointer :: l_pimlt ! ice melting
     logical, pointer :: l_tidy  ! Tidying
     logical, pointer :: l_tidy2 ! Tidying
  end type process_switch

  type(process_switch), save, target :: pswitch

  real(wp) :: contact_efficiency = 0.0001   ! Arbitrary efficiency for contact nucleation
  real(wp) :: immersion_efficiency = 1.0 ! Arbitrary efficiency for immersion/condensation freezing

  real(wp) :: max_mu = 35.0 ! Maximum value of shape parameter
  real(wp) :: max_mu_frac = 0.75 ! Fraction of maximum value of shape parameter at which psd limitation kicks in
  ! (I l_limit_psd=.true.)
  real(wp) :: fix_mu = 2.5 ! Fixed value for shape parameter (1M/2M)

  logical :: l_inhom_revp = .true.  ! Inhomogeneous evaporation of raindrop

  ! Aerosol-cloud switches
  integer :: aerosol_option = 0 ! Determines how many aerosol modes to use
  !  1: 
  !   soluble_modes(:) = (/ 0, 2, 2/)
  !   insoluble_modes(:) = (/ 0, 0/)
  !  2:
  !   soluble_modes(:) = (/ 0, 2, 2/)
  !   insoluble_modes(:) = (/ 0, 2/)

  logical :: l_separate_rain=.false.   ! Use a separate rain category for active aerosol
  logical :: l_process=.false.         ! process aerosol and/or dust
  logical :: l_passivenumbers=.false.  ! Use passive numbers in aerosol recovery
  logical :: l_passivenumbers_ice=.false.  ! Use passive numbers in ice aerosol recovery
  integer :: process_level = 0
  ! 0: no processing:      l_process=l_passivenumbers=.false.
  ! 1: passive processing: l_process=l_passivenumbers=.true.
  ! 2: full processing:    l_process=.true, l_passivenumbers=.false.
  ! 3: passive processing of ice only: l_process=l_passivenumbers_ice=.true., l_passivenumbers=.false.

  ! aerosol process switches
  logical :: l_aacc=.true.
  logical :: l_ased=.true.
  logical :: l_dsaut=.true.
  logical :: l_aact = .true.   ! activation
  logical :: l_aaut = .true.   ! autoconversion
  logical :: l_aacw = .true.   ! accretion cloud by rain
  logical :: l_aevp = .true.   ! evaporation of cloudd
  logical :: l_asedr = .true.  ! sedimentation of rain
  logical :: l_arevp = .true.  ! evaporation of rain
  logical :: l_asedl = .true.  ! sedimentation of cloud
  logical :: l_atidy = .true.  ! tidy process 1
  logical :: l_atidy2 = .true. ! tidy process 2
  logical :: l_dnuc = .true.   ! ice nucleation
  logical :: l_dsub = .true.   ! sublimation of ice
  logical :: l_dsedi = .true.  ! sedimentation of ice
  logical :: l_dseds = .true.  ! sedimentation of snow
  logical :: l_dsedg = .true.  ! sedimentation of graupel
  logical :: l_dssub = .true.  ! sublimation of snow
  logical :: l_dgsub = .true.  ! sublimation of graupel
  logical :: l_dhomc = .true.  ! homogeneous freezing of cloud
  logical :: l_dhomr = .true.  ! homogeneous freezing of rain
  logical :: l_dimlt = .true.  ! melting of ice
  logical :: l_dsmlt = .true.  ! melting of snow
  logical :: l_dgmlt = .true.  ! melting of graupel
  logical :: l_diacw = .true.  ! riming: ice collects cloud
  logical :: l_dsacw = .true.  ! riming: snow collects cloud
  logical :: l_dgacw = .true.  ! riming: graupel collects cloud
  logical :: l_dsacr = .true.  ! riming: snow collects rain
  logical :: l_dgacr = .true.  ! riming: graupel collects rain
  logical :: l_draci = .true.  ! riming: rain collects ice

  logical :: l_raci_g = .true. ! Allow rain collecting ice to go to graupel if rain mass is significant
  ! Goes to snow otherwise

  logical :: l_onlycollect = .false. ! Only do collection processes and sedimentation

  integer :: i_aerosed_method = 1 ! method to use for aerosol sedimentation
  ! 1= Use hydrometeor mass fallspeeds

  ! cloud activation
  integer :: iopt_act=0
  ! 0=fixed_cloud
  ! 1=90% of total aerosol number
  ! 2=Twomey
  ! 3=Abdul-Razzak and Ghan
  ! 5=Shipway (2015)
  integer, parameter :: iopt_shipway_act = 5

  logical :: l_active_inarg2000 = .false.  ! consider activated aerosol in activation calculation

  integer :: iopt_rcrit = 0 ! method used to calculate rcrit (only use 0)

  ! cloud - rain autoconversion
  integer :: iopt_auto = 1

  ! rain - cloud accretion
  integer :: iopt_accr = 1

  ! ice - nucleation
  integer :: iopt_inuc = 1

  logical :: l_itotsg = .false. ! Consider existing snow and graupel numbers when nucleating ice

  logical :: l_passive = .false. ! only do sedimentation and don't apply conversions

  logical :: l_passive3m = .false. ! Don't use 3rd moment to determine distribution parameters

  logical :: l_limit_psd = .true. ! limit distribution parameters to prevent large mean particle sizes

  logical :: l_nudge_to_cooper = .true. ! nudge the addition of ice number concentration back to the Cooper curve

  ! Some checks are made to ensure consistency, this will
  ! override them (but only to be used by Ben and Adrian!)
  logical :: l_override_checks = .false.

  logical :: l_SM_fix_n0 = .true. ! If running in single moment, then fix n0.
  ! If false, then we fix na and nb (c.f. the LEM microphysics)

  logical :: l_kfsm = .false. ! Switch for Kalli's single moment code.

contains

  subroutine set_mphys_switches(in_option, in_aerosol_option)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    use casim_parent_mod, only: casim_parent, parent_um

    implicit none

    character(len=*), parameter :: RoutineName='SET_MPHYS_SWITCHES'

    integer, intent(in) :: in_option, in_aerosol_option

    integer :: iq,iproc,idgproc ! counters
    logical :: l_on
    real :: k1,k2,k3

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (.not. mphys_is_set) then

      if (in_option == 11111 ) then
        l_kfsm = .true.
      end if
      
      if (casim_parent == parent_um .and. l_sed_eulexp) then 
         !AH - this forces the a long sedimentation step, which will prevent sedimentation 
         !     substepping. These settings assumes the microphysics substep is 2 minutes.
         max_step_length = 120.0
         max_sed_length = 120.0
      endif
     
      call derive_logicals()
      option=in_option
      aerosol_option=in_aerosol_option
      ! Set level of aerosol processing
      select case(process_level)
      case default ! 0
        l_process=.false.
        l_passivenumbers=.false.
        l_passivenumbers_ice=.false.
      case(1)
        l_process=.true.
        l_passivenumbers=.true.
        l_passivenumbers_ice=.true.
      case(2)
        l_process=.true.
        l_passivenumbers=.false.
        l_passivenumbers_ice=.false.
      case(3)
        l_process=.true.
        l_passivenumbers=.false.
        l_passivenumbers_ice=.true.
      case(4)
        l_process=.true.
        l_passivenumbers=.true.
        l_passivenumbers_ice=.false.
      end select

      if (l_warm) l_passivenumbers_ice=.false. ! Override if there's no ice

      select case (option)
      case(1)
        option=23000
      case(2)
        option=12000
      case(3)
        option=11000
      case(4)
        option=13000
      case(5)
        option=22000
      case(6)
        option=22220
      case(7)
        option=22222
      end select
      nq_l=option/10000
      nq_r=(option-nq_l*10000)/1000
      if (.not. l_warm) then
        nq_i=(option-nq_l*10000-nq_r*1000)/100
        nq_s=(option-nq_l*10000-nq_r*1000-nq_i*100)/10
        nq_g=(option-nq_l*10000-nq_r*1000-nq_i*100-nq_s*10)
      end if

      ntotalq=2+nq_l+nq_r+nq_i+nq_s+nq_g

      ! Logicals...
      if (nq_l>1) l_2mc=.true.
      if (nq_r>1) l_2mr=.true.
      if (nq_r>2) l_3mr=.true.
      if (nq_i>1) l_2mi=.true.
      if (nq_s>1) l_2ms=.true.
      if (nq_s>2) l_3ms=.true.
      if (nq_g>1) l_2mg=.true.
      if (nq_g>2) l_3mg=.true.

      ! Names
      allocate(hydro_names(ntotalq))

      !-----------------
      ! Allocate indices
      !-----------------
      iq=0
      call allocq(i_qv, iq, hydro_names, 'qv')
      call allocq(i_th, iq, hydro_names, 'th')

      ! first moments
      call allocq(i_ql, iq, hydro_names, 'ql')
      call allocq(i_qr, iq, hydro_names, 'qr')
      if (.not. l_warm) then
        if (nq_i>0) call allocq(i_qi,iq, hydro_names, 'qi')
        if (nq_s>0) call allocq(i_qs,iq, hydro_names, 'qs')
        if (nq_g>0) call allocq(i_qg,iq, hydro_names, 'qg')
      end if

      ! second moments
      i_nstart=iq+1
      if (l_2mc) call allocq(i_nl,iq, hydro_names, 'nl')
      if (l_2mr) call allocq(i_nr,iq, hydro_names, 'nr')
      if (.not. l_warm) then
        if (l_2mi) call allocq(i_ni,iq, hydro_names, 'ni')
        if (l_2ms) call allocq(i_ns,iq, hydro_names, 'ns')
        if (l_2mg) call allocq(i_ng,iq, hydro_names, 'ng')
      end if

      ! third moments
      i_m3start=iq+1
      if (l_3mr) call allocq(i_m3r,iq, hydro_names, 'm3r')
      if (.not. l_warm) then
        if (l_3ms) call allocq(i_m3s,iq, hydro_names, 'm3s')
        if (l_3mg) call allocq(i_m3g,iq, hydro_names, 'm3g')
      end if

      ! Set params

      if (l_kfsm) then
        ! Kalli's single moment in use, so use his new
        ! parameters
        rain_params = rain_params_kf
        ice_params  = ice_params_kf
      else
        ! Use the original parameters
        rain_params = rain_params_orig
        ice_params  = ice_params_orig
      end if

      cloud_params%i_1m=i_ql
      cloud_params%i_2m=i_nl

      rain_params%i_1m=i_qr
      rain_params%i_2m=i_nr
      rain_params%i_3m=i_m3r

      ice_params%i_1m=i_qi
      ice_params%i_2m=i_ni

      snow_params%i_1m=i_qs
      snow_params%i_2m=i_ns
      snow_params%i_3m=i_m3s

      graupel_params%i_1m=i_qg
      graupel_params%i_2m=i_ng
      graupel_params%i_3m=i_m3g

      ! Mphys distribution parameters set in mphys_parameters
      !if (.not. l_kfsm) then
      !  rain_params%fix_mu = fix_mu
      !  graupel_params%fix_mu = fix_mu
      !end if

      !snow_params%fix_mu=fix_mu
      !graupel_params%fix_mu=fix_mu

      ! Set complexity
      if (l_warm) then
        hydro_complexity%nspecies=2
      else
        hydro_complexity%nspecies=5
      end if

      allocate(hydro_complexity%nmoments(hydro_complexity%nspecies))
      hydro_complexity%nmoments(1)=nq_l
      hydro_complexity%nmoments(2)=nq_r
      if (.not. l_warm) then
        hydro_complexity%nmoments(3)=nq_i
        hydro_complexity%nmoments(4)=nq_s
        hydro_complexity%nmoments(5)=nq_g
      end if

      !-----------------------------------------
      ! Set logicals for 1/2/3 moments in params
      !-----------------------------------------
      if (i_ql > 0)  cloud_params%l_1m=.true.
      if (i_qr > 0)  rain_params%l_1m=.true.
      if (i_qi > 0)  ice_params%l_1m=.true.
      if (i_qs > 0)  snow_params%l_1m=.true.
      if (i_qg > 0)  graupel_params%l_1m=.true.

      if (i_nl > 0)  cloud_params%l_2m=.true.
      if (i_nr > 0)  rain_params%l_2m=.true.
      if (i_m3r > 0) rain_params%l_3m=.true.

      if (i_ni > 0)  ice_params%l_2m=.true.
      if (i_ns > 0)  snow_params%l_2m=.true.
      if (i_m3s > 0) snow_params%l_3m=.true.
      if (i_ng > 0)  graupel_params%l_2m=.true.
      if (i_m3g > 0) graupel_params%l_3m=.true.

      !-------------------------------------
      ! Set thresholds for various qfields
      !-------------------------------------

      allocate(thresh_small(0:ntotalq)) ! Note that 0 index will be written to
      k3=p2-p1
      k1=(p3-p2)/k3
      k2=(p1-p3)/k3
      ! if a variable is not used
      thresh_small(i_qv)=qv_small
      thresh_small(i_th)=th_small
      thresh_small(i_ql)=ql_small
      thresh_small(i_qr)=qr_small
      thresh_small(i_nl)=nl_small
      thresh_small(i_nr)=nr_small
      thresh_small(i_m3r)=exp(-k1*log(qr_small)-k2*log(nr_small))
      thresh_small(i_qi)=qi_small
      thresh_small(i_ni)=ni_small
      thresh_small(i_qs)=qs_small
      thresh_small(i_ns)=ns_small
      thresh_small(i_m3s)=exp(-k1*log(qs_small)-k2*log(ns_small))
      thresh_small(i_qg)=qg_small
      thresh_small(i_ng)=ng_small
      thresh_small(i_m3g)=exp(-k1*log(qg_small)-k2*log(ng_small))

      allocate(thresh_tidy(0:ntotalq)) ! Note that 0 index will be written to
      ! if a variable is not used
      thresh_tidy(i_qv)=qv_tidy
      thresh_tidy(i_th)=th_tidy
      thresh_tidy(i_ql)=ql_tidy
      thresh_tidy(i_qr)=qr_tidy
      thresh_tidy(i_nl)=nl_tidy
      thresh_tidy(i_nr)=nr_tidy
      thresh_tidy(i_m3r)=exp(-2.0*k1*log(qr_tidy)-2.0*k2*log(nr_tidy))
      thresh_tidy(i_qi)=qi_tidy
      thresh_tidy(i_qs)=qs_tidy
      thresh_tidy(i_ni)=ni_tidy
      thresh_tidy(i_ns)=ns_tidy
      thresh_tidy(i_m3s)=exp(-2.0*k1*log(qs_tidy)-2.0*k2*log(ns_tidy))
      thresh_tidy(i_qg)=qg_tidy
      thresh_tidy(i_ng)=ng_tidy
      thresh_tidy(i_m3g)=exp(-2.0*k1*log(qs_tidy)-2.0*k2*log(ns_tidy))

      allocate(thresh_sig(0:ntotalq)) ! Note that 0 index will be written to
      thresh_sig(i_ql)=ql_sig
      thresh_sig(i_qr)=qr_sig
      thresh_sig(i_qs)=qs_sig
      thresh_sig(i_qi)=qi_sig
      thresh_sig(i_qg)=qg_sig

      allocate(thresh_large(0:ntotalq)) ! Note that 0 index will be written to
      thresh_large(i_ql)=ql_large
      thresh_large(i_qr)=qr_large
      thresh_large(i_qs)=qs_large
      thresh_large(i_qi)=qi_large
      thresh_large(i_qg)=qg_large
      thresh_large(i_nl)=nl_large
      thresh_large(i_nr)=nr_large
      thresh_large(i_ns)=ns_large
      thresh_large(i_ni)=ni_large
      thresh_large(i_ng)=ng_large
      thresh_large(i_m3r)=exp(-k1*log(qr_large)-k2*log(nr_large))
      thresh_large(i_m3s)=exp(-k1*log(qs_large)-k2*log(ns_large))
      thresh_large(i_m3g)=exp(-k1*log(qg_large)-k2*log(ng_large))

      select case(aerosol_option)
      case default
        ! Use defaults or assign directly through namelists
      case (1)
        soluble_modes(:)=(/ 0, 2, 2/)
        insoluble_modes(:)=(/ 0, 0/)
      case (2)
        soluble_modes(:)=(/ 2, 2, 2/)
        insoluble_modes(:)=(/ 0, 0/)
      case (3)
        soluble_modes(:)=(/ 2, 2, 2/)
        insoluble_modes(:)=(/ 0, 2/)
      case (4)
        soluble_modes(:)=(/ 2, 2, 2/)
        insoluble_modes(:)=(/ 0, 2/)
      end select

      ! NB active_X switches have been made obsolete by
      ! l_process, l_passivenumbers and l_separate_rain switches
      ! so should gradually be removed
      if (l_process) then
        if (l_warm)then
          active_cloud(:)=(/.true., .false./)
          active_ice(:)=(/.false., .false./)
        else
          active_cloud(:)=(/.true., .true./)
          active_ice(:)=(/.true., .true./)
        end if
        if (l_separate_rain) active_rain(:)=(/.true./)
        if (l_passivenumbers) active_number(1)=.true.
        if (l_passivenumbers_ice) active_number(2)=.true.
      end if

      nactivea=count(active_cloud(isol:isol))+count(active_ice(isol:isol))+count(active_rain)+count(active_number(isol:isol))
      nactived=count(active_cloud(iinsol:iinsol))+count(active_ice(iinsol:iinsol))+count(active_number(iinsol:iinsol))
      ! Increase ntotala by sum(soluble_modes) in Jan 2019 to make space for
      ! hygroscopicity properties from UKCA
      ntotala=2*sum(soluble_modes)+sum(insoluble_modes)+count(active_cloud)+&
           count(active_ice)+count(active_rain)+count(active_number)
      allocate(aero_names(ntotala))

      !-----------------
      ! Allocate indices
      !-----------------
      iq=0
      if (soluble_modes(1) > 1) call alloca(i_am1,iq,aero_names, 'AitkenSolMass', aero_index, i_ccn=imass)
      if (soluble_modes(1) > 0) call alloca(i_an1,iq,aero_names, 'AitkenSolNumber', aero_index, i_ccn=inumber)
      if (soluble_modes(2) > 1) call alloca(i_am2,iq,aero_names, 'AccumSolMass', aero_index, i_ccn=imass)
      if (soluble_modes(2) > 0) call alloca(i_an2,iq,aero_names, 'AccumSolNumber', aero_index, i_ccn=inumber)
      if (soluble_modes(3) > 1) call alloca(i_am3,iq,aero_names, 'CoarseSolMass', aero_index, i_ccn=imass)
      if (soluble_modes(3) > 0) call alloca(i_an3,iq,aero_names, 'CoarseSolNumber', aero_index, i_ccn=inumber)
      if (active_cloud(isol)) call alloca(i_am4,iq,aero_names, 'ActiveSolCloud', aero_index)
      if (active_rain(isol)) call alloca(i_am5,iq,aero_names, 'ActiveSolRain', aero_index)
      if (insoluble_modes(2) > 1) call alloca(i_am6,iq,aero_names, 'CoarseInsolMass', aero_index, i_in=imass)
      if (insoluble_modes(2) > 0) call alloca(i_an6,iq,aero_names, 'CoarseInsolNumber', aero_index, i_in=inumber)
      if (active_ice(iinsol)) call alloca(i_am7,iq,aero_names, 'ActiveInsolIce', aero_index)
      if (active_ice(isol)) call alloca(i_am8,iq,aero_names, 'ActiveSolIce', aero_index)
      if (active_cloud(iinsol)) call alloca(i_am9,iq,aero_names, 'ActiveInsolCloud', aero_index)
      if (insoluble_modes(1) > 1) call alloca(i_am10,iq,aero_names, 'AccumInsolMass', aero_index, i_in=imass)
      if (insoluble_modes(1) > 0) call alloca(i_an10,iq,aero_names, 'AccumInsolNumber', aero_index, i_in=inumber)
      if (active_number(isol)) call alloca(i_an11,iq,aero_names, 'ActiveSolNumber', aero_index)
      if (active_number(iinsol)) call alloca(i_an12,iq,aero_names, 'ActiveInsolNumber', aero_index)
      ! The order of these matters, safest to add new ones at the
      ! end. This doesn't do anything except assign numbers to
      ! ak, and the numbers should be 18,19,20
      if (soluble_modes(1) > 0) call alloca(i_ak1,iq,aero_names,'AitkenSolBk',aero_index,i_ccn=ihyg)
      if (soluble_modes(2) > 0) call alloca(i_ak2,iq,aero_names,'AccumSolBk', aero_index,i_ccn=ihyg)
      if (soluble_modes(3) > 0) call alloca(i_ak3,iq,aero_names,'CoarseSolBk',aero_index,i_ccn=ihyg)

      aero_complexity%nspecies=count(soluble_modes > 0)+count(insoluble_modes > 0)+count(active_cloud)&
           +count(active_ice)+count(active_rain)+count(active_number)

      if (soluble_modes(1) > 0) aero_index%i_aitken=1
      if (soluble_modes(2) > 0) aero_index%i_accum=aero_index%i_aitken+1
      if (soluble_modes(3) > 0) aero_index%i_coarse=aero_index%i_accum+1
      if (insoluble_modes(1) > 0) aero_index%i_accum_dust=1
      if (insoluble_modes(2) > 0) aero_index%i_coarse_dust=aero_index%i_accum_dust+1

      if (l_process) then
        allocate(thresh_atidy(0:ntotala)) ! Note that 0 index will be written to
        ! if a variable is not used
        thresh_atidy(i_am1)=aeromass_small
        thresh_atidy(i_an1)=aeronumber_small
        thresh_atidy(i_am2)=aeromass_small
        thresh_atidy(i_an2)=aeronumber_small
        thresh_atidy(i_am3)=aeromass_small
        thresh_atidy(i_an3)=aeronumber_small
        thresh_atidy(i_am4)=aeromass_small
        thresh_atidy(i_am5)=aeromass_small
        thresh_atidy(i_am6)=aeromass_small
        thresh_atidy(i_an6)=aeronumber_small
        thresh_atidy(i_am7)=aeromass_small
        thresh_atidy(i_am8)=aeromass_small
        thresh_atidy(i_am9)=aeromass_small
        thresh_atidy(i_am10)=aeromass_small
        thresh_atidy(i_an10)=aeronumber_small
        thresh_atidy(i_an11)=aeronumber_small
        thresh_atidy(i_an12)=aeronumber_small
      end if

      !--------------------------
      ! Set the process choices - THESE NEED TIDYING WITH LOGICALS
      !--------------------------
      iopt_accr=1
      iopt_auto=1

      !-----------------------------------------------------
      ! Automatically set cloud sedimentation ON with Paul's
      ! new cloud fraction scheme.
      ! Do this here, so it is ahead of allocation of process
      ! indices etc.
      !-----------------------------------------------------

      if (l_prf_cfrac) then
        l_psedl = .true.
        l_psedi = .true.
      end if

      !--------------------------
      ! Allocate process indices
      !--------------------------
      iproc=0
      idgproc=0
      ! AH - Allocate all warm processes even if procs are not used
      call allocp(i_cond, iproc, idgproc, 'pcond')
      call allocp(i_praut, iproc, idgproc, 'praut')
      call allocp(i_pracw, iproc, idgproc, 'pracw')
      call allocp(i_pracr, iproc, idgproc, 'pracr')
      call allocp(i_prevp, iproc, idgproc, 'prevp')
      call allocp(i_psedl, iproc, idgproc, 'psedl')
      call allocp(i_psedr, iproc, idgproc, 'psedr')
      call allocp(i_tidy, iproc, idgproc, 'ptidy')
      call allocp(i_tidy2, iproc, idgproc, 'ptidy2')
      if (.not. l_warm) then
        call allocp(i_inuc, iproc, idgproc, 'pinuc')
        call allocp(i_idep, iproc, idgproc, 'pidep')
        call allocp(i_iacw, iproc, idgproc, 'piacw')
        call allocp(i_saut, iproc, idgproc, 'psaut')
        call allocp(i_sdep, iproc, idgproc, 'psdep')
        call allocp(i_sacw, iproc, idgproc, 'psacw')
        call allocp(i_gdep, iproc, idgproc, 'pgdep')
        call allocp(i_pseds, iproc, idgproc, 'pseds')
        call allocp(i_psedi, iproc, idgproc, 'psedi')
        call allocp(i_psedg, iproc, idgproc, 'psedg')
        call allocp(i_saci, iproc, idgproc, 'psaci')
        call allocp(i_raci, iproc, idgproc, 'praci')
        call allocp(i_sacr, iproc, idgproc, 'psacr')
        call allocp(i_gacr, iproc, idgproc, 'pgacr')
        call allocp(i_gacw, iproc, idgproc, 'pgacw')
        call allocp(i_gaci, iproc, idgproc, 'pgaci')
        call allocp(i_gacs, iproc, idgproc, 'pgacs')
        call allocp(i_iagg, iproc, idgproc, 'piagg')
        call allocp(i_sagg, iproc, idgproc, 'psagg')
        call allocp(i_gagg, iproc, idgproc, 'pgagg')
        call allocp(i_sbrk, iproc, idgproc, 'psbrk')
        call allocp(i_gshd, iproc, idgproc, 'pgshd')
        call allocp(i_ihal, iproc, idgproc, 'pihal')
        call allocp(i_smlt, iproc, idgproc, 'psmlt')
        call allocp(i_gmlt, iproc, idgproc, 'pgmlt')
        call allocp(i_homr, iproc, idgproc, 'phomr')
        call allocp(i_homc, iproc, idgproc, 'phomc')
        call allocp(i_ssub, iproc, idgproc, 'pssub')
        call allocp(i_gsub, iproc, idgproc, 'pgsub')
        call allocp(i_isub, iproc, idgproc, 'pisub')
        call allocp(i_imlt, iproc, idgproc, 'pimlt')
        call allocp(i_iics, iproc, idgproc, 'piics')
        call allocp(i_idps, iproc, idgproc, 'pidps')
      end if
      hydro_complexity%nprocesses=iproc

      ! allocate process indices for aerosols...
      iproc=0
      idgproc=idgproc
      if (l_process) then
        if (l_aact) call allocp(i_aact, iproc, idgproc, 'aact')
        if (l_aaut) call allocp(i_aaut, iproc, idgproc, 'aaut')
        if (l_aacw) call allocp(i_aacw, iproc, idgproc, 'aacw')
        if (l_aevp) call allocp(i_aevp, iproc, idgproc, 'aevp')
        if (l_asedr) call allocp(i_asedr, iproc, idgproc, 'asedr')
        if (l_arevp) call allocp(i_arevp, iproc, idgproc, 'arevp')
        if (l_asedl) call allocp(i_asedl, iproc, idgproc, 'asedl')
        !... additional tidying processes (Need to sort these)
        if (l_atidy) call allocp(i_atidy, iproc, idgproc, 'atidy')
        if (l_atidy2) call allocp(i_atidy2, iproc, idgproc, 'atidy2')
        !... ice related processes
        l_on= .not. l_warm
        if (l_dnuc) call allocp(i_dnuc, iproc, idgproc, 'dnuc', l_onoff=l_on)
        if (l_dsub) call allocp(i_dsub, iproc, idgproc, 'dsub', l_onoff=l_on)
        if (l_dsedi) call allocp(i_dsedi, iproc, idgproc, 'dsedi', l_onoff=l_on)
        if (l_dseds) call allocp(i_dseds, iproc, idgproc, 'dseds', l_onoff=l_on)
        if (l_dsedg) call allocp(i_dsedg, iproc, idgproc, 'dsedg', l_onoff=l_on)
        if (l_dssub) call allocp(i_dssub, iproc, idgproc, 'dssub', l_onoff=l_on)
        if (l_dgsub) call allocp(i_dgsub, iproc, idgproc, 'dgsub', l_onoff=l_on)
        if (l_dhomc) call allocp(i_dhomc, iproc, idgproc, 'dhomc', l_onoff=l_on)
        if (l_dhomr) call allocp(i_dhomr, iproc, idgproc, 'dhomr', l_onoff=l_on)
        if (l_dimlt) call allocp(i_dimlt, iproc, idgproc, 'dimlt', l_onoff=l_on)
        if (l_dsmlt) call allocp(i_dsmlt, iproc, idgproc, 'dsmlt', l_onoff=l_on)
        if (l_dgmlt) call allocp(i_dgmlt, iproc, idgproc, 'dgmlt', l_onoff=l_on)
        if (l_diacw) call allocp(i_diacw, iproc, idgproc, 'diacw', l_onoff=l_on)
        if (l_dsacw) call allocp(i_dsacw, iproc, idgproc, 'dsacw', l_onoff=l_on)
        if (l_dgacw) call allocp(i_dgacw, iproc, idgproc, 'dgacw', l_onoff=l_on)
        if (l_dsacr) call allocp(i_dsacr, iproc, idgproc, 'dsacr', l_onoff=l_on)
        if (l_dgacr) call allocp(i_dgacr, iproc, idgproc, 'dgacr', l_onoff=l_on)
        if (l_draci) call allocp(i_draci, iproc, idgproc, 'draci', l_onoff=l_on)
      end if
      aero_complexity%nprocesses = iproc

      ! Set switches for having ice breakup and droplet shattering
      if (l_sip_icebreakup) then
        l_pgaci = .true.
        l_pgacs = .true.
        l_piics = .true.
        l_no_pgacs_in_sumprocs = .true.
      end if

      if (l_sip_dropletshatter) then
        l_phomr = .true.
        l_pidps = .true.
      end if

      !----------------------------------------------------
      ! Ensure we only do this at the start of the run
      !----------------------------------------------------
      mphys_is_set=.true.
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine set_mphys_switches

  ! Allocate an index to a q variable
  subroutine allocq(i, iq, names, p_name)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ALLOCQ'
    
    integer, intent(out) :: i
    integer, intent(inout) :: iq
    character(10), intent(inout) :: names(:)
    character(*), intent(in) :: p_name

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    i=iq+1
    iq=i
    names(iq)=adjustr(trim(p_name))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine allocq

  subroutine allocp(proc, iproc, idgproc, p_name, l_onoff)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ALLOCP'
    
    type(process_name), intent(inout) :: proc
    integer, intent(inout) :: iproc, idgproc
    character(*), intent(in) :: p_name
    logical, optional, intent(in) :: l_onoff

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    iproc=iproc+1
    idgproc=idgproc+1
    proc%id=iproc
    !    proc%unique_id=iproc
    proc%p_name=adjustr(trim(p_name))
    if (present(l_onoff)) then
      proc%on=l_onoff
    else
      proc%on=.true.
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine allocp

  ! Allocate an index to a aerosol variable
  ! similar to allocq, but also adds information
  ! to aero_index if variable should act as an in
  ! or a ccn.
  subroutine alloca(i, iq, names, p_name, aero_index, i_ccn, i_in)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ALLOCA'
 
    integer, intent(out) :: i
    integer, intent(inout) :: iq
    character(20), intent(inout) :: names(:)
    character(*), intent(in) :: p_name
    type(aerosol_index), intent(inout) :: aero_index
    ! if present and >0 should be set to imass or inumber or ihyg
    integer, optional, intent(in) :: i_ccn, i_in

    integer :: is_ccn, is_in, nin, nccn

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    is_ccn=0
    is_in=0
    if (present(i_ccn)) is_ccn=i_ccn
    if (present(i_in)) is_in=i_in

    i=iq+1
    iq=i
    names(iq)=adjustr(trim(p_name))
    if (is_ccn == imass) then
      ! represents mass of aerosol which can act as ccn
      nccn=count(aero_index%ccn_m > 0)+1
      aero_index%ccn_m(nccn)=iq
      aero_index%nccn=nccn
    else if (is_ccn == inumber) then
      ! represents number of aerosol which can act as ccn
      nccn=count(aero_index%ccn_n > 0)+1
      aero_index%ccn_n(nccn)=iq
      aero_index%nccn=nccn
    else if (is_ccn == ihyg) then
      nccn=count(aero_index%ccn_k > 0)+1
      aero_index%ccn_k(nccn)=iq
      aero_index%nccn=nccn
    end if
    if (is_in == imass) then
      ! represents mass of aerosol which can act as in
      nin=count(aero_index%in_m > 0)+1
      aero_index%in_m(nin)=iq
      aero_index%nin=nin
    else if (is_in == inumber) then
      ! represents number of aerosol which can act as in
      nin=count(aero_index%in_n > 0)+1
      aero_index%in_n(nin)=iq
      aero_index%nin=nin
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine alloca

  ! Routine sets logical switches which depend on or are overridden by other switches
  ! Transfer namelist logicals to derived type for ease of use later
  subroutine derive_logicals()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='DERIVE_LOGICALS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (l_kfsm) then
      ! Kalli's single moment scheme in use, so we need to switch off
      ! all of the process rates which are not required.

      l_psaut = .false.
      l_psdep = .false.
      l_psacw = .false.
      l_pseds = .false.
      l_psaci = .false.
      l_psacr = .false.
      l_pgacs = .false.
      l_psagg = .false.
      l_psbrk = .false.
      l_pihal = .false.
      l_piics = .false.
      l_pidps = .false.
      l_psmlt = .false.
      l_pssub = .false.

      ! But we need to turn on ice and snow sedimentation as these
      ! are on in his runs.

      l_psedl = .true.
      l_psedi = .true.

    end if ! l_kfsm

    pswitch%l_pcond=>l_pcond ! Condensation
    pswitch%l_praut=>l_praut ! Autoconversion cloud -> rain
    pswitch%l_pracw=>l_pracw ! Accretion  cloud -> rain
    pswitch%l_pracr=>l_pracr ! aggregation of rain drops
    pswitch%l_prevp=>l_prevp ! evaporation of rain
    pswitch%l_psedl=>l_psedl ! sedimentation of cloud
    pswitch%l_psedr=>l_psedr ! sedimentation of rain
    pswitch%l_pinuc=>l_pinuc ! ice nucleation
    pswitch%l_pidep=>l_pidep ! ice deposition
    pswitch%l_piacw=>l_piacw ! ice accreting water
    pswitch%l_psaut=>l_psaut ! ice autoconversion ice -> snow
    pswitch%l_psdep=>l_psdep ! vapour deposition onto snow
    pswitch%l_psacw=>l_psacw ! snow accreting water
    pswitch%l_pgdep=>l_pgdep ! vapour deposition onto graupel
    pswitch%l_pseds=>l_pseds ! snow sedimentation
    pswitch%l_psedi=>l_psedi ! ice sedimentation
    pswitch%l_psedg=>l_psedg ! graupel sedimentation
    pswitch%l_psaci=>l_psaci ! snow accreting ice
    pswitch%l_praci=>l_praci ! rain accreting ice
    pswitch%l_psacr=>l_psacr ! snow accreting rain
    pswitch%l_pgacr=>l_pgacr ! graupel accreting rain
    pswitch%l_pgacw=>l_pgacw ! graupel accreting cloud water
    pswitch%l_pgaci=>l_pgaci ! graupel accreting ice
    pswitch%l_pgacs=>l_pgacs ! graupel accreting snow
    pswitch%l_piagg=>l_piagg ! aggregation of ice particles
    pswitch%l_psagg=>l_psagg ! aggregation of snow particles
    pswitch%l_pgagg=>l_pgagg ! aggregation of graupel particles
    pswitch%l_psbrk=>l_psbrk ! break up of snow flakes
    pswitch%l_pgshd=>l_pgshd ! shedding of liquid from graupel
    pswitch%l_pihal=>l_pihal ! hallet mossop
    pswitch%l_piics=>l_piics ! ice-ice collision
    pswitch%l_pidps=>l_pidps ! droplet shattering
    pswitch%l_psmlt=>l_psmlt ! snow melting
    pswitch%l_pgmlt=>l_pgmlt ! graupel melting
    pswitch%l_phomr=>l_phomr ! homogeneous freezing of rain
    pswitch%l_phomc=>l_phomc ! homogeneous freezing of cloud droplets
    pswitch%l_pssub=>l_pssub ! sublimation of snow
    pswitch%l_pgsub=>l_pgsub ! sublimation of graupel
    pswitch%l_pisub=>l_pisub ! sublimation of ice
    pswitch%l_pimlt=>l_pimlt ! ice melting
    pswitch%l_tidy =>l_ptidy ! Tidying
    pswitch%l_tidy2=>l_ptidy2 ! Tidying

    if (.not. l_rain) then
      pswitch%l_praut=.false.
      pswitch%l_pracw=.false.
      pswitch%l_pracr=.false.
      pswitch%l_prevp=.false.
    end if

    if (l_onlycollect) then
      pswitch%l_praut=.false.
      pswitch%l_prevp=.false.
      pswitch%l_psaut=.false.
      pswitch%l_pihal=.false.
      pswitch%l_phomr=.false.
      pswitch%l_pinuc=.false.
      pswitch%l_pidep=.false.
      pswitch%l_psdep=.false.
      pswitch%l_pgdep=.false.
      pswitch%l_psmlt=.false.
      pswitch%l_pgmlt=.false.
      pswitch%l_pimlt=.false.
      pswitch%l_pcond=.false.
    end if

    if (.not. l_sg) then
      pswitch%l_psaut=.false.
      pswitch%l_psacw=.false.
      pswitch%l_psaci=.false.
      pswitch%l_praci=.false.
      pswitch%l_psacr=.false.
      pswitch%l_pgacw=.false.
      pswitch%l_pgacr=.false.
      pswitch%l_psagg=.false.
      pswitch%l_psdep=.false.
      pswitch%l_psmlt=.false.
      pswitch%l_pseds=.false.
      pswitch%l_pihal=.false.
      pswitch%l_piics=.false.
      pswitch%l_pidps=.false.
      l_g=.false.
    end if

    if (.not. l_g) then
      pswitch%l_pgacw=.false.
      pswitch%l_pgacr=.false.
      pswitch%l_pgdep=.false.
      pswitch%l_pgmlt=.false.
      pswitch%l_psedg=.false.
      pswitch%l_pgacs=.false.
      pswitch%l_pgaci=.false.
      pswitch%l_pgshd=.false.
      ! this should not be switched off
      pswitch%l_phomr=.false.
      pswitch%l_piics=.false.
      pswitch%l_pidps=.false.
    end if

    if (.not. l_halletmossop) pswitch%l_pihal=.false.
    if (.not. l_sip_icebreakup) then
      pswitch%l_pgacs=.false.
      pswitch%l_piics=.false.
    end if
    if (.not. l_sip_dropletshatter) pswitch%l_pidps=.false.

    if (.not. (pswitch%l_pidep .or. pswitch%l_psdep .or. pswitch%l_pgdep)) l_idep=.false.
    if (.not. (pswitch%l_pisub .or. pswitch%l_pssub .or. pswitch%l_pgsub)) l_isub=.false.
    if (.not. (pswitch%l_praut .or. pswitch%l_pracw .or. pswitch%l_piacw .or. &
         pswitch%l_psacw .or. pswitch%l_pgacw .or. pswitch%l_phomc)) l_pos1=.false.
    if (.not. (pswitch%l_praci .or. pswitch%l_psaci .or. pswitch%l_pgaci .or. &
         pswitch%l_psaut .or. pswitch%l_pisub .or. pswitch%l_pimlt)) l_pos2=.false.
    if (.not. (pswitch%l_prevp .or. pswitch%l_psacr .or. pswitch%l_pgacr .or. pswitch%l_phomr)) l_pos3=.false.
    if (.not. (pswitch%l_pgacs .or. pswitch%l_psmlt .or. pswitch%l_psacr .or. pswitch%l_pssub)) l_pos4=.false.
    if (.not. (pswitch%l_praut .or. pswitch%l_pracw)) l_pos5=.false.
    if (.not. (pswitch%l_prevp)) l_pos6=.false.

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine derive_logicals
end module mphys_switches
