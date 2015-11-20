MODULE mphys_switches

  USE variable_precision, ONLY: wp
  USE mphys_parameters, ONLY: cloud_params, rain_params, ice_params   &
     , snow_params, graupel_params, naero, p1,p2,p3
  USE thresholds, ONLY: thresh_small, th_small, qv_small   &
     , ql_small, qr_small, nl_small, nr_small, m3r_small &
     , qi_small, qs_small, ni_small, ns_small, m3s_small &
     , qg_small, ng_small, m3g_small &
     , thresh_tidy, th_tidy, qv_tidy &
     , ql_tidy, qr_tidy, nl_tidy, nr_tidy, m3r_tidy  &
     , qi_tidy, qs_tidy, ni_tidy, ns_tidy, m3s_tidy  &
     , qg_tidy, ng_tidy, m3g_tidy &
     , thresh_sig, qs_sig, ql_sig, qr_sig, qi_sig, qg_sig &
     , thresh_large, qs_large, ql_large, qr_large, qi_large &
     , qg_large, ns_large, nl_large, nr_large, ni_large &
     , ng_large, thresh_atidy, aeromass_small, aeronumber_small

  USE process_routines 

  IMPLICIT NONE

  LOGICAL :: mphys_is_set = .FALSE.

  ! Microphysics option are now specified through a 5 digit integer, such that
  ! each digit represents the number of moments of cloud, rain, ice, snow, graupel
  ! respectively.  E.g. 23000 would use double moment cloud and triple moment rain
  INTEGER :: option = 22222

  !-------------------------------------------------
  ! Switches to determine what aerosol we should use
  ! (replaces old aerosol_option)
  !-------------------------------------------------
  ! number of moments to use for soluble modes
  ! (/ aitken, accumulation, coarse /)
  INTEGER :: soluble_modes(3) = (/ 0, 0, 0/)
  ! number of moments to use for insoluble modes
  ! (/ accumulation, coarse /)
  INTEGER :: insoluble_modes(2) = (/ 0, 0/)
  ! do we want to carry activated aersols in cloud/rain?
  ! (/ for soluble aerosol, for insoluble aerosol /)
  LOGICAL :: active_cloud(2) = (/.FALSE., .FALSE./)
  ! do we want to carry activated aersols in ice/snow/graupel?
  ! (/ for soluble aerosol, for insoluble aerosol /)
  LOGICAL :: active_ice(2) = (/.FALSE., .FALSE./)
  ! do we want to carry a separate activated category in rain?
  ! (/ for soluble aerosol /)
  LOGICAL :: active_rain(1) = (/.FALSE./)
  ! do we want to retain information amout activated number
  ! so that we can use this to replenish aerosol (rather than
  ! using the processed hydrometeor numbers)
  ! (/ for soluble aerosol, for insoluble aerosol /)
  LOGICAL :: active_number(2) = (/.FALSE., .FALSE./)


  LOGICAL :: l_warm = .FALSE.  ! Only use warm rain microphysics

  CHARACTER(10), ALLOCATABLE :: hydro_names(:)
  CHARACTER(20), ALLOCATABLE :: aero_names(:)

  ! index to determine if something is a mass or a number
  INTEGER :: inumber=1
  INTEGER :: imass=2
  ! index to determine if something is a soluble(ccn) or insoluble(IN) (NB currently can't be both)
  INTEGER :: isol=1
  INTEGER :: iinsol=2

  ! standard hydrometeor indices
  INTEGER :: i_qv  = 0 ! water vapour
  INTEGER :: i_ql  = 0 ! cloud liquid mass
  INTEGER :: i_nl  = 0 ! cloud liquid number
  INTEGER :: i_qr  = 0 ! rain mass
  INTEGER :: i_nr  = 0 ! rain number
  INTEGER :: i_m3r = 0 ! rain 3rd moment
  INTEGER :: i_qi  = 0 ! ice mass
  INTEGER :: i_ni  = 0 ! ice number
  INTEGER :: i_qs  = 0 ! snow mass
  INTEGER :: i_ns  = 0 ! snow number
  INTEGER :: i_m3s = 0 ! snow 3rd moment
  INTEGER :: i_qg  = 0 ! graupel mass
  INTEGER :: i_ng  = 0 ! graupel number
  INTEGER :: i_m3g = 0 ! graupel 3rd moment
  INTEGER :: i_qh  = 0 ! hail mass
  INTEGER :: i_nh  = 0 ! hail number
  INTEGER :: i_m3h = 0 ! graupel 3rd moment
  INTEGER :: i_th  = 0 ! theta

  ! information about location indices of different moments
  INTEGER :: i_qstart = 1 ! First index in for q variables (always 1)
  INTEGER :: i_nstart     ! First index in for n variables
  INTEGER :: i_m3start    ! First index in for m3 variables
  INTEGER :: i_hstart = 3 ! First index in for hydrometeors (always 3)

  INTEGER :: ntotalq   ! total number of q variables (includes theta and qv)
  INTEGER :: ntotala   ! total number of aerosol variables

  INTEGER :: nactivea  ! total number of active aerosol(soluble) variables
  INTEGER :: nactived  ! total number of active aerosol(insoluble) variables

  ! number of moments for each hydrometeor
  INTEGER :: nq_l = 0  ! cloud liquid
  INTEGER :: nq_r = 0  ! rain
  INTEGER :: nq_i = 0  ! cloud ice
  INTEGER :: nq_s = 0  ! snow
  INTEGER :: nq_g = 0  ! graupel

  ! logicals for moments
  LOGICAL :: l_2mc = .FALSE.
  LOGICAL :: l_2mr = .FALSE.
  LOGICAL :: l_3mr = .FALSE.

  LOGICAL :: l_2mi = .FALSE.
  LOGICAL :: l_2ms = .FALSE.
  LOGICAL :: l_3ms = .FALSE.
  LOGICAL :: l_2mg = .FALSE.
  LOGICAL :: l_3mg = .FALSE.

  ! aerosol indices (soluble)
  INTEGER :: i_am1 = 0    ! aitken aerosol mass
  INTEGER :: i_an1 = 0    ! aitken aerosol number
  INTEGER :: i_am2 = 0    ! accumulation aerosol mass
  INTEGER :: i_an2 = 0    ! accumulation aerosol number
  INTEGER :: i_am3 = 0    ! coarse aerosol mass
  INTEGER :: i_an3 = 0    ! coarse aerosol number
  INTEGER :: i_am4 = 0    ! activated aerosol mass
  INTEGER :: i_am5 = 0    ! activated aerosol mass (in rain)

  ! aerosol indices (insoluble)
  INTEGER :: i_am6 = 0    ! dust mass
  INTEGER :: i_an6 = 0    ! dust number
  INTEGER :: i_am7 = 0    ! activated dust mass

  ! aerosol indices (soluble/insoluble)
  INTEGER :: i_am8 = 0    ! soluble mass in ice (incl. snow/graupel)
  INTEGER :: i_am9 = 0    ! dust mass in water (cloud/rain)

  ! aerosol indices (accumulation mode insoluble)
  INTEGER :: i_am10 = 0    ! dust mass (accum)
  INTEGER :: i_an10 = 0    ! dust number (accum)

  ! aerosol indices (activated number)
  INTEGER :: i_an11 = 0    ! soluble
  INTEGER :: i_an12 = 0    ! insoluble

  TYPE :: complexity
     INTEGER :: nspecies
     INTEGER, POINTER :: nmoments(:)
     INTEGER :: nprocesses = 0
  END TYPE complexity

  INTEGER, PARAMETER :: maxmodes=3  !< maximum ccn/in modes

  TYPE :: aerosol_index              !< contains indexing imformation for aerosol as ccn/in
     INTEGER :: nccn            = 0  !< How many aerosol act as ccn
     INTEGER :: nin             = 0  !< How many aerosol act as in
     INTEGER :: ccn_m(maxmodes) = 0  !< list of ccn mass indices
     INTEGER :: ccn_n(maxmodes) = 0  !< list of ccn number indices
     INTEGER :: in_m(maxmodes)  = 0  !< list of in mass indices
     INTEGER :: in_n(maxmodes)  = 0  !< list of in number indices
     ! Some useful indices for the different modes
     INTEGER :: i_aitken = 0     ! atiken mode
     INTEGER :: i_accum  = 0     ! accumulation mode
     INTEGER :: i_coarse = 0     ! Coarse mode
     INTEGER :: i_accum_dust  = 0     ! accumulation mode dust
     INTEGER :: i_coarse_dust = 0     ! Coarse mode dust
  END TYPE aerosol_index

  TYPE(complexity), SAVE :: hydro_complexity
  TYPE(complexity), SAVE :: aero_complexity

  TYPE(aerosol_index), SAVE :: aero_index

  ! substepping
  REAL(wp) :: max_step_length = 10.0
  REAL(wp) :: max_sed_length = 2.0

  INTEGER :: nsubsteps, nsubseds

  ! process switches
  ! Some of these switches are obsolete or inactive - review these
  ! and then remove this message

  LOGICAL :: l_abelshipway=.FALSE.
  LOGICAL :: l_sed_3mdiff=.FALSE.
  LOGICAL :: l_cons=.FALSE.
  LOGICAL :: l_inuc=.TRUE.

  LOGICAL :: l_sg             = .TRUE.  ! run with snow and graupel
  LOGICAL :: l_g              = .TRUE.  ! run with graupel
  LOGICAL :: l_halletmossop   = .TRUE.

  LOGICAL :: l_harrington     = .FALSE.  ! Use Jerry's method for ice autoconvertion

  LOGICAL :: l_condensation   = .TRUE.  ! condense water
  LOGICAL :: l_evaporation    = .TRUE.  ! evaporate rain
  LOGICAL :: l_rain           = .TRUE.  ! rain sources
  LOGICAL :: l_sed            = .TRUE.  ! sedimentation
  LOGICAL :: l_sedl           = .TRUE.  ! sedimentation
  LOGICAL :: l_boussinesq     = .FALSE. ! set rho=1 everywhere
  INTEGER :: diag_mu_option   = -999    ! select diagnostic mu           &

  LOGICAL :: l_iaut           = .TRUE.  ! autoconversion of ice/snow
  LOGICAL :: l_imelt          = .TRUE.  ! melting of ice/snow/graupel
  LOGICAL :: l_iacw           = .TRUE.  ! ice/snow/graupel accumulation cloud water
  LOGICAL :: l_idep           = .TRUE.  ! deposition of ice/snow/graupel
  LOGICAL :: l_isub           = .TRUE.  ! sublimation of ice/snow/graupel

  LOGICAL :: l_pos1           = .TRUE.   ! switch on positivity check
  LOGICAL :: l_pos2           = .TRUE.   ! switch on positivity check
  LOGICAL :: l_pos3           = .TRUE.   ! switch on positivity check
  LOGICAL :: l_pos4           = .TRUE.   ! switch on positivity check
  LOGICAL :: l_pos5           = .TRUE.   ! switch on positivity check
  LOGICAL :: l_pos6           = .TRUE.   ! switch on positivity check

  LOGICAL, TARGET :: l_pcond   = .TRUE.  ! Condensation
  LOGICAL, TARGET :: l_praut   = .TRUE.  ! Autoconversion cloud -> rain
  LOGICAL, TARGET :: l_pracw   = .TRUE.  ! Accretion  cloud -> rain
  LOGICAL, TARGET :: l_pracr   = .TRUE.  ! aggregation of rain drops
  LOGICAL, TARGET :: l_prevp   = .TRUE.  ! evaporation of rain
  LOGICAL, TARGET :: l_psedl   = .FALSE.  ! sedimentation of cloud
  LOGICAL, TARGET :: l_psedr   = .TRUE.  ! sedimentation of rain
  LOGICAL, TARGET :: l_ptidy   = .TRUE.  ! tidying term 1
  LOGICAL, TARGET :: l_ptidy2  = .TRUE.  ! tidying term 2
  LOGICAL, TARGET :: l_pinuc   = .TRUE.  ! ice nucleation
  LOGICAL, TARGET :: l_pidep   = .TRUE.  ! ice deposition
  LOGICAL, TARGET :: l_piacw   = .TRUE.  ! ice accreting water
  LOGICAL, TARGET :: l_psaut   = .TRUE.  ! ice autoconversion ice -> snow
  LOGICAL, TARGET :: l_psdep   = .TRUE.  ! vapour deposition onto snow
  LOGICAL, TARGET :: l_psacw   = .TRUE.  ! snow accreting water
  LOGICAL, TARGET :: l_pgdep   = .TRUE.  ! vapour deposition onto graupel
  LOGICAL, TARGET :: l_pseds   = .TRUE.  ! snow sedimentation
  LOGICAL, TARGET :: l_psedi   = .FALSE.  ! ice sedimentation
  LOGICAL, TARGET :: l_psedg   = .TRUE.  ! graupel sedimentation
  LOGICAL, TARGET :: l_psaci   = .TRUE.  ! snow accreting ice
  LOGICAL, TARGET :: l_praci   = .TRUE.  ! rain accreting ice
  LOGICAL, TARGET :: l_psacr   = .TRUE.  ! snow accreting rain
  LOGICAL, TARGET :: l_pgacr   = .TRUE.  ! graupel accreting rain
  LOGICAL, TARGET :: l_pgacw   = .TRUE.  ! graupel accreting cloud water
  LOGICAL, TARGET :: l_pgaci   = .TRUE.  ! graupel accreting ice
  LOGICAL, TARGET :: l_pgacs   = .TRUE.  ! graupel accreting snow
  LOGICAL, TARGET :: l_piagg   = .FALSE. ! aggregation of ice particles
  LOGICAL, TARGET :: l_psagg   = .TRUE.  ! aggregation of snow particles
  LOGICAL, TARGET :: l_pgagg   = .FALSE. ! aggregation of graupel particles
  LOGICAL, TARGET :: l_psbrk   = .TRUE.  ! break up of snow flakes
  LOGICAL, TARGET :: l_pgshd   = .TRUE.  ! shedding of liquid from graupel
  LOGICAL, TARGET :: l_pihal   = .TRUE.  ! hallet mossop
  LOGICAL, TARGET :: l_psmlt   = .TRUE.  ! snow melting
  LOGICAL, TARGET :: l_pgmlt   = .TRUE.  ! graupel melting
  LOGICAL, TARGET :: l_phomr   = .TRUE.  ! homogeneous freezing of rain
  LOGICAL, TARGET :: l_phomc   = .TRUE.  ! homogeneous freezing of cloud droplets
  LOGICAL, TARGET :: l_pssub   = .TRUE.  ! sublimation of snow
  LOGICAL, TARGET :: l_pgsub   = .TRUE.  ! sublimation of graupel
  LOGICAL, TARGET :: l_pisub   = .TRUE.  ! sublimation of ice
  LOGICAL, TARGET :: l_pimlt   = .TRUE.  ! ice melting

  LOGICAL :: l_tidy_conserve_E = .true.
  LOGICAL :: l_tidy_conserve_q = .true.

  TYPE :: process_switch
     LOGICAL, POINTER :: l_pcond ! Condensation
     LOGICAL, POINTER :: l_praut ! Autoconversion cloud -> rain
     LOGICAL, POINTER :: l_pracw ! Accretion  cloud -> rain
     LOGICAL, POINTER :: l_pracr ! aggregation of rain drops
     LOGICAL, POINTER :: l_prevp ! evaporation of rain
     LOGICAL, POINTER :: l_psedl ! sedimentation of cloud
     LOGICAL, POINTER :: l_psedr ! sedimentation of rain
     LOGICAL, POINTER :: l_pinuc ! ice nucleation
     LOGICAL, POINTER :: l_pidep ! ice deposition
     LOGICAL, POINTER :: l_piacw ! ice accreting water
     LOGICAL, POINTER :: l_psaut ! ice autoconversion ice -> snow
     LOGICAL, POINTER :: l_psdep ! vapour deposition onto snow
     LOGICAL, POINTER :: l_psacw ! snow accreting water
     LOGICAL, POINTER :: l_pgdep ! vapour deposition onto graupel
     LOGICAL, POINTER :: l_pseds ! snow sedimentation
     LOGICAL, POINTER :: l_psedi ! ice sedimentation
     LOGICAL, POINTER :: l_psedg ! graupel sedimentation
     LOGICAL, POINTER :: l_psaci ! snow accreting ice
     LOGICAL, POINTER :: l_praci ! rain accreting ice
     LOGICAL, POINTER :: l_psacr ! snow accreting rain
     LOGICAL, POINTER :: l_pgacr ! graupel accreting rain
     LOGICAL, POINTER :: l_pgacw ! graupel accreting cloud water
     LOGICAL, POINTER :: l_pgaci ! graupel accreting ice
     LOGICAL, POINTER :: l_pgacs ! graupel accreting snow
     LOGICAL, POINTER :: l_piagg ! aggregation of ice particles
     LOGICAL, POINTER :: l_psagg ! aggregation of snow particles
     LOGICAL, POINTER :: l_pgagg ! aggregation of graupel particles
     LOGICAL, POINTER :: l_psbrk ! break up of snow flakes
     LOGICAL, POINTER :: l_pgshd ! shedding of liquid from graupel
     LOGICAL, POINTER :: l_pihal ! hallet mossop
     LOGICAL, POINTER :: l_psmlt ! snow melting
     LOGICAL, POINTER :: l_pgmlt ! graupel melting
     LOGICAL, POINTER :: l_phomr ! homogeneous freezing of rain
     LOGICAL, POINTER :: l_phomc ! homogeneous freezing of cloud droplets
     LOGICAL, POINTER :: l_pssub ! sublimation of snow
     LOGICAL, POINTER :: l_pgsub ! sublimation of graupel
     LOGICAL, POINTER :: l_pisub ! sublimation of ice
     LOGICAL, POINTER :: l_pimlt ! ice melting
     LOGICAL, POINTER :: l_tidy  ! Tidying
     LOGICAL, POINTER :: l_tidy2 ! Tidying
  END TYPE process_switch

  TYPE(process_switch), SAVE, TARGET :: pswitch


  REAL(wp) :: contact_efficiency = 0.0001   ! Arbitrary efficiency for contact nucleation
  REAL(wp) :: immersion_efficiency = 1.0 ! Arbitrary efficiency for immersion/condensation freezing

  REAL(wp) :: max_mu = 35 ! Maximum value of shape parameter
  REAL(wp) :: max_mu_frac = 0.75 ! Fraction of maximum value of shape parameter at which psd limitation kicks in
  ! (I l_limit_psd=.true.)
  REAL(wp) :: fix_mu = 2.5 ! Fixed value for shape parameter (1M/2M)

  LOGICAL :: l_inhom_revp = .FALSE.  ! Inhomogeneous evaporation of raindrop

  ! Aerosol-cloud switches
  INTEGER :: aerosol_option = 0 ! Determines how many aerosol modes to use
  !  1: 
  !   soluble_modes(:) = (/ 0, 2, 2/)
  !   insoluble_modes(:) = (/ 0, 0/)
  !  2:
  !   soluble_modes(:) = (/ 0, 2, 2/)
  !   insoluble_modes(:) = (/ 0, 2/)


  LOGICAL :: l_separate_rain=.FALSE.   ! Use a separate rain category for active aerosol
  LOGICAL :: l_process=.FALSE.         ! process aerosol and/or dust
  LOGICAL :: l_passivenumbers=.FALSE.  ! Use passive numbers in aerosol recovery
  LOGICAL :: l_passivenumbers_ice=.FALSE.  ! Use passive numbers in ice aerosol recovery
  INTEGER :: process_level = 0
  ! 0: no processing:      l_process=l_passivenumbers=.false.
  ! 1: passive processing: l_process=l_passivenumbers=.true.
  ! 2: full processing:    l_process=.true, l_passivenumbers=.false.
  ! 3: passive processing of ice only: l_process=l_passivenumbers_ice=.true., l_passivenumbers=.false.

  ! aerosol process switches
  LOGICAL :: l_aacc=.TRUE.
  LOGICAL :: l_ased=.TRUE.
  LOGICAL :: l_dsaut=.TRUE.
  LOGICAL :: l_aact = .TRUE.   ! activation
  LOGICAL :: l_aaut = .TRUE.   ! autoconversion
  LOGICAL :: l_aacw = .TRUE.   ! accretion cloud by rain
  LOGICAL :: l_aevp = .TRUE.   ! evaporation of cloudd
  LOGICAL :: l_asedr = .TRUE.  ! sedimentation of rain
  LOGICAL :: l_arevp = .TRUE.  ! evaporation of rain
  LOGICAL :: l_asedl = .TRUE.  ! sedimentation of cloud
  LOGICAL :: l_atidy = .TRUE.  ! tidy process 1
  LOGICAL :: l_atidy2 = .TRUE. ! tidy process 2
  LOGICAL :: l_dnuc = .TRUE.   ! ice nucleation
  LOGICAL :: l_dsub = .TRUE.   ! sublimation of ice
  LOGICAL :: l_dsedi = .TRUE.  ! sedimentation of ice
  LOGICAL :: l_dseds = .TRUE.  ! sedimentation of snow
  LOGICAL :: l_dsedg = .TRUE.  ! sedimentation of graupel
  LOGICAL :: l_dssub = .TRUE.  ! sublimation of snow
  LOGICAL :: l_dgsub = .TRUE.  ! sublimation of graupel
  LOGICAL :: l_dhomc = .TRUE.  ! homogeneous freezing of cloud
  LOGICAL :: l_dhomr = .TRUE.  ! homogeneous freezing of rain
  LOGICAL :: l_dimlt = .TRUE.  ! melting of ice
  LOGICAL :: l_dsmlt = .TRUE.  ! melting of snow
  LOGICAL :: l_dgmlt = .TRUE.  ! melting of graupel
  LOGICAL :: l_diacw = .TRUE.  ! riming: ice collects cloud
  LOGICAL :: l_dsacw = .TRUE.  ! riming: snow collects cloud
  LOGICAL :: l_dgacw = .TRUE.  ! riming: graupel collects cloud
  LOGICAL :: l_dsacr = .TRUE.  ! riming: snow collects rain
  LOGICAL :: l_dgacr = .TRUE.  ! riming: graupel collects rain
  LOGICAL :: l_draci = .TRUE.  ! riming: rain collects ice

  LOGICAL :: l_raci_g = .TRUE. ! Allow rain collecting ice to go to graupel if rain mass is significant
  ! Goes to snow otherwise

  LOGICAL :: l_onlycollect = .FALSE. ! Only do collection processes and sedimentation

  INTEGER :: i_aerosed_method = 1 ! method to use for aerosol sedimentation
  ! 1= Use hydrometeor mass fallspeeds

  ! cloud activation
  INTEGER :: iopt_act=0
  ! 0=fixed_cloud
  ! 1=90% of total aerosol number
  ! 2=Twomey
  ! 3=Abdul-Razzak and Ghan
  LOGICAL :: l_active_inarg2000 = .FALSE.  ! consider activated aerosol in activation calculation

  INTEGER :: iopt_rcrit = 0 ! method used to calculate rcrit (only use 0)

  ! cloud - rain autoconversion
  INTEGER :: iopt_auto = 1

  ! rain - cloud accretion
  INTEGER :: iopt_accr = 1

  ! ice - nucleation
  INTEGER :: iopt_inuc = 1

  LOGICAL :: l_itotsg = .FALSE. ! Consider existing snow and graupel numbers when nucleating ice

  LOGICAL :: l_passive = .FALSE. ! only do sedimentation and don't apply conversions

  LOGICAL :: l_passive3m = .FALSE. ! Don't use 3rd moment to determine distribution parameters

  LOGICAL :: l_limit_psd = .TRUE. ! limit distribution parameters to prevent large mean particle sizes

  ! Some checks are made to ensure consistency, this will
  ! override them (but only to be used by Ben and Adrian!)
  LOGICAL :: l_override_checks = .FALSE.

  LOGICAL :: l_SM_fix_n0 = .TRUE. ! If running in single moment, then fix n0.  
  ! If false, then we fix na and nb (c.f. the LEM microphysics)

CONTAINS

  SUBROUTINE set_mphys_switches(in_option, in_aerosol_option)

    INTEGER, INTENT(IN) :: in_option, in_aerosol_option

    INTEGER :: iq,iproc,idgproc ! counters

    LOGICAL :: l_on

    REAL :: k1,k2,k3

    IF (.NOT. mphys_is_set) THEN

      CALL derive_logicals
      option=in_option
      aerosol_option=in_aerosol_option
      ! Set level of aerosol processing
      SELECT CASE(process_level)
      CASE default ! 0
        l_process=.FALSE.
        l_passivenumbers=.FALSE.
        l_passivenumbers_ice=.FALSE.
      CASE(1)
        l_process=.TRUE.
        l_passivenumbers=.TRUE.
        l_passivenumbers_ice=.TRUE.
      CASE(2)
        l_process=.TRUE.
        l_passivenumbers=.FALSE.
        l_passivenumbers_ice=.FALSE.
      CASE(3)
        l_process=.TRUE.
        l_passivenumbers=.FALSE.
        l_passivenumbers_ice=.TRUE.
      CASE(4)
        l_process=.TRUE.
        l_passivenumbers=.TRUE.
        l_passivenumbers_ice=.FALSE.
      END SELECT

      if (l_warm)l_passivenumbers_ice=.FALSE. ! Override if there's no ice

      SELECT CASE (option)
      CASE(1)
        option=23000
      CASE(2)
        option=12000
      CASE(3)
        option=11000
      CASE(4)
        option=13000
      CASE(5)
        option=22000
      CASE(6)
        option=22220
      CASE(7)
        option=22222
      END SELECT
      nq_l = option/1e4
      nq_r = (option-nq_l*1e4)/1e3
      IF (.NOT. l_warm) THEN
        nq_i = (option-nq_l*1e4-nq_r*1e3)/1e2
        nq_s = (option-nq_l*1e4-nq_r*1e3-nq_i*1e2)/1e1
        nq_g = (option-nq_l*1e4-nq_r*1e3-nq_i*1e2-nq_s*1e1)/1e0
      END IF

      ntotalq= 2 + nq_l + nq_r + nq_i + nq_s + nq_g

      ! Logicals...
      IF (nq_l>1)l_2mc=.TRUE.
      IF (nq_r>1)l_2mr=.TRUE.
      IF (nq_r>2)l_3mr=.TRUE.
      IF (nq_i>1)l_2mi=.TRUE.
      IF (nq_s>1)l_2ms=.TRUE.
      IF (nq_s>2)l_3ms=.TRUE.
      IF (nq_g>1)l_2mg=.TRUE.
      IF (nq_g>2)l_3mg=.TRUE.

      ! Names
      ALLOCATE(hydro_names(ntotalq))

      !-----------------
      ! Allocate indices
      !-----------------
      iq=0
      CALL allocq(i_qv, iq, hydro_names, 'qv')
      CALL allocq(i_th, iq, hydro_names, 'th')

      ! first moments
      CALL allocq(i_ql, iq, hydro_names, 'ql')
      CALL allocq(i_qr, iq, hydro_names, 'qr')
      IF (.NOT. l_warm) THEN
        IF (nq_i>0)CALL allocq(i_qi,iq, hydro_names, 'qi')
        IF (nq_s>0)CALL allocq(i_qs,iq, hydro_names, 'qs')
        IF (nq_g>0)CALL allocq(i_qg,iq, hydro_names, 'qg')
      END IF

      ! second moments
      i_nstart=iq+1
      IF (l_2mc)CALL allocq(i_nl,iq, hydro_names, 'nl')
      IF (l_2mr)CALL allocq(i_nr,iq, hydro_names, 'nr')
      IF (.NOT. l_warm) THEN
        IF (l_2mi)CALL allocq(i_ni,iq, hydro_names, 'ni')
        IF (l_2ms)CALL allocq(i_ns,iq, hydro_names, 'ns')
        IF (l_2mg)CALL allocq(i_ng,iq, hydro_names, 'ng')
      END IF

      ! third moments
      i_m3start=iq+1
      IF (l_3mr)CALL allocq(i_m3r,iq, hydro_names, 'm3r')
      IF (.NOT. l_warm) THEN
        IF (l_3ms)CALL allocq(i_m3s,iq, hydro_names, 'm3s')
        IF (l_3mg)CALL allocq(i_m3g,iq, hydro_names, 'm3g')
      END IF

      ! Set params
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

      rain_params%fix_mu=fix_mu
      snow_params%fix_mu=fix_mu
      graupel_params%fix_mu=fix_mu

      ! Set complexity
      IF (l_warm) THEN
        hydro_complexity%nspecies = 2
      ELSE
        hydro_complexity%nspecies = 5
      END IF

      ALLOCATE(hydro_complexity%nmoments(hydro_complexity%nspecies))
      hydro_complexity%nmoments(1) = nq_l
      hydro_complexity%nmoments(2) = nq_r
      IF (.NOT. l_warm) THEN
        hydro_complexity%nmoments(3) = nq_i
        hydro_complexity%nmoments(4) = nq_s
        hydro_complexity%nmoments(5) = nq_g
      END IF

      !-----------------------------------------
      ! Set logicals for 1/2/3 moments in params
      !-----------------------------------------
      IF (i_ql > 0)  cloud_params%l_1m=.TRUE.
      IF (i_qr > 0)  rain_params%l_1m=.TRUE.
      IF (i_qi > 0)  ice_params%l_1m=.TRUE.
      IF (i_qs > 0)  snow_params%l_1m=.TRUE.
      IF (i_qg > 0)  graupel_params%l_1m=.TRUE.

      IF (i_nl > 0)  cloud_params%l_2m=.TRUE.
      IF (i_nr > 0)  rain_params%l_2m=.TRUE.
      IF (i_m3r > 0) rain_params%l_3m=.TRUE.

      IF (i_ni > 0)  ice_params%l_2m=.TRUE.
      IF (i_ns > 0)  snow_params%l_2m=.TRUE.
      IF (i_m3s > 0) snow_params%l_3m=.TRUE.
      IF (i_ng > 0)  graupel_params%l_2m=.TRUE.
      IF (i_m3g > 0) graupel_params%l_3m=.TRUE.

      !-------------------------------------
      ! Set thresholds for various qfields
      !-------------------------------------

      ALLOCATE(thresh_small(0:ntotalq)) ! Note that 0 index will be written to
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
      thresh_small(i_m3r)=EXP(-k1*LOG(qr_small) -k2*LOG(nr_small))
      thresh_small(i_qi)=qi_small
      thresh_small(i_ni)=ni_small
      thresh_small(i_qs)=qs_small
      thresh_small(i_ns)=ns_small
      thresh_small(i_m3s)=EXP(-k1*LOG(qs_small) -k2*LOG(ns_small))
      thresh_small(i_qg)=qg_small
      thresh_small(i_ng)=ng_small
      thresh_small(i_m3g)=EXP(-k1*LOG(qg_small) -k2*LOG(ng_small))

      ALLOCATE(thresh_tidy(0:ntotalq)) ! Note that 0 index will be written to
      ! if a variable is not used
      thresh_tidy(i_qv)=qv_tidy
      thresh_tidy(i_th)=th_tidy
      thresh_tidy(i_ql)=ql_tidy
      thresh_tidy(i_qr)=qr_tidy
      thresh_tidy(i_nl)=nl_tidy
      thresh_tidy(i_nr)=nr_tidy
      thresh_tidy(i_m3r)=EXP(-2.0*k1*LOG(qr_tidy) -2.0*k2*LOG(nr_tidy))
      thresh_tidy(i_qi)=qi_tidy
      thresh_tidy(i_qs)=qs_tidy
      thresh_tidy(i_ni)=ni_tidy
      thresh_tidy(i_ns)=ns_tidy
      thresh_tidy(i_m3s)=EXP(-2.0*k1*LOG(qs_tidy) -2.0*k2*LOG(ns_tidy))
      thresh_tidy(i_qg)=qg_tidy
      thresh_tidy(i_ng)=ng_tidy
      thresh_tidy(i_m3g)=EXP(-2.0*k1*LOG(qs_tidy) -2.0*k2*LOG(ns_tidy))

      ALLOCATE(thresh_sig(0:ntotalq)) ! Note that 0 index will be written to
      thresh_sig(i_ql)=ql_sig
      thresh_sig(i_qr)=qr_sig
      thresh_sig(i_qs)=qs_sig
      thresh_sig(i_qi)=qi_sig
      thresh_sig(i_qg)=qg_sig

      ALLOCATE(thresh_large(0:ntotalq)) ! Note that 0 index will be written to
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
      thresh_large(i_m3r)=EXP(-k1*LOG(qr_large) -k2*LOG(nr_large))
      thresh_large(i_m3s)=EXP(-k1*LOG(qs_large) -k2*LOG(ns_large))
      thresh_large(i_m3g)=EXP(-k1*LOG(qg_large) -k2*LOG(ng_large))

      SELECT CASE(aerosol_option)
      CASE default
        ! Use defaults or assign directly through namelists
      CASE (1)
        soluble_modes(:) = (/ 0, 2, 2/)
        insoluble_modes(:) = (/ 0, 0/)
      CASE (2)
        soluble_modes(:) = (/ 0, 2, 2/)
        insoluble_modes(:) = (/ 0, 2/)
      END SELECT

      ! NB active_X switches have been made obsolete by
      ! l_process, l_passivenumbers and l_separate_rain switches
      ! so should gradually be removed
      IF (l_process) THEN
         if (l_warm)then
            active_cloud(:) = (/.TRUE., .FALSE./)
            active_ice(:) = (/.FALSE., .FALSE./)
         else
            active_cloud(:) = (/.TRUE., .TRUE./)
            active_ice(:) = (/.TRUE., .TRUE./)
         end if
        IF (l_separate_rain)active_rain(:) = (/.TRUE./)
        IF (l_passivenumbers)active_number(1) = .TRUE.
        IF (l_passivenumbers_ice)active_number(2) = .TRUE.
      END IF

      nactivea = COUNT(active_cloud(isol:isol))      &
         +       COUNT(active_ice(isol:isol))    &
         +       COUNT(active_rain)         &
         +       COUNT(active_number(isol:isol))

      nactived = COUNT(active_cloud(iinsol:iinsol))     &
         +       COUNT(active_ice(iinsol:iinsol))   &
         +       COUNT(active_number(iinsol:iinsol))

      ntotala = SUM(soluble_modes)          &
         +      SUM(insoluble_modes)    &
         +      COUNT(active_cloud)     &
         +      COUNT(active_ice)       &
         +      COUNT(active_rain)      &
         +      COUNT(active_number)

      ALLOCATE(aero_names(ntotala))


      !-----------------
      ! Allocate indices
      !-----------------
      iq=0
      IF (soluble_modes(1) > 1)   CALL alloca(i_am1,iq,aero_names, 'AitkenSolMass', aero_index, i_ccn=imass)
      IF (soluble_modes(1) > 0)   CALL alloca(i_an1,iq,aero_names, 'AitkenSolNumber', aero_index, i_ccn=inumber)
      IF (soluble_modes(2) > 1)   CALL alloca(i_am2,iq,aero_names, 'AccumSolMass', aero_index, i_ccn=imass)
      IF (soluble_modes(2) > 0)   CALL alloca(i_an2,iq,aero_names, 'AccumSolNumber', aero_index, i_ccn=inumber)
      IF (soluble_modes(3) > 1)   CALL alloca(i_am3,iq,aero_names, 'CoarseSolMass', aero_index, i_ccn=imass)
      IF (soluble_modes(3) > 0)   CALL alloca(i_an3,iq,aero_names, 'CoarseSolNumber', aero_index, i_ccn=inumber)
      IF (active_cloud(isol))     CALL alloca(i_am4,iq,aero_names, 'ActiveSolCloud', aero_index)
      IF (active_rain(isol))      CALL alloca(i_am5,iq,aero_names, 'ActiveSolRain', aero_index)
      IF (insoluble_modes(2) > 1) CALL alloca(i_am6,iq,aero_names, 'CoarseInsolMass', aero_index, i_in=imass)
      IF (insoluble_modes(2) > 0) CALL alloca(i_an6,iq,aero_names, 'CoarseInsolNumber', aero_index, i_in=inumber)
      IF (active_ice(iinsol))     CALL alloca(i_am7,iq,aero_names, 'ActiveInsolIce', aero_index)
      IF (active_ice(isol))       CALL alloca(i_am8,iq,aero_names, 'ActiveSolIce', aero_index)
      IF (active_cloud(iinsol))   CALL alloca(i_am9,iq,aero_names, 'ActiveInsolCloud', aero_index)
      IF (insoluble_modes(1) > 1) CALL alloca(i_am10,iq,aero_names, 'AccumInsolMass', aero_index, i_in=imass)
      IF (insoluble_modes(1) > 0) CALL alloca(i_an10,iq,aero_names, 'AccumInsolNumber', aero_index, i_in=inumber)
      IF (active_number(isol))    CALL alloca(i_an11,iq,aero_names, 'ActiveSolNumber', aero_index)
      IF (active_number(iinsol))  CALL alloca(i_an12,iq,aero_names, 'ActiveInsolNumber', aero_index)

      aero_complexity%nspecies = COUNT(soluble_modes > 0)       &
         +                       COUNT(insoluble_modes > 0) &
         +                       COUNT(active_cloud)        &
         +                       COUNT(active_ice)          &
         +                       COUNT(active_rain)         &
         +                       COUNT(active_number)

      IF (soluble_modes(1) > 0)   aero_index%i_aitken=1
      IF (soluble_modes(2) > 0)   aero_index%i_accum=aero_index%i_aitken + 1
      IF (soluble_modes(3) > 0)   aero_index%i_coarse=aero_index%i_accum + 1
      IF (insoluble_modes(1) > 0) aero_index%i_accum_dust=1
      IF (insoluble_modes(2) > 0) aero_index%i_coarse_dust=aero_index%i_accum_dust + 1


      IF (l_process) THEN
        ALLOCATE(thresh_atidy(0:ntotala)) ! Note that 0 index will be written to
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
        thresh_atidy(i_am10)=aeronumber_small
        thresh_atidy(i_an10)=aeronumber_small
        thresh_atidy(i_an11)=aeronumber_small
        thresh_atidy(i_an12)=aeronumber_small
      END IF


      !--------------------------
      ! Set the process choices - THESE NEED TIDYING WITH LOGICALS
      !--------------------------
      iopt_accr = 1
      iopt_auto = 1

      !--------------------------
      ! Allocate process indices
      !--------------------------
      iproc=0
      idgproc=0
      IF (l_pcond) CALL allocp(i_cond, iproc, idgproc, 'pcond')
      IF (l_praut) CALL allocp(i_praut, iproc, idgproc, 'praut')
      IF (l_pracw) CALL allocp(i_pracw, iproc, idgproc, 'pracw')
      IF (l_pracr) CALL allocp(i_pracr, iproc, idgproc, 'pracr')
      IF (l_prevp) CALL allocp(i_prevp, iproc, idgproc, 'prevp')
      IF (l_psedl) CALL allocp(i_psedl, iproc, idgproc, 'psedl')
      IF (l_psedr) CALL allocp(i_psedr, iproc, idgproc, 'psedr')
      IF (l_ptidy) CALL allocp(i_tidy, iproc, idgproc, 'ptidy')
      IF (l_ptidy2) CALL allocp(i_tidy2, iproc, idgproc, 'ptidy2')
      IF (.NOT. l_warm) THEN
        IF (l_pinuc) CALL allocp(i_inuc, iproc, idgproc, 'pinuc')
        IF (l_pidep) CALL allocp(i_idep, iproc, idgproc, 'pidep')
        IF (l_piacw) CALL allocp(i_iacw, iproc, idgproc, 'piacw')
        IF (l_psaut) CALL allocp(i_saut, iproc, idgproc, 'psaut')
        IF (l_psdep) CALL allocp(i_sdep, iproc, idgproc, 'psdep')
        IF (l_psacw) CALL allocp(i_sacw, iproc, idgproc, 'psacw')
        IF (l_pgdep) CALL allocp(i_gdep, iproc, idgproc, 'pgdep')
        IF (l_pseds) CALL allocp(i_pseds, iproc, idgproc, 'pseds')
        IF (l_psedi) CALL allocp(i_psedi, iproc, idgproc, 'psedi')
        IF (l_psedg) CALL allocp(i_psedg, iproc, idgproc, 'psedg')
        IF (l_psaci) CALL allocp(i_saci, iproc, idgproc, 'psaci')
        IF (l_praci) CALL allocp(i_raci, iproc, idgproc, 'praci')
        IF (l_psacr) CALL allocp(i_sacr, iproc, idgproc, 'psacr')
        IF (l_pgacr) CALL allocp(i_gacr, iproc, idgproc, 'pgacr')
        IF (l_pgacw) CALL allocp(i_gacw, iproc, idgproc, 'pgacw')
        IF (l_pgaci) CALL allocp(i_gaci, iproc, idgproc, 'pgaci')
        IF (l_pgacs) CALL allocp(i_gacs, iproc, idgproc, 'pgacs')
        IF (l_piagg) CALL allocp(i_iagg, iproc, idgproc, 'piagg')
        IF (l_psagg) CALL allocp(i_sagg, iproc, idgproc, 'psagg')
        IF (l_pgagg) CALL allocp(i_gagg, iproc, idgproc, 'pgagg')
        IF (l_psbrk) CALL allocp(i_sbrk, iproc, idgproc, 'psbrk')
        IF (l_pgshd) CALL allocp(i_gshd, iproc, idgproc, 'pgshd')
        IF (l_pihal) CALL allocp(i_ihal, iproc, idgproc, 'pihal')
        IF (l_psmlt) CALL allocp(i_smlt, iproc, idgproc, 'psmlt')
        IF (l_pgmlt) CALL allocp(i_gmlt, iproc, idgproc, 'pgmlt')
        IF (l_phomr) CALL allocp(i_homr, iproc, idgproc, 'phomr')
        IF (l_phomc) CALL allocp(i_homc, iproc, idgproc, 'phomc')
        IF (l_pssub) CALL allocp(i_ssub, iproc, idgproc, 'pssub')
        IF (l_pgsub) CALL allocp(i_gsub, iproc, idgproc, 'pgsub')
        IF (l_pisub) CALL allocp(i_isub, iproc, idgproc, 'pisub')
        IF (l_pimlt) CALL allocp(i_imlt, iproc, idgproc, 'pimlt')
      END IF
      hydro_complexity%nprocesses = iproc

      ! allocate process indices for aerosols...
      iproc=0
      idgproc=idgproc
      IF (l_process) THEN
        IF (l_aact) CALL allocp(i_aact, iproc, idgproc, 'aact')
        IF (l_aaut) CALL allocp(i_aaut, iproc, idgproc, 'aaut')
        IF (l_aacw) CALL allocp(i_aacw, iproc, idgproc, 'aacw')
        IF (l_aevp) CALL allocp(i_aevp, iproc, idgproc, 'aevp')
        IF (l_asedr) CALL allocp(i_asedr, iproc, idgproc, 'asedr')
        IF (l_arevp) CALL allocp(i_arevp, iproc, idgproc, 'arevp')
        IF (l_asedl) CALL allocp(i_asedl, iproc, idgproc, 'asedl')
        !... additional tidying processes (Need to sort these)
        IF (l_atidy) CALL allocp(i_atidy, iproc, idgproc, 'atidy')
        IF (l_atidy2) CALL allocp(i_atidy2, iproc, idgproc, 'atidy2')
        !... ice related processes
        l_on= .NOT. l_warm
        IF (l_dnuc) CALL allocp(i_dnuc, iproc, idgproc, 'dnuc', l_onoff=l_on)
        IF (l_dsub) CALL allocp(i_dsub, iproc, idgproc, 'dsub', l_onoff=l_on)
        IF (l_dsedi) CALL allocp(i_dsedi, iproc, idgproc, 'dsedi', l_onoff=l_on)
        IF (l_dseds) CALL allocp(i_dseds, iproc, idgproc, 'dseds', l_onoff=l_on)
        IF (l_dsedg) CALL allocp(i_dsedg, iproc, idgproc, 'dsedg', l_onoff=l_on)
        IF (l_dssub) CALL allocp(i_dssub, iproc, idgproc, 'dssub', l_onoff=l_on)
        IF (l_dgsub) CALL allocp(i_dgsub, iproc, idgproc, 'dgsub', l_onoff=l_on)
        IF (l_dhomc) CALL allocp(i_dhomc, iproc, idgproc, 'dhomc', l_onoff=l_on)
        IF (l_dhomr) CALL allocp(i_dhomr, iproc, idgproc, 'dhomr', l_onoff=l_on)
        IF (l_dimlt) CALL allocp(i_dimlt, iproc, idgproc, 'dimlt', l_onoff=l_on)
        IF (l_dsmlt) CALL allocp(i_dsmlt, iproc, idgproc, 'dsmlt', l_onoff=l_on)
        IF (l_dgmlt) CALL allocp(i_dgmlt, iproc, idgproc, 'dgmlt', l_onoff=l_on)
        IF (l_diacw) CALL allocp(i_diacw, iproc, idgproc, 'diacw', l_onoff=l_on)
        IF (l_dsacw) CALL allocp(i_dsacw, iproc, idgproc, 'dsacw', l_onoff=l_on)
        IF (l_dgacw) CALL allocp(i_dgacw, iproc, idgproc, 'dgacw', l_onoff=l_on)
        IF (l_dsacr) CALL allocp(i_dsacr, iproc, idgproc, 'dsacr', l_onoff=l_on)
        IF (l_dgacr) CALL allocp(i_dgacr, iproc, idgproc, 'dgacr', l_onoff=l_on)
        IF (l_draci) CALL allocp(i_draci, iproc, idgproc, 'draci', l_onoff=l_on)

      END IF
      aero_complexity%nprocesses = iproc

      !----------------------------------------------------
      ! Ensure we only do this at the start of the run
      !----------------------------------------------------
      mphys_is_set=.TRUE.

    END IF

  END SUBROUTINE set_mphys_switches

  SUBROUTINE allocq(i, iq, names, name)

    ! Allocate an index to a q variable

    INTEGER, INTENT(OUT) :: i
    INTEGER, INTENT(INOUT) :: iq

    CHARACTER(10), INTENT(INOUT) :: names(:)
    CHARACTER(*), INTENT(IN) :: name

    i=iq+1
    iq=i
    names(iq) = ADJUSTR(TRIM(name))

  END SUBROUTINE allocq

  SUBROUTINE allocp(proc, iproc, idgproc, name, l_onoff)

    TYPE(process_name) :: proc
    INTEGER :: iproc, idgproc
    CHARACTER(*) :: name
    LOGICAL, OPTIONAL, INTENT(IN) :: l_onoff

    iproc=iproc+1
    idgproc=idgproc+1
    proc%id=iproc
    !    proc%unique_id=iproc
    proc%name=ADJUSTR(TRIM(name))
    IF (PRESENT(l_onoff)) THEN
      proc%on=l_onoff
    ELSE
      proc%on=.TRUE.
    END IF

  END SUBROUTINE allocp

  SUBROUTINE alloca(i, iq, names, name, aero_index   &
     , i_ccn, i_in)

    ! Allocate an index to a aerosol variable
    ! similar to allocq, but also adds information
    ! to aero_index if variable should act as an in
    ! or a ccn.

    INTEGER, INTENT(OUT) :: i
    INTEGER, INTENT(INOUT) :: iq

    CHARACTER(20), INTENT(INOUT) :: names(:)
    CHARACTER(*), INTENT(IN) :: name

    TYPE(aerosol_index), INTENT(INOUT) :: aero_index
    ! if present and >0 should be set to imass or inumber
    INTEGER, OPTIONAL, INTENT(IN) :: i_ccn, i_in

    INTEGER :: is_ccn, is_in, nin, nccn

    is_ccn=0
    is_in=0
    IF (PRESENT(i_ccn))is_ccn=i_ccn
    IF (PRESENT(i_in))is_in=i_in


    i=iq+1
    iq=i
    names(iq) = ADJUSTR(TRIM(name))
    IF (is_ccn == imass) THEN
      ! represents mass of aerosol which can act as ccn
      nccn=COUNT(aero_index%ccn_m > 0) + 1
      aero_index%ccn_m(nccn) = iq
      aero_index%nccn=nccn
    ELSE IF (is_ccn == inumber) THEN
      ! represents number of aerosol which can act as ccn
      nccn=COUNT(aero_index%ccn_n > 0) + 1
      aero_index%ccn_n(nccn) = iq
      aero_index%nccn=nccn
    END IF
    IF (is_in == imass) THEN
      ! represents mass of aerosol which can act as in
      nin=COUNT(aero_index%in_m > 0) + 1
      aero_index%in_m(nin) = iq
      aero_index%nin=nin
    ELSE IF (is_in == inumber) THEN
      ! represents number of aerosol which can act as in
      nin=COUNT(aero_index%in_n > 0) + 1
      aero_index%in_n(nin) = iq
      aero_index%nin=nin
    END IF

  END SUBROUTINE alloca

  SUBROUTINE derive_logicals()

    ! Routine sets logical switches which depend on or are overridden by other switches

    ! Transfer namelist logicals to derived type for ease of use later


    pswitch%l_pcond => l_pcond ! Condensation
    pswitch%l_praut => l_praut ! Autoconversion cloud -> rain
    pswitch%l_pracw => l_pracw ! Accretion  cloud -> rain
    pswitch%l_pracr => l_pracr ! aggregation of rain drops
    pswitch%l_prevp => l_prevp ! evaporation of rain
    pswitch%l_psedl => l_psedl ! sedimentation of cloud
    pswitch%l_psedr => l_psedr ! sedimentation of rain
    pswitch%l_pinuc => l_pinuc ! ice nucleation
    pswitch%l_pidep => l_pidep ! ice deposition
    pswitch%l_piacw => l_piacw ! ice accreting water
    pswitch%l_psaut => l_psaut ! ice autoconversion ice -> snow
    pswitch%l_psdep => l_psdep ! vapour deposition onto snow
    pswitch%l_psacw => l_psacw ! snow accreting water
    pswitch%l_pgdep => l_pgdep ! vapour deposition onto graupel
    pswitch%l_pseds => l_pseds ! snow sedimentation
    pswitch%l_psedi => l_psedi ! ice sedimentation
    pswitch%l_psedg => l_psedg ! graupel sedimentation
    pswitch%l_psaci => l_psaci ! snow accreting ice
    pswitch%l_praci => l_praci ! rain accreting ice
    pswitch%l_psacr => l_psacr ! snow accreting rain
    pswitch%l_pgacr => l_pgacr ! graupel accreting rain
    pswitch%l_pgacw => l_pgacw ! graupel accreting cloud water
    pswitch%l_pgaci => l_pgaci ! graupel accreting ice
    pswitch%l_pgacs => l_pgacs ! graupel accreting snow
    pswitch%l_piagg => l_piagg ! aggregation of ice particles
    pswitch%l_psagg => l_psagg ! aggregation of snow particles
    pswitch%l_pgagg => l_pgagg ! aggregation of graupel particles
    pswitch%l_psbrk => l_psbrk ! break up of snow flakes
    pswitch%l_pgshd => l_pgshd ! shedding of liquid from graupel
    pswitch%l_pihal => l_pihal ! hallet mossop
    pswitch%l_psmlt => l_psmlt ! snow melting
    pswitch%l_pgmlt => l_pgmlt ! graupel melting
    pswitch%l_phomr => l_phomr ! homogeneous freezing of rain
    pswitch%l_phomc => l_phomc ! homogeneous freezing of cloud droplets
    pswitch%l_pssub => l_pssub ! sublimation of snow
    pswitch%l_pgsub => l_pgsub ! sublimation of graupel
    pswitch%l_pisub => l_pisub ! sublimation of ice
    pswitch%l_pimlt => l_pimlt ! ice melting
    pswitch%l_tidy  => l_ptidy ! Tidying
    pswitch%l_tidy2 => l_ptidy2 ! Tidying

    IF (.NOT. l_rain) THEN
      pswitch%l_praut=.FALSE.
      pswitch%l_pracw=.FALSE.
      pswitch%l_pracr=.FALSE.
      pswitch%l_prevp=.FALSE.
    END IF

    IF (l_onlycollect) THEN
      pswitch%l_praut=.FALSE.
      pswitch%l_prevp=.FALSE.
      pswitch%l_psaut=.FALSE.
      pswitch%l_pihal=.FALSE.
      pswitch%l_phomr=.FALSE.
      pswitch%l_pinuc=.FALSE.
      pswitch%l_pidep=.FALSE.
      pswitch%l_psdep=.FALSE.
      pswitch%l_pgdep=.FALSE.
      pswitch%l_psmlt=.FALSE.
      pswitch%l_pgmlt=.FALSE.
      pswitch%l_pimlt=.FALSE.
      pswitch%l_pcond=.FALSE.
    END IF

    IF (.NOT. l_sg) THEN
      pswitch%l_psaut=.FALSE.
      pswitch%l_psacw=.FALSE.
      pswitch%l_psaci=.FALSE.
      pswitch%l_praci=.FALSE.
      pswitch%l_psacr=.FALSE.
      pswitch%l_pgacw=.FALSE.
      pswitch%l_pgacr=.FALSE.
      pswitch%l_psagg=.FALSE.
      pswitch%l_psdep=.FALSE.
      pswitch%l_psmlt=.FALSE.
      pswitch%l_pseds=.FALSE.
      pswitch%l_pihal=.FALSE.
      l_g=.FALSE.
    END IF

    IF (.NOT. l_g) THEN
      pswitch%l_pgacw=.FALSE.
      pswitch%l_pgacr=.FALSE.
      pswitch%l_pgdep=.FALSE.
      pswitch%l_pgmlt=.FALSE.
      pswitch%l_psedg=.FALSE.
    END IF

    IF (.NOT. l_halletmossop) THEN
      pswitch%l_pihal=.FALSE.
    END IF

    IF (.NOT. (pswitch%l_pidep .OR.     &
       pswitch%l_psdep .OR. &
       pswitch%l_pgdep)) l_idep=.FALSE.


    IF (.NOT. (pswitch%l_pisub .OR.     &
       pswitch%l_pssub .OR. &
       pswitch%l_pgsub)) l_isub=.FALSE.

    IF (.NOT. (                   &
       pswitch%l_praut .OR. &
       pswitch%l_pracw .OR. &
       pswitch%l_piacw .OR. &
       pswitch%l_psacw .OR. &
       pswitch%l_pgacw .OR. &
       pswitch%l_phomc      &
       )) l_pos1=.FALSE.

    IF (.NOT. (                   &
       pswitch%l_praci .OR. &
       pswitch%l_psaci .OR. &
       pswitch%l_pgaci .OR. &
       pswitch%l_psaut .OR. &
       pswitch%l_pisub .OR. &
       pswitch%l_pimlt      &
       )) l_pos2=.FALSE.

    IF (.NOT. (                   &
       pswitch%l_prevp .OR. &
       pswitch%l_psacr .OR. &
       pswitch%l_pgacr .OR. &
       pswitch%l_phomr      &
       )) l_pos3=.FALSE.

    IF (.NOT. (                   &
       pswitch%l_pgacs .OR. &
       pswitch%l_psmlt .OR. &
       pswitch%l_psacr .OR. &
       pswitch%l_pssub      &
       )) l_pos4=.FALSE.

    IF (.NOT. (                   &
       pswitch%l_praut .OR. &
       pswitch%l_pracw      &
       )) l_pos5=.FALSE.

    IF (.NOT. (                   &
       pswitch%l_prevp      &
       )) l_pos6=.FALSE.

  END SUBROUTINE derive_logicals

END MODULE mphys_switches
