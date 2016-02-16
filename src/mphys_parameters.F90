module mphys_parameters
  use variable_precision, only: wp, iwp
  use special, only: pi

  implicit none

  ! some numbers
  integer :: nz ! number of levels
  integer :: nq ! number of q fields
  integer :: naero ! number of aerosol fields
  integer, parameter :: nspecies = 5 ! Number of different hyrdrometeors

  integer :: nprocs     ! number of process rates stored
  integer :: naeroprocs ! number of process rates stored

  real(wp) :: parent_dt ! Parent model timestep
  !$OMP THREADPRIVATE(nprocs, naeroprocs)

  type hydro_params
     integer  :: id                 ! identifier unique to species
     real(wp) :: p1, p2, p3         ! moments used
     real(wp) :: sp1, sp2, sp3      ! sedimentation moments
     real(wp) :: c_x, d_x, e_x      ! mass-diameter relation
     real(wp) :: a_x, b_x, f_x, g_x ! fallspeed
     real(wp) :: a2_x, b2_x, f2_x   ! fallspeed Abel-Shipway
     real(wp) :: density            ! density
     integer  :: i_1m, i_2m, i_3m   ! moment indices (set elsewhere)
     logical  :: l_1m, l_2m, l_3m   ! moment logicals (set elsewhere)
     real(wp) :: Dmin, Dmax         ! size limits (for mean particle diameter)
     real(wp) :: fix_mu, fix_N0     ! fixed parameter values for lower order representations (or initialization)
     real(wp) :: maxv               ! Maximum realistic bulk fall speed (m/s)
  end type hydro_params

  type(hydro_params) :: cloud_params = hydro_params(   &
       1,                                   & ! identifier unique to species
       3.0,             0.0,      -999,     & ! moments used
       3.0,             0.0,      -999,     & ! sedimentation moments
       pi*997.0/6.0,     3.0,      0.0,     & ! mass-diameter relation
       3.0e7,           2.0,     0.0,    .5,& ! fallspeed
       0.0,             0.0,      0.0,      & ! fallspeed Abel-Shipway
       997.0,                               & ! density
       -999,         -999,    -999,         & ! moment indices (set elsewhere)
       .false.,      .false.,    .false.,   & ! moment logicals (set elsewhere)
       .01e-6,       50.0e-6,               & ! size limits (m)
       0.0,          50.0e6,                & ! fixed parameter values for lower order representations
       .2                                   & ! Maximum realistic bulk fall speed (m/s)
       )

  type(hydro_params) :: rain_params = hydro_params(   &
       2,                                   & ! identifier unique to species
       3.0,             0.0,      6.0,      & ! moments used
       3.0,             0.0,      6.0,      & ! sedimentation moments
       pi*997.0/6.0,     3.0,      0.0,     & ! mass-diameter relation
       !     4854.1,         1.,    195.,    .5,  & ! fallspeed
       130.0,         .5,    0.0,    .5,    & ! fallspeed
       -446.009, 0.782127, 4085.35,         & ! fallspeed Abel-Shipway
       997.0,                               & ! density
       -999,         -999,    -999,         & ! moment indices (set elsewhere)
       .false.,      .false.,    .false.,   & ! moment logicals (set elsewhere)
       0.1e-6,       0.002,                 & ! size limits (m)
       2.5,          5.0e6,                 & ! fixed parameter values for lower order representations
       10.0                                 & ! Maximum realistic bulk fall speed (m/s)
       )

  type(hydro_params) :: ice_params = hydro_params(   &
       3,                                   & ! identifier unique to species
       3.0,             0.0,      6.0,      & ! moments used
       3.0,             0.0,      6.0,      & ! sedimentation moments
       pi*200.0/6.0,     3.0,      0.0,     & ! mass-diameter relation
       71.34,         .6635,   0.0,    .5,  & ! fallspeed
       0.0,           0.0,     0.0,         & ! fallspeed Abel-Shipway
       200.0,                               & ! density
       -999,         -999,    -999,         & ! moment indices (set elsewhere)
       .false.,      .false.,    .false.,   & ! moment logicals (set elsewhere)
       0.1e-6,       200.0e-6,              & ! size limits (m)
       0.0,          0.05e6,                & ! fixed parameter values for lower order representations
       1.0                                  & ! Maximum realistic bulk fall speed (m/s)
       )

  type(hydro_params) :: snow_params = hydro_params(   &
       4,                                   & ! identifier unique to species
       3.0,             0.0,      6.0,      & ! moments used
       3.0,             0.0,      6.0,      & ! sedimentation moments
       pi*100.0/6.0,     3.0,      0.0,     & ! mass-diameter relation
       4.84,          .25,      0.0,    .5, & ! fallspeed
       0.0,           0.0,     0.0,         & ! fallspeed Abel-Shipway
       100.0,                               & ! density
       -999,         -999,    -999,         & ! moment indices (set elsewhere)
       .false.,      .false.,    .false.,   & ! moment logicals (set elsewhere)
       0.1e-6,       0.005,                 & ! size limits (m)
       2.5,          0.1e6,                 & ! fixed parameter values for lower order representations
       4.0                                  & ! Maximum realistic bulk fall speed (m/s)
       )

  type(hydro_params) :: graupel_params = hydro_params(   &
       5,                                   & ! identifier unique to species
       3.0,             0.0,      6.0,      & ! moments used
       3.0,             0.0,      6.0,      & ! sedimentation moments
       pi*500.0/6.0,     3.0,      0.0,     & ! mass-diameter relation
       253.0,          .734,     0.0,  .422,& ! fallspeed
       !     114.5,          .5,     0.,  .422,   & ! fallspeed
       0.0,           0.0,     0.0,         & ! fallspeed Abel-Shipway
       500.0,                               & ! density
       -999,         -999,    -999,         & ! moment indices (set elsewhere)
       .false.,      .false.,    .false.,   & ! moment logicals (set elsewhere)
       0.1e-6,       .005,                  & ! size limits (m)
       2.5,          0.02e6,                & ! fixed parameter values for lower order representations
       20.0                                 &  ! Maximum realistic bulk fall speed (m/s)
       )


  ! mass-diameter
  real(wp) :: c_r=pi*997.0/6.0
  real(wp) :: d_r=3.0

  ! fall-speed
  !  real(wp) :: a_r=4854.1
  !  real(wp) :: b_r=1.
  !  real(wp) :: f_r=195.
  !  real(wp) :: a_r=841.99667
  !  real(wp) :: b_r=.8
  !  real(wp) :: f_r=0.
  real(wp) :: a_r=130.0
  real(wp) :: b_r=.5
  real(wp) :: f_r=.0

  ! Abel Shipway (2nd term) coefficients
  real(wp) :: a2_r=-446.009
  real(wp) :: b2_r=0.782127
  real(wp) :: f2_r=4085.35

  ! moments (these need to be moved)
  real(wp) :: p1=3.0  ! number
  real(wp) :: p2=0.0  ! mass
  real(wp) :: p3=6.0  ! radar reflectivity
  real(wp) :: sp1=3.0  !
  real(wp) :: sp2=1.5  !
  real(wp) :: sp3=0.0  !

  ! mu initialisations
  real(wp) :: mu_aut=2.0 ! initial mu for autoconversion cloud to rain
  real(wp) :: mu_saut=2.0 ! initial mu for autoconversion ice to snow

  ! ventilation coefficients
  real(wp) :: vent_1=.78
  real(wp) :: vent_2=.31

  ! cloud droplet activation
  ! Twomey law C1s^K1
  real(wp) :: C1=100.0
  real(wp) :: K1=0.4

  ! ice-phase
  real(wp) :: nucleated_ice_radius
  real(wp) :: nucleated_ice_mass

  real(wp) :: DImax=0.000125  !< maximum ice diameter before autoconversion (m)
  real(wp) :: DI2S=0.00033   !< mean size of ice particles autoconverting to snow (m)
  real(wp) :: tau_saut=60.0   !< Timescale for autoconversion of snow (s)

  real(wp) :: DSbrk=0.004     !< Threshold diameter for snow breakup (m)
  real(wp) :: tau_sbrk=60.0   !< Timescale for breakup of snow (s)

  real(wp) :: dN_hallet_mossop=350.0e6 !< Number of ice splinters formed per kg of rimed liquid (/kg)
  real(wp) :: M0_hallet_mossop=1.0E-18 !< Mass of newly formed ice splinter (kg)

  real(wp) :: DR_melt=0.001   !< Mean diameter of rain from melt (m)

  real(wp) :: T_hom_freeze=-38 !< Temperature threshold for homogeneous freezing (C)

  ! aerosol options
  real(wp) :: sigma_arc=1.5    !< width of distribution of aerosol activated in cloud+rain
  real(wp) :: sigma_disg=1.5    !< width of distribution of dust activated in ice, snow + graupel


  ! concentration factors for diagnostic aerosol partitioning
  real(wp) :: alpha_s_si=0.1   ! soluble concentration in snow/ soluble concentration in ice
  real(wp) :: alpha_s_gi=0.01  ! soluble concentration in graupel/ soluble concentration in ice
  real(wp) :: alpha_i_si=0.1   ! insoluble concentration in snow/ insoluble concentration in ice
  real(wp) :: alpha_i_gi=0.01  ! insoluble concentration in graupel/ insoluble concentration in ice
  real(wp) :: alpha_s_ri=1.0   ! soluble concentration in rain/ soluble concentration in droplets
  real(wp) :: alpha_i_ri=1.0   ! insoluble concentration in rain/ insoluble concentration in droplets

  ! Some fixed values used for checking for exact equality.
  real(wp), parameter :: ZERO_REAL_WP=0.0 
end module mphys_parameters
