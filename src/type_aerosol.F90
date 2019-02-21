module type_aerosol
  ! declares derived types used by aerosols

  use variable_precision, only: wp

  implicit none

  type aerosol_type
     integer  :: nmodes            ! number of modes
     real(wp), pointer :: Nd(:)    ! number concentration
     ! for each mode (m-3)
     real(wp), pointer :: rd(:)    ! radius for each mode (m)
     real(wp), pointer :: sigma(:) ! dispersion of
     ! distribution
     real(wp), pointer :: vantHoff(:)   ! Van't Hoff factor
     real(wp), pointer :: massMole(:)   ! mass per mole of aerosol
     real(wp), pointer :: density(:)    ! density of aerosol
     real(wp), pointer :: epsv(:)       ! volume fraction
     real(wp), pointer :: beta(:)  ! parameter describing distribution
     ! of soluble fraction of aerosol
     ! within particle volume
  end type aerosol_type

  type aerosol_phys
     integer  :: nmodes                 ! number of modes
     real(wp), pointer :: N(:)          ! number concentration
     ! for each mode (m-3)
     real(wp), pointer :: M(:)          ! mass concentration
     ! for each mode (m-3)
     real(wp), pointer :: rd(:)         ! mean radius for each mode (m)
     real(wp), pointer :: sigma(:)      ! dispersion of
     ! distribution
     real(wp), pointer :: rpart(:)      ! smallest radius of the distribution
     ! (usually 0 unless using Shipway
     ! activated mass scheme)
  end type aerosol_phys

  type aerosol_chem
     integer  :: nmodes                ! number of modes
     real(wp), pointer :: vantHoff(:)  ! Van't Hoff factor
     real(wp), pointer :: massMole(:)  ! mass per mole of aerosol
     real(wp), pointer :: density(:)   ! density of aerosol
     real(wp), pointer :: epsv(:)      ! volume fraction
     real(wp), pointer :: beta(:)      ! parameter describing distribution
     real(wp), pointer :: bk(:) ! combined parameter for volume-weighted hygroscopicity
     ! of soluble fraction of aerosol
     ! within particle volume
     ! (see Shipway and Abel 2010 for details)
  end type aerosol_chem

  type aerosol_active
     ! Additional information about activated aerosol
     real(wp) :: rcrit=999.0       ! radius of smallest activated particle across all species
     real(wp) :: mact=0.0          ! total mass of activated particles across all species
     real(wp) :: nact=0.0          ! total number of activated particles across all species
     real(wp) :: mact_mean=0.0     ! mean mass of activated particles across all species
     real(wp) :: rd = 0.0          ! mean radius
     real(wp) :: sigma = 0.0       ! dispersion of distribution
     real(wp) :: rcrit1=999.0      ! radius of smallest activated particle in species #1
     real(wp) :: mact1=0.0         ! total mass of activated particles in species #1
     real(wp) :: nact1=0.0         ! total number of activated particles in species #1
     real(wp) :: nratio1=0.0       ! fraction of particles in species #1 containint activated particles
     real(wp) :: mact1_mean=0.0    ! mean mass of activated particles in species #1
     real(wp) :: rcrit2=999.0      ! radius of smallest activated particle in species #2
     real(wp) :: mact2=0.0         ! total mass of activated particles in species #2
     real(wp) :: nact2=0.0         ! total number of activated particles in species #2
     real(wp) :: nratio2=0.0       ! fraction of particles in species #2 containint activated particles
     real(wp) :: mact2_mean=0.0    ! mean mass of activated particles in species #2
     real(wp) :: rcrit3=999.0      ! radius of smallest activated particle in species #3
     real(wp) :: mact3=0.0         ! total mass of activated particles in species #3
     real(wp) :: nact3=0.0         ! total number of activated particles in species #3
     real(wp) :: nratio3=0.0       ! fraction of particles in species #1 containint activated particles
     real(wp) :: mact3_mean=0.0    ! mean mass of activated particles in species #3
  end type aerosol_active
end module type_aerosol
