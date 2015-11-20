MODULE type_aerosol
  ! declares derived types used by aerosols

USE variable_precision, ONLY: wp

IMPLICIT NONE

TYPE aerosol_type
  INTEGER  :: nmodes            ! number of modes
  REAL(wp), POINTER :: Nd(:)    ! number concentration
                                   ! for each mode (m-3)
  REAL(wp), POINTER :: rd(:)    ! radius for each mode (m)
  REAL(wp), POINTER :: sigma(:) ! dispersion of
                                   ! distribution
  REAL(wp), POINTER :: vantHoff(:)   ! Van't Hoff factor
  REAL(wp), POINTER :: massMole(:)   ! mass per mole of aerosol
  REAL(wp), POINTER :: density(:)    ! density of aerosol
  REAL(wp), POINTER :: epsv(:)       ! volume fraction
  REAL(wp), POINTER :: beta(:)  ! parameter describing distribution
                                   ! of soluble fraction of aerosol
                                   ! within particle volume
END TYPE aerosol_type

TYPE aerosol_phys
  INTEGER  :: nmodes                 ! number of modes
  REAL(wp), POINTER :: N(:)          ! number concentration
                                        ! for each mode (m-3)
  REAL(wp), POINTER :: M(:)          ! mass concentration
                                        ! for each mode (m-3)
  REAL(wp), POINTER :: rd(:)         ! mean radius for each mode (m)
  REAL(wp), POINTER :: sigma(:)      ! dispersion of
                                        ! distribution
  REAL(wp), POINTER :: rpart(:)      ! smallest radius of the distribution
                                        ! (usually 0 unless using Shipway
                                        ! activated mass scheme)
END TYPE aerosol_phys

TYPE aerosol_chem
  INTEGER  :: nmodes                ! number of modes
  REAL(wp), POINTER :: vantHoff(:)  ! Van't Hoff factor
  REAL(wp), POINTER :: massMole(:)  ! mass per mole of aerosol
  REAL(wp), POINTER :: density(:)   ! density of aerosol
  REAL(wp), POINTER :: epsv(:)      ! volume fraction
  REAL(wp), POINTER :: beta(:)      ! parameter describing distribution
                                       ! of soluble fraction of aerosol
                                       ! within particle volume
                                       ! (see Shipway and Abel 2010 for details)
END TYPE aerosol_chem

TYPE aerosol_active
     ! Additional information about activated aerosol
  REAL(wp) :: rcrit=999.0       ! radius of smallest activated particle across all species
  REAL(wp) :: mact=0.0          ! total mass of activated particles across all species
  REAL(wp) :: nact=0.0          ! total number of activated particles across all species
  REAL(wp) :: mact_mean=0.0     ! mean mass of activated particles across all species
  REAL(wp) :: rd = 0.0          ! mean radius
  REAL(wp) :: sigma = 0.0       ! dispersion of distribution
  REAL(wp) :: rcrit1=999.0      ! radius of smallest activated particle in species #1
  REAL(wp) :: mact1=0.0         ! total mass of activated particles in species #1
  REAL(wp) :: nact1=0.0         ! total number of activated particles in species #1
  REAL(wp) :: nratio1=0.0       ! fraction of particles in species #1 containint activated particles
  REAL(wp) :: mact1_mean=0.0    ! mean mass of activated particles in species #1
  REAL(wp) :: rcrit2=999.0      ! radius of smallest activated particle in species #2
  REAL(wp) :: mact2=0.0         ! total mass of activated particles in species #2
  REAL(wp) :: nact2=0.0         ! total number of activated particles in species #2
  REAL(wp) :: nratio2=0.0       ! fraction of particles in species #2 containint activated particles
  REAL(wp) :: mact2_mean=0.0    ! mean mass of activated particles in species #2
  REAL(wp) :: rcrit3=999.0      ! radius of smallest activated particle in species #3
  REAL(wp) :: mact3=0.0         ! total mass of activated particles in species #3
  REAL(wp) :: nact3=0.0         ! total number of activated particles in species #3
  REAL(wp) :: nratio3=0.0       ! fraction of particles in species #1 containint activated particles
  REAL(wp) :: mact3_mean=0.0    ! mean mass of activated particles in species #3
END TYPE aerosol_active

END MODULE type_aerosol
