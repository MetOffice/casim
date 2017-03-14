MODULE generic_diagnostic_variables

IMPLICIT NONE
SAVE

! Method:
! To add a new diagnostic item, include a logical flag for that
! diagnostic and a real allocatable array of the appropriate
! dimensions to the list below.

! Each parent model sets the logical flag (outside of CASIM)
! and then CASIM uses that flag to allocate the appropriate space
! as required. 

TYPE diaglist

  ! 2D variable logical flags
  LOGICAL :: l_precip         = .FALSE.
  LOGICAL :: l_surface_rain   = .FALSE.
  LOGICAL :: l_surface_snow   = .FALSE.
  LOGICAl :: l_surface_graup  = .FALSE.

  ! 3D variable logical flags
  ! ... to be added ...

  ! 2D variable arrays
  REAL, ALLOCATABLE :: precip(:,:)
  REAL, ALLOCATABLE :: SurfaceRainR(:,:)
  REAL, ALLOCATABLE :: SurfaceSnowR(:,:)
  REAL, ALLOCATABLE :: SurfaceGraupR(:,:)

  ! 3D variable arrays
  ! ... to be added ...

END TYPE diaglist

TYPE (diaglist) :: casdiags

CONTAINS

SUBROUTINE allocate_diagnostic_space(is, ie, js, je, ks, ke)

IMPLICIT NONE

INTEGER, INTENT(IN) :: is ! Start of i array
INTEGER, INTENT(IN) :: ie ! End of i array
INTEGER, INTENT(IN) :: js ! Start of j array
INTEGER, INTENT(IN) :: je ! End of j array
INTEGER, INTENT(IN) :: ks ! Start of k array
INTEGER, INTENT(IN) :: ke ! End of k array

IF ( casdiags % l_precip) THEN

  ALLOCATE ( casdiags % precip(is:ie, js:je) )
  casdiags % precip(:,:) = 0.0

END IF ! casdiags % l_precip

IF ( casdiags % l_surface_rain ) THEN

  ALLOCATE ( casdiags % SurfaceRainR(is:ie, js:je) )
  casdiags % SurfaceRainR(:,:) = 0.0

END IF

IF ( casdiags % l_surface_snow ) THEN

  ALLOCATE ( casdiags % SurfaceSnowR(is:ie, js:je) )
  casdiags % SurfaceSnowR(:,:) = 0.0

END IF

IF ( casdiags % l_surface_graup ) THEN

  ALLOCATE ( casdiags % SurfaceGraupR(is:ie, js:je) )
  casdiags % SurfaceGraupR(:,:) = 0.0

END IF


END SUBROUTINE allocate_diagnostic_space

SUBROUTINE deallocate_diagnostic_space()

IMPLICIT NONE

IF ( ALLOCATED ( casdiags %  SurfaceGraupR )) THEN
  DEALLOCATE ( casdiags % SurfaceGraupR )
END IF

IF ( ALLOCATED ( casdiags %  SurfaceSnowR )) THEN
  DEALLOCATE ( casdiags % SurfaceSnowR )
END IF

IF ( ALLOCATED ( casdiags %  SurfaceRainR )) THEN
  DEALLOCATE ( casdiags % SurfaceRainR )
END IF

IF ( ALLOCATED ( casdiags %  precip )) THEN
  DEALLOCATE ( casdiags % precip )
END IF

END SUBROUTINE deallocate_diagnostic_space

END MODULE generic_diagnostic_variables
