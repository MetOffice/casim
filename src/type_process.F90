MODULE type_process
  ! declares derived types used by processes

USE variable_precision, ONLY: wp

IMPLICIT NONE

!! Allocatable arrays are quicker than pointers, but
!! not supported by all compilers
#if ALLOCATABLE_TYPE == 1
TYPE :: process_rate
  REAL(wp), ALLOCATABLE :: source(:)
END TYPE process_rate
#else
TYPE :: process_rate
  REAL(wp), POINTER :: source(:) ! tendency to apply to target q variable (/s)
END TYPE process_rate
#endif

TYPE :: process_name
  INTEGER :: id          ! Id for array indexing
  INTEGER :: unique_id   ! Unique id for diagnostic identification
  CHARACTER(20) :: name  ! Process name
  LOGICAL :: on          ! is the process going to be used, i.e. on=.true.
END TYPE process_name

END MODULE type_process
