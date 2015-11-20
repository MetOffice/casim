MODULE mphys_die
  !
  ! Routine to exit the microphysics following an error there
  !
#if DEF_MODEL==MODEL_UM
USE ereport_mod, ONLY : ereport
#elif DEF_MODEL==MODEL_KiD
USE runtime, ONLY: time
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE throw_mphys_error(itype, routine, info)

INTEGER, OPTIONAL :: itype  ! type of error
    !  1 = incorrect specification of options
    !  0 = unknown
CHARACTER(*) :: routine
CHARACTER(*), OPTIONAL :: info ! error information

INTEGER, PARAMETER :: nstandard_types=3
CHARACTER(100) :: stdinfo(nstandard_types) =     &
       (/ 'incorrect specification of options      ' &
       ,  'Bad values found                        ' &
       ,  'unknown error                           ' &
       /)

CHARACTER(1000) :: str

REAL :: minus_one=-1.

str='Error in microphysics: '

IF (itype <= nstandard_types) str=TRIM(str)//TRIM(stdinfo(itype))

IF (PRESENT(info)) str=TRIM(str)//' Additional information: '//TRIM(info)
#if DEF_MODEL==MODEL_KiD
PRINT*, 'Runtime is:' , time
#endif

#if DEF_MODEL==MODEL_UM
CALL Ereport(routine, itype, str)
#else
PRINT*, routine,':', TRIM(str)
PRINT*, (minus_one)**0.5
STOP
#endif

END SUBROUTINE throw_mphys_error

END MODULE mphys_die
