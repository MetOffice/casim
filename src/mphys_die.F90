! Routine to exit the microphysics following an error there
module mphys_die
#if DEF_MODEL==MODEL_UM
  use ereport_mod, only : ereport
#elif DEF_MODEL==MODEL_KiD
  use runtime, only: time
#endif

  implicit none
  private

  public throw_mphys_error
contains

  subroutine throw_mphys_error(itype, routine, info)
    integer, optional, intent(in) :: itype  ! type of error 1 = incorrect specification of options, 0 = unknown
    character(*), intent(in) :: routine
    character(*), optional, intent(in) :: info ! error information

    integer, parameter :: nstandard_types=3
    character(100) :: stdinfo(nstandard_types) =     &
         (/ 'incorrect specification of options      ' &
         ,  'Bad values found                        ' &
         ,  'unknown error                           ' &
         /)
    character(1000) :: str
    real :: minus_one=-1.

    str='Error in microphysics: '
    if (itype <= nstandard_types) str=trim(str)//trim(stdinfo(itype))
    if (present(info)) str=trim(str)//' Additional information: '//trim(info)
#if DEF_MODEL==MODEL_KiD
    print*, 'Runtime is:' , time
#endif

#if DEF_MODEL==MODEL_UM
    call Ereport(routine, itype, str)
#else
    print*, routine,':', trim(str)
    print*, (minus_one)**0.5
    stop
#endif
  end subroutine throw_mphys_error
end module mphys_die
