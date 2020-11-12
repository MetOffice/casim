! Routine to exit the microphysics following an error there
module mphys_die

  ! Module for the KiD and MONC models - UM has its own version

! #if DEF_MODEL==MODEL_KiD
!   use runtime, only: time
! #endif

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='MPHYS_DIE'

  integer, parameter :: incorrect_opt = 1
  integer, parameter :: bad_values    = 2
  integer, parameter :: warn          = -1
  integer, parameter :: std_msg_len   = 400 ! Length of a standard message

  character(len=std_msg_len) :: std_msg = ''

  public throw_mphys_error, incorrect_opt, bad_values, warn, std_msg, mphys_message
  
contains

  subroutine throw_mphys_error(itype, routine, info)

    ! If modifying the subroutine or argument list, ensure that the 
    ! UM version of mphys_die is also modified to give the same answers

    implicit none

    integer,      intent(in) :: itype  ! type of error 1 = Incorrect specification of options
                                       !               2 = Bad values found
                                       !               3 = Unknown error
                                       !              <0 = Warning, code will continue
    character(*), intent(in) :: routine
    character(*), intent(in) :: info ! error information

    integer, parameter :: nstandard_types=3
    character(100) :: stdinfo(nstandard_types) =     &
         (/ 'Incorrect specification of options      ' &
         ,  'Bad values found                        ' &
         ,  'Unknown error                           ' &
         /)
    character(1000) :: str
    real :: minus_one=-1.

    character(len=*), parameter :: RoutineName='THROW_MPHYS_ERROR'

    if (itype > 0) then

      !------------------------------------------------------------
      ! Produce error message
      !------------------------------------------------------------

      str='Error in CASIM microphysics: '

      if ( itype <= 3 ) then
        str=trim(str)//trim(stdinfo(itype))
      else
        str=trim(str)//trim(stdinfo(3))
      end if

      str=trim(str)//' Additional information: '//trim(info)
! #if DEF_MODEL==MODEL_KiD
!       print*, 'Runtime is:' , time
! #endif
      print*, routine,':', trim(str)
      print*, (minus_one)**0.5
      stop

    else if ( itype < 0 ) then

      !------------------------------------------------------------
      ! Produce warning message
      !------------------------------------------------------------

      str='Warning from CASIM microphysics! '
      str=trim(str)//' Message: '//trim(info)
! #if DEF_MODEL==MODEL_KiD
!       print*, 'Runtime is:' , time
! #endif
      print*, routine,':', trim(str)

    end if

  end subroutine throw_mphys_error


  subroutine mphys_message(routine, msg)

    implicit none

    character(*), intent(in) :: routine ! Routine providing the message
    character(*), intent(in) :: msg ! Message

    character(len=*), parameter :: RoutineName='MPHYS_MESSAGE'

    character(2000) :: str

    str = '| Message from CASIM microphysics | Routine:' 
    str = trim(str)//trim(routine)//' | Message: '
    str = trim(str)//trim(msg)//' |'

    print *, trim(str)

  end subroutine mphys_message

end module mphys_die
