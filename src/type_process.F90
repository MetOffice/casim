! declares derived types used by processes
module type_process
  use variable_precision, only: wp

  implicit none

  !! Allocatable arrays are quicker than pointers, but
  !! not supported by all compilers
#if ALLOCATABLE_TYPE == 1
  type :: process_rate
     real(wp), allocatable :: source(:)
  end type process_rate
#else
  type :: process_rate
     real(wp), pointer :: source(:) ! tendency to apply to target q variable (/s)
  end type process_rate
#endif

  type :: process_name
     integer :: id          ! Id for array indexing
     integer :: unique_id   ! Unique id for diagnostic identification
     character(20) :: name  ! Process name
     logical :: on          ! is the process going to be used, i.e. on=.true.
  end type process_name
end module type_process
