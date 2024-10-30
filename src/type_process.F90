! declares derived types used by processes
module type_process
  use variable_precision, only: wp

  implicit none

  type :: process_rate
     real(wp), allocatable :: column_data(:)
  end type process_rate

  type :: process_name
     integer :: id          ! Id for array indexing
     integer :: unique_id   ! Unique id for diagnostic identification
     character(20) :: p_name  ! Process name
     logical :: on          ! is the process going to be used, i.e. on=.true.
  end type process_name
end module type_process
