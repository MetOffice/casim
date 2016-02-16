module variable_precision
  implicit none

  real :: dummy

  integer, parameter :: defp = kind(dummy) ! default precision

  integer, parameter :: OneByteInt = selected_int_kind(2), TwoByteInt = selected_int_kind(4), &
       FourByteInt = selected_int_kind(9), EightByteInt = selected_int_kind(18)

  integer, parameter :: FourByteReal = selected_real_kind(P =  6, R =  37),EightByteReal = selected_real_kind(P = 13, R =  307)

  integer, parameter :: wp=EightByteReal                  & ! real working precision
       ,iwp=EightByteInt                & ! integer working precision
       ,ncdfp=FourByteReal              & ! netcdf real precision
       ,incdfp=FourByteInt                ! netcdf integer precision

  integer, parameter :: sp=FourByteReal, dp=EightByteReal
end module variable_precision
