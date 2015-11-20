MODULE PRECISION

IMPLICIT NONE

REAL :: dummy

INTEGER, PARAMETER ::    &
     defp = KIND(dummy) ! default precision

INTEGER, PARAMETER ::                       &
       OneByteInt = SELECTED_INT_KIND(2)    &
       ,TwoByteInt = SELECTED_INT_KIND(4)   &
       ,FourByteInt = SELECTED_INT_KIND(9)  &
       ,EightByteInt = SELECTED_INT_KIND(18)

INTEGER, PARAMETER ::                                       &
       FourByteReal = SELECTED_REAL_KIND(P =  6, R =  37)   &
       ,EightByteReal = SELECTED_REAL_KIND(P = 13, R =  307)

INTEGER, PARAMETER ::                   &
!       wp=FourByteReal                  & ! real working precision
       wp=EightByteReal                  & ! real working precision
       ,iwp=EightByteInt                & ! integer working precision
       ,ncdfp=FourByteReal              & ! netcdf real precision
       ,incdfp=FourByteInt                ! netcdf integer precision

INTEGER, PARAMETER ::                   &
       sp=FourByteReal                  &
       ,dp=EightByteReal


END MODULE PRECISION
