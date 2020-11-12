module casim_runtime

use variable_precision, only: wp
use mphys_parameters,   only: zero_real_wp

implicit none
save

integer :: i_here = 0
integer :: j_here = 0
integer :: k_here = 0

! Time in casim in seconds from start. Will be updated with parent model
! timestep as run progresses.
real(wp) :: casim_time = zero_real_wp
real(wp) :: casim_smax = 1.5 ! A reasonable guess at maximum supersaturation
                             ! overwritten by the KiD namelist
real(wp) :: casim_smax_limit_time = 1.0e20 ! A very long time which is 
                                           ! overwritten by the KiD namelist
                                           

end module casim_runtime
