module Preconditioning
  use variable_precision, only: wp
  use passive_fields, only: qws
  use mphys_switches, only: i_qv, i_ql, i_qr, i_qi, i_qs, i_qg , cloud_params, &
       rain_params, ice_params, snow_params, graupel_params, l_cfrac_casim_diag_scheme
  use thresholds, only: thresh_tidy
  implicit none
  private
  character(len=*), parameter, private :: ModuleName='PRECONDITIONING'

  logical, allocatable :: precondition(:,:)

!$OMP THREADPRIVATE(precondition)

  public precondition, preconditioner
contains

  subroutine preconditioner(ixy_inner, qfields)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='PRECONDITIONER'

    integer, intent(in) :: ixy_inner

    real(wp), intent(in) :: qfields(:,:)

    integer :: k
    logical :: l_temp

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k=ubound(precondition,1)-1, 1, -1
      ! Do we have any existing hydrometeor mass?
      l_temp=.false.
      if (cloud_params%l_1m) l_temp=l_temp .or. qfields(k, i_ql) > thresh_tidy(i_ql)
      if (rain_params%l_1m) l_temp=l_temp .or. qfields(k, i_qr) > thresh_tidy(i_qr)
      if (ice_params%l_1m) l_temp=l_temp .or. qfields(k, i_qi) > thresh_tidy(i_qi)
      if (snow_params%l_1m) l_temp=l_temp .or. qfields(k, i_qs) > thresh_tidy(i_qs)
      if (graupel_params%l_1m) l_temp=l_temp .or. qfields(k, i_qg) > thresh_tidy(i_qg)
      ! Do we have supersaturation
      l_temp=l_temp .or. qws(k,ixy_inner) < qfields(k, i_qv)
      ! Do we meet heterogeneous freezing condtion
      ! Need to add this for ice phase...
      ! l_temp = l_temp .or. Si > 0.25
      ! Do we have something above which might fall down
      l_temp=l_temp .or. precondition(k+1,ixy_inner)
      l_temp = l_temp .or. l_cfrac_casim_diag_scheme !DPG - To prevent early quit from
        !CASIM if we are sub-saturated since want to allow cloud scheme to operate even
        !if we are subsaturated.
      ! qsat doesn't work at very low pressures,
      ! so if qsaturation is 0.0 then don't do microphysics
      if (qws(k,ixy_inner) <= 1.0e-6) l_temp=.false.
      if (qfields(k,i_qv) <= 3.0e-6) l_temp=.false.

      ! OK, that's all...
      precondition(k,ixy_inner)=l_temp
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine preconditioner
end module Preconditioning
