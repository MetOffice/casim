MODULE Preconditioning

USE variable_precision, ONLY: wp
USE passive_fields, ONLY: qws
USE mphys_switches, ONLY: i_qv, i_ql, i_qr, i_qi, i_qs, i_qg   &
     , cloud_params, rain_params, ice_params, snow_params, &
     graupel_params
USE thresholds, ONLY: thresh_tidy
IMPLICIT NONE

LOGICAL, ALLOCATABLE :: precondition(:)

CONTAINS

SUBROUTINE preconditioner(qfields)

REAL(wp), INTENT(IN) :: qfields(:,:)

INTEGER :: k
LOGICAL :: l_temp

DO k=UBOUND(precondition,1)-1,1,-1
      ! Do we have any existing hydrometeor mass?
  l_temp=.FALSE.
  IF (cloud_params%l_1m)     &
         l_temp = l_temp .OR. qfields(k, i_ql) > thresh_tidy(i_ql)
  IF (rain_params%l_1m)     &
         l_temp = l_temp    .OR. qfields(k, i_qr) > thresh_tidy(i_qr)
  IF (ice_params%l_1m)     &
         l_temp = l_temp    .OR. qfields(k, i_qi) > thresh_tidy(i_qi)
  IF (snow_params%l_1m)     &
         l_temp = l_temp    .OR. qfields(k, i_qs) > thresh_tidy(i_qs)
  IF (graupel_params%l_1m)     &
         l_temp = l_temp    .OR. qfields(k, i_qg) > thresh_tidy(i_qg)
      ! Do we have supersaturation
  l_temp = l_temp .OR. qws(k) < qfields(k, i_qv)
      ! Do we meet heterogeneous freezing condtion
      ! Need to add this for ice phase...
      ! l_temp = l_temp .or. Si > 0.25
      ! Do we have something above which might fall down
  l_temp =l_temp .OR. precondition(k+1)

      ! qsat doesn't work at very low pressures,
      ! so if qsaturation is 0.0 then don't do microphysics
  IF (qws(k) <= 1.0e-6)l_temp = .FALSE.
  IF (qfields(k,i_qv) <= 3.0e-6)l_temp = .FALSE.

      ! OK, that's all...
  precondition(k) = l_temp
END DO


END SUBROUTINE preconditioner

END MODULE Preconditioning
