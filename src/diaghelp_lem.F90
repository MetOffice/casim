module diaghelp_lem
#if DEF_MODEL==MODEL_LEM
  use variable_precision, only: wp
  use type_process, only: process_name, process_rate
  use process_routines, only: i_cond, i_praut, i_pracw, i_pracr, i_prevp, i_psedr, i_tidy, i_tidy2, i_inuc, &
       i_aact, i_asedr, i_arevp, i_aevp, i_atidy, i_atidy2, i_idep, i_isub, i_smlt, i_gmlt, i_saut, i_gaci, &
       i_saci, i_raci, i_homc, i_ihal
  use com_prametr, only : jjp, kkp, iipep

  implicit none
  integer :: i_here, j_here, k_here
  integer :: n_sub, n_subsed
  integer :: koff
  integer :: diaglevel = 5

  integer, parameter :: maxndgfields = 100
  integer, parameter :: maxindex = 200
  integer :: ndgfields
  real(wp) :: dgfields(0:jjp+1,kkp,maxndgfields,0:iipep+1)
  character(100) :: dgnames(maxndgfields)

  integer, allocatable :: dgindex(:,:), dgindex2(:,:)
#endif
end module diaghelp_lem

