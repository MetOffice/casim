MODULE diaghelp_lem
#if DEF_MODEL==MODEL_LEM
USE variable_precision, ONLY: wp
USE type_process, ONLY: process_name, process_rate
USE process_routines, ONLY:    &
     i_cond, i_praut, i_pracw, i_pracr, i_prevp, i_psedr &
     , i_tidy, i_tidy2, i_inuc                             &
     , i_aact, i_asedr, i_arevp, i_aevp, i_atidy, i_atidy2 &
     , i_idep, i_isub, i_smlt, i_gmlt &
     , i_saut, i_gaci, i_saci, i_raci, i_homc, i_ihal
USE com_prametr, ONLY : jjp, kkp, iipep
IMPLICIT NONE

INTEGER :: i_here, j_here, k_here

INTEGER :: n_sub, n_subsed

INTEGER :: koff

INTEGER :: diaglevel = 5


INTEGER, PARAMETER :: maxndgfields = 100
INTEGER, PARAMETER :: maxindex = 200
INTEGER :: ndgfields
REAL(wp) :: dgfields(0:jjp+1,kkp,maxndgfields,0:iipep+1)
CHARACTER(100) :: dgnames(maxndgfields)

INTEGER, ALLOCATABLE :: dgindex(:,:), dgindex2(:,:)


#endif
END MODULE diaghelp_lem

