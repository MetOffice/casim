MODULE process_routines

USE variable_precision, ONLY: wp

USE type_process, ONLY: process_name, process_rate
USE mphys_parameters, ONLY: nprocs, naeroprocs, hydro_params, &
   ZERO_REAL_WP

IMPLICIT NONE

  ! NB These ids are now overwritten in mphys_switches to only allocate those
  ! that are needed given namelist switches

TYPE(process_name) :: i_cond  = process_name(0, 1, 'pcond', on=.FALSE.)
TYPE(process_name) :: i_praut = process_name(0, 2, 'praut', on=.FALSE.)
TYPE(process_name) :: i_pracw = process_name(0, 3, 'pracw', on=.FALSE.)
TYPE(process_name) :: i_pracr = process_name(0, 4, 'pracr', on=.FALSE.)
TYPE(process_name) :: i_prevp = process_name(0, 5, 'prevp', on=.FALSE.)
TYPE(process_name) :: i_psedl = process_name(0, 6, 'psedl', on=.FALSE.)
TYPE(process_name) :: i_psedr = process_name(0, 7, 'psedr', on=.FALSE.)
TYPE(process_name) :: i_tidy  = process_name(0, 8, 'ptidy', on=.FALSE.)
TYPE(process_name) :: i_tidy2 = process_name(0, 9, 'ptidy2', on=.FALSE.)
TYPE(process_name) :: i_inuc  = process_name(0, 10, 'pinuc', on=.FALSE.)
TYPE(process_name) :: i_idep  = process_name(0, 11, 'pidep', on=.FALSE.)
TYPE(process_name) :: i_iacw  = process_name(0, 12, 'piacw', on=.FALSE.)
TYPE(process_name) :: i_saut  = process_name(0, 13, 'psaut', on=.FALSE.)
TYPE(process_name) :: i_sdep  = process_name(0, 14, 'psdep', on=.FALSE.)
TYPE(process_name) :: i_sacw  = process_name(0, 15, 'psacw', on=.FALSE.)
TYPE(process_name) :: i_gdep  = process_name(0, 16, 'pgdep', on=.FALSE.)
TYPE(process_name) :: i_pseds = process_name(0, 17, 'pseds', on=.FALSE.)
TYPE(process_name) :: i_psedi = process_name(0, 18, 'psedi', on=.FALSE.)
TYPE(process_name) :: i_psedg = process_name(0, 19, 'psedg', on=.FALSE.)
TYPE(process_name) :: i_saci  = process_name(0, 20, 'psaci', on=.FALSE.)
TYPE(process_name) :: i_raci  = process_name(0, 21, 'praci', on=.FALSE.)
TYPE(process_name) :: i_sacr  = process_name(0, 22, 'psacr', on=.FALSE.)
TYPE(process_name) :: i_gacr  = process_name(0, 23, 'pgacr', on=.FALSE.)
TYPE(process_name) :: i_gacw  = process_name(0, 24, 'pgacw', on=.FALSE.)
TYPE(process_name) :: i_gaci  = process_name(0, 25, 'pgaci', on=.FALSE.)
TYPE(process_name) :: i_gacs  = process_name(0, 26, 'pgacs', on=.FALSE.)
TYPE(process_name) :: i_iagg  = process_name(0, 27, 'piagg', on=.FALSE.)
TYPE(process_name) :: i_sagg  = process_name(0, 28, 'psagg', on=.FALSE.)
TYPE(process_name) :: i_gagg  = process_name(0, 29, 'pgagg', on=.FALSE.)
TYPE(process_name) :: i_sbrk  = process_name(0, 30, 'psbrk', on=.FALSE.)
TYPE(process_name) :: i_gshd  = process_name(0, 31, 'pgshd', on=.FALSE.)
TYPE(process_name) :: i_ihal  = process_name(0, 32, 'pihal', on=.FALSE.)
TYPE(process_name) :: i_smlt  = process_name(0, 33, 'psmlt', on=.FALSE.)
TYPE(process_name) :: i_gmlt  = process_name(0, 34, 'pgmlt', on=.FALSE.)
TYPE(process_name) :: i_homr  = process_name(0, 35, 'phomr', on=.FALSE.)
TYPE(process_name) :: i_homc  = process_name(0, 36, 'phomc', on=.FALSE.)
TYPE(process_name) :: i_ssub  = process_name(0, 37, 'pssub', on=.FALSE.)
TYPE(process_name) :: i_gsub  = process_name(0, 38, 'pgsub', on=.FALSE.)
TYPE(process_name) :: i_isub  = process_name(0, 39, 'pisub', on=.FALSE.)
TYPE(process_name) :: i_imlt  = process_name(0, 40, 'pimlt', on=.FALSE.)
! aerosol processes
TYPE(process_name)  :: i_aact  = process_name(0, 101, 'aact', on=.FALSE.)
TYPE(process_name)  :: i_aaut  = process_name(0, 102, 'aaut', on=.FALSE.)
TYPE(process_name)  :: i_aacw  = process_name(0, 103, 'aacw', on=.FALSE.)
TYPE(process_name)  :: i_aevp  = process_name(0, 104, 'aevp', on=.FALSE.)
TYPE(process_name)  :: i_asedr = process_name(0, 105, 'asedr', on=.FALSE.)
TYPE(process_name)  :: i_arevp = process_name(0, 106, 'arevp', on=.FALSE.)
TYPE(process_name)  :: i_asedl = process_name(0, 107, 'asedl', on=.FALSE.)
  !... additional tidying processes (Need to sort out location for these)
TYPE(process_name)  :: i_atidy = process_name(0, 108, 'atidy', on=.FALSE.)
TYPE(process_name)  :: i_atidy2 = process_name(0, 109, 'atidy2', on=.FALSE.)
  !... ice related processes
TYPE(process_name)  :: i_dnuc  = process_name(0, 110, 'dnuc', on=.FALSE.)
TYPE(process_name)  :: i_dsub  = process_name(0, 111, 'dsub', on=.FALSE.)
TYPE(process_name)  :: i_dsedi = process_name(0, 112, 'dsedi', on=.FALSE.)
TYPE(process_name)  :: i_dseds = process_name(0, 113, 'dseds', on=.FALSE.)
TYPE(process_name)  :: i_dsedg = process_name(0, 114, 'dsedg', on=.FALSE.)
TYPE(process_name)  :: i_dssub  = process_name(0, 115, 'dssub', on=.FALSE.)
TYPE(process_name)  :: i_dgsub  = process_name(0, 116, 'dgsub', on=.FALSE.)
TYPE(process_name)  :: i_dhomc  = process_name(0, 117, 'dhomc', on=.FALSE.)
TYPE(process_name)  :: i_dhomr  = process_name(0, 118, 'dhomr', on=.FALSE.)
TYPE(process_name)  :: i_dimlt  = process_name(0, 119, 'dimlt', on=.FALSE.)
TYPE(process_name)  :: i_dsmlt  = process_name(0, 120, 'dsmlt', on=.FALSE.)
TYPE(process_name)  :: i_dgmlt  = process_name(0, 121, 'dgmlt', on=.FALSE.)
TYPE(process_name)  :: i_diacw  = process_name(0, 122, 'diacw', on=.FALSE.)
TYPE(process_name)  :: i_dsacw  = process_name(0, 123, 'dsacw', on=.FALSE.)
TYPE(process_name)  :: i_dgacw  = process_name(0, 124, 'dgacw', on=.FALSE.)
TYPE(process_name)  :: i_dsacr  = process_name(0, 125, 'dsacr', on=.FALSE.)
TYPE(process_name)  :: i_dgacr  = process_name(0, 126, 'dgacr', on=.FALSE.)
TYPE(process_name)  :: i_draci  = process_name(0, 127, 'draci', on=.FALSE.)

CONTAINS

SUBROUTINE allocate_procs(procs, nz, nprocs, ntotalq)

    ! Allocate space to store the microphysical process rates

TYPE(process_rate), INTENT(INOUT) :: procs(:,:)
INTEGER, INTENT(IN) :: nz      ! number of height levels
INTEGER, INTENT(IN) :: nprocs  ! number of physical processes
INTEGER, INTENT(IN) :: ntotalq ! number of q or aerosol fields
INTEGER :: iproc, k

DO iproc=1,nprocs
  DO k=1,nz
    ALLOCATE(procs(k,iproc)%source(ntotalq))
  END DO
END DO

CALL zero_procs(procs)

END SUBROUTINE allocate_procs

SUBROUTINE zero_procs(procs)

TYPE(process_rate), INTENT(INOUT) :: procs(:,:)
INTEGER :: iproc, k

INTEGER :: lb1, lb2, ub1, ub2

lb1=LBOUND(procs,1)
ub1=UBOUND(procs,1)
lb2=LBOUND(procs,2)
ub2=UBOUND(procs,2)
DO iproc=lb2,ub2
  DO k=lb1,ub1
    procs(k,iproc)%source(:)=0.0
  END DO
END DO

END SUBROUTINE zero_procs

SUBROUTINE deallocate_procs(procs)

TYPE(process_rate), INTENT(INOUT) :: procs(:,:)
INTEGER :: k, iproc

DO iproc=LBOUND(procs,2),UBOUND(procs,2)
  DO k=LBOUND(procs,1),UBOUND(procs,1)
!! Allocatable arrays are quicker/safer than pointers, but
!! not supported by all compilers
#if ALLOCATABLE_TYPE == 1
    IF (ALLOCATED(procs(k,iproc)%source))DEALLOCATE(procs(k,iproc)%source)
#else
    IF (ASSOCIATED(procs(k,iproc)%source))DEALLOCATE(procs(k,iproc)%source)
#endif
  END DO
END DO

END SUBROUTINE deallocate_procs

END MODULE process_routines
