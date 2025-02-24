module process_routines
  use variable_precision, only: wp
  use type_process, only: process_name, process_rate
  use mphys_parameters, only: hydro_params, ZERO_REAL_WP

  implicit none

  character(len=*), parameter, private :: ModuleName='PROCESS_ROUTINES'

  ! NB These ids are now overwritten in mphys_switches to only allocate those
  ! that are needed given namelist switches
  type(process_name) :: i_cond  = process_name(0, 1, 'pcond', on=.false.)
  type(process_name) :: i_praut = process_name(0, 2, 'praut', on=.false.)
  type(process_name) :: i_pracw = process_name(0, 3, 'pracw', on=.false.)
  type(process_name) :: i_pracr = process_name(0, 4, 'pracr', on=.false.)
  type(process_name) :: i_prevp = process_name(0, 5, 'prevp', on=.false.)
  type(process_name) :: i_psedl = process_name(0, 6, 'psedl', on=.false.)
  type(process_name) :: i_psedr = process_name(0, 7, 'psedr', on=.false.)
  type(process_name) :: i_tidy  = process_name(0, 8, 'ptidy', on=.false.)
  type(process_name) :: i_tidy2 = process_name(0, 9, 'ptidy2', on=.false.)
  type(process_name) :: i_inuc  = process_name(0, 10, 'pinuc', on=.false.)
  type(process_name) :: i_idep  = process_name(0, 11, 'pidep', on=.false.)
  type(process_name) :: i_iacw  = process_name(0, 12, 'piacw', on=.false.)
  type(process_name) :: i_saut  = process_name(0, 13, 'psaut', on=.false.)
  type(process_name) :: i_sdep  = process_name(0, 14, 'psdep', on=.false.)
  type(process_name) :: i_sacw  = process_name(0, 15, 'psacw', on=.false.)
  type(process_name) :: i_gdep  = process_name(0, 16, 'pgdep', on=.false.)
  type(process_name) :: i_pseds = process_name(0, 17, 'pseds', on=.false.)
  type(process_name) :: i_psedi = process_name(0, 18, 'psedi', on=.false.)
  type(process_name) :: i_psedg = process_name(0, 19, 'psedg', on=.false.)
  type(process_name) :: i_saci  = process_name(0, 20, 'psaci', on=.false.)
  type(process_name) :: i_raci  = process_name(0, 21, 'praci', on=.false.)
  type(process_name) :: i_sacr  = process_name(0, 22, 'psacr', on=.false.)
  type(process_name) :: i_gacr  = process_name(0, 23, 'pgacr', on=.false.)
  type(process_name) :: i_gacw  = process_name(0, 24, 'pgacw', on=.false.)
  type(process_name) :: i_gaci  = process_name(0, 25, 'pgaci', on=.false.)
  type(process_name) :: i_gacs  = process_name(0, 26, 'pgacs', on=.false.)
  type(process_name) :: i_iagg  = process_name(0, 27, 'piagg', on=.false.)
  type(process_name) :: i_sagg  = process_name(0, 28, 'psagg', on=.false.)
  type(process_name) :: i_gagg  = process_name(0, 29, 'pgagg', on=.false.)
  type(process_name) :: i_sbrk  = process_name(0, 30, 'psbrk', on=.false.)
  type(process_name) :: i_gshd  = process_name(0, 31, 'pgshd', on=.false.)
  type(process_name) :: i_ihal  = process_name(0, 32, 'pihal', on=.false.)
  type(process_name) :: i_smlt  = process_name(0, 33, 'psmlt', on=.false.)
  type(process_name) :: i_gmlt  = process_name(0, 34, 'pgmlt', on=.false.)
  type(process_name) :: i_homr  = process_name(0, 35, 'phomr', on=.false.)
  type(process_name) :: i_homc  = process_name(0, 36, 'phomc', on=.false.)
  type(process_name) :: i_ssub  = process_name(0, 37, 'pssub', on=.false.)
  type(process_name) :: i_gsub  = process_name(0, 38, 'pgsub', on=.false.)
  type(process_name) :: i_isub  = process_name(0, 39, 'pisub', on=.false.)
  type(process_name) :: i_imlt  = process_name(0, 40, 'pimlt', on=.false.)
  type(process_name) :: i_iics  = process_name(0, 41, 'piics', on=.false.)
  type(process_name) :: i_idps  = process_name(0, 42, 'pidps', on=.false.)
  ! aerosol processes
  type(process_name)  :: i_aact  = process_name(0, 101, 'aact', on=.false.)
  type(process_name)  :: i_aaut  = process_name(0, 102, 'aaut', on=.false.)
  type(process_name)  :: i_aacw  = process_name(0, 103, 'aacw', on=.false.)
  type(process_name)  :: i_aevp  = process_name(0, 104, 'aevp', on=.false.)
  type(process_name)  :: i_asedr = process_name(0, 105, 'asedr', on=.false.)
  type(process_name)  :: i_arevp = process_name(0, 106, 'arevp', on=.false.)
  type(process_name)  :: i_asedl = process_name(0, 107, 'asedl', on=.false.)
  !... additional tidying processes (Need to sort out location for these)
  type(process_name)  :: i_atidy = process_name(0, 108, 'atidy', on=.false.)
  type(process_name)  :: i_atidy2 = process_name(0, 109, 'atidy2', on=.false.)
  !... ice related processes
  type(process_name)  :: i_dnuc  = process_name(0, 110, 'dnuc', on=.false.)
  type(process_name)  :: i_dsub  = process_name(0, 111, 'dsub', on=.false.)
  type(process_name)  :: i_dsedi = process_name(0, 112, 'dsedi', on=.false.)
  type(process_name)  :: i_dseds = process_name(0, 113, 'dseds', on=.false.)
  type(process_name)  :: i_dsedg = process_name(0, 114, 'dsedg', on=.false.)
  type(process_name)  :: i_dssub  = process_name(0, 115, 'dssub', on=.false.)
  type(process_name)  :: i_dgsub  = process_name(0, 116, 'dgsub', on=.false.)
  type(process_name)  :: i_dhomc  = process_name(0, 117, 'dhomc', on=.false.)
  type(process_name)  :: i_dhomr  = process_name(0, 118, 'dhomr', on=.false.)
  type(process_name)  :: i_dimlt  = process_name(0, 119, 'dimlt', on=.false.)
  type(process_name)  :: i_dsmlt  = process_name(0, 120, 'dsmlt', on=.false.)
  type(process_name)  :: i_dgmlt  = process_name(0, 121, 'dgmlt', on=.false.)
  type(process_name)  :: i_diacw  = process_name(0, 122, 'diacw', on=.false.)
  type(process_name)  :: i_dsacw  = process_name(0, 123, 'dsacw', on=.false.)
  type(process_name)  :: i_dgacw  = process_name(0, 124, 'dgacw', on=.false.)
  type(process_name)  :: i_dsacr  = process_name(0, 125, 'dsacr', on=.false.)
  type(process_name)  :: i_dgacr  = process_name(0, 126, 'dgacr', on=.false.)
  type(process_name)  :: i_draci  = process_name(0, 127, 'draci', on=.false.)

contains

  ! Allocate space to store the microphysical process rates
  subroutine allocate_procs(nxy_inner, procs, nz, nprocs, ntotalq)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ALLOCATE_PROCS'

    integer, intent(in) :: nxy_inner
    type(process_rate), intent(inout) :: procs(:,:,:)
    integer, intent(in) :: nz      ! number of height levels
    integer, intent(in) :: nprocs  ! number of physical processes
    integer, intent(in) :: ntotalq ! number of q or aerosol fields

    integer :: iproc, iq, ixy

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do ixy=1, nxy_inner
     do iproc=1, nprocs
       do iq=1,ntotalq
         allocate(procs(iq,iproc,ixy)%column_data(nz))
       end do
       end do
     end do
    
    do ixy=1, nxy_inner
       call zero_procs(procs(:,:,ixy))
    end do

    !  call zero_procs(procs)

     IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

   end subroutine allocate_procs

  subroutine zero_procs(procs, iprocs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='ZERO_PROCS'

    type(process_rate), intent(inout) :: procs(:,:)
    type(process_name), intent(in), optional :: iprocs(:)

    integer :: iproc, nproc, iq, i
    integer :: lb1, lb2, ub1, ub2

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    lb1=lbound(procs,1)
    ub1=ubound(procs,1)
    
    if(present(iprocs)) then 
       nproc = size(iprocs)
       do i = 1, nproc
          iproc = iprocs(i)%id
          do iq=lb1,ub1
             procs(iq,iproc)%column_data(:)=0.0
          enddo
       enddo
    else
       lb2=lbound(procs,2)
       ub2=ubound(procs,2)
       do iproc=lb2, ub2
          do iq=lb1,ub1
             procs(iq,iproc)%column_data(:)=0.0
          end do
       end do
    endif
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine zero_procs
 

  subroutine deallocate_procs(nxy_inner, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='DEALLOCATE_PROCS'

    integer, intent(in) :: nxy_inner
    type(process_rate), intent(inout) :: procs(:,:,:)

    integer :: iproc, iq, ixy

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do ixy=1, nxy_inner
    do iproc=lbound(procs,2), ubound(procs,2)
      do iq=lbound(procs,1), ubound(procs,1)
             if (allocated(procs(iq,iproc,ixy)%column_data)) deallocate(procs(iq,iproc,ixy)%column_data)
          end do
      end do
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine deallocate_procs
end module process_routines
