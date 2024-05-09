module sum_process
  use variable_precision, only: wp
! use mphys_die, only: throw_mphys_error, incorrect_opt, std_msg
  use type_process, only: process_name, process_rate
  use mphys_switches, only: i_th, i_ql, i_qr, i_qs, i_qi, i_qg, l_warm, ntotalq, ntotala
  use passive_fields, only: rexner
! use passive_fields, only: rho
  use mphys_constants, only: cp, Lv, Ls
  use mphys_parameters, only: ZERO_REAL_WP
! use mphys_parameters, only: hydro_params, snow_params, rain_params, graupel_params
! use m3_incs, only: m3_inc_type2
  use process_routines, only: process_name
  
  implicit none
  private

  ! allocated and deallocated in micro_main (initialise and finalise_micromain)
  real(wp), allocatable :: tend_temp(:,:) ! Temporary storage for accumulated tendendies
  real(wp), allocatable :: aerosol_tend_temp(:,:) ! Temporary storage for accumulated aerosol tendendies

!$OMP THREADPRIVATE(tend_temp, aerosol_tend_temp)

  character(len=*), parameter, private :: ModuleName='SUM_PROCESS'

  public sum_aprocs, sum_procs,  tend_temp, aerosol_tend_temp
contains

  subroutine sum_aprocs(dst, nz, procs, tend, iprocs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SUM_APROCS'

    real(wp), intent(in) :: dst  ! step length (s)
    integer, intent(in) :: nz ! number of points in a column
    type(process_rate), intent(in) :: procs(:,:)
    type(process_name), intent(in) :: iprocs(:)
    real(wp), intent(inout) :: tend(:,:)

    !real(wp), allocatable :: tend_temp(:,:) ! Temporary storage for accumulated tendendies

    integer :: k, iq, iproc, i
    integer :: nproc

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !allocate(tend_temp(lbound(tend,1):ubound(tend,1), lbound(tend,2):ubound(tend,2)))
    aerosol_tend_temp=ZERO_REAL_WP

    nproc=size(iprocs)

    do i=1, nproc
      if (iprocs(i)%on) then
        iproc=iprocs(i)%id
        do iq=1, ntotala
           do k=1,nz
         ! if (.not. all(procs(k,iproc)%source(:)==ZERO_REAL_WP)) then
              aerosol_tend_temp(k, iq)=aerosol_tend_temp(k, iq) + &
                   procs(iq,iproc)%column_data(k)*dst
            end do
         ! end if
        end do
      end if
    end do

    ! Add on tendencies to those already passed in.
    tend=tend+aerosol_tend_temp
    !deallocate(tend_temp)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine sum_aprocs

  subroutine sum_procs(ixy_inner, dst, nz, procs, tend, iprocs, l_thermalexchange, i_thirdmoment, qfields, l_passive )

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SUM_PROCS'

    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dst  ! step length (s)
    integer, intent(in) :: nz ! number of points in a column
    type(process_rate), intent(in) :: procs(:,:)
    type(process_name), intent(in) :: iprocs(:)
    real(wp), intent(inout) :: tend(:,:)
    logical, intent(in), optional :: l_thermalexchange  ! Calculate the thermal exchange terms
    integer, intent(in), optional :: i_thirdmoment  ! Calculate the tendency of the third moment
    real(wp), intent(in), optional :: qfields(:,:) ! Required for debugging or with i_thirdmoment
    logical, intent(in), optional :: l_passive ! If true don't apply final tendency (testing diagnostics with pure sedimentation)

    !real(wp), allocatable :: tend_temp(:,:) ! Temporary storage for accumulated tendendies
!   type(hydro_params) :: params
!   real(wp) :: dm1,dm2,dm3,m1,m2,m3
    integer :: k, iq, iproc, i
    integer :: nproc

    logical :: do_thermal
!   logical :: do_third  ! Currently not plumbed in, so explicitly calculated for each process
    logical :: do_update ! update the tendency
!   integer :: third_type

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do_thermal=.false.
    if (present(l_thermalexchange)) do_thermal=l_thermalexchange

    ! do_third=.false.
    ! third_type = 0 ! Set up a default value
    ! if (present(i_thirdmoment)) then
    !   do_third=.true.
    !   third_type=i_thirdmoment
    ! end if

    do_update=.true.
    if (present(l_passive)) do_update= .not. l_passive

    ! this allocation is based on the shape of tend, which is (k, proc) originally.
    ! switch tend_temp to be the same. Should this allocation be done here? 
    !allocate(tend_temp(lbound(tend,1):ubound(tend,1), lbound(tend,2):ubound(tend,2)))
    tend_temp=ZERO_REAL_WP

    nproc=size(iprocs)

    do i=1, nproc
       iproc=iprocs(i)%id
        do iq=1, ntotalq
          do k = 1, nz
             tend_temp(k, iq)=tend_temp(k, iq)+procs(iq,iproc)%column_data(k)*dst
          enddo
       enddo
    enddo
  
    ! if (do_third) then
    !    ! calculate increment to third moment based on collected increments from
    !    ! q and n NB This overwrites any previously calculated values
    !    ! Rain
    !    params=rain_params
    !    if (params%l_3m) then
    !       do k = 1, nz
    !          m1=qfields(k, params%i_1m)*rho(k)/params%c_x
    !          m2=qfields(k, params%i_2m)
    !          m3=qfields(k, params%i_3m)
    !          dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
    !          dm2=tend_temp(k, params%i_2m)
    !          if (dm1 < -.99*m1 .or. dm2 < -.99*m2) then
    !             dm1=-m1
    !             dm2=-m2
    !             dm3=-m3
    !          else
    !             if (m3> 0.0 .and. m1 > 0.0 .and. m2 > 0.0 .and. (abs(dm1) > 0.0 .or. abs(dm2) > 0.0)) then
    !                select case (third_type)
    !                case (2)
    !                   call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
    !                case default
    !                   write(std_msg, '(A)') 'rain i_thirdmoment incorrectly set'
    !                   call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
    !                        std_msg )
    !                end select
    !                tend_temp(k, params%i_3m)=dm3
    !             end if
    !          end if
    !       enddo
    !    end if
       
    !    ! Snow
    !    params=snow_params
    !    if (params%l_3m) then
    !       do k = 1, nz
    !          m1=qfields(k, params%i_1m)*rho(k)/params%c_x
    !          m2=qfields(k, params%i_2m)
    !          m3=qfields(k, params%i_3m)
    !          dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
    !          dm2=tend_temp(k, params%i_2m)
    !          if (dm1 < -.99*m1 .or. dm2 < -.99*m2) then
    !             dm1=-m1
    !             dm2=-m2
    !             dm3=-m3
    !          else
    !             if (m3> 0.0 .and. m1 > 0.0 .and. m2 > 0.0 .and. (abs(dm1) > 0.0 .or. abs(dm2) > 0.0)) then
    !                select case (third_type)
    !                case (2)
    !                   call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
    !                case default
    !                   write(std_msg, '(A)') 'snow i_thirdmoment incorrectly set'
    !                   call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
    !                        std_msg)
    !                end select
    !                tend_temp(k, params%i_3m)=dm3
    !             end if
    !          end if
    !       enddo
    !    end if
       
    !    ! Graupel
    !    params=graupel_params
    !    if (params%l_3m) then
    !       do k = 1, nz
    !          m1=qfields(k, params%i_1m)*rho(k)/params%c_x
    !          m2=qfields(k, params%i_2m)
    !          m3=qfields(k, params%i_3m)
    !          dm1=tend_temp(k,params%i_1m)*rho(k)/params%c_x
    !          dm2=tend_temp(k,params%i_2m)
    !          if (dm1 < -.99*m1 .or. dm2 < -.99*m2) then
    !             dm1=-m1
    !             dm2=-m2
    !             dm3=-m3
    !          else
    !             if (m3 > 0.0 .and. m1 > 0.0 .and. m2 > 0.0 .and. (abs(dm1) > 0.0 .or. abs(dm2) > 0.0)) then
    !                select case (third_type)
    !                case (2)
    !                   call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
    !                case default
    !                   write(std_msg, '(A)') 'graupel i_thirdmoment incorrectly set'
    !                   call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
    !                        std_msg)
    !                end select
    !                tend_temp(k,params%i_3m)=dm3
    !             end if
    !          end if
    !       enddo
    !    end if

    ! endif ! endif for l_dothird
    
    ! Calculate the thermal exchange values
    ! (this overwrites anything that was already stored in the theta tendency)
    if (do_thermal) then
       do k=1,nz
          tend_temp(k, i_th)=(tend_temp(k, i_ql)+tend_temp(k,i_qr))*Lv/cp * rexner(k,ixy_inner)
       enddo
       if (.not. l_warm) then 
          do k=1, nz
             tend_temp(k, i_th)=tend_temp(k,i_th)+ &
               (tend_temp(k, i_qi)+tend_temp(k, i_qs)+tend_temp(k,i_qg))*Ls/cp *rexner(k,ixy_inner)
          end do
       endif
    end if

    if (do_update) then
      tend = tend+tend_temp
    endif

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine sum_procs
end module sum_process
