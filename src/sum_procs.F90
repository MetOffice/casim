module sum_process
  use variable_precision, only: wp
  use mphys_die, only: throw_mphys_error, incorrect_opt, std_msg
  use type_process, only: process_name, process_rate
  use mphys_switches, only: i_th, i_qv, i_ql, i_qr, i_qs, i_qi, i_qg, l_warm, i_m3r, ntotalq, ntotala, &
       i_an11, i_am4, i_nl, i_am9, aero_index, i_nr, i_ng, i_ns, i_am7, i_am8
  use passive_fields, only: rexner, rho
  use mphys_constants, only: cp, Lv, Ls
  use mphys_parameters, only: hydro_params, snow_params, rain_params, graupel_params, parent_dt, ZERO_REAL_WP
  use m3_incs, only: m3_inc_type2

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='SUM_PROCESS'

  public sum_aprocs, sum_procs
contains

  subroutine sum_aprocs(dst, procs, tend, iprocs, names, afields)

    implicit none

    character(len=*), parameter :: RoutineName='SUM_APROCS'

    real(wp), intent(in) :: dst  ! step length (s)
    type(process_rate), intent(in) :: procs(:,:)
    type(process_name), intent(in) :: iprocs(:)
    real(wp), intent(inout) :: tend(:,:)
    character(20), intent(in) :: names(:)
    real(wp), intent(in), optional :: afields(:,:) ! Required for debugging

    real(wp), allocatable :: tend_temp(:,:) ! Temporary storage for accumulated tendendies

    integer :: k, iq, iproc, i, nz, idgproc
    integer :: nproc, index, idg

    character(2) :: char
    character(100) :: name

    allocate(tend_temp(lbound(tend,1):ubound(tend,1), lbound(tend,2):ubound(tend,2)))
    tend_temp=ZERO_REAL_WP

    nproc=size(iprocs)
    nz=size(procs,1)

    do i=1, nproc
      if (iprocs(i)%on) then
        iproc=iprocs(i)%id
        do k=1,nz
          if (.not. all(procs(k,iproc)%source(:)==ZERO_REAL_WP)) then
            do iq=1, ntotala
              tend_temp(k, iq)=tend_temp(k, iq) + procs(k,iproc)%source(iq)*dst
            end do
          end if
        end do
      end if
    end do

    ! Add on tendencies to those already passed in.
    tend=tend+tend_temp
    deallocate(tend_temp)
  end subroutine sum_aprocs

  subroutine sum_procs(dst, procs, tend, iprocs, names, l_thermalexchange, i_thirdmoment, qfields, l_passive )

    implicit none

    character(len=*), parameter :: RoutineName='SUM_PROCS'

    real(wp), intent(in) :: dst  ! step length (s)
    type(process_rate), intent(in) :: procs(:,:)
    type(process_name), intent(in) :: iprocs(:)
    real(wp), intent(inout) :: tend(:,:)
    logical, intent(in), optional :: l_thermalexchange  ! Calculate the thermal exchange terms
    integer, intent(in), optional :: i_thirdmoment  ! Calculate the tendency of the third moment
    real(wp), intent(in), optional :: qfields(:,:) ! Required for debugging or with i_thirdmoment
    logical, intent(in), optional :: l_passive ! If true don't apply final tendency (testing diagnostics with pure sedimentation)
    character(10), intent(in) :: names(:)

    real(wp), allocatable :: tend_temp(:,:) ! Temporary storage for accumulated tendendies
    type(hydro_params) :: params
    real(wp) :: dm1,dm2,dm3,m1,m2,m3
    integer :: k, iq, iproc, isource, nsource, i, nz, idgproc
    integer :: nproc, index, idg

    character(2) :: char
    character(100) :: name

    logical :: do_thermal
    logical :: do_third  ! Currently not plumbed in, so explicitly calculated for each process
    logical :: do_update ! update the tendency
    integer :: third_type
    real :: qh ! total hydrometeor

    do_thermal=.false.
    if (present(l_thermalexchange)) do_thermal=l_thermalexchange

    do_third=.false.
    if (present(i_thirdmoment)) then
      do_third=.true.
      third_type=i_thirdmoment
    end if

    do_update=.true.
    if (present(l_passive)) do_update= .not. l_passive

    allocate(tend_temp(lbound(tend,2):ubound(tend,2), lbound(tend,1):ubound(tend,1)))
    tend_temp=ZERO_REAL_WP

    nproc=size(iprocs)
    nz=size(procs,1)

    do i=1, nproc
      if (iprocs(i)%on) then
        iproc=iprocs(i)%id
        do k=1, nz
          do iq=1, ntotalq
            tend_temp(iq, k)=tend_temp(iq, k)+procs(k,iproc)%source(iq)*dst
          end do

          if (do_third) then
            ! calculate increment to third moment based on collected increments from
            ! q and n NB This overwrites any previously calculated values
            ! Rain
            params=rain_params
            if (params%l_3m) then
              m1=qfields(k, params%i_1m)*rho(k)/params%c_x
              m2=qfields(k, params%i_2m)
              m3=qfields(k, params%i_3m)
              dm1=tend_temp(params%i_1m, k)*rho(k)/params%c_x
              dm2=tend_temp(params%i_2m, k)
              if (dm1 < -.99*m1 .or. dm2 < -.99*m2) then
                dm1=-m1
                dm2=-m2
                dm3=-m3
              else
                if (m3> 0.0 .and. m1 > 0.0 .and. m2 > 0.0 .and. (abs(dm1) > 0.0 .or. abs(dm2) > 0.0)) then
                  select case (i_thirdmoment)
                  case (2)
                    call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
                  case default
                    write(std_msg, '(A)') 'rain i_thirdmoment incorrectly set'
                    call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                                           std_msg )
                  end select
                  tend_temp(params%i_3m, k)=dm3
                end if
              end if
            end if
            ! Snow
            params=snow_params
            if (params%l_3m) then
              m1=qfields(k, params%i_1m)*rho(k)/params%c_x
              m2=qfields(k, params%i_2m)
              m3=qfields(k, params%i_3m)
              dm1=tend_temp(params%i_1m, k)*rho(k)/params%c_x
              dm2=tend_temp(params%i_2m, k)
              if (dm1 < -.99*m1 .or. dm2 < -.99*m2) then
                dm1=-m1
                dm2=-m2
                dm3=-m3
              else
                if (m3> 0.0 .and. m1 > 0.0 .and. m2 > 0.0 .and. (abs(dm1) > 0.0 .or. abs(dm2) > 0.0)) then
                  select case (i_thirdmoment)
                  case (2)
                    call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
                  case default
                    write(std_msg, '(A)') 'snow i_thirdmoment incorrectly set'
                    call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                                           std_msg)
                  end select
                  tend_temp(params%i_3m, k)=dm3
                end if
              end if
            end if
            ! Graupel
            params=graupel_params
            if (params%l_3m) then
              m1=qfields(k, params%i_1m)*rho(k)/params%c_x
              m2=qfields(k, params%i_2m)
              m3=qfields(k, params%i_3m)
              dm1=tend_temp(params%i_1m, k)*rho(k)/params%c_x
              dm2=tend_temp(params%i_2m, k)
              if (dm1 < -.99*m1 .or. dm2 < -.99*m2) then
                dm1=-m1
                dm2=-m2
                dm3=-m3
              else
                if (m3 > 0.0 .and. m1 > 0.0 .and. m2 > 0.0 .and. (abs(dm1) > 0.0 .or. abs(dm2) > 0.0)) then
                  select case (i_thirdmoment)
                  case (2)
                    call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
                  case default
                    write(std_msg, '(A)') 'graupel i_thirdmoment incorrectly set'
                    call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                                           std_msg)
                  end select
                  tend_temp(params%i_3m, k)=dm3
                end if
              end if
            end if

          end if

        end do
      end if
    end do

    ! Calculate the thermal exchange values
    ! (this overwrites anything that was already stored in the theta tendency)
    if (do_thermal) then
      do k=1, nz
        tend_temp(i_th, k)=(tend_temp(i_ql, k)+tend_temp(i_qr, k))*Lv/cp * rexner(k)
        if (.not. l_warm) then
          tend_temp(i_th, k)=tend_temp(i_th, k)+(tend_temp(i_qi, k)+tend_temp(i_qs, k)+tend_temp(i_qg, k))*Ls/cp *rexner(k)
        end if
      end do
    end if

    ! Add on tendencies to those already passed in.
    if (do_update) then
      do k=lbound(tend, 2), ubound(tend, 2)
        tend(:,k)=tend(:,k)+tend_temp(k,:)
      end do
    end if    
    deallocate(tend_temp)
  end subroutine sum_procs
end module sum_process
