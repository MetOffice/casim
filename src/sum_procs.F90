module sum_process
  use variable_precision, only: wp
  use mphys_die, only: throw_mphys_error
  use type_process, only: process_name, process_rate
  use mphys_switches, only: i_th, i_qv, i_ql, i_qr, i_qs, i_qi, i_qg, l_warm, i_m3r, ntotalq, ntotala, &
       i_an11, i_am4, i_nl, i_am9, aero_index, i_nr, i_ng, i_ns, i_am7, i_am8, nsubsteps
  use passive_fields, only: rexner, rho
  use mphys_constants, only: cp, Lv, Ls
  use mphys_parameters, only: hydro_params, snow_params, rain_params, graupel_params, parent_dt, ZERO_REAL_WP
  use m3_incs, only: m3_inc_type2

#if DEF_MODEL==MODEL_KiD
  use parameters, only: diaglevel, nx
  use diagnostics, only: save_dg, i_dgtime, i_here, j_here, n_sub, n_subsed
  use runtime, only: time
#elif DEF_MODEL==MODEL_LEM_DIAG
  use diaghelp_lem, only: dgfields, i_here, j_here, koff , diaglevel, dgindex, dgindex2, n_sub, n_subsed
  use com_params, only: time
  use extra_dgs
#elif DEF_MODEL==MODEL_UM
  use mphys_casim_diagnostics, only: ProcessRates, nProcessDiags, ProcessKeys, ProcessQs, PhaseChanges
  use diaghelp_um, only: i_here, j_here, l_debug_um, debug_i, debug_j, debug_pe, n_sub, n_subsed
  use UM_ParCore, only: mype
  use timestep_mod, only: timestep_number
#elif  DEF_MODEL==MODEL_MONC
  use diaghelp_monc, only: i_here, j_here
#endif

  implicit none
  private

#if DEF_MODEL==MODEL_UM
  integer :: idg_proc

  ! Temporary for debugging
  integer :: max_prints(10)=0
#endif

  public sum_aprocs, sum_procs
contains

  subroutine sum_aprocs(dst, procs, tend, iprocs, names, afields)
    real(wp), intent(in) :: dst  ! step length (s)
    type(process_rate), intent(in) :: procs(:,:)
    type(process_name), intent(in) :: iprocs(:)
    real(wp), intent(inout) :: tend(:,:)
    character(20), intent(in) :: names(:)
    real(wp), intent(in), optional :: afields(:,:) ! Required for debugging

    real(wp), allocatable :: tend_temp(:,:) ! Temporary storage for accumulated tendendies

    integer :: k, iq, iproc, i, nz, idgproc
    integer :: nproc, index, idg

#if DEF_MODEL==MODEL_UM
    integer :: diaglevel
#endif

    character(2) :: char
    character(100) :: name

#if DEF_MODEL==MODEL_UM
    diaglevel=0
#endif

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

#if DEF_MODEL==MODEL_LEM_DIAG
              if (diaglevel > 4 .and. procs(k,iproc)%source(iq) /= ZERO_REAL_WP) then ! Slows model significantly so
                ! Only really used for development
                ! This only works without substepping.
                name=trim(iprocs(i)%name)//'_'//trim(adjustl(names(iq)))

                idg=req_dgproc(trim(name))
                if (idg > 0 .and. j_here<=jenddg .and. j_here>=jstartdg) then
                  dgprocs(j_here-jstartdg+1, k, i_here, idg)=procs(k,iproc)%source(iq)
                end if
              end if
#elif DEF_MODEL==MODEL_KiD
              if (diaglevel > 4 .and. procs(k,iproc)%source(iq) /= ZERO_REAL_WP) then ! Slows model significantly so
                name=trim(iprocs(i)%name)//'_'//trim(adjustl(names(iq)))
                if (nx==1) then
                  call save_dg(k, procs(k, iproc)%source(iq), name, i_dgtime)
                else
                  call save_dg(k, i_here, procs(k, iproc)%source(iq), name, i_dgtime)
                end if
              end if
#endif
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

#if DEF_MODEL==MODEL_UM
    integer :: diaglevel
#endif

    character(2) :: char
    character(100) :: name

    logical :: do_thermal
    logical :: do_third  ! Currently not plumbed in, so explicitly calculated for each process
    logical :: do_update ! update the tendency
    integer :: third_type
    real :: qh ! total hydrometeor

#if DEF_MODEL==MODEL_UM
    diaglevel = 0
#endif

    do_thermal=.false.
    if (present(l_thermalexchange)) do_thermal=l_thermalexchange

    do_third=.false.
    if (present(i_thirdmoment)) then
      do_third=.true.
      third_type=i_thirdmoment
    end if

    do_update=.true.
    if (present(l_passive)) do_update= .not. l_passive

    allocate(tend_temp(lbound(tend,1):ubound(tend,1), lbound(tend,2):ubound(tend,2)))
    tend_temp=ZERO_REAL_WP

    nproc=size(iprocs)
    nz=size(procs,1)

    do i=1, nproc
      if (iprocs(i)%on) then
        iproc=iprocs(i)%id
        do k=1, nz
          do iq=1, ntotalq
            tend_temp(k, iq)=tend_temp(k, iq)+procs(k,iproc)%source(iq)*dst
#if DEF_MODEL==MODEL_LEM_DIAG
            if (diaglevel > 4 .and. procs(k,iproc)%source(iq) /= ZERO_REAL_WP) then 
              ! Slows model significantly so only really used for development. This only works without substepping.
              name=trim(iprocs(i)%name)//'_'//trim(adjustl(names(iq)))
              idg=req_dgproc(trim(name))
              if (idg > 0 .and. j_here <= jenddg .and. j_here >= jstartdg) then
                dgprocs(j_here-jstartdg+1, k, i_here, idg)=procs(k,iproc)%source(iq)
              end if
            end if
#elif DEF_MODEL==MODEL_KiD
            if (diaglevel > 4 .and. procs(k,iproc)%source(iq) /= ZERO_REAL_WP) then
              ! Slows model significantly so only really used for development. This only works without substepping.
              name=trim(iprocs(i)%name)//'_'//trim(adjustl(names(iq)))
              if (nx == 1) then
                call save_dg(k, procs(k,iproc)%source(iq), name, i_dgtime)
              else
                call save_dg(k, i_here, procs(k,iproc)%source(iq), name, i_dgtime)
              end if
            end if
#endif
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
              dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
              dm2=tend_temp(k, params%i_2m)
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
                    call throw_mphys_error(1, 'sum_procs', 'i_thirdmoment incorrectly set')
                  end select
                  tend_temp(k, params%i_3m)=dm3
                end if
              end if
            end if
            ! Snow
            params=snow_params
            if (params%l_3m) then
              m1=qfields(k, params%i_1m)*rho(k)/params%c_x
              m2=qfields(k, params%i_2m)
              m3=qfields(k, params%i_3m)
              dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
              dm2=tend_temp(k, params%i_2m)
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
                    call throw_mphys_error(1, 'sum_procs', 'i_thirdmoment incorrectly set')
                  end select
                  tend_temp(k, params%i_3m)=dm3
                end if
              end if
            end if
            ! Graupel
            params=graupel_params
            if (params%l_3m) then
              m1=qfields(k, params%i_1m)*rho(k)/params%c_x
              m2=qfields(k, params%i_2m)
              m3=qfields(k, params%i_3m)
              dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
              dm2=tend_temp(k, params%i_2m)
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
                    call throw_mphys_error(1, 'sum_procs', 'i_thirdmoment incorrectly set')
                  end select
                  tend_temp(k, params%i_3m)=dm3
                end if
              end if
            end if
#if DEF_MODEL==MODEL_KiD
            params=rain_params
            if (params%l_3m) then
              write(name, '(i1.1)')
              name =  trim(name)//'_'//trim(adjustl(names(params%i_3m)))
              if (nx == 1) then
                call save_dg(k, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
              else
                call save_dg(k, i_here, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
              end if
            end if
            params=snow_params
            if (params%l_3m) then
              write(name, '(i1.1)')
              name =  trim(name)//'_'//trim(adjustl(names(params%i_3m)))
              if (nx == 1) then
                call save_dg(k, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
              else
                call save_dg(k, i_here, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
              end if
            end if
            params=graupel_params
            if (params%l_3m) then
              write(name, '(i1.1)')
              name =  trim(name)//'_'//trim(adjustl(names(params%i_3m)))
              if (nx == 1) then
                call save_dg(k, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
              else
                call save_dg(k, i_here, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
              end if
            end if
#endif
          end if

        end do
      end if
    end do
!!$
!!$IF (PRESENT(qfields)) THEN
!!$       ! Limit condensate values !!!! PRAGMATIC HACK !!!!
!!$  DO k=1,nz
!!$    IF (qfields(k, rain_params%i_1m) > 6.0e-3) THEN
!!$      IF (tend_temp(k,i_qr) > 0.0) THEN
!!$        tend_temp(k,i_qr) = 0.0
!!$        tend_temp(k,i_nr) = 0.0
!!$      END IF
!!$    END IF
!!$    IF (.NOT. l_warm) THEN
!!$      IF (qfields(k, snow_params%i_1m) > 6.0e-3) THEN
!!$        IF (tend_temp(k,i_qs) > 0.0) THEN
!!$          tend_temp(k,i_qs) = 0.0
!!$          tend_temp(k,i_ns) = 0.0
!!$        END IF
!!$      END IF
!!$      IF (qfields(k, graupel_params%i_1m) > 6.0e-3) THEN
!!$        IF (tend_temp(k,i_qg) > 0.0) THEN
!!$          tend_temp(k,i_qg) = 0.0
!!$          tend_temp(k,i_ng) = 0.0
!!$        END IF
!!$      END IF
!!$    END IF
!!$!       if (qfields(k, rain_params%i_1m) > 10.e-3 .and. &
!!$!          qfields(k, cloud_params%i_1m) > 1.5e-3)then
!!$!          if (tend_temp(k,i_ql) > 0.0)then
!!$!             tend_temp(k,i_ql) = 0.0
!!$!             tend_temp(k,i_nl) = 0.0
!!$!          end if
!!$!       end if
!!$  END DO
!!$END IF

    ! Calculate the thermal exchange values
    ! (this overwrites anything that was already stored in the theta tendency)
    if (do_thermal) then
      do k=1, nz
        tend_temp(k,i_th)=(tend_temp(k,i_ql)+tend_temp(k,i_qr))*Lv/cp * rexner(k)
        if (.not. l_warm) then
          tend_temp(k,i_th)=tend_temp(k,i_th)+(tend_temp(k,i_qi)+tend_temp(k,i_qs)+tend_temp(k,i_qg))*Ls/cp *rexner(k)
        end if
      end do
    end if

    ! Add on tendencies to those already passed in.
    if (do_update) tend=tend+tend_temp
    deallocate(tend_temp)
  end subroutine sum_procs
end module sum_process
