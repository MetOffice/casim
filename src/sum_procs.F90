MODULE sum_process

USE variable_precision, ONLY: wp
USE mphys_die, ONLY: throw_mphys_error
USE type_process, ONLY: process_name, process_rate
USE mphys_switches, ONLY: i_th, i_qv, i_ql, i_qr, i_qs, i_qi, i_qg, l_warm, i_m3r, ntotalq, ntotala, &
       i_an11, i_am4, i_nl, i_am9, aero_index, i_nr, i_ng, i_ns, i_am7, i_am8, nsubsteps
USE passive_fields, ONLY: rexner, rho
USE mphys_constants, ONLY: cp, Lv, Ls
USE mphys_parameters, ONLY: hydro_params, snow_params, rain_params, graupel_params, parent_dt, ZERO_REAL_WP
USE m3_incs, ONLY: m3_inc_type2

#if DEF_MODEL==MODEL_KiD
USE parameters, ONLY: diaglevel, nx
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, j_here,   &
       n_sub, n_subsed
USE runtime, ONLY: time
#elif DEF_MODEL==MODEL_LEM_DIAG
USE diaghelp_lem, ONLY: dgfields, i_here, j_here, koff   &
     , diaglevel, dgindex, dgindex2, &
       n_sub, n_subsed
USE com_params, ONLY: time
USE extra_dgs
#elif DEF_MODEL==MODEL_UM
USE mphys_casim_diagnostics, ONLY: ProcessRates, nProcessDiags, ProcessKeys, ProcessQs, &
       PhaseChanges
USE diaghelp_um, ONLY: i_here, j_here, l_debug_um, debug_i, debug_j, debug_pe, &
       n_sub, n_subsed
USE UM_ParCore, ONLY: mype
USE timestep_mod, ONLY: timestep_number
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here
#endif

IMPLICIT NONE


#if DEF_MODEL==MODEL_UM
INTEGER :: idg_proc

  ! Temporary for debugging
INTEGER :: max_prints(10)=0
#endif

CONTAINS

SUBROUTINE sum_aprocs(dst, procs, tend, iprocs, names, afields)

REAL(wp), INTENT(IN) :: dst  ! step length (s)
TYPE(process_rate), INTENT(IN) :: procs(:,:)
TYPE(process_name), INTENT(IN) :: iprocs(:)
REAL(wp), INTENT(INOUT) :: tend(:,:)
CHARACTER(20) :: names(:)

REAL(wp), INTENT(IN), OPTIONAL :: afields(:,:) ! Required for debugging

REAL(wp), ALLOCATABLE :: tend_temp(:,:) ! Temporary storage for accumulated tendendies

INTEGER :: k, iq, iproc, i, nz, idgproc
INTEGER :: nproc, index, idg

#if DEF_MODEL==MODEL_UM
INTEGER :: diaglevel
#endif

CHARACTER(2) :: char
CHARACTER(100) :: name

#if DEF_MODEL==MODEL_UM
diaglevel = 0
#endif

ALLOCATE(tend_temp(LBOUND(tend,1):UBOUND(tend,1), LBOUND(tend,2):UBOUND(tend,2)))
tend_temp=ZERO_REAL_WP

nproc=SIZE(iprocs)
nz=SIZE(procs,1)

DO i=1, nproc
  IF (iprocs(i)%on) THEN
    iproc=iprocs(i)%id
    DO k=1,nz
      IF (.NOT. ALL(procs(k,iproc)%source(:)==ZERO_REAL_WP)) THEN
        DO iq=1, ntotala

          tend_temp(k, iq) = tend_temp(k, iq) + procs(k,iproc)%source(iq)*dst

#if DEF_MODEL==MODEL_LEM_DIAG
          IF (diaglevel > 4 .AND. procs(k,iproc)%source(iq) /= ZERO_REAL_WP) THEN ! Slows model significantly so
                ! Only really used for development
                ! This only works without substepping.
            name =  TRIM(iprocs(i)%name)//'_'//TRIM(ADJUSTL(names(iq)))

            idg = req_dgproc(TRIM(name))
            IF (idg > 0 .AND. j_here<=jenddg .AND. j_here>=jstartdg) THEN
              dgprocs(j_here-jstartdg+1, k, i_here, idg) = procs(k,iproc)%source(iq)
            END IF

          END IF
#elif DEF_MODEL==MODEL_KiD
          IF (diaglevel > 4 .AND. procs(k,iproc)%source(iq) /= ZERO_REAL_WP) THEN ! Slows model significantly so
            name =  TRIM(iprocs(i)%name)//'_'//TRIM(ADJUSTL(names(iq)))
            IF (nx==1) THEN
              CALL save_dg(k, procs(k, iproc)%source(iq), name, i_dgtime)
            ELSE
              CALL save_dg(k, i_here, procs(k, iproc)%source(iq), name, i_dgtime)
            END IF
          END IF
#endif


        END DO
      END IF

    END DO
  END IF
END DO

    ! Add on tendencies to those already passed in.
tend = tend + tend_temp

DEALLOCATE(tend_temp)

END SUBROUTINE sum_aprocs

SUBROUTINE sum_procs(dst, procs, tend, iprocs, names,         &
     l_thermalexchange, i_thirdmoment, qfields, l_passive )

REAL(wp), INTENT(IN) :: dst  ! step length (s)
TYPE(process_rate), INTENT(IN) :: procs(:,:)
TYPE(process_name), INTENT(IN) :: iprocs(:)
REAL(wp), INTENT(INOUT) :: tend(:,:)
LOGICAL, INTENT(IN), OPTIONAL :: l_thermalexchange  ! Calculate the thermal exchange terms
INTEGER, INTENT(IN), OPTIONAL :: i_thirdmoment  ! Calculate the tendency of the third moment
REAL(wp), INTENT(IN), OPTIONAL :: qfields(:,:) ! Required for debugging or with i_thirdmoment
LOGICAL, INTENT(IN), OPTIONAL :: l_passive ! if true don't apply final tendency
                                                   ! (for testing diagnostics with pure sedimentation)
CHARACTER(10) :: names(:)

REAL(wp), ALLOCATABLE :: tend_temp(:,:) ! Temporary storage for accumulated tendendies

TYPE(hydro_params) :: params
REAL(wp) :: dm1,dm2,dm3,m1,m2,m3

INTEGER :: k, iq, iproc, isource, nsource, i, nz, idgproc
INTEGER :: nproc, index, idg


#if DEF_MODEL==MODEL_UM
INTEGER :: diaglevel
#endif

CHARACTER(2) :: char
CHARACTER(100) :: name

LOGICAL :: do_thermal
LOGICAL :: do_third  ! Currently not plumbed in, so explicitly calculated for each process
LOGICAL :: do_update ! update the tendency

INTEGER :: third_type

REAL :: qh ! total hydrometeor

#if DEF_MODEL==MODEL_UM
diaglevel = 0
#endif

do_thermal=.FALSE.
IF (PRESENT(l_thermalexchange))do_thermal=l_thermalexchange

do_third=.FALSE.
IF (PRESENT(i_thirdmoment)) THEN
  do_third=.TRUE.
  third_type=i_thirdmoment
END IF

do_update=.TRUE.
IF (PRESENT(l_passive))do_update= .NOT. l_passive

ALLOCATE(tend_temp(LBOUND(tend,1):UBOUND(tend,1), LBOUND(tend,2):UBOUND(tend,2)))
tend_temp=ZERO_REAL_WP

nproc=SIZE(iprocs)
nz=SIZE(procs,1)


DO i=1, nproc
  IF (iprocs(i)%on) THEN
    iproc=iprocs(i)%id
    DO k=1,nz

      IF (.NOT. ALL(procs(k,iproc)%source(:)==ZERO_REAL_WP)) THEN

        DO iq=1, ntotalq
          tend_temp(k, iq) = tend_temp(k, iq) + procs(k,iproc)%source(iq)*dst
#if DEF_MODEL==MODEL_LEM_DIAG
          IF (diaglevel > 4 .AND. procs(k,iproc)%source(iq) /= ZERO_REAL_WP) THEN ! Slows model significantly so
                ! Only really used for development
                ! This only works without substepping.
            name =  TRIM(iprocs(i)%name)//'_'//TRIM(ADJUSTL(names(iq)))

            idg = req_dgproc(TRIM(name))
            IF (idg > 0 .AND. j_here<=jenddg .AND. j_here>=jstartdg) THEN
              dgprocs(j_here-jstartdg+1, k, i_here, idg) = procs(k,iproc)%source(iq)
            END IF
          END IF
#elif DEF_MODEL==MODEL_KiD
          IF (diaglevel > 4 .AND. procs(k,iproc)%source(iq) /= ZERO_REAL_WP) THEN ! Slows model significantly so
            name =  TRIM(iprocs(i)%name)//'_'//TRIM(ADJUSTL(names(iq)))
            IF (nx==1) THEN
              CALL save_dg(k, procs(k,iproc)%source(iq),    &
                   name, i_dgtime)
            ELSE
              CALL save_dg(k, i_here, procs(k,iproc)%source(iq),     &
                   name, i_dgtime)
            END IF
          END IF
#endif

        END DO

      END IF

      IF (do_third) THEN
            ! calculate increment to third moment based on collected increments from
            ! q and n NB This overwrites any previously calculated values
            ! Rain
        params=rain_params
        IF (params%l_3m) THEN
          m1=qfields(k, params%i_1m)*rho(k)/params%c_x
          m2=qfields(k, params%i_2m)
          m3=qfields(k, params%i_3m)
          dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
          dm2=tend_temp(k, params%i_2m)
          IF (dm1 < -.99*m1 .OR. dm2 < -.99*m2) THEN
            dm1 = -m1
            dm2 = -m2
            dm3 = -m3
          ELSE
            IF (m3>0.0 .AND. m1 >0.0 .AND. m2 > 0.0 .AND. (ABS(dm1) >0.0 .OR. ABS(dm2) > 0.0)) THEN
              SELECT CASE (i_thirdmoment)
                CASE (2)
                  CALL m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
                CASE default
                  CALL throw_mphys_error(1, 'sum_procs', 'i_thirdmoment incorrectly set')
              END SELECT
              tend_temp(k, params%i_3m) = dm3
            END IF
          END IF
        END IF
            ! Snow
        params=snow_params
        IF (params%l_3m) THEN
          m1=qfields(k, params%i_1m)*rho(k)/params%c_x
          m2=qfields(k, params%i_2m)
          m3=qfields(k, params%i_3m)
          dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
          dm2=tend_temp(k, params%i_2m)
          IF (dm1 < -.99*m1 .OR. dm2 < -.99*m2) THEN
            dm1 = -m1
            dm2 = -m2
            dm3 = -m3
          ELSE
            IF (m3>0.0 .AND. m1 >0.0 .AND. m2 > 0.0 .AND. (ABS(dm1) >0.0 .OR. ABS(dm2) > 0.0)) THEN
              SELECT CASE (i_thirdmoment)
                CASE (2)
                  CALL m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
                CASE default
                  CALL throw_mphys_error(1, 'sum_procs', 'i_thirdmoment incorrectly set')
              END SELECT
              tend_temp(k, params%i_3m) = dm3
            END IF
          END IF
        END IF
            ! Graupel
        params=graupel_params
        IF (params%l_3m) THEN
          m1=qfields(k, params%i_1m)*rho(k)/params%c_x
          m2=qfields(k, params%i_2m)
          m3=qfields(k, params%i_3m)
          dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
          dm2=tend_temp(k, params%i_2m)
          IF (dm1 < -.99*m1 .OR. dm2 < -.99*m2) THEN
            dm1 = -m1
            dm2 = -m2
            dm3 = -m3
          ELSE
            IF (m3>0.0 .AND. m1 >0.0 .AND. m2 > 0.0 .AND. (ABS(dm1) >0.0 .OR. ABS(dm2) > 0.0)) THEN
              SELECT CASE (i_thirdmoment)
                CASE (2)
                  CALL m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
                CASE default
                  CALL throw_mphys_error(1, 'sum_procs', 'i_thirdmoment incorrectly set')
              END SELECT
              tend_temp(k, params%i_3m) = dm3
            END IF
          END IF
        END IF
#if DEF_MODEL==MODEL_KiD
        params=rain_params
        IF (params%l_3m) THEN
          WRITE(name, '(i1.1)')
          name =  TRIM(name)//'_'//TRIM(ADJUSTL(names(params%i_3m)))
          IF (nx==1) THEN
            CALL save_dg(k, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
          ELSE
            CALL save_dg(k, i_here, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
          END IF
        END IF
        params=snow_params
        IF (params%l_3m) THEN
          WRITE(name, '(i1.1)')
          name =  TRIM(name)//'_'//TRIM(ADJUSTL(names(params%i_3m)))
          IF (nx==1) THEN
            CALL save_dg(k, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
          ELSE
            CALL save_dg(k, i_here, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
          END IF
        END IF
        params=graupel_params
        IF (params%l_3m) THEN
          WRITE(name, '(i1.1)')
          name =  TRIM(name)//'_'//TRIM(ADJUSTL(names(params%i_3m)))
          IF (nx==1) THEN
            CALL save_dg(k, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
          ELSE
            CALL save_dg(k, i_here, tend_temp(k, params%i_3m), 'd3m_type2'//name, i_dgtime)
          END IF
        END IF
#endif
      END IF

    END DO
  END IF
END DO
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
IF (do_thermal) THEN
  DO k=1,nz
    tend_temp(k,i_th) = (tend_temp(k,i_ql) + tend_temp(k,i_qr))*Lv/cp * rexner(k)
    IF (.NOT. l_warm) THEN
      tend_temp(k,i_th) = tend_temp(k,i_th) +      &
             (tend_temp(k,i_qi) + tend_temp(k,i_qs) + tend_temp(k,i_qg))*Ls/cp *rexner(k)
    END IF
  END DO
END IF

    ! Add on tendencies to those already passed in.
IF (do_update)tend = tend + tend_temp

DEALLOCATE(tend_temp)

END SUBROUTINE sum_procs

END MODULE sum_process
