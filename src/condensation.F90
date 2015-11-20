MODULE condensation

USE mphys_die, ONLY: throw_mphys_error
USE variable_precision, ONLY: wp
USE passive_fields, ONLY: rho, pressure, w, exner, qws
USE mphys_switches, ONLY: i_qv, i_ql, i_nl, i_th,      &
     i_qr, i_qi, i_qs, i_qg, hydro_complexity, l_warm, &
     i_am4, i_am1, i_an1, i_am2, i_an2, i_am3, i_an3,  &
     i_am6, i_an6, i_am9, i_an11, i_an12,  &
     cloud_params, l_process, l_passivenumbers,l_passivenumbers_ice, aero_index, &
     active_cloud, active_number

USE process_routines, ONLY: process_rate, i_cond, i_aact
USE mphys_constants, ONLY: Lv, cp
USE qsat_funs, ONLY: qsaturation, dqwsatdt
USE thresholds, ONLY: ql_small, w_small, nl_small, ss_small,   &
     thresh_tidy
USE activation, ONLY: activate
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
USE special, ONLY: pi

USE which_mode_to_use

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, k_here
USE runtime, ONLY: time
#elif  DEF_MODEL==MODEL_LEM
USE diaghelp_lem, ONLY: k_here, i_here, j_here
USE com_params, ONLY: time
#elif  DEF_MODEL==MODEL_UM
USE casim_switches, ONLY: l_cfrac_casim_diag_scheme
USE UM_ParCore, ONLY: mype
USE diaghelp_um, ONLY: i_here, j_here
USE cloud_frac_scheme, ONLY: cloud_frac_casim_mphys
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here
#endif


IMPLICIT NONE

LOGICAL :: l_notransfer=.TRUE.  ! don't transfer aerosol from one mode to another.

CONTAINS

SUBROUTINE condevp(dt, k, qfields, aerofields, procs, aerophys, aerochem,   &
     aeroact, dustphys, dustchem, dustliq, aerosol_procs, rhcrit_lev)

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN), TARGET :: qfields(:,:), aerofields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! aerosol fields
TYPE(aerosol_phys), INTENT(IN) :: aerophys(:)
TYPE(aerosol_chem), INTENT(IN) :: aerochem(:)
TYPE(aerosol_active), INTENT(IN) :: aeroact(:)

TYPE(aerosol_phys), INTENT(IN) :: dustphys(:)
TYPE(aerosol_chem), INTENT(IN) :: dustchem(:)
TYPE(aerosol_active), INTENT(IN) :: dustliq(:)

    ! optional aerosol fields to be processed
TYPE(process_rate), INTENT(INOUT), OPTIONAL, TARGET :: aerosol_procs(:,:)

REAL(wp) :: dmass, dnumber, dmac, dmad, dnumber_a, dnumber_d
REAL(wp) :: dmac1, dmac2, dnac1, dnac2
REAL(wp), ALLOCATABLE :: dnccn_all(:),dmac_all(:)

REAL(wp) :: th
REAL(wp) :: qv
REAL(wp) :: cloud_mass
REAL(wp) :: cloud_number
TYPE(process_rate), POINTER :: this_proc
TYPE(process_rate), POINTER :: aero_proc

REAL(wp) :: qs, dqsdt, qsatfac

REAL(wp) :: tau   ! timescale for adjustment of condensate
REAL(wp) :: w_act ! vertical velocity to use for activation

REAL(wp) :: rmean ! mean aerosol diameter

! local variables for diagnostics cloud scheme (if needed)
REAL(wp) :: cloud_mass_new, abs_liquid_t

REAL(wp), INTENT(IN) :: rhcrit_lev

LOGICAL :: l_docloud  ! do we want to do the calculation of cond/evap

tau=dt ! adjust instantaneously

    ! Initializations
dmac = 0.0
dmad = 0.0
dnumber = 0.0
dnumber_a = 0.0
dnumber_d = 0.0

ALLOCATE(dnccn_all(aero_index%nccn))
ALLOCATE(dmac_all(aero_index%nccn))
dnccn_all= 0.0
dmac_all = 0.0

    ! Set pointers for convenience
cloud_mass = qfields(k, i_ql)
IF (cloud_params%l_2m) cloud_number = qfields(k, i_nl)

th = qfields(k, i_th)
#if  DEF_MODEL==MODEL_UM
IF (l_cfrac_casim_diag_scheme) THEN
  qv = qfields(k, i_qv)+cloud_mass
ELSE
  qv = qfields(k, i_qv)
END IF
IF (l_cfrac_casim_diag_scheme) THEN
       ! work out saturation vapour pressure/mixing ratio based on
       ! liquid water temperature
  abs_liquid_T = (th*exner(k)) - ((lv * cloud_mass ) / cp )
  qs = qsaturation(abs_liquid_T, pressure(k)/100.0)
  PRINT *, 't_l, cloud, th', abs_liquid_T, cloud_mass, qs, k
ELSE
  qs = qsaturation(th*exner(k), pressure(k)/100.0)
END IF
#else
qv = qfields(k, i_qv)
qs = qsaturation(th*exner(k), pressure(k)/100.0)
#endif

l_docloud=.TRUE.
IF (qs==0.0) THEN
  l_docloud=.FALSE.
END IF

IF ((qv/qs > 1.0 - ss_small .OR. cloud_mass > 0.0) .AND. l_docloud) THEN

#if  DEF_MODEL==MODEL_UM
  IF (l_cfrac_casim_diag_scheme) THEN

          !Call Smith scheme before setting up microphysics vars, to work out
          ! cloud fraction, which is used to derive in-cloud mass and number
          !
          !IMPORTANT - qv is total water at this stage!
    CALL cloud_frac_casim_mphys(k, pressure(k), abs_liquid_T, rhcrit_lev,         &
         qs, qv, cloud_mass, qfields(k,i_qr), cloud_mass_new )

    dmass = MAX (-cloud_mass, (cloud_mass_new - cloud_mass))/dt

          ! qv is total water, convert back to vapour to be consistent with updtaing
          ! and rest of code
          !qv = qfields(k, i_qv)

  ELSE
    dqsdt = dqwsatdt(qs, th*exner(k))
    qsatfac = 1.0/(1.0 + Lv/Cp*dqsdt)
    dmass = MAX ( -cloud_mass, (qv - qs) * qsatfac )/dt
  END IF
#else
  dqsdt = dqwsatdt(qs, th*exner(k))
  qsatfac = 1.0/(1.0 + Lv/Cp*dqsdt)
  dmass = MAX ( -cloud_mass, (qv - qs) * qsatfac )/dt
#endif


#if DEF_MODEL==MODEL_KiD
  CALL save_dg(k_here, i_here, qv, 'qv_in_cond', i_dgtime)
  CALL save_dg(k_here, i_here, qs, 'qs_in_cond', i_dgtime)
  CALL save_dg(k_here, i_here, qv/(qs+1e-20), 'rh_in_cond', i_dgtime)
  CALL save_dg(k_here, i_here, dmass, 'dmass_in_cond', i_dgtime)
#endif

  IF (dmass > 0.0_wp) THEN ! condensation

    IF (dmass*dt + cloud_mass > ql_small) THEN ! is it worth bothering with?

      IF (cloud_params%l_2m) THEN

            ! If significant cloud formed then assume minimum velocity of 0.1m/s
        w_act = MAX(w(k), 0.1_wp)

        CALL activate(tau, cloud_mass, cloud_number, w_act, rho(k), dnumber, dmac, &
             th*exner(k), pressure(k), aerophys(k), aerochem(k), aeroact(k),   &
             dustphys(k), dustchem(k), dustliq(k),           &
             dnccn_all, dmac_all, dnumber_d, dmad)

        dnumber_a = dnumber
      END IF
    ELSE
      dmass = 0.0 ! not worth doing anything
    END IF
  ELSE  ! evaporation
    IF (cloud_mass > thresh_tidy(i_ql)) THEN ! anything significant to remove or just noise?
      IF (dmass*dt + cloud_mass < ql_small) THEN  ! Remove all cloud

            ! Remove small quantities.
        dmass=-cloud_mass/dt
            ! liberate all number and aerosol
        IF (cloud_params%l_2m) THEN
          dnumber=-cloud_number/dt

              !============================
              ! aerosol processing
              !============================
          IF (l_process) THEN
            dmac = -aeroact(k)%mact1/dt
            dmad = -dustliq(k)%mact1/dt

            IF (l_passivenumbers) THEN
              dnumber_a = -aeroact(k)%nact1/dt
            ELSE
              dnumber_a = dnumber
            END IF
            IF (l_passivenumbers_ice) THEN
              dnumber_d = -dustliq(k)%nact1/dt
            ELSE
              dnumber_d = dnumber
            END IF



            IF (aero_index%i_accum >0 .AND. aero_index%i_coarse >0) THEN
                  ! We have both accumulation and coarse modes
              IF (dnumber_a*dmac<= 0) THEN
                    !print*, 'COND - sort out tidy routines...',k,i_here,aerofields(k,i_am4),  dmac, dnumber_a &
                    !     , aerofields(k,i_an11), aeroact(k)%nact1, aeroact(k)%nact, aeroact(k)%nact2
                dnumber_a = dmac/1.0e-18/dt
              END IF
              CALL which_mode(dmac, dnumber_a,                                 &
                   aerophys(k)%rd(aero_index%i_accum), aerophys(k)%rd(aero_index%i_coarse), &
                   aerochem(k)%density(aero_index%i_accum),     &
                   dmac1, dmac2, dnac1, dnac2)

              dmac_all(aero_index%i_accum) = dmac1  ! put it back into accumulation mode
              dnccn_all(aero_index%i_accum)= dnac1
              dmac_all(aero_index%i_coarse) = dmac2  ! put it back into coarse mode
              dnccn_all(aero_index%i_coarse)= dnac2
            END IF
          END IF
        END IF
      ELSE ! Still some cloud will be left behind
        dnumber  = 0.0 ! we assume no change in number during evap
        dnccn_all= 0.0 ! we assume no change in number during evap
        dmac     = 0.0 ! No aerosol processing required
        dmac_all = 0.0 ! No aerosol processing required
        dnumber_a= 0.0 ! No aerosol processing required
        dnumber_d= 0.0 ! No aerosol processing required
      END IF
    ELSE  ! Nothing significant here to remove - the tidying routines will deal with this
      dmass    = 0.0 ! no need to do anything since this is now just numerical noise
      dnumber  = 0.0 ! we assume no change in number during evap
      dnccn_all= 0.0 ! we assume no change in number during evap
      dmac     = 0.0 ! No aerosol processing required
      dmac_all = 0.0 ! No aerosol processing required
      dnumber_a= 0.0 ! No aerosol processing required
      dnumber_d= 0.0 ! No aerosol processing required
    END IF
  END IF

  IF (dmass /= 0.0_wp) THEN
    this_proc => procs(k, i_cond%id)

    this_proc%source(i_qv) = -dmass
    this_proc%source(i_ql) = dmass

    IF (cloud_params%l_2m) THEN
      this_proc%source(i_nl) = dnumber
#if DEF_MODEL==MODEL_KiD
      CALL save_dg(k_here, dnumber, 'dnumber_in_cond', i_dgtime )
      IF (aero_index%i_aitken > 0)CALL save_dg(k_here,  dnccn_all(aero_index%i_aitken), 'dnccn_aitken', i_dgtime)
      IF (aero_index%i_accum > 0)CALL save_dg(k_here,  dnccn_all(aero_index%i_accum), 'dnccn_accum', i_dgtime)
      IF (aero_index%i_coarse > 0)CALL save_dg(k_here,  dnccn_all(aero_index%i_coarse), 'dnccn_coarse', i_dgtime)
#endif
    END IF

    NULLIFY(this_proc)

        !============================
        ! aerosol processing
        !============================
    IF (l_process) THEN
      aero_proc => aerosol_procs(k, i_aact%id)

      aero_proc%source(i_am4)   = dmac
      IF (l_passivenumbers)aero_proc%source(i_an11) = dnumber_a
      IF (l_passivenumbers_ice)aero_proc%source(i_an12) = dnumber_d

      IF (aero_index%i_aitken > 0) THEN
        aero_proc%source(i_am1) =     &
               -dmac_all(aero_index%i_aitken)
        aero_proc%source(i_an1) =     &
               -dnccn_all(aero_index%i_aitken)
      END IF
      IF (aero_index%i_accum > 0) THEN
        aero_proc%source(i_am2) =     &
               -dmac_all(aero_index%i_accum)
        aero_proc%source(i_an2) =     &
               -dnccn_all(aero_index%i_accum)

      END IF
      IF (aero_index%i_coarse > 0) THEN
        aero_proc%source(i_am3) =     &
               -dmac_all(aero_index%i_coarse)
        aero_proc%source(i_an3) =     &
               -dnccn_all(aero_index%i_coarse)
      END IF

      IF (.NOT. l_warm .AND. dmad /=0.0) THEN
            ! We may have some dust in the liquid...
        aero_proc%source(i_am9) = dmad
        aero_proc%source(i_am6) = -dmad ! < USING COARSE
        aero_proc%source(i_an6) = -dnumber_d ! < USING COARSE
      END IF

      NULLIFY(aero_proc)
    END IF
  END IF

END IF

DEALLOCATE(dnccn_all)
DEALLOCATE(dmac_all)

END SUBROUTINE condevp

END MODULE condensation
