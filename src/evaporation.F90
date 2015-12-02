MODULE evaporation

USE variable_precision, ONLY: wp, iwp
USE special, ONLY: pi, GammaFunc

USE passive_fields, ONLY: rho, qws, TdegK, exner

USE mphys_switches, ONLY:    &
     i_qr, i_nr, i_m3r, i_qv, i_ql &
     , hydro_complexity, l_2mr, l_3mr &
     , i_am2, i_an2, i_am3, i_an3, i_am4, i_am5        &
     , l_aevp, l_process, active_cloud, active_rain, isol, aero_index  &
     , l_separate_rain, i_am6, i_an6, i_am9, l_warm &
     , i_qi, i_qs, i_qg, i_an11, i_an12, l_passivenumbers, l_passivenumbers_ice&
     , l_inhom_revp
USE mphys_constants, ONLY: rhow, visair, rho0, Lv   &
     ,ka, Dv, Rv, cp
USE mphys_parameters, ONLY: c_r, vent_1, vent_2, a_r, b_r, f_r   &
     ,p1, p2, p3, hydro_params, rain_params
USE process_routines, ONLY: process_rate, i_prevp, i_arevp
USE thresholds, ONLY: ql_small, qr_small, ss_small, qr_tidy

USE m3_incs, ONLY: m3_inc_type2
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active

USE ventfac, ONLY: ventilation

USE which_mode_to_use

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime
#elif DEF_MODEL==MODEL_UM
USE diaghelp_um, ONLY: i_here, j_here
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE revp(dt, k, qfields, aerofields, aerophys, aerochem, aeroact, dustliq, procs, aerosol_procs, l_sigevap)

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN) :: qfields(:,:), aerofields(:,:)
TYPE(aerosol_phys), INTENT(IN) :: aerophys(:)
TYPE(aerosol_chem), INTENT(IN) :: aerochem(:)
TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustliq(:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)
LOGICAL, INTENT(OUT) :: l_sigevap ! Determines if there is significant evaporation

REAL(wp) :: dmass, dnumber, dnumber_a, dnumber_d
REAL(wp) :: m1, m2, m3, dm1, dm3

REAL(wp) :: n0, lam, mu
REAL(wp) :: V_r, AB

REAL(wp) :: rain_mass
REAL(wp) :: rain_number
REAL(wp) :: rain_m3
REAL(wp) :: qv
TYPE(process_rate), POINTER :: this_proc
TYPE(process_rate), POINTER :: aero_proc

LOGICAL :: l_rain_test ! conditional test on rain

REAL(wp) :: dmac, dmac1, dmac2, dnac1, dnac2, dmacd

l_sigevap=.FALSE.

qv  = qfields(k, i_qv)
rain_mass   = qfields(k, i_qr)
IF (l_2mr)rain_number = qfields(k, i_nr)
IF (l_3mr)rain_m3 = qfields(k, i_m3r)

IF (qv/qws(k) < 1.0-ss_small .AND. qfields(k, i_ql) == 0.0 .AND. rain_mass > qr_tidy) THEN

  this_proc => procs(k, i_prevp%id)

  m1=rain_mass/c_r
  IF (l_2mr)m2=rain_number
  IF (l_3mr)m3=rain_m3

  l_rain_test=.FALSE.
  IF (l_2mr) l_rain_test=rain_number>0
  IF (rain_mass > qr_small .AND. (.NOT. l_2mr .OR. l_rain_test)) THEN

    n0 = dist_n0(k,rain_params%id)
    mu = dist_mu(k,rain_params%id)
    lam = dist_lambda(k,rain_params%id)

    CALL ventilation(k, V_r, n0, lam, mu, rain_params)

    AB = 1.0/(Lv**2/(Rv*ka)*rho(k)*TdegK(k)**(-2)    &
           + 1.0/(Dv*qws(k)))

    dmass = (1.0-qv/qws(k))*V_r*AB

  ELSE

    dmass=rain_mass/dt

  END IF

  dm1 = dmass/c_r
  IF (l_2mr) THEN
    dnumber=0.0
    IF (l_inhom_revp)dnumber = dm1*m2/m1
  END IF
  IF (l_3mr)dm3 = dm1*m3/m1

  IF (l_2mr) THEN
    IF (dnumber*dt > rain_number .OR. dmass*dt >= rain_mass-qr_tidy) THEN
      dmass=rain_mass/dt
      dnumber=rain_number/dt
      IF (l_3mr)dm3=m3/dt
    END IF
  END IF

  this_proc%source(i_qr) = -dmass
  this_proc%source(i_qv) = dmass

  IF (dmass*dt/rain_mass > .8)l_sigevap=.TRUE.

  IF (l_2mr) THEN
    this_proc%source(i_nr) = -dnumber
  END IF
  IF (l_3mr) THEN
    this_proc%source(i_m3r) = -dm3
  END IF

  NULLIFY(this_proc)

      !============================
      ! aerosol processing
      !============================
  IF (l_process .AND. ABS(dnumber) >0) THEN
    aero_proc => aerosol_procs(k, i_arevp%id)

    dmac=dnumber*aeroact(k)%nratio2*aeroact(k)%mact2_mean
    IF (l_separate_rain) THEN
      aero_proc%source(i_am5) = -dmac
    ELSE
      aero_proc%source(i_am4) = -dmac
    END IF

        ! Return aerosol
    IF (aero_index%i_accum >0 .AND. aero_index%i_coarse >0) THEN
          ! Coarse and accumulation mode being used. Which one to return to?
      CALL which_mode(dmac, dnumber*aeroact(k)%nratio2,                        &
           aerophys(k)%rd(aero_index%i_accum), aerophys(k)%rd(aero_index%i_coarse), &
           aerochem(k)%density(aero_index%i_accum),     &
           dmac1, dmac2, dnac1, dnac2)
      aero_proc%source(i_am2) = dmac1
      aero_proc%source(i_an2) = dnac1
      aero_proc%source(i_am3) = dmac2
      aero_proc%source(i_an3) = dnac2
    ELSE
      IF (aero_index%i_accum >0) THEN
        aero_proc%source(i_am2) = dmac
        aero_proc%source(i_an2) = dnumber
      END IF
      IF (aero_index%i_coarse >0) THEN
        aero_proc%source(i_am3) = dmac
        aero_proc%source(i_an3) = dnumber
      END IF
    END IF

    dmacd=dnumber*dustliq(k)%nratio2*dustliq(k)%mact2_mean
    IF (.NOT. l_warm .AND. dmacd /=0.0) THEN
      aero_proc%source(i_am9) = -dmacd
      aero_proc%source(i_am6) = dmacd
      aero_proc%source(i_an6) = dnumber*dustliq(k)%nratio2
    END IF

    IF (l_passivenumbers) THEN
      dnumber_a = -dnumber*aeroact(k)%nratio2
      aero_proc%source(i_an11) = dnumber_a
    END IF

    IF (l_passivenumbers_ice) THEN
      dnumber_d = -dnumber*dustliq(k)%nratio2
      aero_proc%source(i_an12) = dnumber_d
    END IF

    NULLIFY(aero_proc)
  END IF
END IF

END SUBROUTINE revp


END MODULE evaporation
