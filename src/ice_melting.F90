MODULE ice_melting

USE variable_precision, ONLY: wp
USE process_routines, ONLY: process_rate, process_name,   &
     i_imlt, i_smlt, i_gmlt, i_sacw, i_sacr, i_gacw, i_gacr, i_gshd, &
     i_dimlt, i_dsmlt, i_dgmlt
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
USE passive_fields, ONLY: TdegC, qws0, rho

USE mphys_parameters, ONLY: ice_params, snow_params, graupel_params,   &
     rain_params, DR_melt, hydro_params, ZERO_REAL_WP

USE mphys_switches, ONLY: i_qv, i_am4, i_am7, i_am8, i_am9, l_process
USE mphys_constants, ONLY: Lv, Lf, Ka, Cwater, Cice, cp, Dv
USE thresholds, ONLY: thresh_tidy

USE m3_incs, ONLY: m3_inc_type2, m3_inc_type3, m3_inc_type4
USE ventfac, ONLY: ventilation
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE melting(dt, k, params, qfields, procs, aeroice, dustact,   &
     aerosol_procs)

    !< Subroutine to calculate rate of melting of ice species
    !<
    !< OPTIMISATION POSSIBILITIES: Shouldn't have to recalculate all 3m quantities
    !<                             If just rescaling mass conversion for dry mode
    !

REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
TYPE(hydro_params), INTENT(IN) :: params
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! aerosol fields
TYPE(aerosol_active), INTENT(IN) :: aeroice(:), dustact(:)

    ! optional aerosol fields to be processed
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

TYPE(process_name) :: iproc, iaproc ! processes selected depending on
    ! which species we're modifying

TYPE(process_name) :: i_acw, i_acr ! accretion processes

REAL(wp) :: qv

REAL(wp) :: dmass, dnumber, dm1, dm2, dm3, dm3_r

REAL(wp) :: number, mass, m1, m2, m3

TYPE(process_rate), POINTER :: this_proc, aero_proc

REAL(wp) :: n0, lam, mu, V_x

REAL(wp) :: acc_correction

LOGICAL :: l_meltall ! do we melt everything?

REAL(wp) :: dmac, dmad

l_meltall=.FALSE.
mass = qfields(k, params%i_1m)

IF (mass > thresh_tidy(params%i_1m) .AND. TdegC(k) > 0.0) THEN

  if (params%l_2m)number = qfields(k, params%i_2m)

  IF (params%id == ice_params%id) THEN ! instantaneous removal
    iaproc=i_dimlt
    iproc=i_imlt

    this_proc => procs(k, iproc%id)

    dmass=mass/dt
    
    if (params%l_2m)then
      dnumber=number/dt
      this_proc%source(ice_params%i_2m) = -dnumber
      IF (rain_params%l_2m) THEN
        this_proc%source(rain_params%i_2m) = dnumber
      END IF
    end if

    this_proc%source(ice_params%i_1m) = -dmass
    this_proc%source(rain_params%i_1m) = dmass

    IF (rain_params%l_3m) THEN
      m1=qfields(k, rain_params%i_1m)/rain_params%c_x
      m2=qfields(k, rain_params%i_2m)
      m3=qfields(k, rain_params%i_3m)

      dm1=dt*dmass/rain_params%c_x
      dm2=dt*dnumber

      IF (m1 > 0.0) THEN
        CALL m3_inc_type2(m1, m2, m3, rain_params%p1,      &
             rain_params%p2, rain_params%p3, dm1, dm2, dm3_r)
      ELSE
        CALL m3_inc_type3(rain_params%p1, rain_params%p2, rain_params%p3,      &
             dm1, dm2, dm3_r, rain_params%fix_mu)
      END IF
      dm3_r=dm3_r/dt
      this_proc%source(rain_params%i_3m) = dm3_r
    END IF



  ELSE

    IF (params%id==snow_params%id) THEN
      i_acw = i_sacw
      i_acr = i_sacr
      iproc = i_smlt
      iaproc=i_dsmlt
    ELSE IF (params%id==graupel_params%id) THEN
      i_acw = i_gacw
      i_acr = i_gacr
      iproc = i_gmlt
      iaproc=i_dgmlt
    END IF

    acc_correction = 0.0

    IF (i_acw%on)acc_correction = procs(k, i_acw%id)%source(params%i_1m)
    IF (i_acr%on)acc_correction = acc_correction +     &
           procs(k, i_acr%id)%source(params%i_1m)
    IF (params%id==graupel_params%id .AND. i_gshd%on) THEN
      acc_correction = acc_correction +     &
             procs(k, i_gshd%id)%source(params%i_1m)
    END IF

    qv = qfields(k, i_qv)

    m1=mass/params%c_x
    IF (params%l_2m) number = qfields(k, params%i_2m)
    IF (params%l_3m) m3 = qfields(k, params%i_3m)

    n0 = dist_n0(k,params%id)
    mu = dist_mu(k,params%id)
    lam = dist_lambda(k,params%id)

    CALL ventilation(k, V_x, n0, lam, mu, params)

    dmass = (1.0/(rho(k)*Lf))*(Ka*TdegC(k) + Lv*Dv*rho(k)*(qv - qws0(k))) * V_x&
           + (Cwater*TdegC(k)/Lf)*acc_correction

    dmass = MAX(dmass, ZERO_REAL_WP) ! ensure positive
        !

    if (dmass == ZERO_REAL_WP)return  ! No need to do anything

    dmass = MIN(dmass, mass/dt) ! ensure we don't remove too much
    IF (dmass*dt > 0.95*mass) THEN ! we're pretty much removing everything
      l_meltall=.TRUE.
      dmass=mass/dt
    END IF

        !--------------------------------------------------
        ! Apply spontaneous rain breakup if drops are large
        !--------------------------------------------------
        !< RAIN BREAKUP TO BE ADDED

    this_proc => procs(k, iproc%id)
    this_proc%source(params%i_1m) = -dmass
    this_proc%source(rain_params%i_1m) = dmass

    IF (params%l_2m) THEN
      dnumber = dmass*number/mass
      this_proc%source(params%i_2m) = -dnumber
      this_proc%source(rain_params%i_2m) = dnumber
    END IF

    IF (params%l_3m) THEN

      IF (l_meltall) THEN
        dm3=-m3/dt
      ELSE

        dm1=-dt*dmass/params%c_x
        dm2=-dt*dnumber
        m2=number

        CALL m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
        dm3=dm3/dt

      END IF

      this_proc%source(params%i_3m) = dm3
    END IF

    IF (rain_params%l_3m) THEN
      IF (params%l_3m) THEN
        CALL m3_inc_type4(dm3, rain_params%c_x, params%c_x, params%p3, dm3_r)
      ELSE
        m1=qfields(k, rain_params%i_1m)/rain_params%c_x
        m2=qfields(k, rain_params%i_2m)
        m3=qfields(k, rain_params%i_3m)

        dm1=dt*dmass/rain_params%c_x
        dm2=dt*dnumber
        CALL m3_inc_type2(m1, m2, m3, rain_params%p1,     &
             rain_params%p2, rain_params%p3, dm1, dm2, dm3_r &
             , rain_params%fix_mu)
        dm3_r=dm3_r/dt
      END IF
      this_proc%source(rain_params%i_3m) = dm3_r
    END IF

  END IF
      !----------------------
      ! Aerosol processing...
      !----------------------

  IF (l_process) THEN

    aero_proc => aerosol_procs(k, iaproc%id)

    IF (params%id == ice_params%id) THEN
      dmac = dnumber*aeroice(k)%nratio1*aeroice(k)%mact1_mean
      dmad = dnumber*dustact(k)%nratio1*dustact(k)%mact1_mean
    ELSE IF (params%id == snow_params%id) THEN
      dmac = dnumber*aeroice(k)%nratio2*aeroice(k)%mact2_mean
      dmad = dnumber*dustact(k)%nratio2*dustact(k)%mact2_mean
    ELSE IF (params%id == graupel_params%id) THEN
      dmac = dnumber*aeroice(k)%nratio3*aeroice(k)%mact3_mean
      dmad = dnumber*dustact(k)%nratio3*dustact(k)%mact3_mean
    END IF

    aero_proc%source(i_am8) = -dmac
    aero_proc%source(i_am4) = dmac

    aero_proc%source(i_am9) = dmad
    aero_proc%source(i_am7) = -dmad

    NULLIFY(aero_proc)

  END IF


  NULLIFY(this_proc)

END IF

END SUBROUTINE melting


END MODULE ice_melting
