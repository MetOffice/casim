module homogeneous
  use variable_precision, only: wp
  use passive_fields, only: rho, pressure, w, exner
  use mphys_switches, only: i_qv, i_ql, i_qi, i_ni, i_th , hydro_complexity, i_am6, i_an2, l_2mi, l_2ms, l_2mg &
       , i_ns, i_ng, iopt_inuc, i_am7, i_an6, i_am9 , i_m3r, i_m3g, i_qr, i_qg, i_nr, i_ng, i_nl          &
       , isol, i_am4, i_am8, active_ice, l_process
  use process_routines, only: process_rate, i_homr, i_homc, i_dhomc, i_dhomr
  use mphys_parameters, only: rain_params, graupel_params, cloud_params, ice_params, T_hom_freeze
  use mphys_constants, only: Lf, cp
  use qsat_funs, only: qsaturation, qisaturation
  use thresholds, only: thresh_small, thresh_tidy
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use special, only: GammaFunc, pi
  use m3_incs, only: m3_inc_type2, m3_inc_type4

  implicit none

  character(len=*), parameter, private :: ModuleName='HOMOGENEOUS'

contains

  !> Calculates immersion freezing of rain drops
  !> See Bigg 1953
  subroutine ihom_rain(dt, k, qfields, aeroact, dustliq, procs, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    real(wp), intent(in), target :: qfields(:,:)
    type(aerosol_active), intent(in) :: aeroact(:), dustliq(:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    real(wp) :: dmass, dnumber, dmac, coef, dmadl
    real(wp) :: dm1, dm2, dm3, dm3_g, m1, m2, m3
    real(wp) :: n0, lam, mu

    real(wp) :: th
    real(wp) :: qv, qr, nr
    type(process_rate), pointer :: this_proc
    type(process_rate), pointer :: aero_proc

    real(wp) :: Tc

    real(wp), parameter :: A_bigg = 0.66, B_bigg = 100.0

    logical :: l_condition, l_freezeall
    logical :: l_ziegler=.true.

    character(len=*), parameter :: RoutineName='IHOM_RAIN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    th = qfields(k, i_th)
    Tc = th*exner(k) - 273.15
    qr = qfields(k, i_qr)
    if (rain_params%l_2m)then
      nr = qfields(k, i_nr)
    else
      nr = 1000. ! PRAGMATIC SM HACK
    end if

    l_condition=(Tc < -4.0 .and. qr > thresh_small(i_qr))
    l_freezeall=.false.
    if (l_condition) then
      this_proc => procs(k, i_homr%id)

      n0 = dist_n0(k,rain_params%id)
      mu = dist_mu(k,rain_params%id)
      lam = dist_lambda(k,rain_params%id)

      if (l_ziegler) then
        dnumber = B_bigg*(exp(-A_bigg*Tc)-1.0)*rho(k)*qr/rain_params%density
        dnumber = min(dnumber, nr/dt)
        dmass = (qr/nr)*dnumber
      else
        coef = B_bigg*(pi/6.0)*(exp(-A_bigg*Tc)-1.0)/rho(k)               &
             * n0 /(lam*lam*lam)/GammaFunc(1.0 + mu)

        dmass = coef * rain_params%c_x * lam**(-rain_params%d_x)          &
             * GammaFunc(4.0 + mu + rain_params%d_x)

        dnumber = coef                                                    &
             * GammaFunc(4.0 + mu)

      end if

      ! PRAGMATIC HACK - FIX ME
      ! If most of the drops are frozen, do all of them
      if (dmass*dt >0.95*qr .or. dnumber*dt > 0.95*nr) then
        dmass=qr/dt
        dnumber=nr/dt
        l_freezeall=.true.
      end if

      this_proc%source(i_qr) = -dmass
      this_proc%source(i_qg) = dmass

      if (rain_params%l_2m) then
        this_proc%source(i_nr) = -dnumber
      end if
      if (graupel_params%l_2m) then
        this_proc%source(i_ng) = dnumber
      end if

      if (rain_params%l_3m) then
        if (l_freezeall) then
          dm3=-qfields(k,i_m3r)/dt
        else
          m1=qr/rain_params%c_x
          m2=qfields(k,i_nr)
          m3=qfields(k,i_m3r)

          dm1=-dt*dmass/rain_params%c_x
          dm2=-dt*dnumber

          call m3_inc_type2(m1, m2, m3, rain_params%p1, rain_params%p2, rain_params%p3, dm1, dm2, dm3)
          dm3=dm3/dt
        end if
        this_proc%source(i_m3r) = dm3
      end if

      if (graupel_params%l_3m) then
        if (rain_params%l_3m) then
          call m3_inc_type4(dm3, graupel_params%c_x, rain_params%c_x, rain_params%p3, dm3_g)
        else
          m1=qfields(k,i_qg)/graupel_params%c_x
          m2=qfields(k,i_ng)
          m3=qfields(k,i_m3g)

          dm1=-dt*dmass/graupel_params%c_x
          dm2=-dt*dnumber
          call m3_inc_type2(m1, m2, m3, graupel_params%p1, graupel_params%p2, graupel_params%p3, dm1, dm2, dm3_g)
          dm3_g=dm3_g/dt
        end if
        this_proc%source(i_m3g) = dm3_g
      end if

      if (l_process) then
        aero_proc => aerosol_procs(k, i_dhomr%id)
        dmac = dnumber*aeroact(k)%mact2_mean

        aero_proc%source(i_am8) = dmac
        aero_proc%source(i_am4) = -dmac

        ! Dust already in the liquid phase
        dmadl = dnumber*dustliq(k)%mact2_mean*dustliq(k)%nratio2
        if (dmadl /=0.0) then
          aero_proc%source(i_am9) = -dmadl
          aero_proc%source(i_am7) = dmadl
        end if
        nullify(aero_proc)
      end if
      nullify(this_proc)
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ihom_rain

  !> Calculates homogeneous freezing of cloud drops
  !> See Wisener 1972
  subroutine ihom_droplets(dt, k, qfields, aeroact, dustliq, procs, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    real(wp), intent(in), target :: qfields(:,:)
    ! aerosol fields
    type(aerosol_active), intent(in) :: aeroact(:), dustliq(:)

    type(process_rate), intent(inout), target :: procs(:,:)
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    real(wp) :: dmass, dnumber, dmac, coef, dmadl
    real(wp) :: dm1, dm2, dm3, m1, m2, m3
    real(wp) :: n0, lam, mu

    real(wp) :: th
    real(wp) :: qv, ql
    type(process_rate), pointer :: this_proc
    type(process_rate), pointer :: aero_proc

    real(wp) :: Tc

    logical :: l_condition

    character(len=*), parameter :: RoutineName='IHOM_DROPLETS'


    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    th = qfields(k, i_th)
    Tc = th*exner(k) - 273.15
    ql = qfields(k, i_ql)

    l_condition=(Tc < T_hom_freeze .and. ql > thresh_tidy(i_ql))
    if (l_condition) then
      this_proc => procs(k, i_homc%id)
      dmass=min(ql, cp*(T_hom_freeze - Tc)/Lf)/dt
      dnumber=dmass*qfields(k, i_nl)/ql

      this_proc%source(i_ql)=-dmass
      this_proc%source(i_qi)=dmass

      if (cloud_params%l_2m) then
        this_proc%source(i_nl)=-dnumber
      end if
      if (ice_params%l_2m) then
        this_proc%source(i_ni)=dnumber
      end if

      if (l_process) then
        aero_proc=>aerosol_procs(k, i_dhomc%id)
        dmac=dnumber*aeroact(k)%mact1_mean

        aero_proc%source(i_am8)=dmac
        aero_proc%source(i_am4)=-dmac

        ! Dust already in the liquid phase
        dmadl=dnumber*dustliq(k)%mact1_mean*dustliq(k)%nratio1
        if (dmadl /=0.0) then
          aero_proc%source(i_am9)=-dmadl
          aero_proc%source(i_am7)=dmadl
        end if
        nullify(aero_proc)
      end if
      nullify(this_proc)
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ihom_droplets
end module homogeneous
