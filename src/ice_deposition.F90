module ice_deposition
  use variable_precision, only: wp, iwp
  use passive_fields, only: rho, pressure, w, exner, TdegK
  use mphys_switches, only: i_qv, i_qi, i_ni, i_th   &
       , hydro_complexity, i_am6, i_an2, l_2mi, l_2ms, l_2mg &
       , i_ns, i_ng, iopt_inuc, i_am7, i_an6                   &
       , l_process, l_passivenumbers_ice, l_passivenumbers, active_number &
       , active_ice, iinsol, i_an12 &
       , i_qr, i_ql, i_qs, i_qg, i_am8, i_am2, i_an11
  use type_process, only: process_name
  use process_routines, only: process_rate, i_idep,    &
       i_dsub, i_sdep, i_gdep, i_dssub, i_dgsub &
       , i_isub, i_ssub, i_gsub, i_iacw, i_raci, i_sacw, i_sacr  &
       , i_gacw, i_gacr
  use mphys_parameters, only: hydro_params, ice_params, rain_params, cloud_params
  use mphys_constants, only: Ls, cp,  Lv, Lf, ka, Dv, Rv
  use qsat_funs, only: qsaturation, qisaturation
  use thresholds, only: thresh_small
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active

  use distributions, only: dist_lambda, dist_mu, dist_n0
  use ventfac, only: ventilation
  use special, only: pi, Gammafunc
  use m3_incs, only: m3_inc_type2

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='ICE_DEPOSITION'

  logical :: l_latenteffects = .false.

  public idep
contains

  !< Subroutine to determine the deposition/sublimation onto/from
  !< ice, snow and graupel.  There is no source/sink for number
  !< when undergoing deposition, but there is a sink when sublimating.
  subroutine idep(dt, k, params, qfields, procs, dustact, aeroice, aerosol_procs)

    implicit none

    character(len=*), parameter :: RoutineName='IDEP'

    real(wp), intent(IN) :: dt
    integer, intent(IN) :: k
    type(hydro_params), intent(IN) :: params
    real(wp), intent(IN) :: qfields(:,:)
    type(process_rate), intent(INOUT), target :: procs(:,:)

    ! aerosol fields
    type(aerosol_active), intent(IN) :: dustact(:), aeroice(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(INOUT), target :: aerosol_procs(:,:)

    type(process_name) :: iproc, iaproc  ! processes selected depending on
    ! which species we're depositing on.

    type(process_name) :: i_acw, i_acr ! collection processes with cloud and rain
    real(wp) :: dmass, dnumber, dmad, dnumber_a, dnumber_d

    real(wp) :: th
    real(wp) :: qv, qh
    real(wp) :: number, mass, m1,m2, m3, dm1,dm2,dm3
    type(process_rate), pointer :: this_proc
    type(process_rate), pointer :: aero_proc

    real(wp) :: qs, qis, Si, Sw, limit
    real(wp) :: n0, lam, mu
    real(wp) :: V_x, AB
    real(wp) :: frac ! fraction of ice below a threshold

    logical :: l_suball

    l_suball=.false. ! do we want to sublimate everything?

    mass=qfields(k, params%i_1m)

    qv=qfields(k, i_qv)
    th=qfields(k, i_th)

    qis=qisaturation(th*exner(k), pressure(k)/100.0)

    if (qv>qis) then
      select case (params%id)
      case (3_iwp) !ice
        iproc=i_idep
        iaproc=i_dsub
        i_acw=i_iacw
        i_acr=i_raci
      case (4_iwp) !snow
        iproc=i_sdep
        iaproc=i_dssub
        i_acw=i_sacw
        i_acr=i_sacr
      case (5_iwp) !graupel
        iproc=i_gdep
        iaproc=i_dgsub
        i_acw=i_gacw
        i_acr=i_gacr
      end select
    else
      select case (params%id)
      case (3_iwp) !ice
        iproc=i_isub
        iaproc=i_dsub
        i_acw=i_iacw
        i_acr=i_raci
      case (4_iwp) !snow
        iproc=i_ssub
        iaproc=i_dssub
        i_acw=i_sacw
        i_acr=i_sacr
      case (5_iwp) !graupel
        iproc=i_gsub
        iaproc=i_dgsub
        i_acw=i_gacw
        i_acr=i_gacr
      end select
    end if

    if (mass > thresh_small(params%i_1m)) then ! if no existing ice, we don't grow/deplete it.

      this_proc=>procs(k, iproc%id)
      if (params%l_2m) number=qfields(k, params%i_2m)
      if (params%l_3m) m3=qfields(k, params%i_3m)

      n0=dist_n0(k,params%id)
      mu=dist_mu(k,params%id)
      lam=dist_lambda(k,params%id)

      call ventilation(k, V_x, n0, lam, mu, params)

      AB=1.0/(Ls*Ls/(Rv*ka*TdegK(k)*TdegK(k))*rho(k)+1.0/(Dv*qis))
      dmass=(qv/qis-1.0)*V_x*AB

      ! Include latent heat effects of collection of rain and cloud
      ! as done in Milbrandt & Yau (2005)
      if (l_latenteffects) then
        dmass=dmass - Lf*Ls/(Rv*ka*TdegK(k)*TdegK(k))         &
             *(procs(k, i_acw%id)%source(cloud_params%i_1m) &
             + procs(k, i_acr%id)%source(rain_params%i_1m))
      end if

      ! Check we haven't become subsaturated and limit if we have (dep only)
      ! NB doesn't account for simultaneous ice/snow growth - checked elsewhere
      if (dmass > 0.0) dmass=min((qv-qis)/dt,dmass)
      ! Check we don't remove too much (sub only)
      if (dmass < 0.0) dmass=max(-mass/dt,dmass)

      this_proc%source(i_qv)=-dmass
      this_proc%source(params%i_1m)=dmass

      if (params%l_2m) then
        dnumber=0.0
        if (dmass < 0.0) dnumber=dmass*number/mass
      end if

      if (-dmass*dt >0.98*mass .or. (params%l_2m .and. -dnumber*dt > 0.98*number)) then
        l_suball=.true.
        dmass=-mass/dt
        dnumber=-number/dt
      end if

      if (params%l_2m) this_proc%source(params%i_2m)=dnumber

      if (dmass < 0.0 .and. l_process) then ! Only process aerosol if sublimating
        aero_proc=>aerosol_procs(k, iaproc%id)
        if (iaproc%id==i_dsub%id) then
          dmad=dnumber*dustact(k)%mact1_mean*dustact(k)%nratio1
          dnumber_d=dnumber*dustact(k)%nratio1
        else if (iaproc%id==i_dssub%id) then
          dmad=dnumber*dustact(k)%mact2_mean*dustact(k)%nratio2
          dnumber_d=dnumber*dustact(k)%nratio2
        else if (iaproc%id==i_dgsub%id) then
          dmad=dnumber*dustact(k)%mact3_mean*dustact(k)%nratio3
          dnumber_d=dnumber*dustact(k)%nratio3
        end if

        aero_proc%source(i_am7)=dmad
        aero_proc%source(i_am6)=-dmad       ! <WARNING: putting back in coarse mode

        if (iaproc%id==i_dsub%id) then
          dmad=dnumber*aeroice(k)%mact1_mean*aeroice(k)%nratio1
          dnumber_a=dnumber*aeroice(k)%nratio1
        else if (iaproc%id==i_dssub%id) then
          dmad=dnumber*aeroice(k)%mact2_mean*aeroice(k)%nratio2
          dnumber_a=dnumber*aeroice(k)%nratio2
        else if (iaproc%id==i_dgsub%id) then
          dmad=dnumber*aeroice(k)%mact3_mean*aeroice(k)%nratio3
          dnumber_a=dnumber*aeroice(k)%nratio3
        end if

        aero_proc%source(i_am8)=dmad
        aero_proc%source(i_am2)=-dmad    ! <WARNING: putting back in accumulation mode

        if (l_passivenumbers_ice) then
          aero_proc%source(i_an12)=dnumber_d
        end if
        aero_proc%source(i_an6)=-dnumber_d  ! <WARNING: putting back in coarse mode

        if (l_passivenumbers) then
          aero_proc%source(i_an11)=dnumber_a
        end if
        aero_proc%source(i_an2)=-dnumber_a  ! <WARNING: putting back in accumulation mode
        nullify(aero_proc)
      end if
      nullify(this_proc)
    end if
  end subroutine idep
end module ice_deposition
