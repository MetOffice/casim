module snow_autoconversion
  use variable_precision, only: wp, iwp
  use passive_fields, only: rho, exner, TdegK, pressure
  use mphys_switches, only: iopt_auto,   &
       i_qi, i_qs, i_ni, i_ns, i_m3s, hydro_complexity, &
       l_dsaut, i_qv, i_th, &
       l_harrington
  use mphys_constants, only: rhow, fixed_ice_number,   &
       Ls, cp,  Lv,ka, Dv, Rv
  use mphys_parameters, only: mu_saut, snow_params, ice_params,   &
       DImax, tau_saut, DI2S
  use process_routines, only: process_rate, i_saut
  use qsat_funs, only: qsaturation, qisaturation
  use thresholds, only: thresh_small
  use special, only: pi
  use m3_incs, only: m3_inc_type2, m3_inc_type3

  use lookup, only: gfunc
  use distributions, only: dist_lambda, dist_mu, dist_n0

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='SNOW_AUTOCONVERSION'

  public saut
contains

  subroutine saut(dt, k, qfields, aerofields, procs, aerosol_procs)

    implicit none

    character(len=*), parameter :: RoutineName='SAUT'

    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    real(wp), intent(in) :: qfields(:,:), aerofields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    real(wp) :: dmass, dnumber
    real(wp) :: m1, m2, m3, dm1,dm2,dm3
    real(wp) :: ice_n0, ice_lam, ice_mu
    real(wp) :: ice_mass
    real(wp) :: ice_number
    real(wp) :: snow_mass
    real(wp) :: snow_number
    real(wp) :: snow_m3
    real(wp) :: th
    real(wp) :: qv
    type(process_rate), pointer :: this_proc
    real(wp) :: p1, p2, p3
    real(wp) :: lami_min , AB, qis
    integer :: i

    logical :: l_condition ! logical condition to switch on process

    l_condition=.true.
    lami_min=0.0

    ice_mass=qfields(k, i_qi)
    if (ice_params%l_2m) then
      ice_number=qfields(k, i_ni)
    else
      ice_number=fixed_ice_number
    end if

    if (ice_mass < thresh_small(i_qi) .or. ice_number < thresh_small(i_ni)) then
      l_condition=.false.
    else
      ice_n0=dist_n0(k,ice_params%id)
      ice_mu=dist_mu(k,ice_params%id)
      ice_lam=dist_lambda(k,ice_params%id)

      lami_min=(1.0 + ice_mu)/DImax

    end if

    if (l_condition) then
      if (l_harrington) then
        qv=qfields(k, i_qv)
        th=qfields(k, i_th)
        qis=qisaturation(th*exner(k), pressure(k)/100.0)
        l_condition=qv > qis
      else
        l_condition=ice_lam < lami_min ! LEM autconversion
      end if
    end if

    if (l_condition) then
      this_proc=>procs(k, i_saut%id)

      p1=snow_params%p1
      p2=snow_params%p2
      p3=snow_params%p3

      snow_mass=qfields(k, i_qs)
      if (snow_params%l_2m) snow_number = qfields(k, i_ns)
      if (snow_params%l_3m) snow_m3 = qfields(k, i_m3s)

      if (l_harrington) then
        !< AB This is used elsewhere, so we should do it more efficiently.
        AB=1.0/(Lv*Lv/(Rv*ka*TdegK(k)*TdegK(k))*rho(k)+1.0/(Dv*qis))
        dnumber=4.0/DImax/ice_params%density*(qv-qis)*rho(k)*ice_number*exp(-ice_lam*DImax)*Dv/AB
        dnumber=min(dnumber,0.9*ice_number/dt)
        dmass=pi/6.0*ice_params%density*DImax**3*dnumber
        dmass=min(dmass, 0.7*ice_mass/dt)
      else ! LEM version
        dmass=((lami_min/ice_lam)**ice_params%d_x - 1.0)*ice_mass/tau_saut
        if (ice_mass < dmass*dt) then
          dmass=.5*ice_mass/dt
        end if
        dmass=min(dmass, .5*ice_mass/dt)
        dnumber=dmass/(ice_params%c_x*DI2S**ice_params%d_x)
      end if

      if (ice_number < dnumber*dt) then
        dnumber=dmass/ice_mass * ice_number
      end if

      if (dmass*dt > thresh_small(i_qs)) then
        if (snow_params%l_3m) then
          dm1=dt*dmass/snow_params%c_x
          dm2=dt*dnumber
          call m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu_saut)
          dm3=dm3/dt
        end if

        this_proc%source(i_qi)=-dmass
        this_proc%source(i_qs)=dmass

        if (ice_params%l_2m) then
          this_proc%source(i_ni)=-dnumber
        end if
        if (snow_params%l_2m) then
          this_proc%source(i_ns)=dnumber
        end if
        if (snow_params%l_3m) then
          this_proc%source(i_m3s)=dm3
        end if
      end if
      !==============================
      ! No aerosol processing needed
      !==============================
      nullify(this_proc)
    end if
  end subroutine saut
end module snow_autoconversion
