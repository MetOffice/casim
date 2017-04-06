module ice_melting
  use variable_precision, only: wp
  use process_routines, only: process_rate, process_name, i_imlt, i_smlt, i_gmlt, i_sacw, i_sacr, i_gacw, i_gacr, i_gshd, &
       i_dimlt, i_dsmlt, i_dgmlt
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use passive_fields, only: TdegC, qws0, rho
  use mphys_parameters, only: ice_params, snow_params, graupel_params, rain_params, DR_melt, hydro_params, ZERO_REAL_WP
  use mphys_switches, only: i_qv, i_am4, i_am7, i_am8, i_am9, l_process
  use mphys_constants, only: Lv, Lf, Ka, Cwater, Cice, cp, Dv
  use thresholds, only: thresh_tidy
  use m3_incs, only: m3_inc_type2, m3_inc_type3, m3_inc_type4
  use ventfac, only: ventilation
  use distributions, only: dist_lambda, dist_mu, dist_n0

  implicit none
contains

  !> Subroutine to calculate rate of melting of ice species
  !>
  !> OPTIMISATION POSSIBILITIES: Shouldn't have to recalculate all 3m quantities
  !>                             If just rescaling mass conversion for dry mode
  subroutine melting(dt, k, params, qfields, procs, aeroice, dustact, aerosol_procs)
    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    type(hydro_params), intent(in) :: params
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! aerosol fields
    type(aerosol_active), intent(in) :: aeroice(:), dustact(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    type(process_name) :: iproc, iaproc ! processes selected depending on
    ! which species we're modifying
    type(process_name) :: i_acw, i_acr ! accretion processes
    real(wp) :: qv
    real(wp) :: dmass, dnumber, dm1, dm2, dm3, dm3_r
    real(wp) :: number, mass, m1, m2, m3
    type(process_rate), pointer :: this_proc, aero_proc
    real(wp) :: n0, lam, mu, V_x
    real(wp) :: acc_correction
    logical :: l_meltall ! do we melt everything?
    real(wp) :: dmac, dmad

    l_meltall=.false.
    mass=qfields(k, params%i_1m)
    if (mass > thresh_tidy(params%i_1m) .and. TdegC(k) > 0.0) then
      if (params%l_2m) number=qfields(k, params%i_2m)
      if (params%id == ice_params%id) then ! instantaneous removal
        iaproc=i_dimlt
        iproc=i_imlt
        this_proc=>procs(k, iproc%id)
        dmass=mass/dt

        if (params%l_2m)then
          dnumber=number/dt
          this_proc%source(ice_params%i_2m)=-dnumber
          if (rain_params%l_2m) then
            this_proc%source(rain_params%i_2m)=dnumber
          end if
        end if

        this_proc%source(ice_params%i_1m)=-dmass
        this_proc%source(rain_params%i_1m)=dmass

        if (rain_params%l_3m) then
          m1=qfields(k, rain_params%i_1m)/rain_params%c_x
          m2=qfields(k, rain_params%i_2m)
          m3=qfields(k, rain_params%i_3m)
          dm1=dt*dmass/rain_params%c_x
          dm2=dt*dnumber

          if (m1 > 0.0) then
            call m3_inc_type2(m1, m2, m3, rain_params%p1,      &
                 rain_params%p2, rain_params%p3, dm1, dm2, dm3_r)
          else
            call m3_inc_type3(rain_params%p1, rain_params%p2, rain_params%p3,      &
                 dm1, dm2, dm3_r, rain_params%fix_mu)
          end if
          dm3_r=dm3_r/dt
          this_proc%source(rain_params%i_3m) = dm3_r
        end if
      else
        if (params%id==snow_params%id) then
          i_acw=i_sacw
          i_acr=i_sacr
          iproc=i_smlt
          iaproc=i_dsmlt
        else if (params%id==graupel_params%id) then
          i_acw=i_gacw
          i_acr=i_gacr
          iproc=i_gmlt
          iaproc=i_dgmlt
        end if
        acc_correction=0.0
        if (i_acw%on)acc_correction=procs(k, i_acw%id)%source(params%i_1m)
        if (i_acr%on)acc_correction=acc_correction + procs(k, i_acr%id)%source(params%i_1m)
        if (params%id==graupel_params%id .and. i_gshd%on) then
          acc_correction=acc_correction + procs(k, i_gshd%id)%source(params%i_1m)
        end if

        qv=qfields(k, i_qv)

        m1=mass/params%c_x
        if (params%l_2m) number=qfields(k, params%i_2m)
        if (params%l_3m) m3=qfields(k, params%i_3m)

        n0=dist_n0(k,params%id)
        mu=dist_mu(k,params%id)
        lam=dist_lambda(k,params%id)

        call ventilation(k, V_x, n0, lam, mu, params)

        dmass=(1.0/(rho(k)*Lf))*(Ka*TdegC(k) + Lv*Dv*rho(k)*(qv - qws0(k))) * V_x&
             + (Cwater*TdegC(k)/Lf)*acc_correction

        dmass=max(dmass, ZERO_REAL_WP) ! ensure positive

        if (dmass == ZERO_REAL_WP) return  ! No need to do anything

        dmass=min(dmass, mass/dt) ! ensure we don't remove too much
        if (dmass*dt > 0.95*mass) then ! we're pretty much removing everything
          l_meltall=.true.
          dmass=mass/dt
        end if

        !--------------------------------------------------
        ! Apply spontaneous rain breakup if drops are large
        !--------------------------------------------------
        !< RAIN BREAKUP TO BE ADDED

        this_proc=>procs(k, iproc%id)
        this_proc%source(params%i_1m)=-dmass
        this_proc%source(rain_params%i_1m)=dmass

        if (params%l_2m) then
          dnumber=dmass*number/mass
          this_proc%source(params%i_2m)=-dnumber
          this_proc%source(rain_params%i_2m)=dnumber
        end if
        if (params%l_3m) then
          if (l_meltall) then
            dm3=-m3/dt
          else
            dm1=-dt*dmass/params%c_x
            dm2=-dt*dnumber
            m2=number
            call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
            dm3=dm3/dt
          end if
          this_proc%source(params%i_3m)=dm3
        end if

        if (rain_params%l_3m) then
          if (params%l_3m) then
            call m3_inc_type4(dm3, rain_params%c_x, params%c_x, params%p3, dm3_r)
          else
            m1=qfields(k, rain_params%i_1m)/rain_params%c_x
            m2=qfields(k, rain_params%i_2m)
            m3=qfields(k, rain_params%i_3m)

            dm1=dt*dmass/rain_params%c_x
            dm2=dt*dnumber
            call m3_inc_type2(m1, m2, m3, rain_params%p1, rain_params%p2, rain_params%p3, dm1, dm2, dm3_r, rain_params%fix_mu)
            dm3_r=dm3_r/dt
          end if
          this_proc%source(rain_params%i_3m) = dm3_r
        end if
      end if
      !----------------------
      ! Aerosol processing...
      !----------------------

      if (l_process) then
        aero_proc=>aerosol_procs(k, iaproc%id)

        if (params%id == ice_params%id) then
          dmac=dnumber*aeroice(k)%nratio1*aeroice(k)%mact1_mean
          dmad=dnumber*dustact(k)%nratio1*dustact(k)%mact1_mean
        else if (params%id == snow_params%id) then
          dmac=dnumber*aeroice(k)%nratio2*aeroice(k)%mact2_mean
          dmad=dnumber*dustact(k)%nratio2*dustact(k)%mact2_mean
        else if (params%id == graupel_params%id) then
          dmac=dnumber*aeroice(k)%nratio3*aeroice(k)%mact3_mean
          dmad=dnumber*dustact(k)%nratio3*dustact(k)%mact3_mean
        end if

        aero_proc%source(i_am8)=-dmac
        aero_proc%source(i_am4)=dmac
        aero_proc%source(i_am9)=dmad
        aero_proc%source(i_am7)=-dmad
        nullify(aero_proc)
      end if
      nullify(this_proc)
    end if
  end subroutine melting
end module ice_melting
