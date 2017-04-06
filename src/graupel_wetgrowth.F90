module graupel_wetgrowth
  use variable_precision, only: wp
  use process_routines, only: process_rate, process_name,   &
       i_gshd, i_gacw, i_gacr, i_gaci, i_gacs
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use passive_fields, only: TdegC,  qws0, rho

  use mphys_parameters, only: ice_params, snow_params, graupel_params,   &
       rain_params, T_hom_freeze, DR_melt
  use mphys_switches, only: i_qv
  use mphys_constants, only: Lv, Lf, Ka, Cwater, Cice
  use thresholds, only: thresh_small

  use m3_incs, only: m3_inc_type2
  use ventfac, only: ventilation
  use distributions, only: dist_lambda, dist_mu, dist_n0

  implicit none
  private

  public wetgrowth
contains
  subroutine wetgrowth(dt, k, qfields, procs, aerophys, aerochem, aeroact , aerosol_procs)
    !< Subroutine to determine if all liquid accreted by graupel can be
    !< frozen or if there will be some shedding
    !<
    !< CODE TIDYING: Should move efficiencies into parameters

    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! aerosol fields
    type(aerosol_phys), intent(in) :: aerophys(:)
    type(aerosol_chem), intent(in) :: aerochem(:)
    type(aerosol_active), intent(in) :: aeroact(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), optional :: aerosol_procs(:,:)

    type(process_name) :: iproc ! processes selected depending on
    ! which species we're modifying

    real(wp) :: qv
    real(wp) :: dmass, dnumber, dm1, dm2, dm3
    real(wp) :: number, mass, m1, m2, m3
    real(wp) :: gacw, gacr, gaci, gacs  ! graupel accretion process rates
    real(wp) :: dnumber_s, dnumber_g  ! number conversion rate from snow/graupel
    real(wp) :: dmass_s, dmass_g      ! mass conversion rate from snow/graupel

    type(process_rate), pointer :: this_proc

    real(wp) :: n0, lam, mu, V_x
    real(wp) :: pgacw, pgacr, pgaci, pgacs
    real(wp) :: Eff, sdryfac, idryfac
    real(wp) :: pgaci_dry, pgacs_dry, pgdry
    real(wp) :: pgwet       !< Amount of liquid that graupel can freeze withouth shedding
    real(wp) :: pgacsum     !< Sum of all graupel accretion terms

    mass=qfields(k, graupel_params%i_1m)

    if (mass > thresh_small(graupel_params%i_1m)) then
      pgacw=procs(k, i_gacw%id)%source(graupel_params%i_1m)
      pgacr=procs(k, i_gacr%id)%source(graupel_params%i_1m)
      pgaci=procs(k, i_gaci%id)%source(graupel_params%i_1m)
      pgacs=procs(k, i_gacs%id)%source(graupel_params%i_1m)

      ! Factors for converting wet collection efficiencies to dry ones
      Eff=1.0
      sdryfac=min(1.0_wp, 0.2*exp(0.08*TdegC(k))/Eff)
      Eff=1.0
      idryfac=min(1.0_wp, 0.2*exp(0.08*TdegC(k))/Eff)

      pgaci_dry=idryfac*pgaci
      pgacs_dry=sdryfac*pgacs

      pgdry=pgacw+pgacr+pgaci_dry+pgacs_dry

      pgacsum=pgacw+pgacr+pgaci+pgacs

      if (pgacsum > thresh_small(graupel_params%i_1m)) then
        qv=qfields(k, i_qv)
        mass=qfields(k, graupel_params%i_1m)
        if (graupel_params%l_2m) number=qfields(k, graupel_params%i_2m)

        n0=dist_n0(k,graupel_params%id)
        mu=dist_mu(k,graupel_params%id)
        lam=dist_lambda(k,graupel_params%id)

        call ventilation(k, V_x, n0, lam, mu, graupel_params)

        pgwet=(910.0/graupel_params%density)**0.625*(Lv*(qws0(k)-qv)-Ka*TdegC(k)/rho(k))/(Lf+Cwater*TdegC(k))*V_x
        pgwet=pgwet+(pgaci+pgacs)*(1.0-Cice*TdegC(k)/(Lf+Cwater*TdegC(k)))

        if (pgdry < pgwet .or. TdegC(k) < T_hom_freeze) then ! Dry growth, so use recalculated gaci, gacs
          if (pgaci + pgacs > 0.0) then
            this_proc=>procs(k, i_gacs%id)
            dmass=pgacs_dry
            this_proc%source(graupel_params%i_1m)=dmass
            this_proc%source(snow_params%i_1m)=-dmass
            if (snow_params%l_2m) then
              dnumber=this_proc%source(snow_params%i_2m)*sdryfac
              this_proc%source(snow_params%i_2m)=dnumber
            end if
            this_proc=>procs(k, i_gaci%id)
            dmass=pgaci_dry
            this_proc%source(graupel_params%i_1m)=dmass
            this_proc%source(ice_params%i_1m)=-dmass
            if (ice_params%l_2m) then
              dnumber=this_proc%source(ice_params%i_2m)*idryfac
              this_proc%source(ice_params%i_2m)=dnumber
            end if
          end if
        else ! Wet growth mode so recalculate gacr, gshd
          this_proc=>procs(k, i_gacr%id)
          if (abs(procs(k, i_gacr%id)%source(graupel_params%i_1m)) > thresh_small(graupel_params%i_1m)) then
            dmass=pgacr+pgwet-pgacsum
            if (dmass > 0.0) then
              this_proc%source(graupel_params%i_1m)=dmass
              this_proc%source(rain_params%i_1m)=-dmass
            else
              iproc=i_gshd
              this_proc=>procs(k, iproc%id)
              dmass=min(pgacw, -1.0*dmass)
              this_proc%source(rain_params%i_1m)=dmass
              if (rain_params%l_2m) then
                dnumber=dmass/(rain_params%c_x*DR_melt**3)
                this_proc%source(rain_params%i_2m)=dnumber
              end if
            end if
          end if
        end if
      end if

      !----------------------
      ! Aerosol processing...
      !----------------------
      nullify(this_proc)
    end if
  end subroutine wetgrowth
end module graupel_wetgrowth
