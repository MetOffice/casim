module ice_multiplication
  use variable_precision, only: wp, iwp
  use process_routines, only: process_rate, process_name, i_gacw, i_sacw, i_ihal
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use passive_fields, only: TdegC
  use mphys_parameters, only: ice_params, snow_params, graupel_params, dN_hallet_mossop, M0_hallet_mossop
  use thresholds, only: thresh_small
  use m3_incs, only: m3_inc_type2
  use distributions, only: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
  use diagnostics, only: save_dg, i_dgtime, i_here, k_here
  use runtime, only: time
#endif

  implicit none
contains
  !> Subroutine to determine the ice splintering by Hallet-Mossop
  !> This effect requires prior calculation of the accretion rate of
  !> graupel and snow.
  !> This is a source of ice number and mass and a sink of liquid
  !> (but this is done via the accretion processes already so is
  !> represented here as a sink of snow/graupel)
  !> For triple moment species there is a corresponding change in the
  !> 3rd moment assuming shape parameter is not changed
  !>
  !> AEROSOL: All aerosol sinks/sources are assumed to come from soluble modes
  !
  !> OPTIMISATION POSSIBILITIES:
  subroutine hallet_mossop(dt, k, qfields, procs, aerophys, aerochem, aeroact , aerosol_procs)
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

    real(wp) :: dnumber, dm1, dm2, dm3
    real(wp) :: number, mass, m1, m2, m3
    real(wp) :: gacw, sacw  ! accretion process rates
    real(wp) :: dnumber_s, dnumber_g  ! number conversion rate from snow/graupel
    real(wp) :: dmass_s, dmass_g      ! mass conversion rate from snow/graupel
    real(wp) :: n0, lam, mu
    real(wp) :: Eff  !< splintering efficiency

    type(process_rate), pointer :: this_proc
    type(process_name) :: iproc ! processes selected depending on which species we're modifying

    if (.not. ice_params%l_2m) return

    Eff=1.0 - abs(TdegC(k) + 5.0)/2.5 ! linear increase between -2.5/-7.5 and -5C

    if (Eff > 0.0) then
      sacw=0.0
      gacw=0.0
      if (snow_params%i_1m > 0)   sacw=procs(k, i_sacw%id)%source(snow_params%i_1m)
      if (graupel_params%i_1m > 0) gacw=procs(k, i_gacw%id)%source(graupel_params%i_1m)

      if ((sacw + gacw)*dt > thresh_small(snow_params%i_1m)) then
        iproc=i_ihal
        this_proc=>procs(k, iproc%id)

        dnumber_g=dN_hallet_mossop * Eff * (gacw) ! Number of splinters from graupel
        dnumber_s=dN_hallet_mossop * Eff * (sacw) ! Number of splinters from snow

        dnumber_g=min(dnumber_g, 0.5*gacw/M0_hallet_mossop) ! don't remove more than 50% of rimed liquid
        dnumber_s=min(dnumber_s, 0.5*sacw/M0_hallet_mossop) ! don't remove more than 50% of rimed liquid

        dmass_g=dnumber_g * M0_hallet_mossop
        dmass_s=dnumber_s * M0_hallet_mossop

        !-------------------
        ! Sources for ice...
        !-------------------
        this_proc%source(ice_params%i_1m)=dmass_g+dmass_s
        this_proc%source(ice_params%i_2m)=dnumber_g+dnumber_s

        !-------------------
        ! Sinks for snow...
        !-------------------
        if (sacw > 0.0) then
          this_proc%source(snow_params%i_1m)=-dmass_s
          this_proc%source(snow_params%i_2m)=0.0
        end if

        !---------------------
        ! Sinks for graupel...
        !---------------------
        if (gacw > 0.0) then
          this_proc%source(graupel_params%i_1m)=-dmass_g
          this_proc%source(graupel_params%i_2m)=0.0
        end if

        !----------------------
        ! Aerosol processing...
        !----------------------
      end if
      nullify(this_proc)
    end if
  end subroutine hallet_mossop
end module ice_multiplication
