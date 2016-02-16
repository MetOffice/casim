module accretion
  use variable_precision, only: wp
  use passive_fields, only: rho
  use mphys_switches, only: i_ql, i_qr, i_nl, i_nr, i_m3r, l_2mc, l_2mr, &
       l_aacc, i_am4, i_am5, l_process, active_rain, isol, l_preventsmall
  use mphys_constants, only: fixed_cloud_number
  use mphys_parameters, only: hydro_params
  use process_routines, only: process_rate, i_pracw, i_aacw
  use thresholds, only: ql_small, qr_small
  use sweepout_rate, only: sweepout
  use distributions, only: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
  use diagnostics, only: save_dg, i_dgtime, k_here, i_here
#endif
  
  implicit none

  private

  public racw
contains
  !---------------------------------------------------------------------------
  !> @author
  !> Ben Shipway
  !
  !> @brief
  !> This subroutine calculates increments due to the accretion of
  !> cloud water by rain
  !--------------------------------------------------------------------------- !
  subroutine racw(dt, k, qfields, aerofields, procs, params, aerosol_procs)
    real(wp), intent(in) :: dt !< microphysics time increment (s)
    integer,  intent(in)  :: k  !< k index
    real(wp), intent(in) :: qfields(:,:)     !< hydrometeor fields
    real(wp), intent(in) :: aerofields(:,:)  !< aerosol fields
    type(process_rate), intent(inout), target :: procs(:,:)         !< hydrometeor process rates
    type(process_rate), intent(inout), target :: aerosol_procs(:,:) !< aerosol process rates
    type(hydro_params), intent(in) :: params !< parameters describing hydrometor size distribution/fallspeeds etc.

    real(wp) :: dmass, dnumber, damass
    real(wp) :: m1, m2, dm1, dm2

    real(wp) :: cloud_mass
    real(wp) :: cloud_number
    real(wp) :: rain_mass
    real(wp) :: rain_number

    type(process_rate), pointer :: this_proc
    type(process_rate), pointer :: aero_proc

    real(wp) :: mu, n0, lam
    logical :: l_kk_acw=.true.

    cloud_mass=qfields(k, i_ql)
    if (l_2mc) then
      cloud_number=qfields(k, i_nl)
    else
      cloud_number=fixed_cloud_number
    end if
    rain_mass=qfields(k, i_qr)
    if (l_2mr) rain_number=qfields(k, i_nr)

    if (cloud_mass > ql_small .and. rain_mass > qr_small) then
      this_proc=>procs(k, i_pracw%id)
      if (l_kk_acw) then
        dmass=min(0.9*cloud_mass, 67.0*(cloud_mass*rain_mass)**1.15)
      else
        n0=dist_n0(k,params%id)
        mu=dist_mu(k,params%id)
        lam=dist_lambda(k,params%id)
        dmass=sweepout(n0, lam, mu, params, rho(k))*cloud_mass
      end if

      if (l_preventsmall .and. dmass < qr_small) dmass=0.0
      if (l_2mc) dnumber=dmass/(cloud_mass/cloud_number)

      this_proc%source(i_ql)=-dmass
      this_proc%source(i_qr)=dmass
      if (l_2mc) then
        this_proc%source(i_nl)=-dnumber
      end if

      if (l_aacc .and. l_process) then
        aero_proc=>aerosol_procs(k, i_aacw%id)
        if (active_rain(isol)) then
          damass=dmass/cloud_mass*aerofields(k,i_am4)
          aero_proc%source(i_am4)=-damass
          aero_proc%source(i_am5)=damass
        end if
        nullify(aero_proc)
      end if
      nullify(this_proc)
    end if
  end subroutine racw
end module accretion
