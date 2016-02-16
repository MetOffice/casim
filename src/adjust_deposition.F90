module adjust_deposition
  ! As in Harrington et al. 1995 move some of the depositional
  ! growth on ice into the snow category

  use variable_precision, only: wp, iwp
  use passive_fields, only: rho
  use mphys_switches, only: i_qi, i_qs, i_ns, i_m3s
  use mphys_parameters, only: DImax, snow_params, ice_params
  use distributions, only: dist_lambda, dist_mu
  use process_routines, only: process_rate, i_idep, i_sdep, i_saut
  use m3_incs, only: m3_inc_type2

  implicit none
  private

  public adjust_dep
contains

  subroutine adjust_dep(dt, k, procs, qfields)
    ! only grow ice which is not autoconverted to
    ! snow, c.f. Harrington et al (1995)
    ! This assumes that mu_ice==0, so the fraction becomes
    ! P(mu+2, lambda*DImax) (see Abramowitz & Stegun 6.5.13)

    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    type(process_rate), intent(inout), target :: procs(:,:)
    real(wp), intent(in) :: qfields(:,:)

    type(process_rate), pointer :: ice_dep, snow_dep, ice_aut
    real(wp) :: lam, frac, dmass
    integer :: i_pqi, i_pqai, i_pqs, i_pns, i_pm3s
    real(wp) :: m1,m2,m3,dm1,dm2,dm3

    ice_dep=>procs(k, i_idep%id)
    snow_dep=>procs(k, i_sdep%id)
    ice_aut=>procs(k, i_saut%id)

    if (ice_aut%source(ice_params%i_1m) > 0) then
      lam=dist_lambda(k,ice_params%id)
      frac=1.0-exp(-lam*DImax)*(1.0+lam*DImax)
      dmass=frac*ice_dep%source(ice_params%i_1m)

      ice_dep%source(ice_params%i_1m)=ice_dep%source(ice_params%i_1m)-dmass
      snow_dep%source(snow_params%i_1m)=snow_dep%source(snow_params%i_1m)+dmass
    end if

    nullify(ice_dep)
    nullify(ice_aut)
    nullify(snow_dep)
  end subroutine adjust_dep
end module adjust_deposition
