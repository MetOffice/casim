module breakup
  use variable_precision, only: wp, iwp
  use process_routines, only: process_rate, process_name, i_sbrk
  use mphys_parameters, only: hydro_params, DSbrk, tau_sbrk
  use thresholds, only: qr_small, thresh_small
  use m3_incs, only: m3_inc_type2
  use distributions, only: dist_lambda, dist_mu, dist_n0

#if DEF_MODEL==MODEL_KiD
  use diagnostics, only: save_dg, i_dgtime
#endif

  implicit none
  private

  public ice_breakup
contains
  !< Subroutine to determine the breakup of large particles
  !< This code is specified for just snow, but could be used for
  !< other species (e.g. rain)
  !< For triple moment species there is a corresponding change in the
  !< 3rd moment assuming shape parameter is not changed
  !< NB: Aerosol mass is not modified by this process
  !
  !< OPTIMISATION POSSIBILITIES: strip out shape parameters
  subroutine ice_breakup(dt, k, params, qfields, procs)
    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    type(hydro_params), intent(in) :: params
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    type(process_name) :: iproc ! processes selected depending on which species we're modifying
    real(wp) :: dnumber, dm1, dm2, dm3
    real(wp) :: number, mass, m1, m2, m3
    type(process_rate), pointer :: this_proc
    real(wp) :: n0, lam, mu
    real(wp) :: Dm ! Mass-weighted mean diameter

    select case (params%id)
    case (4_iwp) !snow
      iproc=i_sbrk
    end select

    mass=qfields(k, params%i_1m)

    if (mass > thresh_small(params%i_1m) .and. params%l_2m) then ! if no existing ice, we don't bother
      this_proc=>procs(k, iproc%id)
      number=qfields(k, params%i_2m)
      if (params%l_3m) m3=qfields(k, params%i_3m)
      n0=dist_n0(k,params%id)
      mu=dist_mu(k,params%id)
      lam=dist_lambda(k,params%id)
      Dm=(1.0 + params%d_x + mu)/lam

      if (Dm > DSbrk) then ! Mean size exceeds threshold
        dnumber=(Dm/DSbrk - 1.0)**params%d_x * number / tau_sbrk
        this_proc%source(params%i_2m)=dnumber
        nullify(this_proc)
      end if
    end if
  end subroutine ice_breakup
end module breakup
