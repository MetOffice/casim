module breakup
  use variable_precision, only: wp, iwp
  use process_routines, only: process_rate, process_name, i_sbrk
  use mphys_parameters, only: hydro_params, DSbrk, tau_sbrk
  use thresholds, only: thresh_small
  use distributions, only: dist_lambda, dist_mu

  implicit none

  character(len=*), parameter, private :: ModuleName='BREAKUP'

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
  subroutine ice_breakup(nz, l_Tcold, params, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: nz
    logical, intent(in) :: l_Tcold(:)
    type(hydro_params), intent(in) :: params
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! Local variables
    type(process_name) :: iproc ! processes selected depending on which species we're modifying
    real(wp) :: dnumber
    real(wp) :: num, mass
    real(wp) :: lam, mu
    real(wp) :: Dm ! Mass-weighted mean diameter
    
    integer :: k
    
    character(len=*), parameter :: RoutineName='ICE_BREAKUP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k = 1, nz
       if (l_Tcold(k)) then

          select case (params%id)
          case (4_iwp) !snow
             iproc=i_sbrk
          end select

          mass=qfields(k, params%i_1m)

          if (mass > thresh_small(params%i_1m) .and. params%l_2m) then ! if no existing ice, we don't bother
             num=qfields(k, params%i_2m)
             mu=dist_mu(k,params%id)
             lam=dist_lambda(k,params%id)
             Dm=(1.0 + params%d_x + mu)/lam
             
             if (Dm > DSbrk) then ! Mean size exceeds threshold
                dnumber=(Dm/DSbrk - 1.0)**params%d_x * num / tau_sbrk
                procs(params%i_2m, iproc%id)%column_data(k)=dnumber
             end if
          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ice_breakup
end module breakup
