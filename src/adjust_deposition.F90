module adjust_deposition
  ! As in Harrington et al. 1995 move some of the depositional
  ! growth on ice into the snow category

  use variable_precision, only: wp
  use mphys_parameters, only: DImax, snow_params, ice_params
  use distributions, only: dist_lambda
  use process_routines, only: process_rate, i_idep, i_sdep, i_saut

  implicit none

  character(len=*), parameter, private :: ModuleName='ADJUST_DEPOSITION'

  private

  public adjust_dep
contains

  subroutine adjust_dep(nz, l_Tcold, procs)
    ! only grow ice which is not autoconverted to
    ! snow, c.f. Harrington et al (1995)
    ! This assumes that mu_ice==0, so the fraction becomes
    ! P(mu+2, lambda*DImax) (see Abramowitz & Stegun 6.5.13)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments

    integer, intent(in) :: nz
    logical, intent(in) :: l_Tcold(:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! local variables

    !type(process_rate), pointer :: ice_dep, snow_dep, ice_aut
    real(wp) :: lam, frac, dmass

    integer :: k

    character(len=*), parameter :: RoutineName='ADJUST_DEP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k = 1, nz
       if (l_Tcold(k)) then
          if (procs(ice_params%i_1m, i_saut%id)%column_data(k) > 0) then
             lam=dist_lambda(k,ice_params%id)
             frac=1.0-exp(-lam*DImax)*(1.0+lam*DImax)
             dmass=frac*procs(ice_params%i_1m, i_idep%id)%column_data(k)

             procs(ice_params%i_1m, i_idep%id)%column_data(k)= &
                  procs(ice_params%i_1m, i_idep%id)%column_data(k)-dmass
             procs(snow_params%i_1m,i_sdep%id)%column_data(k)= &
                  procs(snow_params%i_1m,i_sdep%id)%column_data(k)+dmass
          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine adjust_dep
end module adjust_deposition
