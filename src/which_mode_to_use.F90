! Routine decides which modes to put back
! re-evaporated (or potentially any other) aerosol
module which_mode_to_use
  use variable_precision, only: wp
  use mphys_switches, only: l_aeroproc_midway
  use lognormal_funcs, only: MNtoRm
  use mphys_constants, only: pi

  implicit none

  private

  character(len=*), parameter, private :: ModuleName='WHICH_MODE_TO_USE'

  integer, parameter :: imethod = 2 ! method to use
  integer, parameter :: iold_method = 1
  integer, parameter :: isimple_method = 2 ! simple method
  real(wp), parameter :: r_thresh_fixed = 0.5e-6 ![m] for simple method
  !set threshold to arithmetically half way between accum and coarse sizes

  real(wp) :: max_accumulation_mean_radius = 0.25e-6
  real(wp) :: min_coarse_mean_radius = 1.0e-6

  public which_mode
contains

  ! Subroutine calculates how much of the total increments
  ! to mass and number should be sent to each of two modes.
  subroutine which_mode(dm, dn, r1_in, r2_in, density, sigma, dm1, dm2, dn1, dn2)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='WHICH_MODE'

    real(wp), intent(in) :: dm  ! total mass increment
    real(wp), intent(in) :: dn  ! total number increment

    real(wp), intent(in) :: r1_in  ! mean radius of mode 1
    real(wp), intent(in) :: r2_in  ! mean radius of mode 2

    real(wp), intent(in) :: density ! density of aerosol, assumed the same 
                                    ! across all modes

    real(wp), intent(in) :: sigma ! sigma of the aerosol that is being 
                                  ! evaporated/added

    real(wp), intent(out) :: dm1  ! mass increment to mode 1
    real(wp), intent(out) :: dn1  ! number increment to mode 1
    real(wp), intent(out) :: dm2  ! mass increment to mode 2
    real(wp), intent(out) :: dn2  ! number increment to mode 2

    ! local variables
    real(wp) :: rm    ! mean radius of increment
    real(wp) :: gamma_var ! convenience variable

    real(wp) :: r1, r2     ! r1 and r2 (possibly modified) 
    real(wp) :: r1_3, r2_3 ! r1**3 and r2**3

    real(wp) :: ftpi ! 4/3*pi
    real(wp) :: rftpi ! 1./(4/3*pi)

    real(wp) :: r_thresh

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dm1=0.0
    dm2=0.0
    dn1=0.0
    dn2=0.0

    if (dm*dn > 0.0) then 
      ! dm and dn should be positive and of the same sign
      r1=min(max_accumulation_mean_radius, r1_in)
      r2=max(min_coarse_mean_radius, r2_in)

      if (l_aeroproc_midway) then
        r_thresh = (r1+r2)/2.0
      else
        r_thresh = r2
      end if

      ftpi = pi *4./3.
      rftpi=1./ftpi

      if (l_aeroproc_midway) then
        rm = MNtoRm(dm,dn,density,sigma) !DPG_bug_changes
      else
         ! Keep original code to preserve answers when switched off
        rm = (rftpi*dm/dn/density)**(1.0/3.0)
      end if

      select case(imethod)
      case default
        if (rm .ge. r_thresh) then
          dm1=0.0
          dn1=0.0
          dm2=dm
          dn2=dn
        else if (rm < r1 .OR. l_aeroproc_midway) then
          dm1=dm
          dn1=dn
          dm2=0.0
          dn2=0.0
        else
          r1_3=r1*r1*r1
          r2_3=r2*r2*r2
          gamma_var=(rm*rm*rm-r1_3)/(r2_3-r1_3)
          dn2=dn*(gamma_var)
          dn1=dn - dn2
          if (l_aeroproc_midway) then
            dm2=ftpi*density*r2_3*dn2*exp(4.5*log(sigma)**2) !DPG_bug_changes
          else
            dm2=FTPI*density*r2_3*dn2
          end if

          dm1=dm-dm2
        end if
      case(isimple_method)
        if (rm .ge. r_thresh_fixed) then
          dm1=0.0
          dn1=0.0
          dm2=dm
          dn2=dn
        else
          dm1=dm
          dn1=dn
          dm2=0.0
          dn2=0.0
        end if
      end select
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine which_mode
end module which_mode_to_use
