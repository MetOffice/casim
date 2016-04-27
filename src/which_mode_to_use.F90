! Routine decides which modes to put back
! re-evaporated (or potentially any other) aerosol
module which_mode_to_use
  use variable_precision, only: wp
#if DEF_MODEL==MODEL_UM
  use ukca_option_mod, only: l_ukca
  implicit none
#else
  implicit none
  logical :: l_ukca = .false.
#endif
  private

  integer :: imethod = 1 ! method to use
  integer, parameter :: iukca_method = 2 ! method to use if UKCA is used

  real(wp) :: max_accumulation_mean_radius = 0.25e-6
  real(wp) :: min_coarse_mean_radius = 1.0e-6

  public which_mode
contains

  ! Subroutine calculates how much of the total increments
  ! to mass and number should be sent to each of two modes.
  subroutine which_mode(dm, dn, r1_in, r2_in, density, dm1, dm2, dn1, dn2)
    real(wp), intent(in) :: dm  ! total mass increment
    real(wp), intent(in) :: dn  ! total number increment

    real(wp), intent(in) :: r1_in  ! mean radius of mode 1
    real(wp), intent(in) :: r2_in  ! mean radius of mode 2

    real(wp), intent(in) :: density ! density of aerosol, assumed the same across all modes

    real(wp), intent(out) :: dm1  ! mass increment to mode 1
    real(wp), intent(out) :: dn1  ! number increment to mode 1
    real(wp), intent(out) :: dm2  ! mass increment to mode 2
    real(wp), intent(out) :: dn2  ! number increment to mode 2

    ! local variables
    real(wp) :: rm    ! mean radius of increment
    real(wp) :: gamma ! convenience variable

    real(wp) :: r1, r2     ! r1 and r2 (possibly modified) 
    real(wp) :: r1_3, r2_3 ! r1**3 and r2**3

    real(wp) :: FTPI ! 4/3*pi
    real(wp) :: RFTPI ! 1./(4/3*pi)


    dm1=0.0
    dm2=0.0
    dn1=0.0
    dn2=0.0

    if (l_ukca) imethod=iukca_method
    if (dm*dn > 0.0) then 
      ! dm and dn should be positive and of the same sign
      r1=min(max_accumulation_mean_radius, r1)
      r2=max(min_coarse_mean_radius, r2)

      ftpi=3.14159*4./3.
      rftpi=1./ftpi
      rm=(rftpi*dm/dn/density)**(1.0/3.0)

      select case(imethod)
      case (iukca_method)
        ! Call ukca subroutine, or add UKCA code here?
      case default
        if (rm > r2) then
          dm1=0.0
          dn1=0.0
          dm2=dm
          dn2=dn
        else if (rm < r1) then
          dm1=dm
          dn1=dn
          dm2=0.0
          dn2=0.0
        else
          r1_3=r1*r1*r1
          r2_3=r2*r2*r2
          gamma=(rm*rm*rm-r1_3)/(r2_3-r1_3)
          dn2=dn*(gamma)
          dn1=dn - dn2
          dm2=FTPI*density*r2_3*dn2
          dm1=dm-dm2
        end if
      end select
    end if
  end subroutine which_mode
end module which_mode_to_use
