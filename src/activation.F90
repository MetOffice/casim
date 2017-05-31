module activation
  use variable_precision, only: wp
  use mphys_constants, only: fixed_cloud_number
  use mphys_parameters, only: C1, K1
  use mphys_switches, only: iopt_act, aero_index,   &
       i_am2, i_an2, i_am4, i_am1, i_an1, i_am3, i_an3, i_am5, i_am6, i_an6, i_am7, l_warm
  use aerosol_routines, only: aerosol_active, aerosol_phys, aerosol_chem, &
       abdulRazzakGhan2000, invert_partial_moment, upperpartial_moment_logn, &
       invert_partial_moment_approx, invert_partial_moment_betterapprox
  use special, only: erfinv, erf, pi
  use thresholds, only: w_small, nl_tidy, ccn_tidy

  implicit none

  character(len=*), parameter, private :: ModuleName='ACTIVATION'

  private

  public activate
contains

  subroutine activate(dt, cloud_mass, cloud_number, w, rho, dnumber, dmac, T, p, &
       aerophys, aerochem, aeroact, dustphys, dustchem, dustliq, dnccn_all, dmac_all, &
       dnumber_d, dmass_d)

    implicit none

    character(len=*), parameter :: RoutineName='ACTIVATE'

    real(wp), intent(in) :: dt
    real(wp), intent(in) :: cloud_mass, cloud_number, w, rho, T, p
    real(wp), intent(out) ::  dnumber, dmac
    type(aerosol_phys), intent(in) :: aerophys
    type(aerosol_chem), intent(in) :: aerochem
    type(aerosol_active), intent(in) :: aeroact
    type(aerosol_phys), intent(in) :: dustphys
    type(aerosol_chem), intent(in) :: dustchem
    type(aerosol_active), intent(in) :: dustliq
    real(wp), intent(out) :: dnccn_all(:),dmac_all(:)
    real(wp), intent(out) :: dnumber_d, dmass_d ! activated dust number and mass
    real(wp) :: active, nccn, Smax, rcrit, nccn_active
    real(wp) :: Nd, rm, sigma, density

    real :: work, dmac_i, diff
    integer :: imode
    logical :: l_useactive

    l_useactive=.false.
    dmac=0.0
    dnumber_d=0.0
    dmass_d=0.0
    dnumber=0.0
    dmac=0.0
    dmac_all=0.0
    dnccn_all=0.0

    select case(iopt_act)
    case default
      ! fixed number
      active=fixed_cloud_number
    case(1)
      ! activate 100% aerosol
      active=sum(aerophys%N(:))
    case(2)
      ! simple Twomey law Cs^k expressed as
      ! a function of w (Rogers and Yau 1989)
      active=0.88*C1**(2.0/(K1+2.0))*(7.0E-2*(w*100.0)**1.5)**(K1/(K1+2.0))*1.0e6/rho
    case(3)
      ! Use scheme of Abdul-Razzak and Ghan
      if (w > w_small .and. sum(aerophys%N(:)) > ccn_tidy) then
        call AbdulRazzakGhan2000(w, p, T, aerophys, aerochem, dnccn_all, Smax, aeroact, &
             nccn_active, l_useactive)
        active=sum(dnccn_all(:))

        if (Smax > 0.02 .and. .not. l_warm) then
          ! need better model than this. Could do partitioning in the same way as CCN
          dnumber_d=.01*dustphys%N(1)
          dmass_d=dnumber_d*dustphys%M(1)/dustphys%N(1)
        end if
      else
        active=0.0
      end if
    end select

    select case(iopt_act)
    case default
      ! fixed number, so no need to calculate aerosol changes
      dnumber=max(0.0_wp,(active-cloud_number))
    case (1:4)
      if (active > nl_tidy) then
        if (l_useactive) then
          dnumber=active
        else
          dnumber=max(0.0_wp,(active-cloud_number))
          ! Rescale to ensure total removal of aerosol number=creation of cloud number
          dnccn_all = dnccn_all*(dnumber/(sum(dnccn_all) + tiny(dnumber)))
        end if
        ! Need to make this consistent with all aerosol_options
        do imode = 1, aero_index%nccn
          Nd=aerophys%N(imode)
          if (Nd > ccn_tidy) then
            rm=aerophys%rd(imode)
            sigma=aerophys%sigma(imode)
            density=aerochem%density(imode)
            
            rcrit=invert_partial_moment_betterapprox(dnccn_all(imode), 0.0_wp, Nd, rm, sigma)
            
            dmac_all(imode)=(4.0*pi*density/3.0)*(upperpartial_moment_logn(Nd, rm, sigma, 3.0_wp, rcrit))
            dmac_all(imode)=min(dmac_all(imode),0.999*aerophys%M(imode)) ! Don't remove more than 99.9%
            dmac=dmac+dmac_all(imode)
          end if
        end do
      end if
    end select

    ! Convert to rates rather than increments
    dmac=dmac/dt
    dmac_all=dmac_all/dt
    dnccn_all=dnccn_all/dt
    dnumber=dnumber/dt
    dnumber_d=dnumber_d/dt
    dmass_d=dmass_d/dt

  end subroutine activate
end module activation
