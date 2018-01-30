module activation
  use variable_precision, only: wp
  use mphys_constants, only: fixed_cloud_number
  use mphys_parameters, only: C1, K1
  use mphys_switches, only: iopt_act, aero_index,   &
       i_am2, i_an2, i_am4, i_am1, i_an1, i_am3, i_an3, &
       i_am5, i_am6, i_an6, i_am7, l_warm, iopt_shipway_act
  use aerosol_routines, only: aerosol_active, aerosol_phys, aerosol_chem, &
       abdulRazzakGhan2000, invert_partial_moment, upperpartial_moment_logn, &
       invert_partial_moment_approx, invert_partial_moment_betterapprox, &
       AbdulRazzakGhan2000_dust
  use special, only: erfinv, pi
  use thresholds, only: w_small, nl_tidy, ccn_tidy
  Use shipway_parameters, only: max_nmodes, nmodes, Ndi, &
     rdi, sigmad, bi, betai, use_mode, nd_min
  use shipway_constants, only: Mw, rhow
  use shipway_activation_mod, only: solve_nccn_household, solve_nccn_brent

  implicit none

  character(len=*), parameter, private :: ModuleName='ACTIVATION'

  private

  public activate

  ! Variables used in the call to the Shipway (2015) activation scheme
  
  real(wp) :: ent_fraction=0.
  integer  :: order=2
  integer  :: niter=8
  real(wp) :: smax0=0.001
  real(wp) :: alpha_c=0.05 !kinetic parameter
  real(wp) :: nccni(max_nmodes) ! maximumg number of modes that can be used (usually=3)

contains

  subroutine activate(dt, cloud_mass, cloud_number, w, rho, dnumber, dmac, T, p, &
       aerophys, aerochem, aeroact, dustphys, dustchem, dustliq, dnccn_all, dmac_all, &
       dnumber_d, dmass_d, dnccnd_all, dmad_all)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments

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
    real(wp), intent(out) :: dnccnd_all(:),dmad_all(:)
    real(wp), intent(out) :: dnumber_d, dmass_d ! activated dust number and mass

    ! Local Variables

    real(wp) :: active, nccn, Smax, rcrit, nccn_active, dactive,nccn_dactive
    real(wp) :: Nd, rm, sigma, density

    real :: work, dmac_i, diff
    integer :: imode
    logical :: l_useactive

    integer, parameter :: solve_household = 1
    integer, parameter :: solve_brent = 2
    
    integer :: solve_select

    character(len=*), parameter :: RoutineName='ACTIVATE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_useactive=.false.
    dmac=0.0
    dnumber_d=0.0
    dmass_d=0.0
    dnumber=0.0
    dmac=0.0
    dmac_all=0.0
    dnccn_all=0.0
    solve_select = solve_brent

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
          dactive = 0.01*dustphys%N(1)
          !dnumber_d=.01*dustphys%N(1)
          dmass_d=dactive*dustphys%M(1)/dustphys%N(1)
        end if
      else
        active=0.0
        dactive=0.0
      end if
    case(4)
      ! Use scheme of Abdul-Razzak and Ghan (including for insoluble aerosol by
      ! assuming small amount of soluble material on it)
      if (w > w_small .and. sum(aerophys%N(:)) > ccn_tidy) then
        if (l_warm) then
          call AbdulRazzakGhan2000(w, p, T, aerophys, aerochem, dnccn_all, Smax, aeroact, &
             nccn_active, l_useactive)
          active=sum(dnccn_all(:))
        end if
        if (.not. l_warm) then
           call AbdulRazzakGhan2000_dust(w, p, T, aerophys, aerochem, dnccn_all, Smax, aeroact, &
                nccn_active, nccn_dactive, dustphys, dustchem, dustliq, dnccnd_all, l_useactive)
           active   = sum(dnccn_all(:))
           dactive  = dnccnd_all(aero_index%i_coarse_dust) !SUM(dnccnd_all(:)) i_accum_dust currently not used!
        end if
      else
        active=0.0
        dactive=0.0
      end if

    case(iopt_shipway_act)
      ! Use scheme of Shipway 2015
      ! This is a bit clunk and could be harmonized
      nmodes=aero_index%nccn
      do imode=1,aero_index%nccn
        bi(imode) =aerochem%vantHoff(imode)*aerochem%epsv(imode)* &
           aerochem%density(imode)*Mw/(rhow*aerochem%massMole(imode)) 
        Ndi(imode)=aerophys%N(imode)
        rdi(imode)=aerophys%rd(imode)
        sigmad(imode)=aerophys%sigma(imode)
        betai(imode)=aerochem%beta(imode)
        ! only use the mode if there's significant number
        use_mode(imode) = aerophys%N(imode) > Nd_min 
      end do
      
      if (any(use_mode) .and. cloud_number < sum(Ndi))then

        select case (solve_select)
          case (solve_household)
            call solve_nccn_household( order, niter, smax0, w, T, p, alpha_c, &
                                       ent_fraction, smax, active, nccni     )
          case (solve_brent)
            call solve_nccn_brent(smax0, w, T, p, alpha_c, ent_fraction,      &
                                  smax, active, nccni)
        end select

        dnccn_all(1:aero_index%nccn) = nccni(1:aero_index%nccn)
      else
        active=0.0
      end if
    end select

    select case(iopt_act)
    case default
      ! fixed number, so no need to calculate aerosol changes
      dnumber=max(0.0_wp,(active-cloud_number))
    case (1:5)
      if (active > nl_tidy) then
        if (l_useactive) then
          dnumber=active
          dnumber_d=dactive
        else
          dnumber=max(0.0_wp,(active+dactive-cloud_number))
          ! Rescale to ensure total removal of aerosol number=creation of cloud number
          dnccn_all = dnccn_all*(dnumber/(sum(dnccn_all) + sum(dnccnd_all) + tiny(dnumber)))
          dnccnd_all = dnccnd_all*(dnumber/(sum(dnccn_all) + sum(dnccnd_all) + tiny(dnumber)))
          dnumber = sum(dnccn_all)
          dnumber_d = sum(dnccnd_all)
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
        ! for the insoluble mode
        if (iopt_act == 4) then
          do imode = 1, aero_index%nin
            Nd=dustphys%N(imode)
            if (Nd > ccn_tidy) then
              rm=dustphys%rd(imode)
              sigma=dustphys%sigma(imode)
              density=dustchem%density(imode)

              rcrit=invert_partial_moment_betterapprox(dnccnd_all(imode), 0.0_wp, Nd, rm, sigma)

              dmad_all(imode)=(4.0*pi*density/3.0)*(upperpartial_moment_logn(Nd, rm, sigma, 3.0_wp, rcrit))
              dmad_all(imode)=min(dmac_all(imode),0.999*dustphys%M(imode)) ! Don't remove more than 99.9%
              dmass_d=dmass_d+dmad_all(imode)
            end if
          end do
        end if
      end if
    end select

    ! Convert to rates rather than increments
    dmac=dmac/dt
    dmac_all=dmac_all/dt
    dnccn_all=dnccn_all/dt
    dnumber=dnumber/dt
    dmass_d = dmass_d/dt
    dmad_all = dmad_all/dt
    dnccnd_all=dnccnd_all/dt
    dnumber_d=dnumber_d/dt

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine activate
end module activation
