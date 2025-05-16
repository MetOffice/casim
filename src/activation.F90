module activation
  use variable_precision, only: wp
  use mphys_constants, only: fixed_cloud_number
  use mphys_parameters, only: C1, K1, cloud_params, zero_real_wp
  use mphys_switches, only: iopt_act, iopt_inuc, aero_index,   &
       l_warm,              &
       iopt_shipway_act,l_ukca_casim,      &
       activate_in_cloud
  use aerosol_routines, only: aerosol_active, aerosol_phys, aerosol_chem, &
       abdulRazzakGhan2000, upperpartial_moment_logn, &
       invert_partial_moment_betterapprox, &
       AbdulRazzakGhan2000_dust
  use special, only: pi,GammaFunc
  use thresholds, only: w_small, nl_tidy, ni_tidy, ccn_tidy, ql_small
  Use shipway_parameters, only: max_nmodes, nmodes, Ndi, &
     rdi, sigmad, bi, betai, use_mode, nd_min
  use shipway_constants, only: Mw, rhow, eps, Rd, Dv, Lv, cp, &
      Dv_mean, alpha_c, zetasa, Ru
  use qsat_funs, only: qsaturation,dqwsatdt
  use shipway_activation_mod, only: solve_nccn_household, solve_nccn_brent
  use casim_stph, only: l_rp2_casim, fixed_cloud_number_rp

  implicit none

  character(len=*), parameter, private :: ModuleName='ACTIVATION'

  private

  public activate

  ! Variables used in the call to the Shipway (2015) activation scheme

  real(wp) :: ent_fraction=0.
  integer  :: order=2
  integer  :: niter=8
  real(wp) :: smax0=0.001
  !real(wp) :: alpha_c=0.05 !kinetic parameter
  real(wp) :: nccni(max_nmodes) ! maximum number of modes that can be used (usually=3)

contains

  subroutine activate(dt, cloud_mass, cloud_number, w, rho, dnumber, dmac, T, p,  &
       cfliq,cfliq_old, aerophys, aerochem, aeroact, dustphys, dustchem,   &
       dustliq, dnccn_all, dmac_all, dnumber_d, dmass_d, dnccnd_all,dmad_all,  &
       smax,ait_cdnc,accum_cdnc, tot_cdnc,activated_arg,activated_cloud)
!PRF NB extra cfliq argument missing from activate call in condensation for non-um


    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    use thresholds, only: cfliq_small
    implicit none

    ! Subroutine arguments

    real(wp), intent(in) :: dt
    real(wp), intent(in) :: cloud_mass, cloud_number, w, rho, T, p, &
         cfliq, cfliq_old
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
    real(wp), intent(out) :: smax,ait_cdnc, accum_cdnc,tot_cdnc
    real(wp), intent(out) :: activated_cloud, activated_arg
    ! Local Variables

    real(wp) :: cloud_number_work, cloud_number_work_old, cloud_mass_work
    real(wp) :: cloud_radius_work, smax_cloud, smax_act
    real(wp) :: active, rcrit, nccn_active, dactive,nccn_dactive
    real(wp) :: Nd, rm, sigma, density
    real(wp) :: dnccn_cloud(aero_index%nccn)
    real(wp) :: cf_liquid, cf_thresh, cf_liquid_old

    integer :: imode
    logical :: l_useactive

    real(wp) :: LvT, alpha, lam, tau, qs, &
         Dv_here, erfarg
    real(wp) :: Ak, bigGthermal, bigGdiffusion, gammaL,     &
         gammaR, gammastar, bigG
    real(wp) :: s0i(max_nmodes), sigmas(max_nmodes)
    real(wp) :: m1, m2, j1
    real(wp) :: kwdqsdz, dqsdt
    real(wp), parameter :: smax_act_min = 0.02

    integer, parameter :: solve_household = 1
    integer, parameter :: solve_brent = 2

    real(wp) :: dv_flag=0

    integer :: solve_select

    character(len=*), parameter :: RoutineName='ACTIVATE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Apply RP scheme
    if ( l_rp2_casim ) then
        fixed_cloud_number = fixed_cloud_number_rp
    endif

    cf_thresh = cfliq_small ! should be a small number, shouldn't matter too much
    if (cfliq .gt. cfliq_small) then  !only doing liquid cloud fraction at the moment
      cf_liquid=cfliq
    else
      cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
    endif
    if (cfliq_old .gt. cfliq_small) then  !only doing liquid cloud fraction at the moment
      cf_liquid_old=cfliq_old
    else
      cf_liquid_old=cfliq_small !nonzero value - maybe move cf test higher up
    endif

!make in-cloud number conc
! Should this be the new cloud fraction or the old cloud fraction?
! I think it should be the cloud fraction after advection but before
! the new condensation/cloud scheme call, to be consistent with the CDNC used
! here.

    cloud_number_work_old = cloud_number / cf_liquid_old
    cloud_number_work = cloud_number / cf_liquid
    ! threshold cloud number conc, so that cloud number
    ! is larger 1 droplet per m^3. This ensures stability of the
    ! tau calculation.
    if (cloud_number .lt. 1.0e-6) then
       cloud_number_work_old = 1.0e-6/cf_liquid_old
       cloud_number_work = 1.0e-6/cf_liquid
    end if

! This is cloud_mass_preap2 which is the OLD cloud mass.
    cloud_mass_work = cloud_mass / cf_liquid_old
    l_useactive=.false.
    dmac=0.0
    dnumber_d=0.0
    dmass_d=0.0
    dnumber=0.0
    dmac_all=0.0
    dnccn_all=0.0
    dmad_all=0.0
    dnccnd_all=0.0
    dnccn_cloud(:)=0.0
    smax = 0.0
    ait_cdnc=0.0
    accum_cdnc=0.0
    tot_cdnc=0.0
    activated_arg=0.0
    activated_cloud=0.0
    solve_select = solve_brent
    gammastar=0.0
    bigG=0
    if (int(dv_flag)==1)then
       Dv_here=Dv
    else if(dv_flag <= 0)then
       Dv_here=Dv_mean(T,1.0_wp)
    else
      Dv_here=dv_flag
    end if
    m2=cloud_number_work_old
    m1=cloud_mass_work/cloud_params%c_x
    j1=1.0/(cloud_params%p1-cloud_params%p2)
    qs=(0.029/0.018)*qsaturation(T, p/100.0)
    dqsdt = dqwsatdt (qs,T)

    LvT = Lv -(4217.4-1870.0)*(T-273.15)
    if (cloud_mass_work > ql_small .and. cf_liquid > cf_thresh) then
      ! The following calculations are based on Gordon et al, 2020,
      ! "Development of aerosol activation in the double-moment
      ! Unified Model and evaluation with CLARIFY measurements",
      ! Atmos. Chem. Phys., 20, 10997-11024,
      ! https://doi.org/10.5194/acp-20-10997-2020, 2020
      lam=((GammaFunc(1.0+cloud_params%fix_mu+cloud_params%p1) / &
           GammaFunc(1.0+cloud_params%fix_mu+cloud_params%p2))*(m2/m1))**(j1)
      cloud_radius_work = 0.5*GammaFunc(2.0+cloud_params%fix_mu)/ &
           (lam*GammaFunc(1.0+cloud_params%fix_mu))
      bigGdiffusion = rhow*Ru*T/(p*qs*Dv_here*0.018)
      bigGthermal = (LvT*rhow/(0.024*T))*(LvT*0.018/(Ru*T) -1)
      bigG = 1.0/(bigGdiffusion+bigGthermal)
      gammaR= 0.018*LvT*LvT/(cp*Ru*T*T)
      gammaL = 0.029/(qs*0.018)
      gammastar = 4*pi*rhow*(gammaL+gammaR)/rho
      tau = 1.0/(gammastar*bigG*cloud_number_work_old*rho*cloud_radius_work)
    else
      tau =1000.0  ! This ensures activation will always happen. Activation is
      ! only called when water is condensing
    end if
    if (tau >= 1000.0 .or. tau <= 0.0) then
      tau=1000.0
    end if

    ! This agrees well with the version in Ghan et al
    alpha = 9.8*(LvT/(eps*cp*T)-1.0)/(T*Rd)*(1-ent_fraction)

    smax_cloud= alpha*w*tau

    kwdqsdz = w*5.3e5*15*dqsdt*0.006 ! constant*w*time-threshold*dqsat/dT*dT/dz

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
      ! setup Shipway parameters to ensure calc_nccn works correctly
      nmodes=aero_index%nccn
      Ak = 2.*Mw*zetasa/(Ru*T*rhow)

      do imode=1,aero_index%nccn
        if(l_ukca_casim) then
          bi(imode) =aerochem%bk(imode)*aerochem%epsv(imode)
        else
          bi(imode) =aerochem%vantHoff(imode)*aerochem%epsv(imode)* &
             aerochem%density(imode)*Mw/(rhow*aerochem%massMole(imode))
        endif
        Ndi(imode)=aerophys%N(imode)
        rdi(imode)=aerophys%rd(imode)
        sigmad(imode)=aerophys%sigma(imode)
        betai(imode)=0.5      !aerochem%beta(imode) This is set to 0.5 for
                              !    Shipway not for ARG
        !print *, 'imode, bi', bi(imode), imode, Ak,betai(imode),rdi(imode)
        if (rdi(imode) > epsilon(1.0_wp)) then
           s0i(imode) = rdi(imode)**(-(1.0+betai(imode))) * &
                sqrt(4.0*Ak**3.0/(27.0*bi(imode)))
        else
           s0i(imode) = 0.0
        endif
        sigmas(imode) = sigmad(imode)**(1.+betai(imode))
        ! only use the mode if there's significant number
        use_mode(imode) = aerophys%N(imode) > Nd_min
      end do ! loop over modes

      if (w > w_small .and. sum(aerophys%N(:)) > ccn_tidy) then
         ! The following derivations of nccn is based on an
         ! adaptation of Abdul-Razzak and Ghan scheme, which is
         ! described in Gordon et al, 2020,
         ! "Development of aerosol activation in the double-moment
         ! Unified Model and evaluation with CLARIFY measurements",
         ! Atmos. Chem. Phys., 20, 10997-11024,
         ! https://doi.org/10.5194/acp-20-10997-2020, 2020
        call AbdulRazzakGhan2000(w, p, T, aerophys, aerochem, dnccn_all, Smax_act, aeroact, &
               nccn_active, l_useactive)

        activated_arg = sum(dnccn_all(:))

        if (Smax_act > smax_act_min .and. .not. l_warm) then

          if (iopt_inuc < 4) then
            ! For lower-order ice nucleation options, need to initialise
            ! dactive and dmass_d
            dactive = zero_real_wp
            dmass_d = zero_real_wp
          else
            ! For higher-order ice nucleation schemes, dactive and dmass_d
            ! can be based on dustphys
            dactive = 0.01*dustphys%N(1)

            if ( dustphys%N(1) > ni_tidy ) then
              dmass_d = dactive * dustphys%M(1) / dustphys%N(1)
            else
              ! Prevent divide by zero generating nonsense.
              dmass_d = zero_real_wp
            end if ! dustphys%N(1) > ni_tidy

          end if ! iopt_inuc

        end if ! Smax_act > 0.02 and .not. l_warm

        if(activate_in_cloud ==2 .and. cf_liquid > cf_thresh .and. &
            smax_cloud <= smax_act) then
          ! In this case use a weighted sum of ARG-activation outside cloud and equilibrium inside cloud
          ! should really use the equilibrium*old-cloud-frac+arg*(new-cloud-frac-old-cloud-frac)
          activated_cloud=0
          dnccn_cloud(:)=0
          do imode=1,aero_index%nccn
             if (use_mode(imode)) then
                  ! In ARG2000, error_func=1.0-erf(2.0*log(s_cr(i)/smax)/(3.0*sqrt(2.0)*log(phys%sigma(i))))
                  ! so insert the factor of 2/3 for consistency
                erfarg=2.0*log(smax_cloud/s0i(imode))/(3.0*sqrt(2.)*log(sigmas(imode)))
                dnccn_cloud(imode) = 0.5*Ndi(imode)*(1.0+erf(erfarg))
                activated_cloud=activated_cloud+dnccn_cloud(imode)
             end if
          end do
          ! New CDNC is a sum of CDNC created in old and new cloud
          active = activated_cloud*cf_liquid_old+(cf_liquid-cf_liquid_old)*activated_arg
          active = active/cf_liquid ! This is just because we multiply by CF
          dnccn_all(:) = (dnccn_cloud(:)*cf_liquid_old+(cf_liquid-cf_liquid_old)*dnccn_all(:))/cf_liquid
          ! at the end of the subroutine, and we need to avoid multiplying by
          ! CF twice
        else if(activate_in_cloud==1 .or. cf_liquid <= cf_thresh        &
          .or. smax_cloud > smax_act ) then
          ! In this case use Abdul-Razzak & Ghan
          smax = smax_act
          active=sum(dnccn_all(:))
        else
          dnccn_all(:)=0.0
          active=0.0
        end if ! activate_in_cloud == 2 etc

      else
        active=0.0
        dactive=0.0
      end if

    case(4)
      ! Use scheme of Abdul-Razzak and Ghan (including for insoluble aerosol by
      ! assuming small amount of soluble material on it)
      if (w > w_small .and. (sum(aerophys%N(:))+sum(dustphys%N(:))) > ccn_tidy) then
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
        if(l_ukca_casim) then
          bi(imode) =aerochem%bk(imode)*aerochem%epsv(imode)
        else
          bi(imode) =aerochem%vantHoff(imode)*aerochem%epsv(imode)* &
             aerochem%density(imode)*Mw/(rhow*aerochem%massMole(imode))
        endif
        Ndi(imode)=aerophys%N(imode)
        rdi(imode)=aerophys%rd(imode)
        sigmad(imode)=aerophys%sigma(imode)
        betai(imode)=aerochem%beta(imode)
        ! only use the mode if there's significant number
        use_mode(imode) = aerophys%N(imode) > Nd_min
      end do

      if (any(use_mode) .and. cloud_number_work < sum(Ndi))then
        if (tau > kwdqsdz .or. activate_in_cloud==1) then
          select case (solve_select)
            case (solve_household)
              call solve_nccn_household( order, niter, smax0, w, T, p, alpha_c, &
                                         ent_fraction, smax, active, nccni     )
            case (solve_brent)
              call solve_nccn_brent(w, T, p, alpha_c, ent_fraction,      &
                                  smax, active, nccni)
          end select
        else
           LvT = Lv -(4217.4-1870.0)*(T-273.15)
           alpha = 9.8*(LvT/(eps*cp*T)-1.0)/(T*Rd)*(1-ent_fraction)
           smax= alpha*w*tau
           !call calc_nccn(smax,active,nccni)
        end if

        dnccn_all(1:aero_index%nccn) = nccni(1:aero_index%nccn)
      else
        active=0.0
      end if
    end select

    if(activate_in_cloud==2) then
      smax = smax_act
      ait_cdnc = tau
      accum_cdnc = smax_cloud
    else if(activate_in_cloud==0) then
      ait_cdnc = tau
      accum_cdnc = dnccn_all(2)
    else
      ait_cdnc = aerophys%N(1)
      accum_cdnc = aerophys%N(2)
    end if
    tot_cdnc = active

    select case(iopt_act)
    case default
      ! fixed number, so no need to calculate aerosol changes
      dnumber=max(0.0_wp,(active-cloud_number_work))
    case (1:5)
      if (active > nl_tidy) then
        if (l_useactive) then
          dnumber=active
          dnumber_d=dactive
        else
          dnumber=max(0.0_wp,(active+dactive-cloud_number_work))
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
          dmass_d=zero_real_wp
          dmad_all(:)=zero_real_wp
          do imode = 1, aero_index%nin
            Nd=dustphys%N(imode)
            if (Nd > ccn_tidy .and. dnumber_d > ccn_tidy) then
              rm=dustphys%rd(imode)
              sigma=dustphys%sigma(imode)
              density=dustchem%density(imode)

              rcrit=invert_partial_moment_betterapprox(dnccnd_all(imode), 0.0_wp, Nd, rm, sigma)

              dmad_all(imode)=(4.0*pi*density/3.0)*(upperpartial_moment_logn(Nd, rm, sigma, 3.0_wp, rcrit))
              dmad_all(imode)=min(dmad_all(imode),0.999*dustphys%M(imode)) ! Don't remove more than 99.9%
              dmass_d=dmass_d+dmad_all(imode)
            end if
          end do
        end if
      end if
    end select

    ! Convert to rates rather than increments .... and back to gridbox means
    dmac=dmac/dt * cf_liquid
    dmac_all=dmac_all/dt * cf_liquid
    dnccn_all=dnccn_all/dt * cf_liquid
    dnumber=dnumber/dt * cf_liquid
    dmass_d = dmass_d/dt * cf_liquid
    dmad_all = dmad_all/dt * cf_liquid
    dnccnd_all=dnccnd_all/dt * cf_liquid
    dnumber_d=dnumber_d/dt * cf_liquid

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine activate
end module activation
