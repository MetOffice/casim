module ice_accretion
use variable_precision, only: wp, iwp
use mphys_die, only: throw_mphys_error, incorrect_opt, std_msg
use passive_fields, only: rho, TdegC, TdegK
use mphys_switches, only: l_process,   &
     cloud_params, rain_params,  ice_params, snow_params, &
     i_am4, i_am7, i_am8, i_am9, l_prf_cfrac, &
     l_g, l_raci_g, l_prf_cfrac, &
     i_cfi, i_cfs, i_cfg, i_cfl, i_cfr, mpof, l_srg
use type_process, only: process_name
use process_routines, only: process_rate,   &
     i_iacw, i_sacw, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs, &
     i_diacw, i_dsacw, i_dgacw, i_dsacr, i_dgacr, i_draci
use mphys_parameters, only: hydro_params
use thresholds, only: thresh_small, cfliq_small, thresh_sig 
use aerosol_routines, only: aerosol_active

use distributions, only: dist_lambda, dist_mu, dist_n0
use sweepout_rate, only: sweepout, binary_collection
!  use sweepout_rate, only: binary_collection_1M2M, sweepout_1M2M
  use casim_stph, only: l_rp2_casim, mpof_casim_rp

implicit none
private

character(len=*), parameter, private :: ModuleName='ICE_ACCRETION'

public iacc
contains

subroutine iacc(ixy_inner, dt, nz, l_Tcold, params_X, params_Y, params_Z, qfields, cffields, procs, &
                l_sigevap, aeroact, dustliq, aerosol_procs, params_snow)
    !
    !< CODE TIDYING: Move efficiencies into parameters

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: ixy_inner
real(wp), intent(in) :: dt
integer, intent(in) :: nz
logical, intent(in) :: l_Tcold(:) 
type(hydro_params), intent(in) :: params_X !< parameters for species which does the collecting
type(hydro_params), intent(in) :: params_Y !< parameters for species which is collected
type(hydro_params), intent(in) :: params_Z !< parameters for species to which resulting amalgamation is sent
type(hydro_params), intent(in), optional :: params_snow 
                                              !< snow parameters to which resulting 
                                              !< amalgamation may be sent when rain collecting snow, where 
                                              !< rain mass determines whether graupel (params_Z) or snow,
                                              !< and snow collecting rain T determines whether graupel 
                                              !< (params_Z) or snow.
real(wp), intent(in), target :: qfields(:,:)
real(wp), intent(in) :: cffields(:,:)
type(process_rate), intent(inout), target :: procs(:,:)

    ! aerosol fields
type(aerosol_active), intent(in) :: aeroact(:)
type(aerosol_active), intent(in) :: dustliq(:)

    ! optional aerosol fields to be processed
type(process_rate), intent(inout), target :: aerosol_procs(:,:)

logical, intent(in) :: l_sigevap(:) ! logical to determine significant evaporation

! Local variables
type(process_name) :: iproc, iaproc  ! processes selected depending on
! which species we're depositing on.

real(wp) :: dmass_Y, dnumber_Y, dmac, dmad
real(wp) :: dmass_X, dmass_Z, dnumber_X, dnumber_Z

real(wp) :: number_X, mass_X
real(wp) :: mass_Y, number_Y

real(wp) :: n0_X, lam_X, mu_X
real(wp) :: n0_Y, lam_Y, mu_Y

real(wp) :: cf_X, cf_Y, overlap_cf

real(wp) :: Eff ! collection efficiencies need to re-evaluate these and put them in properly to mphys_parameters

logical :: l_condition, l_alternate_Z
logical :: l_slow ! fall speed is slow compared to interactive species
logical :: l_aero ! If true then this process will modify aerosol

integer :: k

! local variables for the indexing of the amalgation hydrometeor of X and Y. These 
! are required because when rain collects ice or snow collects rain, the resulting 
! hydrometeor depends on rain mass and T, respectively
integer :: Zid ! id of the hydrometeor
integer :: Zi_1m ! index for mass of hydrometeor
logical :: l_Zi_2m ! logical to determine if double moment resulting hydrometeor
integer :: Zi_2m ! index for number of hydrometeor

character(len=*), parameter :: RoutineName='IACC'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise increments
dmass_X = 0.0_wp
dnumber_X = 0.0_wp

! Apply RP scheme
if ( l_rp2_casim ) then
    mpof = mpof_casim_rp
endif

!setup cloud fractions
do k = 1, nz
  if (l_Tcold(k)) then   
    if (l_prf_cfrac) then
      select case (params_X%id)
      case(1) 
        cf_X=max(cffields(k,i_cfl), cfliq_small)
      case(2)
        cf_X=max(cffields(k,i_cfr), cfliq_small)
      case(3) 
        cf_X=max(cffields(k,i_cfi), cfliq_small)
      case(4)
        cf_X=max(cffields(k,i_cfs), cfliq_small)
      case(5) 
        cf_X=max(cffields(k,i_cfg), cfliq_small)
      end select
      select case (params_Y%id)
      case(1) 
        cf_Y=max(cffields(k,i_cfl), cfliq_small)
      case(2)
        cf_Y=max(cffields(k,i_cfr), cfliq_small)
      case(3) 
        cf_Y=max(cffields(k,i_cfi), cfliq_small)
      case(4)
        cf_Y=max(cffields(k,i_cfs), cfliq_small)
      case(5) 
        cf_Y=max(cffields(k,i_cfg), cfliq_small)
      end select
    else
      cf_X=1.0
      cf_Y=1.0
    endif

    mass_Y=qfields(k, params_Y%i_1m) / cf_Y  !only do cloud species
    mass_X=qfields(k, params_X%i_1m) / cf_X  !only do cloud species

    !prf
    
    l_condition=mass_X * cf_X > thresh_small(params_X%i_1m) .and. & 
                mass_Y * cf_Y > thresh_small(params_Y%i_1m)
    
    ! Check for rain in collector or collected and then check for significant evap
    if ((params_X%id .eq. rain_params%id .or. params_Y%id .eq. rain_params%id) &
         .and. l_sigevap(k))  l_condition =.false.

    
    if (l_condition) then
      ! initialize variables which may not have been set
      number_X=0.0
      dnumber_Y=0.0

      Zid = params_Z%id 
      Zi_1m = params_Z%i_1m 
      if (params_Z%l_2m) then 
        l_Zi_2m = params_Z%l_2m
        Zi_2m = params_Z%i_2m
      endif

      ! determine resulting hydrometeor type from rain collecting ice
      !
      if (params_X%id .eq. rain_params%id .and.  params_Y%id .eq. ice_params%id) then
        if (mass_X > thresh_sig(params_X%i_1m) .and. l_g .and. l_raci_g) then
          ! resulting hydrometeor is graupel
          Zid = params_Z%id
          Zi_1m = params_Z%i_1m
          if (params_Z%l_2m) then
            l_Zi_2m = params_Z%l_2m
            Zi_2m = params_Z%i_2m
          endif
        else ! resulting hydrometeor is snow
          if (PRESENT(params_snow)) then
            Zid = params_snow%id
            Zi_1m = params_snow%i_1m
            if (params_snow%l_2m) then
              l_Zi_2m = params_snow%l_2m
              Zi_2m = params_snow%i_2m
            endif
          else
            write(std_msg, '(A)') 'Resultant of iacc, praci is snow but snow_params '//&
                         'opt arg not passed'
            call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                          std_msg)
          endif
        endif
      endif
      ! determine resulting hydrometeor type from snow collecting rain
      !
      if (params_X%id .eq. snow_params%id .and.  params_Y%id .eq. rain_params%id) then
        if (TdegK(k,ixy_inner) < 268.15 .and. l_g .and. l_srg) then
          ! resulting hydrometeor is graupel
          Zid = params_Z%id 
          Zi_1m = params_Z%i_1m 
          if (params_Z%l_2m) then 
            l_Zi_2m = params_Z%l_2m
            Zi_2m = params_Z%i_2m
          endif
        else ! resulting hydrometeor is snow 
          if (PRESENT(params_snow)) then 
            Zid = params_snow%id 
            Zi_1m = params_snow%i_1m 
            if (params_snow%l_2m) then 
               l_Zi_2m = params_snow%l_2m
               Zi_2m = params_snow%i_2m
            endif
          else
            write(std_msg, '(A)') 'Resultant of iacc, pracr is snow but snow_params '//&
                                  'opt arg not passed'
            call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, std_msg)
          endif
        endif
      endif

      l_alternate_Z=params_X%id /=  Zid ! params_Z%id
      l_aero=.false.

      select case (params_Y%id)
      case(1_iwp) ! cloud liquid is collected
        l_slow=.true. ! Y catagory falls slowly compared to X
        select case (params_X%id)
        case (3_iwp) !ice collects
          iproc=i_iacw
          Eff=0.0 !1.0  !PRF turn off ice collecting cloud drops - too small
          iaproc=i_diacw
          l_aero=.true.
        case (4_iwp) !snow collects
          iproc=i_sacw
          Eff=0.5 !make assumption that collecting area is approx half of circle - similar to operational.
          iaproc=i_dsacw
          l_aero=.true.
        case (5_iwp) !graupel collects
          iproc=i_gacw
          Eff=1.0
          iaproc=i_dgacw
          l_aero=.true.
        end select
      case(2_iwp) ! rain is collected
        l_slow=.false. ! Y catagory falls slowly compared to X
        select case (params_X%id)
        case (4_iwp) !snow collects
          iproc=i_sacr
          Eff=1.0
          iaproc=i_dsacr
          l_aero=.true.
        case (5_iwp) !graupel collects
          iproc=i_gacr
          Eff=1.0
          iaproc=i_dgacr
          l_aero=.true.
        end select
      case(3_iwp)! cloud ice is collected
        l_slow=.true. ! Y catagory falls slowly compared to X
        select case (params_X%id)
        case (2_iwp) !rain collects
          iproc=i_raci
          Eff=1.0          ! Rain collecting ice efficiency to make graupel =~0.5? Lew and Pruppacher 1983
          iaproc=i_draci
          l_aero=.true.
        case (4_iwp) !snow collects
          iproc=i_saci
          !Eff=1.0
          Eff = min(1.0_wp, 0.2*exp(0.08*TdegC(k,ixy_inner)))
        case (5_iwp) !graupel collects
          iproc=i_gaci
          !Eff=1.0
          Eff = min(1.0_wp, 0.2*exp(0.08*TdegC(k,ixy_inner)))
        end select
      case(4_iwp) ! snow is collected
        l_slow=.false. ! Y catagory falls slowly compared to X
        select case (params_X%id)
        case (5_iwp) !graupel collects
          iproc=i_gacs
          !Eff=1.0
          Eff = min(1.0_wp, 0.2*exp(0.08*TdegC(k,ixy_inner)))
        end select
      end select

      if (params_X%l_2m) number_X=qfields(k, params_X%i_2m)  / cf_X

      n0_X=dist_n0(k,params_X%id)
      mu_X=dist_mu(k,params_X%id)
      lam_X=dist_lambda(k,params_X%id)

      if (l_slow) then  ! collected species is approximated to have zero fallspeed
        !if (l_gamma_online) then
        dmass_Y=-Eff*sweepout(n0_X, lam_X, mu_X, params_X, rho(k,ixy_inner))*mass_Y
       !else
       !   dmass_Y=-Eff*sweepout_1M2M(n0_X, lam_X, params_X, rho(k))*mass_Y
       !endif
        dmass_Y=max(dmass_Y, -mass_Y/dt)

        if (params_Y%l_2m) then
          number_Y=qfields(k, params_Y%i_2m) / cf_Y
          dnumber_Y=dmass_Y*number_Y/mass_Y
        end if

        if (l_alternate_Z) then ! We move resulting collision to another species.
          if (params_Y%l_2m) then
            number_Y=qfields(k, params_Y%i_2m)  / cf_Y
          else
            ! we need an else in here
          end if

          !if (l_gamma_online) then
          dmass_X=-Eff*sweepout(n0_X, lam_X, mu_X, params_X, rho(k,ixy_inner), mass_weight=.true.) * number_Y
          !else
          !   dmass_X=-Eff*sweepout_1M2M(n0_X, lam_X, params_X, rho(k), mass_weight=.true.) * number_Y
          !endif

          dmass_X=max(dmass_X, -mass_X/dt)
          dnumber_X=dnumber_Y !dmass_Y*number_X/mass_Y
          dmass_Z=-1.0*(dmass_X + dmass_Y)
          dnumber_Z=-1.0*dnumber_X
        end if
      else  ! both species have significant fall velocity
        if (params_Y%l_2m) then
          number_Y=qfields(k, params_Y%i_2m) / cf_Y
        end if

        n0_Y=dist_n0(k,params_Y%id)
        mu_Y=dist_mu(k,params_Y%id)
        lam_Y=dist_lambda(k,params_Y%id)

        !if (l_gamma_online) then
        dmass_Y=-Eff*binary_collection(n0_X, lam_X, mu_X, n0_Y, lam_Y, mu_Y,&
                   params_X, params_Y, rho(k,ixy_inner), mass_weight=.true.)
        !else
        !   dmass_Y=-Eff*binary_collection_1M2M(n0_X, lam_X, n0_Y, lam_Y, params_X, params_Y, rho(k), mass_weight=.true.)
        !endif
        dmass_Y=max(dmass_Y, -mass_Y/dt)

        if (params_Y%l_2m) then
           !if (l_gamma_online) then
          dnumber_Y=-Eff*binary_collection(n0_X, lam_X, mu_X, n0_Y, lam_Y, &
                                mu_Y, params_X, params_Y, rho(k,ixy_inner))
           !else
           !   dnumber_Y=-Eff*binary_collection_1M2M(n0_X, lam_X, n0_Y, lam_Y, params_X, params_Y, rho(k))
           !endif
        end if

        if (l_alternate_Z) then ! We move resulting collision to another species.
           !if (l_gamma_online) then 
          dmass_X=-Eff*binary_collection(n0_Y, lam_Y, mu_Y, n0_X, lam_X,   &
            mu_X, params_Y, params_X, rho(k,ixy_inner), mass_weight=.true.)
         !else
         !   dmass_X=-Eff*binary_collection_1M2M(n0_Y, lam_Y, n0_X, lam_X, params_Y, params_X, rho(k), mass_weight=.true.)
         !endif
          dmass_X=max(dmass_X, -mass_X/dt)
          dnumber_X=dnumber_Y
          dmass_Z=-1.0*(dmass_X + dmass_Y)
          dnumber_Z=-1.0*dnumber_X
        end if
      end if

      ! PRAGMATIC HACK - FIX ME (ALTHOUGH I QUITE LIKE IT)
      ! If most of the collected particles are to be removed then remove all of them
      if (-dmass_Y*dt >0.95*mass_Y .or. (params_Y%l_2m .and. -dnumber_Y*dt > 0.95*number_Y)) then
        dmass_Y=-mass_Y/dt
        dnumber_Y=-number_Y/dt
      end if
      ! If most of the collecting particles are to be removed then remove all of them
      if (l_alternate_Z .and. (-dmass_X*dt >0.95*mass_X .or.                                          &
           (params_X%l_2m .and. -dnumber_X*dt > 0.95*number_X))) then
        dmass_X=-mass_X/dt
        dnumber_X=-number_X/dt
        dmass_Z = -1.0*( dmass_X + dmass_Y )
        dnumber_Z = -1.0*dnumber_X
      end if

      if (.not. l_alternate_Z) then
        dmass_X=-dmass_Y
        dnumber_X=0.0
        dnumber_Z=0.0
        dmass_Z=0.0
      end if


!convert back to grid mean
     !use mixed-phase overlap function
      if (((params_Y%id == cloud_params%id) .OR. (params_X%id == cloud_params%id)) &
         .OR. ((params_Y%id == rain_params%id) .OR. (params_X%id == rain_params%id)) ) then
        overlap_cf=min(1.0,max(0.0,mpof*min(cf_X, cf_Y) +  max(0.0,(1.0-mpof)*(cf_X+cf_Y-1.0))))
      else
        overlap_cf = min(cf_Y, cf_X)
      endif
       !!overlap_cf=1.0
     
      dmass_Y = dmass_Y * overlap_cf
      dmass_X = dmass_X * overlap_cf
      if (params_Y%l_2m) dnumber_Y = dnumber_Y * overlap_cf
      if (params_X%l_2m) dnumber_X = dnumber_X * overlap_cf
      if (l_alternate_Z) then
        dmass_Z = dmass_Z * overlap_cf
        if (params_Z%l_2m) dnumber_Z = dnumber_Z * overlap_cf
      endif

      procs(params_Y%i_1m, iproc%id)%column_data(k)=dmass_Y
      if (params_Y%l_2m) then
        procs(params_Y%i_2m, iproc%id)%column_data(k)=dnumber_Y 
      end if

      procs(params_X%i_1m, iproc%id)%column_data(k)=dmass_X  
      if (params_X%l_2m) then
        procs(params_X%i_2m, iproc%id)%column_data(k)=dnumber_X  
      end if

      if (l_alternate_Z) then
        procs(Zi_1m, iproc%id)%column_data(k)=dmass_Z
        if (l_Zi_2m) then
          procs(Zi_2m, iproc%id)%column_data(k)=dnumber_Z
        end if
      end if

      !----------------------
      ! Aerosol processing...
      !----------------------

      if (l_process .and. l_aero) then
        ! We note that all processes result in a source of frozen water
        ! i.e. params_Z is either ice, snow or graupel
     
        if (params_Y%id == cloud_params%id) then
          dmac=abs(dnumber_Y)*aeroact(k)%mact1_mean*aeroact(k)%nratio1
          dmad=abs(dnumber_Y)*dustliq(k)%mact1_mean*dustliq(k)%nratio1
        else if (params_Y%id == rain_params%id) then
          dmac=abs(dnumber_Y)*aeroact(k)%mact2_mean*aeroact(k)%nratio2
          dmad=abs(dnumber_Y)*dustliq(k)%mact2_mean*dustliq(k)%nratio2
        else if (params_X%id == cloud_params%id) then ! This is never the case !
          dmac=abs(dnumber_X)*aeroact(k)%mact1_mean*aeroact(k)%nratio1
          dmad=abs(dnumber_X)*dustliq(k)%mact1_mean*dustliq(k)%nratio1
        else if (params_X%id == rain_params%id) then
          dmac=abs(dnumber_X)*aeroact(k)%mact2_mean*aeroact(k)%nratio2
          dmad=abs(dnumber_X)*dustliq(k)%mact2_mean*dustliq(k)%nratio2
        end if

        aerosol_procs(i_am8, iaproc%id)%column_data(k)=dmac
        aerosol_procs(i_am4, iaproc%id)%column_data(k)=-dmac
        aerosol_procs(i_am7, iaproc%id)%column_data(k)=dmad
        aerosol_procs(i_am9, iaproc%id)%column_data(k)=-dmad
      end if
    end if
  end if
enddo

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine iacc
end module ice_accretion
