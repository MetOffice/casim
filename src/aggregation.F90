module aggregation
  use variable_precision, only: wp, iwp
  use passive_fields, only: rho, TdegC
  use mphys_switches, only: i_qr, i_nr, l_2mr, i_ns
! use mphys_switches, only: i_m3r, l_3mr
  use process_routines, only: process_rate,  process_name, i_pracr, i_iagg, i_sagg, i_gagg
  use mphys_parameters, only: d_r, hydro_params
! use mphys_parameters, only: p1, p2, p3, rain_params
  use mphys_constants, only: rhow, rho0
  use thresholds, only: qr_small, thresh_small, nr_small
! use m3_incs, only: m3_inc_type2
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use special, only: pi
  use gauss_casim_micro, only: gaussfunclookup, gaussfunclookup_2d
! use gauss_casim_micro, only: gauss_casim_func

  implicit none

  character(len=*), parameter, private :: ModuleName='AGGREGATION'

  private

  public racr, ice_aggregation
contains

  subroutine racr(ixy_inner, dt, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt           !! dt NEEDED for 3rd moment code
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! Local variables
    real(wp) :: dnumber, Dr, Eff
    real(wp) :: rain_mass
    real(wp) :: rain_number
!! variables below NEEDED for 3rd moment code
!   real(wp) :: m1, m2, m3, dm1, dm2, dm3
!   real(wp) :: rain_m3
    logical :: l_beheng=.true.

    integer :: k ! k index for looping over column

    character(len=*), parameter :: RoutineName='RACR'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dnumber=0.0
    if (l_2mr) then
       do k = 1, ubound(qfields,1)
          rain_mass=qfields(k, i_qr)
          rain_number=qfields(k, i_nr)
          ! if (l_3mr) rain_m3=qfields(k, i_m3r)

          if (rain_mass > qr_small .and. rain_number>nr_small) then
             if (l_beheng) then
                ! This commented code block imposes a limit on the size of 
                ! agg. This is based on a pragmatic choice that drops large
                ! than 600 microns will breakup. This code effectively negates
                ! excessive numerical size sorting in kinematic cases
                 Dr=(.75/pi)*(rain_mass/rain_number/rhow)**(1.0/d_r)
                 if ( Dr < 600.0e-6) then
                    ! Modified from original
                    Eff=.5
                 else
                    Eff=0.0
                 end if
                !
                ! For RA3 testing Beheng 1994 Atmos. Res. (equ 10) was used (below)
                ! Beheng uses cgs units, ie. g m-3, CASIM uses SI. Hence, the dnumber
                ! equation below is multiplied by rho to convert from cgs to SI 
                ! (also the 10e3 factor from Beheng is not included)
                !Eff=1.0 !Beheng 1994 Atmos. Res.
                dnumber=Eff * 8.0 * rain_number * rain_mass * rho(k,ixy_inner)   ! #/kg/s
                 !CFL limitation 
                dnumber=rain_number * min( dnumber*dt/rain_number, 0.5  )/dt                            
             end if
             procs(i_nr, i_pracr%id)%column_data(k)=-dnumber

            ! if (l_3mr) then
            !    m1=rain_mass/rain_params%c_x
            !    m2=rain_number
            !    m3=rain_m3
            !
            !    dm1=0.0
            !    dm2=-Dt*dnumber
            !    call m3_inc_type2(m1, m2, m3, p1, p2, p3, dm1, dm2, dm3)
            !    dm3=dm3/dt
            !    procs(i_m3r, i_pracr%id)%column_data(k) = dm3
            ! end if
          end if
       enddo
    end if



    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine racr

  !< Subroutine to determine the aggregation of
  !< ice, snow and graupel.  Aggregation causes a sink in number
  !< and for triple moment species a corresponding change in the
  !< 3rd moment assuming shape parameter is not changed
  !< NB: Aerosol mass is not modified by this process
  !
  !< CODE TIDYING: Move efficiencies into parameters
  subroutine ice_aggregation(ixy_inner, dt, nz, l_Tcold, params, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    USE special, only: Gammafunc

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    logical, intent(in) :: l_Tcold(:)   
    type(hydro_params), intent(in) :: params
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! Local variables
    type(process_name) :: iproc ! processes selected depending on which species we're modifying
    real(wp) :: dnumber
    real(wp) :: Eff ! collection efficiencies need to re-evaluate these and put them in properly to mphys_parameters
    real(wp) :: mass, gaussterm
    real(wp) :: n0, lam, mu

    integer :: k

    character(len=*), parameter :: RoutineName='ICE_AGGREGATION'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k = 1, nz
       if (l_Tcold(k)) then 
          Eff=min(1.0_wp, 0.2*exp(0.08*TdegC(k,ixy_inner)))
          select case (params%id)
          case (3_iwp) !ice
             iproc=i_iagg
          case (4_iwp) !snow
             iproc=i_sagg
             Eff=MIN(1.0_wp, 0.1*exp(0.08*TdegC(k,ixy_inner)))
          case (5_iwp) !graupel
             iproc=i_gagg
          end select
          
          mass=qfields(k, params%i_1m)

          if (mass > thresh_small(params%i_1m) .and. params%l_2m) then ! if no significant ice, we don't bother

             n0=dist_n0(k,params%id)
             mu=dist_mu(k,params%id)
             lam=dist_lambda(k,params%id)

             n0=n0*(lam**(mu+1.0))/(GammaFunc(1.0+mu))

             if (params%l_3m) then
                !gaussterm = gauss_casim_func(mu, params%b_x)
                call gaussfunclookup_2d(params%id, gaussterm, mu, params%b_x)
             else
                call gaussfunclookup(params%id, gaussterm)
             end if
             
             dnumber=-1.0*gaussterm * pi*0.125*                                &
                 (rho0/rho(k,ixy_inner))**params%g_x*params%a_x*n0*n0*         &
                 Eff*lam**(-(4.0 + 2.0*mu + params%b_x))
      
      
             dnumber=-1.0 * qfields(k, i_ns) *                                 &
                 min( -dnumber*dt/qfields(k, i_ns), 0.5  )/dt            

      
             procs(params%i_2m, iproc%id)%column_data(k)=dnumber
          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ice_aggregation
end module aggregation
