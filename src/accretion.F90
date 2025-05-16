module accretion
  use variable_precision, only: wp
  use passive_fields, only: rho
! use mphys_switches, only: i_m3r, l_3mr
  use mphys_switches, only: i_ql, i_qr, i_nl, l_2mc, &
       l_aacc, i_am4, i_am5, l_process, active_rain, isol, l_preventsmall, &
       l_prf_cfrac, i_cfl, i_cfr, l_kk00
  use mphys_constants, only: fixed_cloud_number
  use mphys_parameters, only: hydro_params
! use mphys_parameters, only: p1, p2, p3, rain_params
  use process_routines, only: process_rate, i_pracw, i_aacw
  use thresholds, only: ql_small, qr_small, cfliq_small
  use sweepout_rate, only: sweepout
  use distributions, only: dist_lambda, dist_mu, dist_n0
! use m3_incs, only: m3_inc_type2
  use casim_stph, only: l_rp2_casim, fixed_cloud_number_rp

  implicit none

  character(len=*), parameter, private :: ModuleName='ACCRETION'

  private

  public racw
contains
  !---------------------------------------------------------------------------
  !> @author
  !> Ben Shipway
  !
  !> @brief
  !> This subroutine calculates increments due to the accretion of
  !> cloud water by rain
  !--------------------------------------------------------------------------- !
  subroutine racw(ixy_inner, dt, qfields, cffields, aerofields, procs, params, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments

    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt !< microphysics time increment (s)   
                         !! dt NEEDED for 3rd moment code
    real(wp), intent(in) :: qfields(:,:)     !< hydrometeor fields
    real(wp), intent(in) :: cffields(:,:)     ! < cloud fractions
    real(wp), intent(in) :: aerofields(:,:)  !< aerosol fields
    type(process_rate), intent(inout), target :: procs(:,:)         !< hydrometeor process rates
    type(process_rate), intent(inout), target :: aerosol_procs(:,:) !< aerosol process rates
    type(hydro_params), intent(in) :: params !< parameters describing hydrometor size distribution/fallspeeds etc.

    ! Local Variables

    real(wp) :: dmass, dnumber, damass
!   real(wp) :: m1, m2, m3, dm1, dm2, dm3

    real(wp) :: cloud_mass
    real(wp) :: cloud_number
    real(wp) :: rain_mass
!   real(wp) :: rain_number
!   real(wp) :: rain_m3
    real(wp) :: cf_liquid, cf_rain


    real(wp) :: mu, n0, lam
    logical :: l_kk_acw=.true.

    integer :: k ! local index for k

    character(len=*), parameter :: RoutineName='RACW'

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

    do k = 1, ubound(qfields,1)
       if (l_prf_cfrac) then
          if (cffields(k,i_cfl) .gt. cfliq_small) then
             ! only doing liquid cloud fraction at the moment
             cf_liquid=cffields(k,i_cfl)
          else
             cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
          endif
          if (cffields(k,i_cfr) .gt. cfliq_small) then
             ! only doing liquid cloud fraction at the moment
             cf_rain=cffields(k,i_cfr)
          else
             cf_rain=cfliq_small !nonzero value - maybe move cf test higher up
          endif
       else
          cf_liquid=1.0_wp
          cf_rain=1.0_wp
       endif
       
       cloud_mass = qfields(k, i_ql) / cf_liquid
       rain_mass = qfields(k, i_qr) / cf_rain
       
       if (l_2mc ) then
          cloud_number=qfields(k, i_nl) / cf_liquid
       else
          cloud_number=fixed_cloud_number / cf_liquid !check that fixed_cloud_number is grid mean
       end if
       
    
       ! if (l_3mr) rain_m3 = qfields(k, i_m3r)
       
       if (cloud_mass*cf_liquid > ql_small .and. rain_mass*cf_rain > qr_small) then
          if (l_kk_acw) then
             !        dmass=min(0.9*cloud_mass, 67.0*(cloud_mass*rain_mass)**1.15)
             if (l_kk00) then
                ! Use KK accretion parametrisation but limit to 90% of cloud mass removal
                dmass = MIN(0.9*cloud_mass, 67.0*(cloud_mass*rain_mass)**1.15)
             else
                ! Use Kogan(2013) accretion parametrisation but limit to 90% 
                ! of cloud mass removal
                dmass = min(0.9*cloud_mass, 8.53*(cloud_mass**1.05)*(rain_mass)**0.98)
             endif
             
          else
             n0=dist_n0(k,params%id)
             mu=dist_mu(k,params%id)
             lam=dist_lambda(k,params%id)
             dmass=sweepout(n0, lam, mu, params, rho(k,ixy_inner))*cloud_mass
          end if
          
          if (l_preventsmall .and. dmass < qr_small) dmass=0.0
          if (l_2mc) dnumber=dmass/(cloud_mass/cloud_number)
          
          
          if (l_prf_cfrac) then
             ! convert back to grid mean
             dmass=dmass*min(cf_liquid, cf_rain)
             dnumber=dnumber*min(cf_liquid, cf_rain)
             cloud_mass=cloud_mass*cf_liquid
          end if
          
          
          procs(i_ql, i_pracw%id)%column_data(k)=-dmass
          procs(i_qr, i_pracw%id)%column_data(k)=dmass
          if (l_2mc) then
             procs(i_nl, i_pracw%id)%column_data(k)=-dnumber
          end if
       
      ! if (l_3mr) then
      !    m1=rain_mass/rain_params%c_x
      !    m2=rain_number
      !    m3=rain_m3
             
      !    dm1=dt*dmass/rain_params%c_x
      !    dm2=0
      !    call m3_inc_type2(m1, m2, m3, p1, p2, p3, dm1, dm2, dm3)
      !    dm3=dm3/dt
      !    procs(i_m3r, i_pracw%id)%column_data(k) = dm3
      ! end if

    
       end if
    enddo

    if (l_aacc .and. l_process) then
       if (active_rain(isol)) then
          do k = 1, ubound(qfields,1)
             dmass = procs(i_qr, i_pracw%id)%column_data(k)
             damass=dmass/cloud_mass*aerofields(k,i_am4)
             aerosol_procs(i_am4, i_aacw%id)%column_data(k)=-damass
             aerosol_procs(i_am5, i_aacw%id)%column_data(k)=damass
          enddo
       end if
    end if
          
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine racw
end module accretion
