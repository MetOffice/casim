module accretion
  use variable_precision, only: wp
  use passive_fields, only: rho, cfliq
  use mphys_switches, only: i_ql, i_qr, i_nl, i_nr, i_m3r, l_2mc, l_2mr, l_3mr, &
       l_aacc, i_am4, i_am5, l_process, active_rain, isol, l_preventsmall, &
       l_prf_cfrac
  use mphys_constants, only: fixed_cloud_number
  use mphys_parameters, only: hydro_params, p1, p2, p3, rain_params
  use process_routines, only: process_rate, i_pracw, i_aacw
  use thresholds, only: ql_small, qr_small, cfliq_small
  use sweepout_rate, only: sweepout
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use m3_incs, only: m3_inc_type2

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
  subroutine racw(dt, qfields, aerofields, procs, params, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments

    real(wp), intent(in) :: dt !< microphysics time increment (s)
    real(wp), intent(in) :: qfields(:,:)     !< hydrometeor fields
    real(wp), intent(in) :: aerofields(:,:)  !< aerosol fields
    type(process_rate), intent(inout), target :: procs(:,:)         !< hydrometeor process rates
    type(process_rate), intent(inout), target :: aerosol_procs(:,:) !< aerosol process rates
    type(hydro_params), intent(in) :: params !< parameters describing hydrometor size distribution/fallspeeds etc.

    ! Local Variables

    real(wp) :: dmass, dnumber, damass
    real(wp) :: m1, m2, m3, dm1, dm2, dm3

    real(wp) :: cloud_mass
    real(wp) :: cloud_number
    real(wp) :: rain_mass
    real(wp) :: rain_number
    real(wp) :: rain_m3
    real(wp) :: cf_liquid


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

    do k = 1, ubound(qfields,1)
    if (l_prf_cfrac) then
      if (cfliq(k) .gt. cfliq_small) then
        ! only doing liquid cloud fraction at the moment
        cf_liquid=cfliq(k)
      else
        cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
      endif

      cloud_mass = qfields(k, i_ql) / cf_liquid

    else
      cloud_mass = qfields(k, i_ql)

    end if

    if (l_2mc .and. l_prf_cfrac) then
      cloud_number=qfields(k, i_nl) / cf_liquid
    else if (l_2mc) then
      cloud_number = qfields(k, i_nl)
    else
      cloud_number=fixed_cloud_number
    end if
    rain_mass=qfields(k, i_qr)
    if (l_2mr) rain_number=qfields(k, i_nr)
    if (l_3mr) rain_m3 = qfields(k, i_m3r)

    if (cloud_mass > ql_small .and. rain_mass > qr_small) then
      if (l_kk_acw) then
        dmass=min(0.9*cloud_mass, 67.0*(cloud_mass*rain_mass)**1.15)
      else
        n0=dist_n0(k,params%id)
        mu=dist_mu(k,params%id)
        lam=dist_lambda(k,params%id)
        dmass=sweepout(n0, lam, mu, params, rho(k))*cloud_mass
      end if

      if (l_preventsmall .and. dmass < qr_small) dmass=0.0
      if (l_2mc) dnumber=dmass/(cloud_mass/cloud_number)


      if (l_prf_cfrac) then
        ! convert back to grid mean
        dmass=dmass*cf_liquid
        dnumber=dnumber*cf_liquid
        cloud_mass=cloud_mass*cf_liquid
      end if


          procs(i_ql, i_pracw%id)%column_data(k)=-dmass
          procs(i_qr, i_pracw%id)%column_data(k)=dmass
          if (l_2mc) then
             procs(i_nl, i_pracw%id)%column_data(k)=-dnumber
          end if
       
!!$          if (l_3mr) then
!!$             m1=rain_mass/rain_params%c_x
!!$             m2=rain_number
!!$             m3=rain_m3
!!$             
!!$             dm1=dt*dmass/rain_params%c_x
!!$             dm2=0
!!$             call m3_inc_type2(m1, m2, m3, p1, p2, p3, dm1, dm2, dm3)
!!$             dm3=dm3/dt
!!$             procs(k, i_pracw%id)%source(i_m3r) = dm3
!!$          end if

          if (l_aacc .and. l_process) then
             if (active_rain(isol)) then
                damass=dmass/cloud_mass*aerofields(k,i_am4)
                aerosol_procs(i_am4, i_aacw%id)%column_data(k)=-damass
                aerosol_procs(i_am5, i_aacw%id)%column_data(k)=damass
             end if
          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine racw
end module accretion
