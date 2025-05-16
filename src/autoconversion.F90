module autoconversion
  use variable_precision, only: wp
  use passive_fields, only: rho
  use mphys_switches, only: i_ql, i_qr, i_nl, i_nr, l_2mc, &
       l_2mr, l_aaut, i_am4, i_am5, cloud_params, rain_params, l_process, &
       l_separate_rain, l_preventsmall, l_prf_cfrac, i_cfl, l_kk00
! use mphys_switches, only: m3r, l_3mr
  use mphys_constants, only: fixed_cloud_number
  use mphys_parameters, only: rain_params
! use mphys_parameters, only: mu_aut
  use process_routines, only: process_rate, i_praut, i_aaut
  use thresholds, only: ql_small, nl_small, qr_small, cfliq_small
! use m3_incs, only: m3_inc_type3
  use casim_stph, only: l_rp2_casim, fixed_cloud_number_rp

  implicit none

  character(len=*), parameter, private :: ModuleName='AUTOCONVERSION'

  private

  public raut
contains

  subroutine raut(ixy_inner, dt, qfields, cffields, aerofields, procs, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: qfields(:,:), aerofields(:,:),    cffields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    ! Local variables
    real(wp) :: dmass, dnumber1, dnumber2, damass
!   real(wp) :: dm1,dm2,dm3
    real(wp) :: cloud_mass
    real(wp) :: cloud_number
!   real(wp) :: p1, p2, p3
!   real(wp) :: k1, k2, k3
    real(wp) :: mu_qc ! < cloud shape parameter (currently only used diagnostically here)
    real(wp) :: cf_liquid

    integer :: k
    character(len=*), parameter :: RoutineName='RAUT'

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
    end if

   do k = 1, ubound(qfields,1)
   if (l_prf_cfrac) then

      if (cffields(k,i_cfl) .gt. cfliq_small) then  !only doing liquid cloud fraction at the moment
        cf_liquid=cffields(k,i_cfl)
      else
        cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
      endif

    else ! l_prf_cfrac

      cf_liquid = 1.0_wp

    end if

    cloud_mass=qfields(k, i_ql)/cf_liquid

    if (l_2mc) then
      cloud_number=qfields(k, i_nl)/cf_liquid
    else
      cloud_number=fixed_cloud_number
    end if

    if (cloud_mass *cf_liquid> ql_small .and. &
         cloud_number * cf_liquid > nl_small) then
      if (l_kk00) then
         dmass = 1350.*cloud_mass**2.47*  &
              (cloud_number/1.e6*rho(k,ixy_inner))**(-1.79)
      else
         ! new method, k13 scheme
         dmass = 7.98e10*cloud_mass**4.22*  &
              (cloud_number/1.e6*rho(k,ixy_inner))**(-3.01)
      endif

      dmass=min(.25*cloud_mass/dt, dmass)
      if (l_preventsmall .and. dmass < qr_small) dmass=0.0
      if (l_2mc) dnumber1=dmass/(cloud_mass/cloud_number)
      mu_qc=min(15.0_wp, (1000.0E6/cloud_number + 2.0))
      ! The following line is correct for l_kk00 = .true., i.e. use KK2000 for 
      ! autoconversion and accretion, which is the default switch in 
      ! mphys_switches. However, need to confirm this is 
      ! correct for Kogan (k13), i.e. l_kk00 = .false., since if the 50.0E-6
      ! is the notional diameter threshold for autoconversion it needs to 
      ! be changed 
      if (l_2mr) dnumber2=dmass/(rain_params%c_x*(mu_qc/3.0)*(50.0E-6)**3)

          ! AH - found that at the cloud edges rain number produced 
          !      by autoconversion can be larger than cloud number removed. This
          !      is not physically possible, so limit the dnumber2 to be the 
          !      same as dnumber1
          if (dnumber2 > dnumber1) Then
             dnumber2 = dnumber1
          endif

          ! if (l_3mr) then
          !    dm1=dt*dmass/rain_params%c_x
          !    dm2=dt*dnumber2
          !    p1=rain_params%p1
          !    p2=rain_params%p2
          !    p3=rain_params%p3
          !    call m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu_aut)
          !    dm3=dm3/dt
          ! end if

!convert back to grid mean
      dmass=dmass*cf_liquid
      dnumber1=dnumber1*cf_liquid
      dnumber2=dnumber2*cf_liquid
      cloud_mass=cloud_mass*cf_liquid !for aerosol processing below

          procs(i_ql, i_praut%id)%column_data(k)=-dmass
          procs(i_qr, i_praut%id)%column_data(k)=dmass

          if (cloud_params%l_2m) then
             procs(i_nl, i_praut%id)%column_data(k)=-dnumber1
          end if
          if (rain_params%l_2m) then
             procs(i_nr, i_praut%id)%column_data(k)=dnumber2
          end if
          ! if (rain_params%l_3m) then
          !    procs(i_m3r, i_praut%id)%column_data(k)=dm3
          ! end if
          
          if (l_separate_rain) then
             if (l_aaut .and. l_process) then
                ! Standard Single soluble mode, 2 activated species
                damass=dmass/cloud_mass*aerofields(k,i_am4)
                aerosol_procs(i_am4, i_aaut%id)%column_data(k)=-damass
                aerosol_procs(i_am5, i_aaut%id)%column_data(k)=damass
             end if
          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine raut
end module autoconversion
