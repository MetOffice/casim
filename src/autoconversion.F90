module autoconversion
  use variable_precision, only: wp
  use passive_fields, only: rho
  use mphys_switches, only: i_ql, i_qr, i_nl, i_nr, i_m3r, hydro_complexity, l_2mc, &
       l_2mr, l_3mr, l_aaut, i_am4, i_am5, cloud_params, rain_params, l_process, l_separate_rain, l_preventsmall
  use mphys_constants, only: rhow, fixed_cloud_number
  use mphys_parameters, only: mu_aut, rain_params
  use process_routines, only: process_rate, i_praut, i_aaut
  use thresholds, only: ql_small, nl_small, qr_small, qr_sig
  use special, only: pi
  use m3_incs, only: m3_inc_type2, m3_inc_type3

  implicit none

  character(len=*), parameter, private :: ModuleName='AUTOCONVERSION'

  private

  public raut
contains

  subroutine raut(dt, k, qfields, aerofields, procs, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    real(wp), intent(in) :: dt
    integer,  intent(in) :: k
    real(wp), intent(in) :: qfields(:,:), aerofields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    ! Local variables
    real(wp) :: dmass, dnumber1, dnumber2, damass
    real(wp) :: m1, m2, m3, dm1,dm2,dm3
    real(wp) :: n0, lam, mu
    real(wp) :: cloud_mass
    real(wp) :: cloud_number
    real(wp) :: rain_mass
    real(wp) :: rain_number
    real(wp) :: rain_m3
    type(process_rate), pointer :: this_proc
    real(wp) :: p1, p2, p3
    real(wp) :: k1, k2, k3
    real(wp) :: mu_qc ! < cloud shape parameter (currently only used diagnostically here)

    integer :: i
    character(len=*), parameter :: RoutineName='RAUT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    cloud_mass=qfields(k, i_ql)
    if (l_2mc) then
      cloud_number=qfields(k, i_nl)
    else
      cloud_number=fixed_cloud_number
    end if

    rain_mass=qfields(k, i_qr)
    if (l_2mr) rain_number=qfields(k, i_nr)
    if (l_3mr) rain_m3=qfields(k, i_m3r)

    if (cloud_mass > ql_small .and. cloud_number > nl_small) then
      this_proc=>procs(k, i_praut%id)
      dmass=1350.0*cloud_mass**2.47*(cloud_number/1.0e6*rho(k))**(-1.79)
      dmass=min(.25*cloud_mass/dt, dmass)
      if (l_preventsmall .and. dmass < qr_small) dmass=0.0
      if (l_2mc) dnumber1=dmass/(cloud_mass/cloud_number)
      mu_qc=min(15.0_wp, (1000.0E6/cloud_number + 2.0))
      if (l_2mr) dnumber2=dmass/(rain_params%c_x*(mu_qc/3.0)*(50.0E-6)**3)

      if (l_3mr) then
        dm1=dt*dmass/rain_params%c_x
        dm2=dt*dnumber2
        p1=rain_params%p1
        p2=rain_params%p2
        p3=rain_params%p3
        m1=rain_mass/rain_params%c_x
        m2=rain_number
        m3=rain_m3
        call m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu_aut)
        dm3=dm3/dt
      end if
      this_proc%source(i_ql)=-dmass
      this_proc%source(i_qr)=dmass

      if (cloud_params%l_2m) then
        this_proc%source(i_nl)=-dnumber1
      end if
      if (rain_params%l_2m) then
        this_proc%source(i_nr)=dnumber2
      end if
      if (rain_params%l_3m) then
        this_proc%source(i_m3r)=dm3
      end if

      if (l_separate_rain) then
        if (l_aaut .and. l_process) then
          ! Standard Single soluble mode, 2 activated species
          damass=dmass/cloud_mass*aerofields(k,i_am4)
          aerosol_procs(k, i_aaut%id)%source(i_am4)=-damass
          aerosol_procs(k, i_aaut%id)%source(i_am5)=damass
        end if
      end if
      nullify(this_proc)
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine raut
end module autoconversion
