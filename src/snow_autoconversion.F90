module snow_autoconversion
  use variable_precision, only: wp
  use passive_fields, only: rho, exner, TdegK, pressure
  use mphys_switches, only:  &
       i_qi, i_qs, i_ni, i_ns, &
       i_qv, i_th, &
       l_harrington
! use mphys_switches, only:  i_m3s
  use mphys_constants, only: fixed_ice_number,   &
       Lv,ka, Dv, Rv
! use mphys_parameters, only: mu_saut
  use mphys_parameters, only: snow_params, ice_params,   &
       DImax, tau_saut, DI2S
  use process_routines, only: process_rate, i_saut
  use qsat_funs, only: qisaturation
  use thresholds, only: thresh_small
  use special, only: pi
! use m3_incs, only: m3_inc_type3

  use distributions, only: dist_lambda, dist_mU

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='SNOW_AUTOCONVERSION'

  public saut
contains

  subroutine saut(ixy_inner, dt, nz, l_Tcold, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SAUT'

    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    logical, intent(in) :: l_Tcold(:) 
    real(wp), intent(in) :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    real(wp) :: dmass, dnumber
!   real(wp) :: dm1,dm2,dm3
    real(wp) :: ice_lam, ice_mu
    real(wp) :: ice_mass
    real(wp) :: ice_number
    real(wp) :: th
    real(wp) :: qv
    real(wp) :: lami_min , AB, qis
    integer :: k

    logical :: l_condition ! logical condition to switch on process

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k = 1, nz
      if (l_Tcold(k)) then
          
        ice_mass=qfields(k, i_qi)
        if (ice_params%l_2m) then
          ice_number=qfields(k, i_ni)
        else
          ice_number=fixed_ice_number
        end if
          
        if (ice_mass >= thresh_small(i_qi) .and. ice_number >= thresh_small(i_ni)) then
          ice_mu=dist_mu(k,ice_params%id)
          ice_lam=dist_lambda(k,ice_params%id)
             
          lami_min=(1.0 + ice_mu)/DImax
             
          if (l_harrington) then
            qv=qfields(k, i_qv)
            th=qfields(k, i_th)
            qis=qisaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)
            l_condition=qv > qis
          else
            l_condition=ice_lam < lami_min ! LEM autconversion
          end if
          
          if (l_condition) then
             
            if (l_harrington) then
              !< AB This is used elsewhere, so we should do it more efficiently.
              AB=1.0/(Lv*Lv/(Rv*ka*TdegK(k,ixy_inner)*TdegK(k,ixy_inner))*   &
                 rho(k,ixy_inner)+1.0/(Dv*qis))
              dnumber=4.0/DImax/ice_params%density*(qv-qis)*                 &
                      rho(k,ixy_inner)*ice_number*exp(-ice_lam*DImax)*Dv/AB
              dnumber=min(dnumber,0.9*ice_number/dt)
              dmass=pi/6.0*ice_params%density*DImax**3*dnumber
              dmass=min(dmass, 0.7*ice_mass/dt)
            else ! LEM version
              dmass=((lami_min/ice_lam)**ice_params%d_x - 1.0)*ice_mass/tau_saut
              if (ice_mass < dmass*dt) then
                dmass=.5*ice_mass/dt
              end if
              dmass=min(dmass, .5*ice_mass/dt)
              dnumber=dmass/(ice_params%c_x*DI2S**ice_params%d_x)
            end if
             
            if (ice_number < dnumber*dt) then
              dnumber=dmass/ice_mass * ice_number
            end if
             

            if (dmass*dt > thresh_small(i_qs)) then
              ! if (snow_params%l_3m) then
              !   dm1=dt*dmass/snow_params%c_x
              !   dm2=dt*dnumber
              !   call m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu_saut)
              !   dm3=dm3/dt
              ! end if
               
              procs(i_qi, i_saut%id)%column_data(k)=-dmass
              procs(i_qs, i_saut%id)%column_data(k)=dmass
                
              if (ice_params%l_2m) then
                 procs(i_ni, i_saut%id)%column_data(k)=-dnumber
              end if
              if (snow_params%l_2m) then
                procs(i_ns, i_saut%id)%column_data(k)=dnumber
              end if
              ! if (snow_params%l_3m) then
              !   procs(i_m3s, i_saut%id)%column_data(k)=dm3
              ! end if
            end if
            !==============================
            ! No aerosol processing needed
            !==============================
          end if
        end if
      end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine saut
end module snow_autoconversion
