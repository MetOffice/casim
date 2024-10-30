module ice_melting
  use variable_precision, only: wp
  use process_routines, only: process_rate, process_name, i_imlt, i_smlt, i_gmlt, i_sacw, i_sacr, i_gacw, i_gacr, i_gshd, &
       i_dimlt, i_dsmlt, i_dgmlt
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use passive_fields, only: TdegC, qws0, rho
  use mphys_parameters, only: ice_params, snow_params, graupel_params, rain_params, DR_melt, hydro_params, ZERO_REAL_WP
  use mphys_switches, only: i_qv, i_am4, i_am7, i_am8, i_am9, l_process, l_gamma_online
  use mphys_constants, only: Lv, Lf, Ka, Cwater, Cice, cp, Dv
  use thresholds, only: thresh_tidy
  use m3_incs, only: m3_inc_type2, m3_inc_type3, m3_inc_type4
  use ventfac, only: ventilation_1M_2M, ventilation_3M
  use distributions, only: dist_lambda, dist_mu, dist_n0

  implicit none

  character(len=*), parameter, private :: ModuleName='ICE_MELTING'

contains

  !> Subroutine to calculate rate of melting of ice species
  !>
  !> OPTIMISATION POSSIBILITIES: Shouldn't have to recalculate all 3m quantities
  !>                             If just rescaling mass conversion for dry mode
  subroutine melting(ixy_inner, dt, nz, params, qfields, cffields, procs, l_sigevap, aeroice, dustact, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    USE mphys_switches,       ONLY: i_cfs, i_cfg

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    type(hydro_params), intent(in) :: params
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)
    real(wp), intent(in) :: cffields(:,:)


    ! aerosol fields
    type(aerosol_active), intent(in) :: aeroice(:), dustact(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    logical, intent(in) :: l_sigevap(:) ! logical to determine significant evaporation

    type(process_name) :: iproc, iaproc ! processes selected depending on
    ! which species we're modifying
    type(process_name) :: i_acw, i_acr ! accretion processes
    real(wp) :: qv
    real(wp) :: dmass, dnumber
!   real(wp) :: dm1, dm2, dm3, dm3_r, m2, m3
    real(wp) :: num, mass, m1
    real(wp) :: n0, lam, mu, V_x
    real(wp) :: acc_correction
    logical :: l_meltall ! do we melt everything?
    real(wp) :: dmac, dmad
    real(wp) :: cf

    integer :: k
    
    character(len=*), parameter :: RoutineName='MELTING'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!intiliase
    dmac=0.0
    dmad=0.0
    dmass=0.0
    dnumber=0.0
    
    do k = 1, nz
       if (.not. l_sigevap(k)) then 
          l_meltall=.false.
          mass=qfields(k, params%i_1m)
          if (mass > thresh_tidy(params%i_1m) .and. TdegC(k,ixy_inner) > 0.0) then
             if (params%l_2m) num=qfields(k, params%i_2m)
             if (params%id == ice_params%id) then ! instantaneous removal
                iaproc=i_dimlt
                iproc=i_imlt
                dmass=mass/dt
                
                if (params%l_2m)then
                   dnumber=num/dt
                   procs(ice_params%i_2m, iproc%id)%column_data(k)=-dnumber
                   if (rain_params%l_2m) then
                      procs(rain_params%i_2m, iproc%id)%column_data(k)=dnumber
                   end if
                end if
                
                procs(ice_params%i_1m, iproc%id)%column_data(k)=-dmass
                procs(rain_params%i_1m, iproc%id)%column_data(k)=dmass
                
                ! 3-moment code retained for future implementation
                ! if (rain_params%l_3m) then
                !   m1=qfields(k, rain_params%i_1m)/rain_params%c_x
                !   m2=qfields(k, rain_params%i_2m)
                !   m3=qfields(k, rain_params%i_3m)
                !   dm1=dt*dmass/rain_params%c_x
                !   dm2=dt*dnumber
                
                !   if (m1 > 0.0) then
                !     call m3_inc_type2(m1, m2, m3, rain_params%p1,      &
                !          rain_params%p2, rain_params%p3, dm1, dm2, dm3_r)
                !   else
                !     call m3_inc_type3(rain_params%p1, rain_params%p2, rain_params%p3,      &
                !          dm1, dm2, dm3_r, rain_params%fix_mu)
                !   end if
                !   dm3_r=dm3_r/dt
                !   procs(rain_params%i_3m, iproc%id)%column_data(k) = dm3_r
                ! end if
             else
                if (params%id==snow_params%id) then
                   i_acw=i_sacw
                   i_acr=i_sacr
                   iproc=i_smlt
                   iaproc=i_dsmlt
                   cf=cffields(k,i_cfs)
                else if (params%id==graupel_params%id) then
                   i_acw=i_gacw
                   i_acr=i_gacr
                   iproc=i_gmlt
                   iaproc=i_dgmlt
                   cf=cffields(k,i_cfg)
                end if
                acc_correction=0.0
                if (i_acw%on)acc_correction=procs(params%i_1m, i_acw%id)%column_data(k)
                if (i_acr%on)acc_correction=acc_correction + procs(params%i_1m, i_acr%id)%column_data(k)
                if (params%id==graupel_params%id .and. i_gshd%on) then
                   acc_correction=acc_correction + procs(params%i_1m, i_gshd%id)%column_data(k)
                end if
                
                qv=qfields(k, i_qv)
                
                m1=mass/params%c_x
                if (params%l_2m) num=qfields(k, params%i_2m)
                ! 3-moment code retained for future implementation
                ! if (params%l_3m) m3=qfields(k, params%i_3m)
                
                n0=dist_n0(k,params%id)
                mu=dist_mu(k,params%id)
                lam=dist_lambda(k,params%id)
                
                if (l_gamma_online) then 
                   call ventilation_3M(ixy_inner, k, V_x, n0, lam, mu, params)
                else 
                   call ventilation_1M_2M(ixy_inner, k, V_x, n0, lam, mu, params)
                endif
                
                dmass=(1.0/(rho(k,ixy_inner)*Lf))*(Ka*TdegC(k,ixy_inner) + Lv*Dv*rho(k,ixy_inner)*(qv - qws0(k,ixy_inner))) * V_x&
                     + (Cwater*TdegC(k,ixy_inner)/Lf)*acc_correction
                dmass=dmass*cf ! grid mean


                dmass=max(dmass, ZERO_REAL_WP) ! ensure positive
                
                
                dmass=min(dmass, mass/dt) ! ensure we don't remove too much
                if (dmass*dt > 0.95*mass) then ! we're pretty much removing everything
                   l_meltall=.true.
                   dmass=mass/dt
                end if
                
                !--------------------------------------------------
                ! Apply spontaneous rain breakup if drops are large
                !--------------------------------------------------
                !< RAIN BREAKUP TO BE ADDED
                
                procs(params%i_1m, iproc%id)%column_data(k)=-dmass
                procs(rain_params%i_1m, iproc%id)%column_data(k)=dmass
                
                if (params%l_2m) then
                   dnumber=dmass*num/mass
                   procs(params%i_2m, iproc%id)%column_data(k)=-dnumber
                   procs(rain_params%i_2m, iproc%id)%column_data(k)=dnumber
                end if
                ! 3-moment code retained for future implementation
                ! if (params%l_3m) then
                !   if (l_meltall) then
                !     dm3=-m3/dt
                !   else
                !     dm1=-dt*dmass/params%c_x
                !     dm2=-dt*dnumber
                !     m2=num
                !     call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
                !     dm3=dm3/dt
                !   end if
                !   procs(params%i_3m, iproc%id)%column_data(k)=dm3
                ! end if
                
                ! if (rain_params%l_3m) then
                !   if (params%l_3m) then
                !     call m3_inc_type4(dm3, rain_params%c_x, params%c_x, params%p3, dm3_r)
                !   else
                !     m1=qfields(k, rain_params%i_1m)/rain_params%c_x
                !     m2=qfields(k, rain_params%i_2m)
                !     m3=qfields(k, rain_params%i_3m)
                
                !     dm1=dt*dmass/rain_params%c_x
                !     dm2=dt*dnumber
                !     call m3_inc_type2(m1, m2, m3, rain_params%p1, rain_params%p2, rain_params%p3, dm1, dm2, dm3_r, rain_params%fix_mu)
                !     dm3_r=dm3_r/dt
                !   end if
                !   procs(rain_params%i_3m, iproc%id)%column_data(k) = dm3_r
                ! end if
             end if
             !----------------------
             ! Aerosol processing...
             !----------------------
             
             if (l_process) then
                
                if (params%id == ice_params%id) then
                   dmac=dnumber*aeroice(k)%nratio1*aeroice(k)%mact1_mean
                   dmac=min(dmac, aeroice(k)%mact1 /dt)
                   dmad=dnumber*dustact(k)%nratio1*dustact(k)%mact1_mean
                   dmad=min(dmad, dustact(k)%mact1 /dt)
                else if (params%id == snow_params%id) then
                   dmac=dnumber*aeroice(k)%nratio2*aeroice(k)%mact2_mean
                   dmac=min(dmac, aeroice(k)%mact2 /dt)
                   dmad=dnumber*dustact(k)%nratio2*dustact(k)%mact2_mean
                   dmad=min(dmad, dustact(k)%mact2 /dt)
                else if (params%id == graupel_params%id) then
                   dmac=dnumber*aeroice(k)%nratio3*aeroice(k)%mact3_mean
                   dmac=min(dmac, aeroice(k)%mact3/dt)
                   dmad=dnumber*dustact(k)%nratio3*dustact(k)%mact3_mean
                   dmad=min(dmad, dustact(k)%mact3/dt)
                end if
                
                aerosol_procs(i_am8, iaproc%id)%column_data(k)=-dmac
                aerosol_procs(i_am4, iaproc%id)%column_data(k)=dmac
                aerosol_procs(i_am9, iaproc%id)%column_data(k)=dmad
                aerosol_procs(i_am7, iaproc%id)%column_data(k)=-dmad
             end if
          end if
       end if
    enddo
    
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine melting
end module ice_melting
