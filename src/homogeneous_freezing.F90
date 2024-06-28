module homogeneous
  use variable_precision, only: wp
  use passive_fields, only: rho, pressure, w, exner
  use mphys_switches, only: i_qv, i_ql, i_qi, i_ni, i_th , hydro_complexity, i_am6, i_an2, l_2mi, l_2ms, l_2mg &
       , i_ns, i_ng, iopt_inuc, i_am7, i_an6, i_am9 , i_m3r, i_m3g, i_qr, i_qg, i_nr, i_ng, i_nl          &
       , isol, i_am4, i_am8, active_ice, l_process, l_prf_cfrac
  use process_routines, only: process_rate, i_homr, i_homc, i_dhomc, i_dhomr
  use mphys_parameters, only: rain_params, graupel_params, cloud_params, ice_params, T_hom_freeze
  use mphys_constants, only:  Ls, cp, Mw, g, Rd, Ru, Lv, Lf, Dv, Rv
  use qsat_funs, only: qsaturation, qisaturation
  use thresholds, only: thresh_small, thresh_tidy
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use special, only: GammaFunc, pi
  use m3_incs, only: m3_inc_type2, m3_inc_type4

  implicit none

  character(len=*), parameter, private :: ModuleName='HOMOGENEOUS'

contains

  !> Calculates immersion freezing of rain drops
  !> See Bigg 1953
  subroutine ihom_rain(ixy_inner, dt, nz, l_Tcold, qfields, l_sigevap, aeroact, dustliq,  &
                      procs, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    logical, intent(in) :: l_Tcold(:)
    real(wp), intent(in), target :: qfields(:,:)
    type(aerosol_active), intent(in) :: aeroact(:), dustliq(:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    logical, intent(in) :: l_sigevap(:) ! logical to determine significant evaporation

    real(wp) :: dmass, dnumber, dmac, coef, dmadl
!    real(wp) :: dm1, dm2, dm3, dm3_g, m1, m2, m3
    real(wp) :: n0, lam, mu

    real(wp) :: th
    real(wp) :: qr, nr

    real(wp) :: Tc

    real(wp), parameter :: A_bigg = 0.66, B_bigg = 100.0

    logical :: l_condition, l_freezeall
    logical :: l_ziegler=.true.

    integer :: k

    character(len=*), parameter :: RoutineName='IHOM_RAIN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k = 1, nz
       if (l_Tcold(k)) then 
          th = qfields(k, i_th)
          Tc = th*exner(k,ixy_inner) - 273.15
          qr = qfields(k, i_qr)
          if (rain_params%l_2m)then
             nr = qfields(k, i_nr)
          else
             nr = 1000. ! PRAGMATIC SM HACK
          end if
          
          l_condition=(Tc < -4.0 .and. qr > thresh_small(i_qr) .and. .not. l_sigevap(k))
          l_freezeall=.false.
          if (l_condition) then
             
             n0 = dist_n0(k,rain_params%id)
             mu = dist_mu(k,rain_params%id)
             lam = dist_lambda(k,rain_params%id)
             
             if (l_ziegler) then
                dnumber = B_bigg*(exp(-A_bigg*Tc)-1.0)*rho(k,ixy_inner)*qr/rain_params%density
                dnumber = min(dnumber, nr/dt)
                dmass = (qr/nr)*dnumber
             else
                coef = B_bigg*(pi/6.0)*(exp(-A_bigg*Tc)-1.0)/rho(k,ixy_inner)               &
                     * n0 /(lam*lam*lam)/GammaFunc(1.0 + mu)
                
                dmass = coef * rain_params%c_x * lam**(-rain_params%d_x)          &
                     * GammaFunc(4.0 + mu + rain_params%d_x)
                
                dnumber = coef                                                    &
                     * GammaFunc(4.0 + mu)
                
             end if
             
             ! PRAGMATIC HACK - FIX ME
             ! If most of the drops are frozen, do all of them
             if (dmass*dt >0.95*qr .or. dnumber*dt > 0.95*nr) then
                dmass=qr/dt
                dnumber=nr/dt
                l_freezeall=.true.
             end if

             procs(i_qr, i_homr%id)%column_data(k) = -dmass
             procs(i_qg, i_homr%id)%column_data(k) = dmass
             
             if (rain_params%l_2m) then
                procs(i_nr, i_homr%id)%column_data(k) = -dnumber
             end if
             if (graupel_params%l_2m) then
                procs(i_ng, i_homr%id)%column_data(k) = dnumber
             end if
             ! 3-moment code is commented for future implementation
             ! if (rain_params%l_3m) then
             !   if (l_freezeall) then
             !     dm3=-qfields(k,i_m3r)/dt
             !   else
             !     m1=qr/rain_params%c_x
             !     m2=qfields(k,i_nr)
             !     m3=qfields(k,i_m3r)
             
             !     dm1=-dt*dmass/rain_params%c_x
             !     dm2=-dt*dnumber
             
             !     call m3_inc_type2(m1, m2, m3, rain_params%p1, rain_params%p2, rain_params%p3, dm1, dm2, dm3)
             !     dm3=dm3/dt
             !   end if
             !   procs(i_m3r, i_homr%id)%column_data(k) = dm3
             ! end if
             
             ! if (graupel_params%l_3m) then
             !   if (rain_params%l_3m) then
             !     call m3_inc_type4(dm3, graupel_params%c_x, rain_params%c_x, rain_params%p3, dm3_g)
             !   else
             !     m1=qfields(k,i_qg)/graupel_params%c_x
             !     m2=qfields(k,i_ng)
             !     m3=qfields(k,i_m3g)
             
             !     dm1=-dt*dmass/graupel_params%c_x
             !     dm2=-dt*dnumber
             !     call m3_inc_type2(m1, m2, m3, graupel_params%p1, graupel_params%p2, graupel_params%p3, dm1, dm2, dm3_g)
             !     dm3_g=dm3_g/dt
             !   end if
             !   procs(i_m3g, i_homr%id)%column_data(k) = dm3_g
             ! end if

             if (l_process) then
                dmac = dnumber*aeroact(k)%mact2_mean

                aerosol_procs(i_am8, i_dhomr%id)%column_data(k) = dmac
                aerosol_procs(i_am4, i_dhomr%id)%column_data(k) = -dmac
                
                ! Dust already in the liquid phase
                dmadl = dnumber*dustliq(k)%mact2_mean*dustliq(k)%nratio2
                if (dmadl /=0.0) then
                   aerosol_procs(i_am9, i_dhomr%id)%column_data(k) = -dmadl
                   aerosol_procs(i_am7, i_dhomr%id)%column_data(k) = dmadl
                end if
             end if
          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ihom_rain

  !> Calculates homogeneous freezing of cloud drops
  !> See Wisener 1972
  subroutine ihom_droplets(ixy_inner, dt, nz, l_Tcold, qfields, aeroact, dustliq, procs, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    logical, intent(in) :: l_Tcold(:) 
    real(wp), intent(in), target :: qfields(:,:)
    ! aerosol fields
    type(aerosol_active), intent(in) :: aeroact(:), dustliq(:)

    type(process_rate), intent(inout), target :: procs(:,:)
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    real(wp) :: dmass, dnumber, dmac, dmadl
!    real(wp) :: m1, m2, m3

    real(wp) :: th
    real(wp) :: qv, ql

    real(wp) :: Tc

    real(wp) :: Tk, cap, rhoi, Ei, Ew, bm, Ai, B0, Bis, aw, dnumberi, &
                ka, min_homog_ni, dniraw, d0_homog
    logical :: l_use_critical_w = .True.  !ni controlled by w and environmental conditions
    logical :: l_use_ni_limit = .False.   !ni limited to max per timestep
    
    
    logical :: l_condition

    integer :: k

    character(len=*), parameter :: RoutineName='IHOM_DROPLETS'


    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k = 1, nz
       if (l_Tcold(k)) then
          th = qfields(k, i_th)
          Tc = th*exner(k,ixy_inner) - 273.15
          ql = qfields(k, i_ql) 
          
          l_condition=(Tc < T_hom_freeze .and. ql > thresh_tidy(i_ql))
          if (l_condition) then
             dmass=min(ql, cp*(T_hom_freeze - Tc)/Lf)/dt    
             if (cloud_params%l_2m) then    
                dnumber=dmass*(qfields(k, i_nl))/ql   
                dnumberi=dnumber ! this will be overwritten if l_use_critical_W is used   
                
                !limit production to 1/cc in timestep - do this properly with KL but need to adapt for droplets
                if (l_use_ni_limit) then
                   dnumberi=min(dnumber, 1e2)
                   dmass=qfields(k, i_ql)*dnumberi/qfields(k, i_nl)
                endif
                
                if (l_use_critical_w) then
                   !! alternative limiter based on w - this is explicit w. Will overwrite dnumber
                   !! Method is based on derivation for sink of vapour to ice defined on Field et al, 2014, 
                   !! Mixed-phase clouds in a turbulent environment. Part 2: Analytic treatment, 
                   !! https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.2175
                   qv = qfields(k, i_qv)
                   Tk = th*exner(k,ixy_inner)
                   ka=((5.69+0.017*(Tk-273.15))*1e-5) * 418.6 
                   
                   cap=1.0 !assume spheres and radius used
                   rhoi=200.0
                   min_homog_ni=1e2 !kg-1
                   d0_homog=50e-6 !m
                   
                   Ei=qisaturation(Tk,pressure(k,ixy_inner)/100.)*pressure(k,ixy_inner)/ &
                        (qisaturation(Tk,pressure(k,ixy_inner)/100.)+0.62198) !vap press over ice [Pa]
                   Ew=qsaturation(Tk,pressure(k,ixy_inner)/100.)*pressure(k,ixy_inner)/ &
                        (qsaturation(Tk,pressure(k,ixy_inner)/100.)+0.62198) !vap press over liq [Pa]
                   
                   bm=1.0/(qv)+Lv*Lf/(cp*Rv*Tk**2)
                   Ai=1.0/(rhoi*Lf**2/(ka*Rv*Tk**2)+rhoi*Rv*Tk/(Ei*Dv))
                   B0=4.0*pi*cap*rhoi*Ai/rho(k,ixy_inner)
                   Bis=bm*B0*(Ew/Ei-1.0)
                   aw=g/(Rd*Tk)*(Lv*Rd/(cp*Rv*Tk)-1.0)
                   
                   dnumberi=max(w(k,ixy_inner),0.0)*(aw/Bis)/d0_homog  / dt  !convert to a rate
                   dniraw=dnumberi
                   !make max just below limit
                   dnumberi=min(dnumber*0.90, max(min_homog_ni, dnumberi))
                endif
             endif


             procs(i_ql, i_homc%id)%column_data(k)=-dmass  
             procs(i_qi, i_homc%id)%column_data(k)=dmass   
             
             if (cloud_params%l_2m) then
                procs(i_nl, i_homc%id)%column_data(k)=-dnumber   
                if (ice_params%l_2m) then
                   procs(i_ni, i_homc%id)%column_data(k)=dnumberi   
                   ! if l_use_critical_w or l_limit_ni is true then 
                   ! dnumber not necessarily equal to dnumberi 
                end if
             end if
             
             
!!!NB if we keep l_use_critical_w approach need to deal with dnumber=/=dnumberi for processing!!

             if (l_process) then
                dmac=dnumber*aeroact(k)%mact1_mean
                
                aerosol_procs(i_am8, i_dhomc%id)%column_data(k)=dmac
                aerosol_procs(i_am4, i_dhomc%id)%column_data(k)=-dmac
                
                ! Dust already in the liquid phase
                dmadl=dnumber*dustliq(k)%mact1_mean*dustliq(k)%nratio1
                if (dmadl /=0.0) then
                   aerosol_procs(i_am9, i_dhomc%id)%column_data(k)=-dmadl
                   aerosol_procs(i_am7, i_dhomc%id)%column_data(k)=dmadl
                end if
             end if
          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ihom_droplets
end module homogeneous
