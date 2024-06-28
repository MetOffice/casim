module graupel_wetgrowth
  use variable_precision, only: wp
  use process_routines, only: process_rate, process_name,   &
       i_gshd, i_gacw, i_gacr, i_gaci, i_gacs
  use passive_fields, only: TdegC,  qws0, rho

  use mphys_parameters, only: ice_params, snow_params, graupel_params,   &
       rain_params, T_hom_freeze, DR_melt
  use mphys_switches, only: i_qv, l_kfsm, l_gamma_online, l_prf_cfrac, i_cfg
  use mphys_constants, only: Lv, Lf, Ka, Cwater, Cice, Dv
  use thresholds, only: thresh_small, cfliq_small

  use ventfac, only: ventilation_3M, ventilation_1M_2M
  use distributions, only: dist_lambda, dist_mu, dist_n0

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='GRAUPEL_WETGROWTH'

  public wetgrowth
contains
  subroutine wetgrowth(ixy_inner, nz, l_Tcold, qfields, cffields, procs, l_sigevap)
    !< Subroutine to determine if all liquid accreted by graupel can be
    !< frozen or if there will be some shedding
    !<
    !< CODE TIDYING: Should move efficiencies into parameters


    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    integer, intent(in) :: nz
    logical, intent(in) :: l_Tcold(:)
    real(wp), intent(in), target :: qfields(:,:)
    real(wp), intent(in) :: cffields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    logical, intent(in) :: l_sigevap(:) ! logical to determine significant evaporation

    ! Local variables

    type(process_name) :: iproc ! processes selected depending on
    ! which species we're modifying

    real(wp) :: qv
    real(wp) :: dmass, dnumber
    real(wp) :: mass

    real(wp) :: n0, lam, mu, V_x
    real(wp) :: pgacw, pgacr, pgaci, pgacs
    real(wp) :: Eff, sdryfac, idryfac
    real(wp) :: pgaci_dry, pgacs_dry, pgdry
    real(wp) :: pgwet       !< Amount of liquid that graupel can freeze withouth shedding
    real(wp) :: pgacsum     !< Sum of all graupel accretion terms
    real(wp) :: cf_graupel

    integer :: k

    character(len=*), parameter :: RoutineName='WETGROWTH'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k = 1, nz
       if (.not. l_sigevap(k) .and. l_Tcold(k)) then 
          if (l_prf_cfrac) then
            if (cffields(k,i_cfg) .gt. cfliq_small) then
               cf_graupel=cffields(k,i_cfg)
            else
               cf_graupel=cfliq_small !nonzero value - maybe move cf test higher up
            endif
          else
            cf_graupel=1.0
          endif

          mass=qfields(k, graupel_params%i_1m) / cf_graupel

          if (mass > thresh_small(graupel_params%i_1m)) then
             pgacw=procs(graupel_params%i_1m, i_gacw%id)%column_data(k)/cf_graupel !convert to ingraupel rates
             pgacr=procs(graupel_params%i_1m, i_gacr%id)%column_data(k)/cf_graupel !convert to ingraupel rates
             pgaci=procs(graupel_params%i_1m, i_gaci%id)%column_data(k)/cf_graupel !convert to ingraupel rates
             
             if (l_kfsm) then
                pgacs = 0.0 !<KF. Single ice prognostics. No process gacs. Set to zero.>
             else
                pgacs=procs(graupel_params%i_1m, i_gacs%id)%column_data(k) /cf_graupel !convert to ingraupel rates
             end if ! l_kfsm
             
             ! Factors for converting wet collection efficiencies to dry ones
             ! AH - Changed 1.0 to 0.001 to match RA and GA tests. 0.001 prevents a divide by zero but is equivalent 
             !      to a 0.0 eff. 
             Eff=0.001
             sdryfac=1.0_wp ! min(1.0_wp, 0.2*exp(0.08*TdegC(k))/Eff)
             Eff=0.001
             idryfac=1.0_wp !min(1.0_wp, 0.2*exp(0.08*TdegC(k))/Eff)
             
             pgaci_dry=idryfac*pgaci
             pgacs_dry=sdryfac*pgacs
             
             pgdry=pgacw+pgacr+pgaci_dry+pgacs_dry
             
             pgacsum=pgacw+pgacr+pgaci+pgacs
             
             if (pgacsum > thresh_small(graupel_params%i_1m)) then
                qv=qfields(k, i_qv)
                mass=qfields(k, graupel_params%i_1m) / cf_graupel
                
                n0=dist_n0(k,graupel_params%id)
                mu=dist_mu(k,graupel_params%id)
                lam=dist_lambda(k,graupel_params%id)
                
                if (l_gamma_online) then 
                   call ventilation_3M(ixy_inner, k, V_x, n0, lam, mu, graupel_params)
                else 
                   call ventilation_1M_2M(ixy_inner, k, V_x, n0, lam, mu, graupel_params)
                endif

                pgwet=(910.0/graupel_params%density)**0.625*(Lv*Dv*(qws0(k,ixy_inner)-qv)- &
                     Ka*TdegC(k,ixy_inner)/rho(k,ixy_inner))/(Lf+Cwater*TdegC(k,ixy_inner))*V_x *cf_graupel ! grid mean
                pgwet=pgwet+(pgaci+pgacs)*(1.0-Cice*TdegC(k,ixy_inner)/(Lf+Cwater*TdegC(k,ixy_inner)))

                if (pgdry < pgwet .or. TdegC(k,ixy_inner) < T_hom_freeze) then ! Dry growth, so use recalculated gaci, gacs
                   if (pgaci + pgacs > 0.0) then

                      if ( .not. l_kfsm ) then
                         ! KF. Single ice prognostics. No process gacs
                         dmass=pgacs_dry * cf_graupel !convert back to gridbox mean
                         procs(graupel_params%i_1m, i_gacs%id)%column_data(k)=dmass
                         procs(snow_params%i_1m, i_gacs%id)%column_data(k)=-dmass
                         if (snow_params%l_2m) then
                            dnumber=procs(snow_params%i_2m, i_gacs%id)%column_data(k)*sdryfac
                            procs(snow_params%i_2m, i_gacs%id)%column_data(k)=dnumber
                         end if
                      end if ! l_kfsm
                      
                      dmass=pgaci_dry * cf_graupel !convert back to gridbox mean
                      procs(graupel_params%i_1m, i_gaci%id)%column_data(k)=dmass
                      procs(ice_params%i_1m, i_gaci%id)%column_data(k)=-dmass
                      if (ice_params%l_2m) then
                         dnumber=procs(ice_params%i_2m, i_gaci%id)%column_data(k)*idryfac
                         procs(ice_params%i_2m, i_gaci%id)%column_data(k)=dnumber
                      end if
                   end if
                else ! Wet growth mode so recalculate gacr, gshd
                   if (abs(procs(graupel_params%i_1m, i_gacr%id)%column_data(k)) > &
                        thresh_small(graupel_params%i_1m)) then
                      dmass=(pgacr+pgwet-pgacsum) * cf_graupel !convert back to gridbox mean
                      if (dmass > 0.0) then
                         procs(graupel_params%i_1m, i_gacr%id)%column_data(k)=dmass
                         procs(rain_params%i_1m, i_gacr%id)%column_data(k)=-dmass
                      else
                         iproc=i_gshd
                         dmass=min(pgacw * cf_graupel, -1.0*dmass)
                         procs(rain_params%i_1m, iproc%id)%column_data(k)=dmass
                         if (rain_params%l_2m) then
                            dnumber=dmass/(rain_params%c_x*DR_melt**3)
                            procs(rain_params%i_2m, iproc%id)%column_data(k)=dnumber
                         end if
                      end if
                   end if
                end if
             end if
          end if
       endif
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine wetgrowth
end module graupel_wetgrowth
