module ice_deposition
  use variable_precision, only: wp, iwp
  use passive_fields, only: rho, pressure, exner, TdegK
  use mphys_switches, only: i_qv, i_th   &
       , i_am6, i_an2  &
       , i_am7, i_an6                   &
       , l_process, l_passivenumbers_ice, l_passivenumbers &
       , i_an12 &
       , i_am8, i_am2, i_an11, l_gamma_online
  use type_process, only: process_name
  use process_routines, only: process_rate, i_idep,    &
       i_dsub, i_sdep, i_gdep, i_dssub, i_dgsub &
       , i_isub, i_ssub, i_gsub, i_iacw, i_raci, i_sacw, i_sacr  &
       , i_gacw, i_gacr
  use mphys_parameters, only: hydro_params, rain_params, cloud_params
  use mphys_constants, only: Ls, Lf, ka, Dv, Rv
  use qsat_funs, only: qisaturation
  use thresholds, only: thresh_small
  use aerosol_routines, only: aerosol_active

  use distributions, only: dist_lambda, dist_mu, dist_n0
  use ventfac, only: ventilation_1M_2M, ventilation_3M
! removing line below changes answers
  use special, only: pi

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='ICE_DEPOSITION'

  logical :: l_latenteffects = .false.

  public idep
contains

  !< Subroutine to determine the deposition/sublimation onto/from
  !< ice, snow and graupel.  There is no source/sink for number
  !< when undergoing deposition, but there is a sink when sublimating.
  subroutine idep(ixy_inner, dt, nz, l_Tcold, params, qfields, cffields, procs, dustact, aeroice, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    USE mphys_switches,       ONLY: i_cfi, i_cfs, i_cfg


    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(IN) :: dt
    integer, intent(IN) :: nz
    logical, intent(in) :: l_Tcold(:) 
    type(hydro_params), intent(IN) :: params
    real(wp), intent(IN) :: qfields(:,:)
    type(process_rate), intent(INOUT), target :: procs(:,:)
    real(wp), intent(in) :: cffields(:,:)

    ! aerosol fields
    type(aerosol_active), intent(IN) :: dustact(:), aeroice(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(INOUT), target :: aerosol_procs(:,:)


    ! Local Variables
    type(process_name) :: iproc, iaproc  ! processes selected depending on
    ! which species we're depositing on.

    type(process_name) :: i_acw, i_acr ! collection processes with cloud and rain
    real(wp) :: dmass, dnumber, dmad, dnumber_a, dnumber_d, dmac

    real(wp) :: th
    real(wp) :: qv
    real(wp) :: num, mass

    real(wp) :: qis(nz)
    real(wp) :: n0, lam, mu
    real(wp) :: V_x, AB
    real(wp) :: cf

    logical :: l_suball

    integer :: k

    character(len=*), parameter :: RoutineName='IDEP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do k = 1, nz
      if (l_Tcold(k)) then
        th=qfields(k, i_th)
        qis(k) = qisaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)
      end if
    end do

    do k = 1, nz
      if (l_Tcold(k)) then 
        l_suball=.false. ! do we want to sublimate everything?
          
        mass=qfields(k, params%i_1m)
          
        qv=qfields(k, i_qv)
          
        if (qv>qis(k)) then
          select case (params%id)
          case (3_iwp) !ice
            iproc=i_idep
            iaproc=i_dsub
            i_acw=i_iacw
            i_acr=i_raci
            cf=cffields(k,i_cfi)
          case (4_iwp) !snow
            iproc=i_sdep
            iaproc=i_dssub
            i_acw=i_sacw
            i_acr=i_sacr
            cf=cffields(k,i_cfs)
          case (5_iwp) !graupel
            iproc=i_gdep
            iaproc=i_dgsub
            i_acw=i_gacw
            i_acr=i_gacr
            cf=cffields(k,i_cfg)
          end select
        else
          select case (params%id)
          case (3_iwp) !ice
            iproc=i_isub
            iaproc=i_dsub
            i_acw=i_iacw
            i_acr=i_raci
            cf=cffields(k,i_cfi)
          case (4_iwp) !snow
            iproc=i_ssub
            iaproc=i_dssub
            i_acw=i_sacw
            i_acr=i_sacr
            cf=cffields(k,i_cfs)
          case (5_iwp) !graupel
            iproc=i_gsub
            iaproc=i_dgsub
            i_acw=i_gacw
            i_acr=i_gacr
            cf=cffields(k,i_cfg)
          end select
        end if
          
        if (mass > thresh_small(params%i_1m)) then ! if no existing ice, we don't grow/deplete it.
             
          if (params%l_2m) num=qfields(k, params%i_2m)
             
          n0=dist_n0(k,params%id)
          mu=dist_mu(k,params%id)
          lam=dist_lambda(k,params%id)
            
          if (l_gamma_online) then
            call ventilation_3M(ixy_inner, k, V_x, n0, lam, mu, params)
          else
            call ventilation_1M_2M(ixy_inner, k, V_x, n0, lam, mu, params)
          endif
             
          AB=1.0/(Ls*Ls/(Rv*ka*TdegK(k,ixy_inner)*TdegK(k,ixy_inner))*rho(k,ixy_inner)+1.0/(Dv*qis(k)))
          dmass=(qv/qis(k)-1.0)*V_x*AB *cf ! grid mean
             
          ! Include latent heat effects of collection of rain and cloud
          ! as done in Milbrandt & Yau (2005)
          if (l_latenteffects) then
            dmass=dmass - Lf*Ls/(Rv*ka*TdegK(k,ixy_inner)*TdegK(k,ixy_inner))                      &
                  *(procs(cloud_params%i_1m, i_acw%id)%column_data(k)          &
                  + procs(rain_params%i_1m, i_acr%id)%column_data(k))
          end if

          ! Check we haven't become subsaturated and limit if we have (dep only)
          ! NB doesn't account for simultaneous ice/snow growth - checked elsewhere
          if (dmass > 0.0) dmass=min((qv-qis(k))/dt,dmass)
          ! Check we don't remove too much (sub only)
          if (dmass < 0.0) dmass=max(-mass/dt,dmass)
             
          if (params%l_2m) then
            dnumber=0.0
            if (dmass < 0.0) dnumber=dmass*num/mass
          end if
             
          if (-dmass*dt >0.98*mass .or. (params%l_2m .and.                     &
              -dnumber*dt > 0.98*num)) then
            l_suball=.true.
            dmass=-mass/dt
            dnumber=-num/dt
          end if
             
          procs(i_qv, iproc%id)%column_data(k)=-dmass
          procs(params%i_1m, iproc%id)%column_data(k)=dmass
             
          if (params%l_2m) procs(params%i_2m,iproc%id)%column_data(k)=dnumber
             
          if (dmass < 0.0 .and. l_process) then ! Only process aerosol if sublimating
            if (iaproc%id==i_dsub%id) then
              dmad=dnumber*dustact(k)%mact1_mean*dustact(k)%nratio1
              dnumber_d=dnumber*dustact(k)%nratio1
            else if (iaproc%id==i_dssub%id) then
              dmad=dnumber*dustact(k)%mact2_mean*dustact(k)%nratio2
              dnumber_d=dnumber*dustact(k)%nratio2
            else if (iaproc%id==i_dgsub%id) then
              dmad=dnumber*dustact(k)%mact3_mean*dustact(k)%nratio3
              dnumber_d=dnumber*dustact(k)%nratio3
            end if

            !checking that mass change from
            !active insol in ice is negative
            !and larger in magnitude than epsilon
            !for sublimation of ice
            if (dmad < -epsilon(dmad)) then 
              aerosol_procs(i_am7, iaproc%id)%column_data(k)=dmad
              aerosol_procs(i_am6, iaproc%id)%column_data(k)=-dmad
              ! <WARNING: putting back in coarse mode
              if (l_passivenumbers_ice) then
                aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d 
              end if
              aerosol_procs(i_an6, iaproc%id)%column_data(k)=-dnumber_d
              ! <WARNING: putting back in coarse mode
            endif

            if (iaproc%id==i_dsub%id) then
              dmac=dnumber*aeroice(k)%mact1_mean*aeroice(k)%nratio1
              dmac=min(dmac,aeroice(k)%mact1/dt)
              dnumber_a=dnumber*aeroice(k)%nratio1
            else if (iaproc%id==i_dssub%id) then
              dmac=dnumber*aeroice(k)%mact2_mean*aeroice(k)%nratio2
              dmac=min(dmac,aeroice(k)%mact2/dt)
              dnumber_a=dnumber*aeroice(k)%nratio2
            else if (iaproc%id==i_dgsub%id) then
              dmac=dnumber*aeroice(k)%mact3_mean*aeroice(k)%nratio3
              dmac=min(dmac,aeroice(k)%mact3/dt)
              dnumber_a=dnumber*aeroice(k)%nratio3
            end if

            !checking that mass change from
            !active sol in ice is negative
            !and larger in magnitude than epsilon
            !for sublimation of ice
            if (dmac < -epsilon(dmac)) then 
              aerosol_procs(i_am8, iaproc%id)%column_data(k)=dmac
              aerosol_procs(i_am2, iaproc%id)%column_data(k)=-dmac
              ! <WARNING: putting back in accumulation mode                  
              if (l_passivenumbers) then
                aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
              end if
              aerosol_procs(i_an2, iaproc%id)%column_data(k)=-dnumber_a
              ! <WARNING: putting back in accumulation mode
            end if
          end if
        end if
      end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine idep
end module ice_deposition
