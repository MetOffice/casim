module evaporation
  use variable_precision, only: wp
  use passive_fields, only: rho, qws, TdegK
! use mphys_switches, only: i_m3r, l_3mr
  use mphys_switches, only: i_qr, i_nr, i_qv, i_ql, l_2mr, i_am2, &
       i_an2, i_am3, i_an3, i_am4, i_am5, l_process, aero_index, &
       l_separate_rain, i_am6, i_an6, i_am9, l_warm, i_an11, i_an12, l_passivenumbers, &
       l_passivenumbers_ice, l_inhom_revp, l_gamma_online, &
       l_bypass_which_mode, iopt_which_mode
  use mphys_constants, only: Lv,ka, Dv, Rv
  use mphys_parameters, only: c_r, rain_params
  use process_routines, only: process_rate, i_prevp, i_arevp
  use thresholds, only: qr_small, ss_small, qr_tidy, ql_tidy
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use ventfac, only: ventilation_1M_2M, ventilation_3M
  use which_mode_to_use, only : which_mode
  use mphys_die, only: throw_mphys_error, incorrect_opt, std_msg
  
  implicit none

  character(len=*), parameter, private :: ModuleName='EVAPORATION'

  private

  public revp
contains

  subroutine revp(ixy_inner, dt, nz, qfields, cffields, aerophys, aerochem, aeroact, dustliq, procs, aerosol_procs, l_sigevap)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    USE mphys_switches,       ONLY: i_cfr

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: qfields(:,:)
    integer, intent(in) :: nz
    type(aerosol_phys), intent(in) :: aerophys(:)
    type(aerosol_chem), intent(in) :: aerochem(:)
    type(aerosol_active), intent(in) :: aeroact(:), dustliq(:)
    type(process_rate), intent(inout) :: procs(:,:)
    type(process_rate), intent(inout) :: aerosol_procs(:,:)
    logical, intent(out) :: l_sigevap(:) ! Determines if there is significant evaporation
    real(wp), intent(in) :: cffields(:,:)

    ! Local variables
    real(wp) :: dmass, dnumber, dnumber_a, dnumber_d
    real(wp) :: m1, m2, dm1
!   real(wp) :: m3, dm3

    real(wp) :: n0, lam, mu
    real(wp) :: V_r, AB

    real(wp) :: rain_mass
    real(wp) :: rain_number
!   real(wp) :: rain_m3
    real(wp) :: qv
    real(wp) :: cf

    logical :: l_rain_test ! conditional test on rain

    real(wp) :: dmac, dmac1, dmac2, dnac1, dnac2, dmacd
    integer :: k

    character(len=*), parameter :: RoutineName='REVP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    do k = 1, nz
       
       l_sigevap(k)=.false.
       
       qv=qfields(k, i_qv)
       rain_mass=qfields(k, i_qr)
       if (l_2mr)rain_number=qfields(k, i_nr)
       ! if (l_3mr)rain_m3=qfields(k, i_m3r)
       
       if (qv/qws(k,ixy_inner) < 1.0-ss_small .and. qfields(k, i_ql) < ql_tidy .and. rain_mass > qr_tidy) then
          
          m1=rain_mass/c_r
          if (l_2mr) m2=rain_number
          ! if (l_3mr) m3=rain_m3

          l_rain_test=.false.
          if (l_2mr) l_rain_test=rain_number>0
          if (rain_mass > qr_small .and. (.not. l_2mr .or. l_rain_test)) then
             n0=dist_n0(k,rain_params%id)
             mu=dist_mu(k,rain_params%id)
             lam=dist_lambda(k,rain_params%id)
             cf=cffields(k,i_cfr)

             if (l_gamma_online) then 
                call ventilation_3M(ixy_inner, k, V_r, n0, lam, mu, rain_params)
             else
                call ventilation_1M_2M(ixy_inner, k, V_r, n0, lam, mu, rain_params)
             endif

             AB=1.0/(Lv**2/(Rv*ka)*rho(k,ixy_inner)*TdegK(k,ixy_inner)**(-2)+1.0/(Dv*qws(k,ixy_inner)))
             dmass=(1.0-qv/qws(k,ixy_inner))*V_r*AB  *cf !grid mean
             
             dmass=MIN(dmass,rain_mass/dt)
           
          else
             dmass=rain_mass/dt
          end if

          dm1=dmass/c_r
          if (l_2mr) then
             dnumber=0.0
             if (l_inhom_revp) dnumber=dm1*m2/m1
          end if
          ! if (l_3mr) dm3=dm1*m3/m1

          if (l_2mr) then
             if (dnumber*dt > rain_number .or. dmass*dt >= rain_mass-qr_tidy) then
                dmass=rain_mass/dt
                dnumber=rain_number/dt
               ! if (l_3mr) dm3=m3/dt
             end if
          end if

          procs(i_qr, i_prevp%id)%column_data(k)=-dmass
          procs(i_qv, i_prevp%id)%column_data(k)=dmass

          if (dmass*dt/rain_mass > .8) l_sigevap(k)=.true.

          if (l_2mr) then
             procs(i_nr, i_prevp%id)%column_data(k)=-dnumber
          end if
          ! if (l_3mr) then
          !    procs(i_m3r, i_prevp%id)%column_data(k)=-dm3
          ! end if
     
      !============================
      ! aerosol processing
      !============================
          if (l_process .and. abs(dnumber) >0) then

             dmac=dnumber*aeroact(k)%nratio2*aeroact(k)%mact2_mean
             dmac=min(dmac,aeroact(k)%mact2/dt)
             if (l_separate_rain) then
                aerosol_procs(i_am5, i_arevp%id)%column_data(k)=-dmac
             else
                aerosol_procs(i_am4, i_arevp%id)%column_data(k)=-dmac
             end if
             if (l_passivenumbers) then
                dnumber_a=-dnumber*aeroact(k)%nratio2
                aerosol_procs(i_an11, i_arevp%id)%column_data(k)=dnumber_a
             end if

             if (l_passivenumbers_ice) then
                dnumber_d=-dnumber*dustliq(k)%nratio2
                aerosol_procs(i_an12, i_arevp%id)%column_data(k)=dnumber_d
             end if

             ! Return aerosol
             if (aero_index%i_accum >0 .and. aero_index%i_coarse >0) then
                ! Coarse and accumulation mode being used. Which one to return to?
                if (l_bypass_which_mode) then
                   !Don't use which_mode - just transfer all to either
                   !accum or coarse mode
                  if (iopt_which_mode.eq.1) then !All to accum
                       dmac1=dmac
                       dmac2=0.0
                       dnac1=dnumber_a
                       dnac2=0.0
                  else if (iopt_which_mode.eq.2) then !All to coarse
                       dmac1=0.0
                       dmac2=dmac
                       dnac1=0.0
                       dnac2=dnumber_a
                  else
                       write(std_msg, '(A)') "incorrect iopt_which_mode option selected"
                       call throw_mphys_error(incorrect_opt,ModuleName, std_msg)
                  endif
                else
                  call which_mode(dmac, dnumber*aeroact(k)%nratio2, aerophys(k)%rd(aero_index%i_accum), &
                     aerophys(k)%rd(aero_index%i_coarse), aerochem(k)%density(aero_index%i_accum),    &
                     aerophys(k)%sigma(aero_index%i_accum),                                           &
                     dmac1, dmac2, dnac1, dnac2)
                end if !end of bypass whichmode
                
                aerosol_procs(i_am2, i_arevp%id)%column_data(k)=dmac1
                aerosol_procs(i_an2, i_arevp%id)%column_data(k)=dnac1
                aerosol_procs(i_am3, i_arevp%id)%column_data(k)=dmac2
                aerosol_procs(i_an3, i_arevp%id)%column_data(k)=dnac2
             else
                if (aero_index%i_accum >0) then
                   aerosol_procs(i_am2, i_arevp%id)%column_data(k) = dmac
                   aerosol_procs(i_an2, i_arevp%id)%column_data(k) = dnumber *             &
                        aeroact(k)%nratio2
                end if
                if (aero_index%i_coarse >0) then
                   aerosol_procs(i_am3, i_arevp%id)%column_data(k) = dmac
                   aerosol_procs(i_an3, i_arevp%id)%column_data(k) = dnumber *             &
                        aeroact(k)%nratio2
                end if
             end if

             dmacd=dnumber*dustliq(k)%nratio2*dustliq(k)%mact2_mean
             if (.not. l_warm .and. dmacd /=0.0) then
                aerosol_procs(i_am9, i_arevp%id)%column_data(k)=-dmacd
                aerosol_procs(i_am6, i_arevp%id)%column_data(k)=dmacd
                aerosol_procs(i_an6, i_arevp%id)%column_data(k)=dnumber*dustliq(k)%nratio2
             end if

          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine revp
end module evaporation
