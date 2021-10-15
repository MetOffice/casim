module ice_multiplication
  use variable_precision, only: wp, iwp
  use process_routines, only: process_rate, process_name, i_gacw, i_sacw, i_ihal
  use passive_fields, only: TdegC
  use mphys_parameters, only: ice_params, snow_params, graupel_params, dN_hallet_mossop, M0_hallet_mossop
  use thresholds, only: thresh_small, cfliq_small
  use m3_incs, only: m3_inc_type2
  use mphys_switches, only: l_prf_cfrac, i_cfs, i_cfg, i_cfl, mpof

  implicit none

  character(len=*), parameter, private :: ModuleName='ICE_MULTIPLICATION'

contains
  !> Subroutine to determine the ice splintering by Hallet-Mossop
  !> This effect requires prior calculation of the accretion rate of
  !> graupel and snow.
  !> This is a source of ice number and mass and a sink of liquid
  !> (but this is done via the accretion processes already so is
  !> represented here as a sink of snow/graupel)
  !> For triple moment species there is a corresponding change in the
  !> 3rd moment assuming shape parameter is not changed
  !>
  !> AEROSOL: All aerosol sinks/sources are assumed to come from soluble modes
  !
  !> OPTIMISATION POSSIBILITIES:
  subroutine hallet_mossop(dt, nz, cffields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    real(wp), intent(in) :: cffields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! Local variables
    real(wp) :: gacw, sacw  ! accretion process rates
    real(wp) :: dnumber_s, dnumber_g  ! number conversion rate from snow/graupel
    real(wp) :: dmass_s, dmass_g      ! mass conversion rate from snow/graupel
    real(wp) :: Eff  !< splintering efficiency
    real(wp) :: cf_snow, cf_graupel, cf_liquid, overlap_cfsnow, overlap_cfgraupel

    integer :: k

    type(process_name) :: iproc ! processes selected depending on which species we're modifying

    character(len=*), parameter :: RoutineName='HALLET_MOSSOP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (.not. ice_params%l_2m) then
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      return
    end if
    
    do k = 1, nz
       if (TdegC(k) < 0.0_wp) then 
          if (l_prf_cfrac) then
             if (cffields(k,i_cfl) .gt. cfliq_small) then
                cf_liquid=cffields(k,i_cfl)
             else
                cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
             endif
             if (cffields(k,i_cfs) .gt. cfliq_small) then
                cf_snow=cffields(k,i_cfs)
             else
                cf_snow=cfliq_small !nonzero value - maybe move cf test higher up
             endif
             if (cffields(k,i_cfg) .gt. cfliq_small) then
                cf_graupel=cffields(k,i_cfg)
             else
                cf_graupel=cfliq_small !nonzero value - maybe move cf test higher up
             endif
          else
             cf_snow=1.0
             cf_graupel=1.0
             cf_liquid=1.0
          endif

          !use mixed-phase overlap function
          overlap_cfsnow=min(1.0,max(0.0,mpof*min(cf_liquid, cf_snow) +         &
               max(0.0,(1.0-mpof)*(cf_liquid+cf_snow-1.0))))
          overlap_cfgraupel=min(1.0,max(0.0,mpof*min(cf_liquid, cf_graupel) +   &
               max(0.0,(1.0-mpof)*(cf_liquid+cf_graupel-1.0))))

          Eff=1.0 - abs(TdegC(k) + 5.0)/2.5 ! linear increase between -2.5/-7.5 and -5C

          if (Eff > 0.0) then
             sacw=0.0
             gacw=0.0
             !! should use cf_overlap as in ice accretion
             if (snow_params%i_1m > 0) &
                  sacw=procs(snow_params%i_1m, i_sacw%id)%column_data(k)/overlap_cfsnow  !insnow process rate
             if (graupel_params%i_1m > 0) &
                  gacw=procs(graupel_params%i_1m, i_gacw%id)%column_data(k)/overlap_cfgraupel ! ingraupel process rate
             
             if ((sacw*overlap_cfsnow + gacw*overlap_cfgraupel)*dt > thresh_small(snow_params%i_1m)) then
                iproc=i_ihal

                dnumber_g=dN_hallet_mossop * Eff * (gacw) ! Number of splinters from graupel
                dnumber_s=dN_hallet_mossop * Eff * (sacw) ! Number of splinters from snow
                
                dnumber_g=min(dnumber_g, 0.5*gacw/M0_hallet_mossop) ! don't remove more than 50% of rimed liquid
                dnumber_s=min(dnumber_s, 0.5*sacw/M0_hallet_mossop) ! don't remove more than 50% of rimed liquid

                dmass_g=dnumber_g * M0_hallet_mossop * overlap_cfgraupel  ! convert back to grid mean
                dmass_s=dnumber_s * M0_hallet_mossop * overlap_cfsnow ! convert back to grid mean
                
                dnumber_g=dnumber_g * overlap_cfgraupel  ! convert back to grid mean
                dnumber_s=dnumber_s * overlap_cfsnow  ! convert back to grid mean
        

                !-------------------
                ! Sources for ice...
                !-------------------
                procs(ice_params%i_1m, iproc%id)%column_data(k)=dmass_g + dmass_s 
                procs(ice_params%i_2m, iproc%id)%column_data(k)=dnumber_g + dnumber_s
                
                !-------------------
                ! Sinks for snow...
                !-------------------
                if (sacw > 0.0) then
                   procs(snow_params%i_1m, iproc%id)%column_data(k)=-dmass_s
                   procs(snow_params%i_2m, iproc%id)%column_data(k)=0.0
                end if
                
                !---------------------
                ! Sinks for graupel...
                !---------------------
                if (gacw > 0.0) then
                   procs(graupel_params%i_1m, iproc%id)%column_data(k)=-dmass_g
                   procs(graupel_params%i_2m, iproc%id)%column_data(k)=0.0
                end if
                
             end if
          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine hallet_mossop
end module ice_multiplication
