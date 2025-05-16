module ice_multiplication
  use variable_precision, only: wp, iwp
  use process_routines, only: process_rate, process_name, i_gacw, i_sacw,      &
                                 i_ihal, i_idps, i_iics, i_gaci, i_gacs, i_homr
  use passive_fields, only: TdegC, TdegK
  use mphys_parameters, only: ice_params, snow_params, graupel_params,         &
                       dN_hallet_mossop, M0_hallet_mossop, dN_droplet_shatter, &
                       P_droplet_shatter, coef_ice_breakup
  use thresholds, only: thresh_small, cfliq_small
  use m3_incs, only: m3_inc_type2
  use mphys_switches, only: l_prf_cfrac, i_cfs, i_cfg, i_cfl, mpof, i_qg,      &
                            i_ng, i_cfi
  use mphys_constants, only: pi
  use casim_stph, only: l_rp2_casim, mpof_casim_rp

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
  !> Subroutine to determine the droplet shattering
  !> This effect requires prior calculation of the Bigg freezing
  !> rate for raindrops.
  !> This is a source of ice number and mass, and a sink of graupel
  !> number and mass (sink for snow number and mass already done in 
  !> the Bigg's raindrop freezing).
  !>
  !> Subroutine to determine the ice-ice collision
  !> This effect requires prior calculation of collision tendency
  !> rates for graupel-snow accretion and snow collecting snow.
  !> This is a source of snow and ice number concentration, and not 
  !> changing the graupel mass and number.
  !>
  !> AEROSOL: All aerosol sinks/sources are assumed to come from soluble modes
  !
  !> OPTIMISATION POSSIBILITIES:
  subroutine hallet_mossop(ixy_inner, dt, nz, cffields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
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

    ! Apply RP scheme
    if ( l_rp2_casim ) then
        mpof = mpof_casim_rp
    endif

    if (.not. ice_params%l_2m) then
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      return
    end if
    
    do k = 1, nz
       if (TdegC(k,ixy_inner) < 0.0_wp) then 
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

          Eff=1.0 - abs(TdegC(k,ixy_inner) + 5.0)/2.5 ! linear increase between -2.5/-7.5 and -5C

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

!--------------------------------------------------------------------------------
!> Subroutine for droplet shattering (Sullivan et al., 2018)
!--------------------------------------------------------------------------------
  subroutine droplet_shattering(ixy_inner, dt, nz, cffields, qfields, procs)

   USE yomhook, ONLY: lhook, dr_hook
   USE parkind1, ONLY: jprb, jpim

   implicit none

   ! Subroutine arguments
   integer, intent(in) :: ixy_inner
   real(wp), intent(in) :: dt
   integer, intent(in) :: nz
   real(wp), intent(in) :: cffields(:,:)
   real(wp), intent(in) :: qfields(:,:)
   type(process_rate), intent(inout), target :: procs(:,:)

   ! Local variables
   real(wp) :: homr_mass, homr_number  ! rate of homogeneous freezing of rain (Bigg, 1953)
   real(wp) :: dnumber_i, dnumber_g  ! number conversion rate for ice crystal and graupel
   real(wp) :: dmass_i, dmass_g      ! mass conversion rate for ice crystal and graupel
   real(wp) :: prob_DS ! temperature-dependent shattering probability
   real(wp) :: cf_graupel, graupel_mass, graupel_number

   integer :: k

   type(process_name) :: iproc ! processes selected depending on which species we're modifying

   character(len=*), parameter :: RoutineName='DROPLET_SHATTERING'

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
      if (TdegC(k,ixy_inner) < 0.0_wp) then 
         if (l_prf_cfrac) then
            if (cffields(k,i_cfg) .gt. cfliq_small) then
               cf_graupel=cffields(k,i_cfg)
            else
               cf_graupel=cfliq_small !nonzero value - maybe move cf test higher up
            endif
         else
            cf_graupel=1.0
         endif
         
         ! Shattering Probability:
         ! Normal distribution centred at 258 K with standard deviation of 3 K
         ! probablity of droplet shatter P_Droplet_shatter defaults to 0.2
         ! maximum of distribution is 0.13298 (Sullivan (2018))
         prob_DS= (P_droplet_shatter / 0.13298) * (1 / (SQRT(2 * pi) * 3))          &
                 * EXP((-(TdegK(k,ixy_inner) - 258)**2) / (18))
         
         if (prob_DS > 0.0) then
            homr_mass = 0.0   ! mass tendency of raindrop frozen
            homr_number = 0.0   ! number tendency of raindrop frozen

            graupel_mass=qfields(k,i_qg)
            graupel_number=qfields(k,i_ng)
            
            if (graupel_params%i_1m > 0) &
                 homr_mass=procs(graupel_params%i_1m, i_homr%id)%column_data(k)/cf_graupel ! ingraupel process rate (kg / kg-1)

            if (graupel_params%i_2m > 0) &
                 homr_number=procs(graupel_params%i_2m, i_homr%id)%column_data(k)/cf_graupel ! ingraupel process rate (number / kg-1)

            if ((homr_mass*cf_graupel)*dt > thresh_small(graupel_params%i_1m)  &
                .and. (homr_number*cf_graupel)*dt > thresh_small(graupel_params%i_2m)) then 
               
               dnumber_i=(1 + prob_DS * dN_droplet_shatter) * homr_number ! Number of splinters from graupel
               ! No more than 50% of the graupels created from the frozen raindrops
               dmass_i=min(dnumber_i * M0_hallet_mossop,0.5*homr_mass) 

               dnumber_i=dnumber_i * cf_graupel ! Convert back to grid-box mean
               dmass_i=dmass_i * cf_graupel

               dmass_g=-dmass_i
               dnumber_g=dmass_g * graupel_number / graupel_mass 

               if (homr_mass > 0.0) then
                  iproc = i_idps
                  !-------------------
                  ! Sources for ice...
                  !-------------------
                  procs(ice_params%i_1m, iproc%id)%column_data(k)=dmass_i
                  procs(ice_params%i_2m, iproc%id)%column_data(k)=dnumber_i
                  
                  !---------------------
                  ! Sinks for graupel...
                  !---------------------
                  procs(graupel_params%i_1m, iproc%id)%column_data(k)=dmass_g
                  procs(graupel_params%i_2m, iproc%id)%column_data(k)=dnumber_g
               end if
            end if
         end if
      end if
   enddo

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

 end subroutine droplet_shattering

!-----------------------------------------------------------------------------------
!> Subroutine for ice-ice collision

 subroutine ice_collision(ixy_inner, dt, nz, cffields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner    
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    real(wp), intent(in) :: cffields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! Local variables
    real(wp) :: gaci, gacs  ! accretion process rates
    real(wp) :: dnumber_i, dnumber_s ! number tendency for the collided hydrometeors
    real(wp) :: BR_fragments ! Temperature-dependent fragments from ice-ice collision
    real(wp) :: cf_snow, cf_graupel, cf_ice, overlap_cfsg, overlap_cfig

    integer :: k

    type(process_name) :: iproc ! processes selected depending on which species we're modifying

    character(len=*), parameter :: RoutineName='ICE_COLLISION'

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
       if (TdegC(k,ixy_inner) < 0.0_wp) then 
          if (l_prf_cfrac) then
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
             if (cffields(k,i_cfi) .gt. cfliq_small) then
               cf_ice=cffields(k,i_cfi)
            else
               cf_ice=cfliq_small !nonzero value - maybe move cf test higher up
            endif
          else
             cf_snow=1.0
             cf_graupel=1.0
             cf_ice=1.0
          endif
         
          overlap_cfsg = min(cf_snow, cf_graupel)
          overlap_cfig = min(cf_ice, cf_graupel)
           
          ! Number of fragments generated based on Takahashi et al., (1995)
          BR_fragments=coef_ice_breakup * ((TdegK(k,ixy_inner) - 252) ** 1.2)      &
                      * EXP(-(TdegK(k,ixy_inner) - 252)/5)

          if (BR_fragments > 0.0) then
             gacs=0.0
             gaci=0.0
             
             if (snow_params%i_2m > 0) &
                  gacs=-procs(snow_params%i_2m, i_gacs%id)%column_data(k)/overlap_cfsg  ! insnow process rate
             if (graupel_params%i_2m > 0) &
                  gaci=-procs(ice_params%i_2m, i_gaci%id)%column_data(k)/overlap_cfig   ! inice process rate
             
             if ((gacs*cf_snow)*dt > thresh_small(snow_params%i_2m)) then
                iproc=i_iics

                dnumber_s=BR_fragments * (gacs) * overlap_cfsg   ! Number of splinters from graupel and convert back to grid mean
                !-------------------
                ! Sources for snow...
                !-------------------
                procs(snow_params%i_2m, iproc%id)%column_data(k)=dnumber_s                
             end if

             if ((gaci*cf_ice)*dt > thresh_small(ice_params%i_2m)) then
                iproc=i_iics

                dnumber_i=BR_fragments * (gaci) * overlap_cfig   ! Number of splinters from graupel and convert back to grid mean
                !-------------------
                ! Sources for ice...
                !-------------------
                procs(ice_params%i_2m, iproc%id)%column_data(k)=dnumber_i
             endif

          end if
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine ice_collision

end module ice_multiplication
