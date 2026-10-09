module ice_multiplication
  use variable_precision, only: wp, iwp
  use process_routines, only: process_rate, process_name, i_gacw, i_sacw,      &
                                 i_ihal, i_idps, i_iics, i_gaci, i_gacs, i_homr, &
                                 i_imo1, i_imo2, i_iicb_i, i_iicb_s, i_iicb_g,   &
                                 i_gacr, i_sacr, i_raci
  use passive_fields, only: TdegC, TdegK, rho
  use mphys_parameters, only: ice_params, snow_params, graupel_params,         &
                       dN_hallet_mossop, M0_hallet_mossop, dN_droplet_shatter, &
                       P_droplet_shatter, coef_ice_breakup, hydro_params,      &
                       rain_params, pthreshr, pthreshi, pthreshs, pthreshg,    &
                       gam1r, gam2r
  use thresholds, only: thresh_small, cfliq_small
  use m3_incs, only: m3_inc_type2
  use mphys_switches, only: l_prf_cfrac, i_cfs, i_cfg, i_cfl, mpof, i_qg,      &
                            i_ng, i_cfi, i_cfr, i_qr, i_qi, i_qs, i_nr
  use mphys_constants, only: pi, rho0, rhow, rhoi, ttr, Lf, Cwater
  use casim_stph, only: l_rp2_casim, mpof_casim_rp
  use special, only: GammaFunc
  use distributions, only: dist_lambda, dist_mu, dist_n0

  implicit none

  !----------------------------------------------------------------------------
  ! Constants for the Phillips et al. secondary ice production schemes
  ! (sip_phillips_mode1, sip_phillips_mode2, sip_phillips_breakup)
  !----------------------------------------------------------------------------
  real(wp), parameter :: eri = 1.0_wp           !< Collision efficiency
  real(wp), parameter :: gamma_liq = 0.072_wp   !< Surface tension of water (J m-2)
  real(wp), parameter :: decrit = 0.2_wp        !< Critical dimensionless energy
                                                !! for splashing (Phillips et al. 2018)
  real(wp), parameter :: phi_mode2 = 0.35_wp    !< Fraction of Mode 2 splash
                                                !! fragments that freeze
  real(wp), parameter :: phi_phillips = 3.5e-3_wp !< Phillips et al. (2017)
                                                !! graupel-graupel coefficient
  real(wp), parameter :: dtt = 10.0e-6_wp       !< Diameter of tiny Mode 1
                                                !! fragments (m)
  real(wp), parameter :: oneoversix = 1.0_wp/6.0_wp
  real(wp), parameter :: oneoverthree = 1.0_wp/3.0_wp
  real(wp), parameter :: oneovernine = 1.0_wp/9.0_wp

  ! 10-point Gauss-Legendre quadrature used for the collision integrals
  ! (DLMF 3.5).
  ! Positive abscissae and weights of the 10-point Gauss-Legendre rule on
  ! [-1,1] (the rule is symmetric about zero).
  integer, parameter, private :: n_half = 5
  real(wp), parameter, private :: gl_node(n_half) = (/                                  &
       0.14887433898163122_wp, 0.43339539412924720_wp, 0.67940956829902440_wp, &
       0.86506336668898450_wp, 0.97390652851717170_wp /)
  real(wp), parameter, private :: gl_weight(n_half) = (/                                &
       0.29552422471475280_wp, 0.26926671930999650_wp, 0.21908636251598200_wp, &
       0.14945134915058040_wp, 0.06667134430868814_wp /)

  abstract interface
    !> Integrand of one variable, evaluated at a vector of abscissae
    function integrand_1d(x) result(f)
      import :: wp
      real(wp), intent(in) :: x(:)
      real(wp) :: f(size(x))
    end function integrand_1d

    !> Integrand of two variables: scalar outer coordinate x and a vector
    !> of inner abscissae y
    function integrand_2d(x, y) result(f)
      import :: wp
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y(:)
      real(wp) :: f(size(y))
    end function integrand_2d

    !> Inner integration limit as a function of the outer coordinate
    function limit_function(x) result(y)
      import :: wp
      real(wp), intent(in) :: x
      real(wp) :: y
    end function limit_function
  end interface

  ! Collisional-breakup pair identifiers.  CB_XY means category X is the
  ! fragmenting particle and Y its collision partner.
  integer, parameter :: CB_II       = 1
  integer, parameter :: CB_IS       = 2
  integer, parameter :: CB_IG       = 3
  integer, parameter :: CB_SI       = 4
  integer, parameter :: CB_SS       = 5
  integer, parameter :: CB_SG       = 6
  integer, parameter :: CB_GG_SMALL = 7
  integer, parameter :: CB_GG_HAIL  = 8

  ! Working state shared between the SIP routines and their integrands.
  ! Thread-private so that columns can be processed in parallel.
  type(hydro_params) :: params_mod, params_mod1
  real(wp) :: rho_mod, n0r, alpha_r, lambda0r, n0i, alpha_i, lambda0i, f_mode2
  real(wp) :: mrthresh, mrupper, mithresh, miupper, milower
  real(wp) :: mifragthresh, mifragupper
  real(wp) :: t_send, lam_freeze, n0_freeze
  integer :: type1_send
!$OMP THREADPRIVATE(params_mod, params_mod1, rho_mod, n0r, alpha_r, lambda0r)
!$OMP THREADPRIVATE(n0i, alpha_i, lambda0i, f_mode2, mrthresh, mrupper)
!$OMP THREADPRIVATE(mithresh, miupper, milower, mifragthresh, mifragupper)
!$OMP THREADPRIVATE(t_send, lam_freeze, n0_freeze, type1_send)

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


!--------------------------------------------------------------------------------
!> Secondary ice production following Phillips et al. (2017, 2018), with
!> collision integrals evaluated numerically over the gamma size
!> distributions (Sun et al., 2025, ACP, doi:10.5194/acp-25-18549-2025):
!>
!>  sip_phillips_mode1   - fragmentation of freezing raindrops (Mode 1):
!>                         spherical freezing of rain, and collisions of rain
!>                         with less massive ice, snow and graupel
!>  sip_phillips_mode2   - fragmentation in collisions of supercooled rain
!>                         with more massive ice, snow and graupel (Mode 2)
!>  sip_phillips_breakup - ice-ice collisional breakup (Phillips et al. 2017)
!>
!> These are alternatives to droplet_shattering and ice_collision above and
!> are switched on with l_sip_phillips_mode1, l_sip_phillips_mode2 and
!> l_sip_phillips_breakup (all default .false.).
!>
!> Contributed by the University of Manchester (P. J. Connolly, M. Sun,
!> B. Z. Portman, R. L. James) under the Horizon Europe project CERTAINTY,
!> grant agreement 101137680.
!--------------------------------------------------------------------------------
    subroutine sip_phillips_mode1(ixy_inner, dt, nz, cffields, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    real(wp), intent(in) :: cffields(:,:)
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! Local variables
    real(wp) :: gacr, sacr, raci1, raci2  ! accretion process rates
    real(wp) :: dnumber_s, dnumber_g,  &
                dnumber_i ! number conversion rate from snow/graupel
    real(wp) :: dmass_s, dmass_g, &
                dmass_i      ! mass conversion rate from snow/graupel
    real(wp) :: cf_snow, cf_graupel, cf_liquid, cf_ice, cf_rain, overlap_cfsnow, &
         overlap_cfgraupel, overlap_cfice, ice_mass, rain_mass, rain_number, m0

    type(process_name) :: iproc ! processes selected depending on which species we're modifying

    integer :: k
    real(wp) :: arg, dummy3, nfreeze, mass_freeze,nfrag_nucc, mfrag_nucc


    character(len=*), parameter :: RoutineName='SIP_PHILLIPS_MODE1'

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

    ! loop over all levels - exactly like the HM routine above
    do k = 1, nz
       if (TdegC(k,ixy_inner) < 0.0_wp) then
          if (l_prf_cfrac) then
             if (cffields(k,i_cfl) .gt. cfliq_small) then
                cf_liquid=cffields(k,i_cfl)
             else
                cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
             endif
             if (cffields(k,i_cfr) .gt. cfliq_small) then
                cf_rain=cffields(k,i_cfr)
             else
                cf_rain=cfliq_small !nonzero value - maybe move cf test higher up
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
             ! added for ice
             if (cffields(k,i_cfi) .gt. cfliq_small) then
                cf_ice=cffields(k,i_cfi)
             else
                cf_ice=cfliq_small !nonzero value - maybe move cf test higher up
             endif
          else
             cf_rain=1.0
             cf_snow=1.0
             cf_graupel=1.0
             cf_liquid=1.0
             cf_ice=1.0
          endif

          !use mixed-phase overlap function
          ! PJC changed to be for rain and added ice
          overlap_cfsnow=min(1.0,max(0.0,mpof*min(cf_rain, cf_snow) +         &
               max(0.0,(1.0-mpof)*(cf_rain+cf_snow-1.0))))
          overlap_cfgraupel=min(1.0,max(0.0,mpof*min(cf_rain, cf_graupel) +   &
               max(0.0,(1.0-mpof)*(cf_rain+cf_graupel-1.0))))
          overlap_cfice=min(1.0,max(0.0,mpof*min(cf_rain, cf_ice) +   &
               max(0.0,(1.0-mpof)*(cf_rain+cf_ice-1.0))))

          sacr=0.0
          gacr=0.0
          raci1=0.0
          raci2=0.0
          if (snow_params%i_1m > 0) then
            sacr=procs(snow_params%i_1m, i_sacr%id)%column_data(k)/overlap_cfsnow  !insnow process rate
            raci1=procs(snow_params%i_1m, i_raci%id)%column_data(k)/overlap_cfsnow ! ingraupel process rate
          endif
          if (graupel_params%i_1m > 0) then
            gacr=procs(graupel_params%i_1m, i_gacr%id)%column_data(k)/overlap_cfgraupel ! ingraupel process rate
            raci2=procs(graupel_params%i_1m, i_raci%id)%column_data(k)/overlap_cfgraupel
          endif


          rho_mod=rho(k,ixy_inner)
          alpha_r=dist_mu(k,rain_params%id)
          lambda0r=dist_lambda(k,rain_params%id)
          arg = 1.0_wp+alpha_r
          n0r=dist_n0(k,rain_params%id)*lambda0r**(arg) / &
            GammaFunc(arg)

          ! for use within module, ICE first++++++++++++++++++++++++++++++++++++++++++++++
          params_mod = ice_params
          alpha_i=dist_mu(k,params_mod%id)
          lambda0i=dist_lambda(k,params_mod%id)
          arg = 1.0_wp+alpha_i
          n0i=dist_n0(k,params_mod%id)*lambda0i**(arg) / &
            GammaFunc(arg)
          t_send = TdegK(k,ixy_inner)

          ! calculate the increase in ice crystal number
          mrthresh = rain_params%c_x*1.e-6_wp**rain_params%d_x
          milower   = params_mod%c_x*1.e-6_wp**params_mod%d_x
          mrthresh  = max(mrthresh,milower)
          mrupper   = rain_params%c_x*(pthreshr/lambda0r)**rain_params%d_x
          miupper   = params_mod%c_x*(pthreshi/lambda0i)**params_mod%d_x
          rain_mass = qfields(k,i_qr)
          rain_number = qfields(k,i_nr)
          ice_mass = qfields(k,i_qi)
          dnumber_i = 0.0
          dmass_i = 0.0
          if((mrupper.gt.mrthresh).and.(rain_mass.gt.thresh_small(rain_params%i_1m)) .and. &
            (ice_mass.gt.thresh_small(params_mod%i_1m))) then
            dummy3 = gl_quad_2d(dintegral_mode1, limit1_mode1, limit2_mode1, &
                mrthresh, mrupper)

            ! it is a rate
            dnumber_i = rho_mod*max(dummy3,0.0)

          endif
          !END ICE------------------------------------------------------------------------

          ! for use within module, SNOW ++++++++++++++++++++++++++++++++++++++++++++++++++
          params_mod = snow_params
          alpha_i=dist_mu(k,params_mod%id)
          lambda0i=dist_lambda(k,params_mod%id)
          arg = 1.0_wp+alpha_i
          n0i=dist_n0(k,params_mod%id)*lambda0i**(arg) / &
            GammaFunc(arg)

          ! calculate the increase in ice crystal number
          mrthresh = rain_params%c_x*1.e-6_wp**rain_params%d_x
          milower   = params_mod%c_x*1.e-6_wp**params_mod%d_x
          mrthresh  = max(mrthresh,milower)
          mrupper   = rain_params%c_x*(pthreshr/lambda0r)**rain_params%d_x
          miupper   = params_mod%c_x*(pthreshs/lambda0i)**params_mod%d_x
          ice_mass = qfields(k,i_qs)
          dnumber_s = 0.0
          dmass_s = 0.0
          if((mrupper.gt.mrthresh).and.(rain_mass.gt.thresh_small(rain_params%i_1m)) .and. &
            (ice_mass.gt.thresh_small(params_mod%i_1m))) then
            dummy3 = gl_quad_2d(dintegral_mode1, limit1_mode1, limit2_mode1, &
                mrthresh, mrupper)
            ! it is a rate
            dnumber_s = rho_mod*max(dummy3,0.0)

            dummy3 = gl_quad_2d(dintegral_mode1_mass, limit1_mode1, limit2_mode1, &
                mrthresh, mrupper)
            dmass_s = rho_mod*max(dummy3,0.0)
            ! calculate the average mass of new splinters
            m0 = dmass_s /max(dnumber_s,1.0e-7)
            ! less than 50% of rain mass accreted onto snow should form splinters
            dmass_s = min(dmass_s,0.5*sacr)
            ! re-scale
            dnumber_s = dmass_s / max(m0, M0_hallet_mossop)
          endif
          !END SNOW-----------------------------------------------------------------------

          ! for use within module, GRAUPEL++++++++++++++++++++++++++++++++++++++++++++++++
          params_mod = graupel_params
          alpha_i=dist_mu(k,params_mod%id)
          lambda0i=dist_lambda(k,params_mod%id)
          arg = 1.0_wp+alpha_i
          n0i=dist_n0(k,params_mod%id)*lambda0i**(arg) / &
            GammaFunc(arg)

          ! calculate the increase in ice crystal number
          mrthresh = rain_params%c_x*1.e-6_wp**rain_params%d_x
          milower   = params_mod%c_x*1.e-6_wp**params_mod%d_x
          mrthresh  = max(mrthresh,milower)
          mrupper   = rain_params%c_x*(pthreshr/lambda0r)**rain_params%d_x
          miupper   = params_mod%c_x*(pthreshg/lambda0i)**params_mod%d_x
          ice_mass = qfields(k,i_qg)
          dnumber_g = 0.0
          dmass_g = 0.0
          if((mrupper.gt.mrthresh).and.(rain_mass.gt.thresh_small(rain_params%i_1m)) .and. &
            (ice_mass.gt.thresh_small(params_mod%i_1m))) then
            dummy3 = gl_quad_2d(dintegral_mode1, limit1_mode1, limit2_mode1, &
                mrthresh, mrupper)

            ! it is a rate
            dnumber_g = rho_mod*max(dummy3,0.0)

            dummy3 = gl_quad_2d(dintegral_mode1_mass, limit1_mode1, limit2_mode1, &
                mrthresh, mrupper)
            dmass_g = rho_mod*max(dummy3,0.0)
            ! calculate the average mass of new splinters
            m0 = dmass_g /max(dnumber_g,1.0e-7)
            ! less than 50% of rain mass accreted onto graupel should form splinters
            dmass_g = min(dmass_g,0.5*gacr)
            ! re-scale
            dnumber_g = dmass_g / max(m0,M0_hallet_mossop)
          endif
          !END GRAUPEL--------------------------------------------------------------------


          ! now for the rain drops that are freezing++++++++++++++++++++++++++++++++++++++
              mrthresh=rain_params%c_x*1.e-6_wp**rain_params%d_x
          nfrag_nucc=0.0
          mfrag_nucc=0.0
          nfreeze=0.0_wp
          mass_freeze=0.0_wp
          if (graupel_params%l_2m) then
             ! number and mass of rain frozen to form graupel
             nfreeze=procs(graupel_params%i_2m, i_homr%id)%column_data(k)*dt
             mass_freeze=procs(graupel_params%i_1m, i_homr%id)%column_data(k)*dt
          end if

          if((mass_freeze.gt.thresh_small(ice_params%i_1m)).and. &
            (nfreeze.gt.thresh_small(ice_params%i_2m))) then
            arg=1._wp+alpha_r
            lam_freeze = (nfreeze/mass_freeze*gam2r /  &
                gam1r * rain_params%c_x)**(1.0/rain_params%d_x)
            n0_freeze = nfreeze / gam1r*lam_freeze**(arg)
            mrupper  = rain_params%c_x*(pthreshr/lam_freeze)**rain_params%d_x
            if(mrupper.gt.mrthresh) then

                ! multiplication according to mode-1
                nfrag_nucc=max(gl_quad_1d(integral_m1,mrthresh,mrupper)/dt,0.0)
                mfrag_nucc = max(gl_quad_1d(integral_m1m,mrthresh,mrupper)/dt,0.0)
                ! calculate the average mass of new splinters
                m0 = mfrag_nucc /max(nfrag_nucc,1.0e-7)
                ! less than 50% of rain mass frozen should form splinters
                mfrag_nucc = min(mfrag_nucc,0.5*mass_freeze/dt)
                ! re-scale
                nfrag_nucc = mfrag_nucc / max(m0, M0_hallet_mossop)
            endif
          endif
          ! end raindrop freezing---------------------------------------------------------
            if ((mfrag_nucc+raci1*overlap_cfsnow +raci2*overlap_cfgraupel + sacr*overlap_cfsnow + &
                gacr*overlap_cfgraupel)*dt > thresh_small(snow_params%i_1m)) then
                iproc=i_imo1

                dmass_i=dmass_i * overlap_cfice  ! convert back to grid mean
                dmass_g=dmass_g * overlap_cfgraupel  ! convert back to grid mean
                dmass_s=dmass_s * overlap_cfsnow ! convert back to grid mean

                dnumber_i=dnumber_i * overlap_cfice + nfrag_nucc  ! convert back to grid mean
                dnumber_g=dnumber_g * overlap_cfgraupel  ! convert back to grid mean
                dnumber_s=dnumber_s * overlap_cfsnow  ! convert back to grid mean

                !-------------------
                ! Sources for ice...
                !-------------------
                procs(ice_params%i_1m, iproc%id)%column_data(k)=dmass_g + dmass_s + &
                    mfrag_nucc
                procs(ice_params%i_2m, iproc%id)%column_data(k)=  &
                    dnumber_g + dnumber_s + dnumber_i

                !-------------------
                ! Sinks for snow...
                !-------------------
                if (sacr > 0.0) then
                   procs(snow_params%i_1m, iproc%id)%column_data(k)=-dmass_s
                   if (snow_params%l_2m) procs(snow_params%i_2m,iproc%id)%column_data(k)=0.0
                end if

                !---------------------
                ! Sinks for graupel...
                !---------------------
                procs(graupel_params%i_1m, iproc%id)%column_data(k)=-mfrag_nucc
                if (gacr > 0.0) then
                   procs(graupel_params%i_1m, iproc%id)%column_data(k)= &
                    procs(graupel_params%i_1m, iproc%id)%column_data(k)-dmass_g
                   if (graupel_params%l_2m) procs(graupel_params%i_2m,iproc%id)%column_data(k)=0.0
                end if

            endif
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine sip_phillips_mode1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mode 1 fragmentation integral over size distribution                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Paul J. Connolly, The University of Manchester
    !>@brief
    !>calculates the number of fragments and their mass
    function integral_m1(x)
        implicit none
        real(wp), dimension(:), intent(in) :: x
        real(wp), dimension(size(x)) :: integral_m1

        real(wp), dimension(size(x)) :: nfrag, diam, jac
        real(wp) :: n,nt,nb,mb,mt
        integer :: i

        ! x is particle mass.  The reconstructed freezing PSD is expressed
        ! in diameter space, so convert D(m) and include dD/dm.
        diam = (x/rain_params%c_x)**(1.0_wp/rain_params%d_x)
        jac  = diam**(1.0_wp-rain_params%d_x) / &
               (rain_params%c_x*rain_params%d_x)

        do i=1,size(x)
            call calculate_mode1(x(i),0._wp,t_send,n,nt,nb,mb,mt)
            nfrag(i)=n
        enddo

        integral_m1=n0_freeze*exp(-lam_freeze*diam)*diam**alpha_r*jac*nfrag
    end function integral_m1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mode 1 fragmentation integral over size distribution                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Paul J. Connolly, The University of Manchester
    !>@brief
    !>calculates the number of fragments and their mass
    function integral_m1m(x)
        implicit none
        real(wp), dimension(:), intent(in) :: x
        real(wp), dimension(size(x)) :: integral_m1m

        real(wp), dimension(size(x)) :: mfrag, diam, jac
        real(wp) :: n,nt,nb,mb,mt
        integer :: i

        ! x is particle mass.  Convert the diameter-space PSD to mass space.
        diam = (x/rain_params%c_x)**(1.0_wp/rain_params%d_x)
        jac  = diam**(1.0_wp-rain_params%d_x) / &
               (rain_params%c_x*rain_params%d_x)

        do i=1,size(x)
            call calculate_mode1(x(i),0._wp,t_send,n,nt,nb,mb,mt)
            mfrag(i)=(nt*mt+nb*mb)
        enddo

        integral_m1m=n0_freeze*exp(-lam_freeze*diam)*diam**alpha_r*jac*mfrag
    end function integral_m1m
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mode 1 fragmentation                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Paul J. Connolly, The University of Manchester
    !>@brief
    !>calculates the number of fragments and their mass
    subroutine calculate_mode1(min1,min2,t,n,nt,nb,mb,mt)
        implicit none
        real(wp), intent(in) :: min1,min2,t
        real(wp), intent(inout) :: n, nt, nb, mb, mt
        real(wp) :: tc, dthresh, x, beta1,log10zeta, log10nabla, t0, zetab, nablab, tb0, &
            sigma, omega, m,d, fac1

        if((min2>min1).or.(min1<=6.55e-11_wp)) then
            ! the ice is more massive than the drop or drop small, don't do it
            n=0._wp
            nt=0._wp
            nb=0._wp
            mb=0._wp
            mt=0._wp
            return
        endif

        d = (6._wp*min1/(pi*rhow))**(oneoverthree)
        tc=t-ttr
        dthresh = min(d,1.6e-3)
        x = log10(dthresh*1000._wp)

        ! table 3, phillips et al.
        beta1 = 0.
        log10zeta = 2.4268_wp*x*x*x + 3.3274_wp*x*x + 2.0783_wp*x + 1.2927_wp
        log10nabla = 0.1242_wp*x*x*x - 0.2316_wp*x*x - 0.9874_wp*x - 0.0827_wp
        t0 = -1.3999_wp*x*x*x - 5.3285_wp*x*x - 3.9847_wp*x - 15.0332_wp

        ! table 4, phillips et al.
        zetab = -0.4651_wp*x*x*x - 1.1072_wp*x*x - 0.4539_wp*x+0.5137_wp
        nablab = 28.5888*x*x*x + 49.8504_wp*x*x + 22.4873_wp*x + 8.0481_wp
        tb0 = 13.3588_wp*x*x*x + 15.7432_wp*x*x - 2.6545_wp*x - 18.4875_wp

        sigma = min(max((d-50.e-6_wp)/10.e-6_wp,0._wp), 1._wp)
        omega = min(max((-3._wp-tc)/3._wp,0._wp),1._wp)

        n = sigma*omega*(10._wp**log10zeta *(10**log10nabla)**2) / &
            ((tc-t0)**2+(10._wp**log10nabla)**2+beta1*tc)

        ! total number of fragments
        n=n*d/dthresh
        ! number of large fragments
        nb = min(sigma*omega*(zetab*nablab**2/((tc-tb0)**2+nablab**2)),n)
        ! number of small fragments
        nt = n-nb

        m=oneoversix*rhow*pi*d**3

        ! mass of large fragments
        mb=0.4_wp*m

        ! mass of small fragments
        mt=oneoversix*rhoi*pi*dtt**3

        if ((mt*nt+mb*nb) > 0.0_wp) then
            fac1=min(min1/(mt*nt+mb*nb),1._wp)
        else
            fac1=1.0_wp
        endif
        nt = nt *fac1
        nb = nb *fac1
        n=nt+nb

    end subroutine calculate_mode1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for mode 1 ice multiplication                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_mode1(x,y)
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_mode1
        real(wp) :: diamr, mr, vr, n,nt,nb,mb,mt
        real(wp), dimension(size(y)) :: mi, diami, delv, vi
        integer :: i

        mr=x
        mi=y
        diamr=(mr/rain_params%c_x)**(1.0_wp/rain_params%d_x)
        diami=(mi/params_mod%c_x)**(1.0_wp/params_mod%d_x)
        ! fall-speeds
        ! fall-speed of rain
        vr=(rain_params%a_x*diamr**rain_params%b_x*exp(-rain_params%f_x*diamr) + &
         rain_params%a2_x*diamr**rain_params%b2_x*exp(-rain_params%f2_x*diamr)) * &
         (rho0/rho_mod)**rain_params%g_x
        ! fall-speed of ice
        vi=params_mod%a_x*diami**params_mod%b_x*(rho0/rho_mod)**params_mod%g_x
        delv=abs(vr-vi)
        ! last bit is to convert to integral over m
        dintegral_mode1=eri*pi*0.25_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-rain_params%d_x)) /  &
            (rain_params%c_x*rain_params%d_x)* &
            (diami**(1.0_wp-params_mod%d_x)) / (params_mod%c_x*params_mod%d_x)

        do i=1,size(y)
            call calculate_mode1(mr,mi(i),t_send,n,nt,nb,mb,mt)
            dintegral_mode1(i)=dintegral_mode1(i)*n
        enddo

    end function dintegral_mode1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for mode 1 ice multiplication                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_mode1_mass(x,y)
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_mode1_mass
        real(wp) :: diamr, mr, vr, n,nt,nb,mb,mt
        real(wp), dimension(size(y)) :: mi, diami, delv, vi
        integer :: i

        mr=x
        mi=y
        diamr=(mr/rain_params%c_x)**(1.0_wp/rain_params%d_x)
        diami=(mi/params_mod%c_x)**(1.0_wp/params_mod%d_x)
        ! fall-speeds
        ! fall-speed of rain
        vr=(rain_params%a_x*diamr**rain_params%b_x*exp(-rain_params%f_x*diamr) + &
         rain_params%a2_x*diamr**rain_params%b2_x*exp(-rain_params%f2_x*diamr)) * &
         (rho0/rho_mod)**rain_params%g_x
        ! fall-speed of ice
        vi=params_mod%a_x*diami**params_mod%b_x*(rho0/rho_mod)**params_mod%g_x
        delv=abs(vr-vi)
        ! last bit is to convert to integral over m
        dintegral_mode1_mass=eri*pi*0.25_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-rain_params%d_x)) /  &
            (rain_params%c_x*rain_params%d_x)* &
            (diami**(1.0_wp-params_mod%d_x)) / (params_mod%c_x*params_mod%d_x)

        do i=1,size(y)
            call calculate_mode1(mr,mi(i),t_send,n,nt,nb,mb,mt)
            dintegral_mode1_mass(i)=dintegral_mode1_mass(i)*(nt*mt+nb*mb)
        enddo

    end function dintegral_mode1_mass
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine sip_phillips_mode2(ixy_inner, dt, nz, cffields, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    real(wp), intent(in) :: cffields(:,:)
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! Local variables
    real(wp) :: gacr, sacr, raci1,raci2  ! accretion process rates
    real(wp) :: dnumber_s, dnumber_g,  &
                dnumber_i ! number conversion rate from snow/graupel
    real(wp) :: dmass_s, dmass_g, &
                dmass_i      ! mass conversion rate from snow/graupel
    real(wp) :: cf_snow, cf_graupel, cf_liquid, cf_ice, cf_rain, overlap_cfsnow, &
         overlap_cfgraupel, overlap_cfice, ice_mass, rain_mass, m0

    type(process_name) :: iproc ! processes selected depending on which species we're modifying

    integer :: k
    real(wp) :: arg, dummy3


    character(len=*), parameter :: RoutineName='SIP_PHILLIPS_MODE2'

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

    ! loop over all levels - exactly like the HM routine above
    do k = 1, nz
       if (TdegC(k,ixy_inner) < 0.0_wp) then
          if (l_prf_cfrac) then
             if (cffields(k,i_cfl) .gt. cfliq_small) then
                cf_liquid=cffields(k,i_cfl)
             else
                cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
             endif
             if (cffields(k,i_cfr) .gt. cfliq_small) then
                cf_rain=cffields(k,i_cfr)
             else
                cf_rain=cfliq_small !nonzero value - maybe move cf test higher up
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
             ! added for ice
             if (cffields(k,i_cfi) .gt. cfliq_small) then
                cf_ice=cffields(k,i_cfi)
             else
                cf_ice=cfliq_small !nonzero value - maybe move cf test higher up
             endif
          else
             cf_rain=1.0
             cf_snow=1.0
             cf_graupel=1.0
             cf_liquid=1.0
             cf_ice=1.0
          endif

          f_mode2=min(-Cwater*(TdegC(k,ixy_inner))/lf,1.0_wp)
          !use mixed-phase overlap function
          ! PJC changed to be for rain and added ice
          overlap_cfsnow=min(1.0,max(0.0,mpof*min(cf_rain, cf_snow) +         &
               max(0.0,(1.0-mpof)*(cf_rain+cf_snow-1.0))))
          overlap_cfgraupel=min(1.0,max(0.0,mpof*min(cf_rain, cf_graupel) +   &
               max(0.0,(1.0-mpof)*(cf_rain+cf_graupel-1.0))))
          overlap_cfice=min(1.0,max(0.0,mpof*min(cf_rain, cf_ice) +   &
               max(0.0,(1.0-mpof)*(cf_rain+cf_ice-1.0))))

          sacr=0.0
          gacr=0.0
          raci1=0.0
          raci2=0.0
          if (snow_params%i_1m > 0) then
            sacr=procs(snow_params%i_1m, i_sacr%id)%column_data(k)/overlap_cfsnow  !insnow process rate
            raci1=procs(snow_params%i_1m, i_raci%id)%column_data(k)/overlap_cfsnow ! ingraupel process rate
          endif
          if (graupel_params%i_1m > 0) then
            gacr=procs(graupel_params%i_1m, i_gacr%id)%column_data(k)/overlap_cfgraupel ! ingraupel process rate
            raci2=procs(graupel_params%i_1m, i_raci%id)%column_data(k)/overlap_cfgraupel
          endif

          rho_mod=rho(k,ixy_inner)
          alpha_r=dist_mu(k,rain_params%id)
          lambda0r=dist_lambda(k,rain_params%id)
          arg = 1.0_wp+alpha_r
          n0r=dist_n0(k,rain_params%id)*lambda0r**(arg) / &
            GammaFunc(arg)

          ! for use within module, ICE first++++++++++++++++++++++++++++++++++++++++++++++
          params_mod = ice_params
          alpha_i=dist_mu(k,params_mod%id)
          lambda0i=dist_lambda(k,params_mod%id)
          arg = 1.0_wp+alpha_i
          n0i=dist_n0(k,params_mod%id)*lambda0i**(arg) / &
            GammaFunc(arg)
          t_send = TdegK(k,ixy_inner)

          ! calculate the increase in ice crystal number
          mrthresh = rain_params%c_x*150.e-6_wp**rain_params%d_x
          mrupper  = rain_params%c_x*(pthreshr/lambda0r)**rain_params%d_x
          miupper  = params_mod%c_x*(pthreshi/lambda0i)**params_mod%d_x
          mrupper  = min(mrupper,miupper)
          rain_mass = qfields(k,i_qr)
          ice_mass = qfields(k,i_qi)
          dnumber_i = 0.0
          dmass_i = 0.0
          if((mrupper.gt.mrthresh).and.(rain_mass.gt.thresh_small(rain_params%i_1m)) .and. &
            (ice_mass.gt.thresh_small(params_mod%i_1m))) then
            dummy3 = gl_quad_2d(dintegral_mode2, limit1_mode2, limit2_mode2, &
                mrthresh, mrupper)

            ! it is a rate
            dnumber_i = rho_mod*max(dummy3,0.0)

            dmass_i = M0_hallet_mossop*dnumber_i
            ! calculate the average mass of new splinters
            m0 = dmass_i /max(dnumber_i,1.0e-7)
            ! less than 50% of rain mass accreted onto ice should form splinters
            dmass_i = min(dmass_i,(raci1+raci2)*0.5)
            ! re-scale
            dnumber_i = dmass_i / max( m0, M0_hallet_mossop)


          endif
          !END ICE------------------------------------------------------------------------

          ! for use within module, SNOW ++++++++++++++++++++++++++++++++++++++++++++++++++
          params_mod = snow_params
          alpha_i=dist_mu(k,params_mod%id)
          lambda0i=dist_lambda(k,params_mod%id)
          arg = 1.0_wp+alpha_i
          n0i=dist_n0(k,params_mod%id)*lambda0i**(arg) / &
            GammaFunc(arg)

          ! calculate the increase in ice crystal number
          mrupper  = rain_params%c_x*(pthreshr/lambda0r)**rain_params%d_x
          miupper  = params_mod%c_x*(pthreshs/lambda0i)**params_mod%d_x
          mrupper  = min(mrupper,miupper)
          ice_mass = qfields(k,i_qs)
          dnumber_s = 0.0
          dmass_s = 0.0
          if((mrupper.gt.mrthresh).and.(rain_mass.gt.thresh_small(rain_params%i_1m)) .and. &
            (ice_mass.gt.thresh_small(params_mod%i_1m))) then
            dummy3 = gl_quad_2d(dintegral_mode2, limit1_mode2, limit2_mode2, &
                mrthresh, mrupper)
            ! it is a rate
            dnumber_s = rho_mod*max(dummy3,0.0)

            dmass_s = M0_hallet_mossop*dnumber_s

            ! calculate the average mass of new splinters
            m0 = dmass_s /max(dnumber_s,1.0e-7)
            ! less than 50% of rain mass accreted onto snow should form splinters
            dmass_s = min(dmass_s,0.5*sacr)
            ! re-scale
            dnumber_s = dmass_s / max( m0, M0_hallet_mossop)

          endif
          !END SNOW-----------------------------------------------------------------------

          ! for use within module, GRAUPEL++++++++++++++++++++++++++++++++++++++++++++++++
          params_mod = graupel_params
          alpha_i=dist_mu(k,params_mod%id)
          lambda0i=dist_lambda(k,params_mod%id)
          arg = 1.0_wp+alpha_i
          n0i=dist_n0(k,params_mod%id)*lambda0i**(arg) / &
            GammaFunc(arg)

          ! calculate the increase in ice crystal number
          mrupper  = rain_params%c_x*(pthreshr/lambda0r)**rain_params%d_x
          miupper  = params_mod%c_x*(pthreshg/lambda0i)**params_mod%d_x
          mrupper  = min(mrupper,miupper)
          ice_mass = qfields(k,i_qg)
          dnumber_g = 0.0
          dmass_g = 0.0
          if((mrupper.gt.mrthresh).and.(rain_mass.gt.thresh_small(rain_params%i_1m)) .and. &
            (ice_mass.gt.thresh_small(params_mod%i_1m))) then
            dummy3 = gl_quad_2d(dintegral_mode2, limit1_mode2, limit2_mode2, &
                mrthresh, mrupper)

            ! it is a rate
            dnumber_g = rho_mod*max(dummy3,0.0)

            dmass_g = M0_hallet_mossop*dnumber_g
            ! calculate the average mass of new splinters
            m0 = dmass_g /max(dnumber_g,1.0e-7)
            ! less than 50% of rain mass accreted onto graupel should form splinters
            dmass_g = min(dmass_g,0.5*gacr)
            ! re-scale
            dnumber_g = dmass_g / max( m0, M0_hallet_mossop)
          endif
          !END GRAUPEL--------------------------------------------------------------------
            if ((raci1*overlap_cfsnow +raci2*overlap_cfgraupel + sacr*overlap_cfsnow + &
                gacr*overlap_cfgraupel)*dt > thresh_small(snow_params%i_1m)) then
                iproc=i_imo2

                dmass_i=dmass_i * overlap_cfice  ! convert back to grid mean
                dmass_g=dmass_g * overlap_cfgraupel  ! convert back to grid mean
                dmass_s=dmass_s * overlap_cfsnow ! convert back to grid mean

                dnumber_i=dnumber_i * overlap_cfice  ! convert back to grid mean
                dnumber_g=dnumber_g * overlap_cfgraupel  ! convert back to grid mean
                dnumber_s=dnumber_s * overlap_cfsnow  ! convert back to grid mean

                !-------------------
                ! Sources for ice...
                !-------------------
                procs(ice_params%i_1m, iproc%id)%column_data(k)=dmass_g + dmass_s
                procs(ice_params%i_2m, iproc%id)%column_data(k)=  &
                    dnumber_g + dnumber_s + dnumber_i

                !-------------------
                ! Sinks for snow...
                !-------------------
                if (sacr > 0.0) then
                   procs(snow_params%i_1m, iproc%id)%column_data(k)=-dmass_s
                   if (snow_params%l_2m) procs(snow_params%i_2m,iproc%id)%column_data(k)=0.0
                end if

                !---------------------
                ! Sinks for graupel...
                !---------------------
                if (gacr > 0.0) then
                   procs(graupel_params%i_1m, iproc%id)%column_data(k)=-dmass_g
                   if (graupel_params%l_2m) procs(graupel_params%i_2m,iproc%id)%column_data(k)=0.0
                end if

            endif
       end if
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine sip_phillips_mode2


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This evaluates the integrand                                                       !
    ! for mode 2 ice multiplication                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_mode2(x,y)
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_mode2
        real(wp) :: diamr, mr, vr
        real(wp), dimension(size(y)) :: mi, diami, delv, vi, k0, de, nfrag, nfrag_freeze1, &
            nfrag_freeze2

        mr=x
        mi=y
        diamr=(mr/rain_params%c_x)**(1.0_wp/rain_params%d_x)
        diami=(mi/params_mod%c_x)**(1.0_wp/params_mod%d_x)
        ! fall-speeds
        ! fall-speed of rain
        vr=(rain_params%a_x*diamr**rain_params%b_x*exp(-rain_params%f_x*diamr) + &
         rain_params%a2_x*diamr**rain_params%b2_x*exp(-rain_params%f2_x*diamr)) * &
         (rho0/rho_mod)**rain_params%g_x
        ! fall-speed of ice
        vi=params_mod%a_x*diami**params_mod%b_x*(rho0/rho_mod)**params_mod%g_x
        delv=abs(vr-vi)
        ! last bit is to convert to integral over m
        dintegral_mode2=eri*pi*0.25_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alpha_r* &
            exp(-lambda0r*diamr)*n0i*diami**alpha_i*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-rain_params%d_x)) /  &
            (rain_params%c_x*rain_params%d_x)* &
            (diami**(1.0_wp-params_mod%d_x)) / (params_mod%c_x*params_mod%d_x)

        ! cke from equation 6
        k0=0.5_wp*(mr*mi/(mr+mi))*(vr-vi)**2
        ! de parameter
        de=k0/(gamma_liq*pi*diamr**2)
        ! number of fragments in spalsh
        nfrag=3.0_wp*max(de-decrit, 0.0_wp)
        ! number of fragments in splash that freeze due to mode 1
        nfrag_freeze1=nfrag*f_mode2
        ! number of fragments in splash that freeze due to mode 2
        nfrag_freeze2=nfrag*(1.0_wp-f_mode2)*phi_mode2

        dintegral_mode2=dintegral_mode2*nfrag_freeze2

    end function dintegral_mode2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sip_phillips_breakup(ixy_inner, dt, nz, cffields, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    real(wp), intent(in) :: cffields(:,:)
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! The CB_* rates below are grid-mean tendencies after multiplying the
    ! in-cloud collision integrals by rho and the appropriate cloud overlap.
    real(wp) :: dnumber_cb_ii, dnumber_cb_is, dnumber_cb_ig
    real(wp) :: dnumber_cb_si, dnumber_cb_ss, dnumber_cb_sg
    real(wp) :: dnumber_cb_gg_small, dnumber_cb_gg_hail, dnumber_cb_gg
    real(wp) :: dmass_cb_si, dmass_cb_ss, dmass_cb_sg, dmass_cb_gg
    real(wp) :: cf_snow, cf_graupel, cf_ice
    real(wp) :: overlap_cf_is, overlap_cf_ig, overlap_cf_sg
    real(wp) :: raw_rate, fragment_mass, total_number, total_mass
    integer :: k
    type(process_name) :: iproc

    character(len=*), parameter :: RoutineName='SIP_PHILLIPS_BREAKUP'

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

    ! Retain the original UM/CASIM breakup-fragment mass:
    ! mass corresponding to a 100-micron ice particle using the CASIM
    ! ice mass-diameter relation.  The collision bookkeeping follows WRF,
    ! but the assumed fragment size remains the original UM choice.
    fragment_mass=ice_params%c_x*(100.e-6_wp)**ice_params%d_x

    do k = 1, nz
       if (TdegC(k,ixy_inner) < 0.0_wp) then

          if (l_prf_cfrac) then
             cf_snow    = max(cffields(k,i_cfs),cfliq_small)
             cf_graupel = max(cffields(k,i_cfg),cfliq_small)
             cf_ice     = max(cffields(k,i_cfi),cfliq_small)
          else
             cf_snow    = 1.0_wp
             cf_graupel = 1.0_wp
             cf_ice     = 1.0_wp
          endif

          overlap_cf_is=min(1.0_wp,max(0.0_wp,mpof*min(cf_ice,cf_snow) + &
               max(0.0_wp,(1.0_wp-mpof)*(cf_ice+cf_snow-1.0_wp))))
          overlap_cf_ig=min(1.0_wp,max(0.0_wp,mpof*min(cf_ice,cf_graupel) + &
               max(0.0_wp,(1.0_wp-mpof)*(cf_ice+cf_graupel-1.0_wp))))
          overlap_cf_sg=min(1.0_wp,max(0.0_wp,mpof*min(cf_snow,cf_graupel) + &
               max(0.0_wp,(1.0_wp-mpof)*(cf_snow+cf_graupel-1.0_wp))))

          rho_mod=rho(k,ixy_inner)
          t_send=TdegC(k,ixy_inner)

          dnumber_cb_ii=0.0_wp
          dnumber_cb_is=0.0_wp
          dnumber_cb_ig=0.0_wp
          dnumber_cb_si=0.0_wp
          dnumber_cb_ss=0.0_wp
          dnumber_cb_sg=0.0_wp
          dnumber_cb_gg_small=0.0_wp
          dnumber_cb_gg_hail=0.0_wp
          dnumber_cb_gg=0.0_wp

          dmass_cb_si=0.0_wp
          dmass_cb_ss=0.0_wp
          dmass_cb_sg=0.0_wp
          dmass_cb_gg=0.0_wp

          !-------------------------------------------------------------------
          ! CB_GG_SMALL: same behaviour as WRF.  Both particles are graupel
          ! in the 0.5--5 mm range and only one triangular half is integrated.
          !-------------------------------------------------------------------
          if (qfields(k,i_qg) > thresh_small(graupel_params%i_1m)) then
             call evaluate_cb_wrf(k,graupel_params,graupel_params, &
                  pthreshg,pthreshg,CB_GG_SMALL, &
                  500.e-6_wp,5.e-3_wp,500.e-6_wp,5.e-3_wp,.true.,raw_rate)
             dnumber_cb_gg_small=rho_mod*cf_graupel*raw_rate
          endif

          !-------------------------------------------------------------------
          ! CB_GG_HAIL: both graupel particles are >5 mm; again integrate one
          ! triangular half.  This deliberately mirrors the present WRF code.
          !-------------------------------------------------------------------
          if (qfields(k,i_qg) > thresh_small(graupel_params%i_1m)) then
             call evaluate_cb_wrf(k,graupel_params,graupel_params, &
                  pthreshg,pthreshg,CB_GG_HAIL, &
                  5.e-3_wp,-1.0_wp,5.e-3_wp,-1.0_wp,.true.,raw_rate)
             dnumber_cb_gg_hail=rho_mod*cf_graupel*raw_rate
          endif
          dnumber_cb_gg=dnumber_cb_gg_small+dnumber_cb_gg_hail
          dmass_cb_gg=dnumber_cb_gg*fragment_mass

          !-------------------------------------------------------------------
          ! CB_IG: ice is particle 1 and is restricted to 0.5--5 mm.
          ! The graupel partner spans its full numerical PSD range.
          ! Fragment mass remains in the ice category, so only Ni changes.
          !-------------------------------------------------------------------
          if ((qfields(k,i_qi) > thresh_small(ice_params%i_1m)) .and. &
              (qfields(k,i_qg) > thresh_small(graupel_params%i_1m))) then
             call evaluate_cb_wrf(k,ice_params,graupel_params, &
                  pthreshi,pthreshg,CB_IG, &
                  500.e-6_wp,5.e-3_wp,1.e-6_wp,-1.0_wp,.false.,raw_rate)
             dnumber_cb_ig=rho_mod*overlap_cf_ig*raw_rate
          endif

          !-------------------------------------------------------------------
          ! CB_SG: snow is particle 1 (0.5--5 mm), graupel spans its full PSD.
          ! New small-ice mass is taken from snow, matching WRF PIICSG.
          !-------------------------------------------------------------------
          if ((qfields(k,i_qs) > thresh_small(snow_params%i_1m)) .and. &
              (qfields(k,i_qg) > thresh_small(graupel_params%i_1m))) then
             call evaluate_cb_wrf(k,snow_params,graupel_params, &
                  pthreshs,pthreshg,CB_SG, &
                  500.e-6_wp,5.e-3_wp,1.e-6_wp,-1.0_wp,.false.,raw_rate)
             dnumber_cb_sg=rho_mod*overlap_cf_sg*raw_rate
             dmass_cb_sg=dnumber_cb_sg*fragment_mass
          endif

          !-------------------------------------------------------------------
          ! CB_II: particle 1 is ice in the 0.5--5 mm interval and particle 2
          ! spans the full ice PSD.  When both are in 0.5--5 mm the integrand
          ! removes the duplicated ordering exactly as in the WRF routine.
          !-------------------------------------------------------------------
          if (qfields(k,i_qi) > thresh_small(ice_params%i_1m)) then
             call evaluate_cb_wrf(k,ice_params,ice_params, &
                  pthreshi,pthreshi,CB_II, &
                  500.e-6_wp,5.e-3_wp,1.e-6_wp,-1.0_wp,.false.,raw_rate)
             dnumber_cb_ii=rho_mod*cf_ice*raw_rate
          endif

          !-------------------------------------------------------------------
          ! CB_SS: analogous to CB_II, but the small-ice fragment mass is
          ! transferred from snow to ice, matching WRF PIICSS.
          !-------------------------------------------------------------------
          if (qfields(k,i_qs) > thresh_small(snow_params%i_1m)) then
             call evaluate_cb_wrf(k,snow_params,snow_params, &
                  pthreshs,pthreshs,CB_SS, &
                  500.e-6_wp,5.e-3_wp,1.e-6_wp,-1.0_wp,.false.,raw_rate)
             dnumber_cb_ss=rho_mod*cf_snow*raw_rate
             dmass_cb_ss=dnumber_cb_ss*fragment_mass
          endif

          !-------------------------------------------------------------------
          ! CB_IS: ice is particle 1 (0.5--5 mm), snow is the full partner PSD.
          ! Since the fragmenting category is already ice, this is number-only.
          ! When snow is also in its 0.5--5 mm breakup range, the integrand
          ! retains CB_IS only when the ice particle is the smaller of the two.
          !-------------------------------------------------------------------
          if ((qfields(k,i_qi) > thresh_small(ice_params%i_1m)) .and. &
              (qfields(k,i_qs) > thresh_small(snow_params%i_1m))) then
             call evaluate_cb_wrf(k,ice_params,snow_params, &
                  pthreshi,pthreshs,CB_IS, &
                  500.e-6_wp,5.e-3_wp,1.e-6_wp,-1.0_wp,.false.,raw_rate)
             dnumber_cb_is=rho_mod*overlap_cf_is*raw_rate

             ! CB_SI: snow is particle 1 (0.5--5 mm), ice is the partner.
             ! In the IS/SI overlap, the integrand retains CB_SI only when
             ! snow is the smaller particle, so each physical collision is
             ! represented by only one of CB_IS or CB_SI.
             call evaluate_cb_wrf(k,snow_params,ice_params, &
                  pthreshs,pthreshi,CB_SI, &
                  500.e-6_wp,5.e-3_wp,1.e-6_wp,-1.0_wp,.false.,raw_rate)
             dnumber_cb_si=rho_mod*overlap_cf_is*raw_rate
             dmass_cb_si=dnumber_cb_si*fragment_mass
          endif

          total_number=dnumber_cb_gg+dnumber_cb_ig+dnumber_cb_sg+ &
               dnumber_cb_ii+dnumber_cb_ss+dnumber_cb_is+dnumber_cb_si
          total_mass=dmass_cb_gg+dmass_cb_sg+dmass_cb_ss+dmass_cb_si

          if (fragment_mass*total_number*dt > thresh_small(ice_params%i_1m)) then

             ! Ice-fragmenting collisions: CB_II + CB_IS + CB_IG.
             ! Fragment material remains in ice, so only Ni changes.
             iproc=i_iicb_i
             procs(ice_params%i_2m,iproc%id)%column_data(k)= &
                  dnumber_cb_ii+dnumber_cb_is+dnumber_cb_ig

             ! Snow-fragmenting collisions: CB_SI + CB_SS + CB_SG.
             ! Fragment number goes to ice; fragment mass moves snow -> ice.
             iproc=i_iicb_s
             procs(ice_params%i_2m,iproc%id)%column_data(k)= &
                  dnumber_cb_si+dnumber_cb_ss+dnumber_cb_sg
             procs(ice_params%i_1m,iproc%id)%column_data(k)= &
                  dmass_cb_si+dmass_cb_ss+dmass_cb_sg
             procs(snow_params%i_1m,iproc%id)%column_data(k)= &
                  -(dmass_cb_si+dmass_cb_ss+dmass_cb_sg)
             if (snow_params%l_2m) procs(snow_params%i_2m,iproc%id)%column_data(k)=0.0_wp

             ! Graupel-fragmenting collisions: CB_GG_SMALL + CB_GG_HAIL.
             ! Fragment number goes to ice; fragment mass moves graupel -> ice.
             iproc=i_iicb_g
             procs(ice_params%i_2m,iproc%id)%column_data(k)=dnumber_cb_gg
             procs(ice_params%i_1m,iproc%id)%column_data(k)=dmass_cb_gg
             procs(graupel_params%i_1m,iproc%id)%column_data(k)=-dmass_cb_gg
             if (graupel_params%l_2m) procs(graupel_params%i_2m,iproc%id)%column_data(k)=0.0_wp

          endif
       endif
    enddo

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine sip_phillips_breakup

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Evaluate one WRF-style collisional-breakup contribution.
    !
    ! d1_low/d1_high define the physical diameter range of particle 1.
    ! d2_low/d2_high define particle 2.  A negative upper bound means use the
    ! species-specific numerical PSD cutoff pthresh/lambda.
    !
    ! triangular=.true. gives y>=x and is used only for the two WRF graupel-
    ! graupel regimes, where both species and mass-diameter relations are the same.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine evaluate_cb_wrf(k,p1,p2,pth1,pth2,pair_id, &
         d1_low,d1_high,d2_low,d2_high,triangular,rate)
        implicit none
        integer, intent(in) :: k,pair_id
        type(hydro_params), intent(in) :: p1,p2
        real(wp), intent(in) :: pth1,pth2,d1_low,d1_high,d2_low,d2_high
        logical, intent(in) :: triangular
        real(wp), intent(out) :: rate
        real(wp) :: arg,dmax1,dmax2,d1u,d2u,dummy3

        rate=0.0_wp
        params_mod=p1
        params_mod1=p2
        type1_send=pair_id

        alpha_r=dist_mu(k,p1%id)
        lambda0r=dist_lambda(k,p1%id)
        alpha_i=dist_mu(k,p2%id)
        lambda0i=dist_lambda(k,p2%id)

        if ((lambda0r <= 0.0_wp) .or. (lambda0i <= 0.0_wp)) return
        if ((dist_n0(k,p1%id) <= 0.0_wp) .or. &
            (dist_n0(k,p2%id) <= 0.0_wp)) return

        arg=1.0_wp+alpha_r
        n0r=dist_n0(k,p1%id)*lambda0r**arg/GammaFunc(arg)
        arg=1.0_wp+alpha_i
        n0i=dist_n0(k,p2%id)*lambda0i**arg/GammaFunc(arg)

        dmax1=pth1/lambda0r
        dmax2=pth2/lambda0i

        if (d1_high > 0.0_wp) then
           d1u=min(d1_high,dmax1)
        else
           d1u=dmax1
        endif
        if (d2_high > 0.0_wp) then
           d2u=min(d2_high,dmax2)
        else
           d2u=dmax2
        endif

        if ((d1u <= d1_low) .or. (d2u <= d2_low)) return

        mrthresh=p1%c_x*d1_low**p1%d_x
        mrupper=p1%c_x*d1u**p1%d_x
        mithresh=p2%c_x*d2_low**p2%d_x
        miupper=p2%c_x*d2u**p2%d_x

        ! Mass-space limits corresponding to the 0.5--5 mm breakup range
        ! for particle 2.  These are used only to identify the overlap
        ! between the directional CB_IS and CB_SI integrations.  Keep the
        ! physical size limits in one place rather than hard-coding diameter
        ! tests in dintegral_collisional_breakup.
        mifragthresh=p2%c_x*(500.e-6_wp)**p2%d_x
        mifragupper=min(p2%c_x*(5000.e-6_wp)**p2%d_x,miupper)

        if (triangular) then
           dummy3 = gl_quad_2d(dintegral_collisional_breakup, limit1_coll_x, limit2_coll, &
                mrthresh, mrupper)
        else
           dummy3 = gl_quad_2d(dintegral_collisional_breakup, limit1_coll, limit2_coll, &
                mrthresh, mrupper)
        endif

        rate=max(dummy3,0.0_wp)
    end subroutine evaluate_cb_wrf
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Integrand for collisional breakup following Phillips et al. (2017).
    !
    ! This follows the agreed WRF/CASIM directional bookkeeping:
    !   * CB_XY means X is the fragmenting parent and Y is the collision partner.
    !   * CB_IS and CB_SI partition the ice-snow collision space rather than
    !     double-counting the region where both particles are breakup-eligible.
    !   * CB_GG is split into SMALL and HAIL regimes.
    !   * dsmall=min(D1,D2) is used in the Phillips ice/snow expressions.
    !   * CB_II and CB_SS remove duplicate orderings only when both particles
    !     are within the 0.5--5 mm particle-1 interval.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dintegral_collisional_breakup(x,y)
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral_collisional_breakup
        real(wp) :: m1, diam1, v1, a0, t0, cpar, nmax, gamma
        real(wp), dimension(size(y)) :: m2, diam2, v2, delv, k0, nfrag
        real(wp), dimension(size(y)) :: dsmall, alpha, apar
        real(wp), parameter :: frimes=0.1_wp

        m1=x
        m2=y

        diam1=(m1/params_mod%c_x)**(1.0_wp/params_mod%d_x)
        diam2=(m2/params_mod1%c_x)**(1.0_wp/params_mod1%d_x)
        dsmall=min(diam1,diam2)

        v1=params_mod%a_x*diam1**params_mod%b_x * &
             (rho0/rho_mod)**params_mod%g_x
        v2=params_mod1%a_x*diam2**params_mod1%b_x * &
             (rho0/rho_mod)**params_mod1%g_x
        delv=abs(v1-v2)

        ! In-cloud collision rate using number-mixing-ratio PSDs.  rho_mod is
        ! applied outside the quadrature to obtain # kg-air^-1 s^-1.
        dintegral_collisional_breakup=eri*pi*0.25_wp*(diam1+diam2)**2* &
             delv*n0r*diam1**alpha_r*exp(-lambda0r*diam1)* &
             n0i*diam2**alpha_i*exp(-lambda0i*diam2)* &
             diam1**(1.0_wp-params_mod%d_x)/(params_mod%c_x*params_mod%d_x)* &
             diam2**(1.0_wp-params_mod1%d_x)/(params_mod1%c_x*params_mod1%d_x)

        alpha=pi*dsmall**2
        apar=0.0_wp
        cpar=0.0_wp
        nmax=0.0_wp
        gamma=1.0_wp

        select case (type1_send)

        case (CB_GG_SMALL)
           a0=3.78e4_wp*(1.0_wp+0.0079_wp/diam1**1.5_wp)
           t0=-15.0_wp
           apar=a0*oneoverthree+max(2.0_wp*a0*oneoverthree- &
                a0*oneovernine*abs(t_send-t0),0.0_wp)
           cpar=6.30e6_wp*phi_phillips
           nmax=100.0_wp
           gamma=0.30_wp

        case (CB_GG_HAIL)
           a0=4.35e5_wp
           t0=-15.0_wp
           apar=a0*oneoverthree+max(2.0_wp*a0*oneoverthree- &
                a0*oneovernine*abs(t_send-t0),0.0_wp)
           cpar=3.31e5_wp
           nmax=1000.0_wp
           gamma=0.54_wp

        case (CB_II,CB_IS,CB_IG,CB_SI,CB_SS,CB_SG)
           if ((t_send >= -17.0_wp) .and. (t_send <= -12.0_wp)) then
              apar=1.41e6_wp*(1.0_wp+100.0_wp*frimes**2)* &
                   (1.0_wp+3.98e-5_wp/dsmall**1.5_wp)
              cpar=3.09e6_wp*frimes
              nmax=100.0_wp
              gamma=0.50_wp-0.25_wp*frimes
           elseif ((t_send < -17.0_wp) .or. &
                  ((t_send > -12.0_wp) .and. (t_send <= -9.0_wp))) then
              apar=1.58e7_wp*(1.0_wp+100.0_wp*frimes**2)* &
                   (1.0_wp+1.33e-4_wp/dsmall**1.5_wp)
              cpar=7.08e6_wp*frimes
              nmax=100.0_wp
              gamma=0.50_wp-0.25_wp*frimes
           else
              dintegral_collisional_breakup=0.0_wp
              return
           endif

        case default
           dintegral_collisional_breakup=0.0_wp
           return
        end select

        k0=0.5_wp*(m1*m2/(m1+m2))*delv**2

        where (alpha*apar > 0.0_wp)
           nfrag=min(alpha*apar* &
                (1.0_wp-exp(-(cpar*k0/(alpha*apar))**gamma)),nmax)
        elsewhere
           nfrag=0.0_wp
        end where

        ! WRF duplicate-pair removal for II/SS.  Particle 1 is restricted to
        ! 0.5--5 mm; if particle 2 is also in that interval, retain one ordering.
        if ((type1_send == CB_II) .or. (type1_send == CB_SS)) then
           where ((m2 >= mrthresh) .and. (m2 <= mrupper) .and. (m2 < m1))
              dintegral_collisional_breakup=0.0_wp
           end where

        elseif (type1_send == CB_IS) then
           ! Particle 1 = ice, particle 2 = snow.  If snow is also in its
           ! breakup range, assign the collision to CB_IS only when ice is
           ! the smaller particle.  The partner eligibility test is done
           ! entirely in mass space using the precomputed 0.5--5 mm bounds;
           ! diameters are compared only to decide which particle is smaller.
           where ((m2 >= mifragthresh) .and. (m2 <= mifragupper) .and. &
                  (diam2 < diam1))
              dintegral_collisional_breakup=0.0_wp
           end where

        elseif (type1_send == CB_SI) then
           ! Particle 1 = snow, particle 2 = ice.  In the overlap, retain
           ! CB_SI only when snow is strictly the smaller particle.  The
           ! equality case is assigned to CB_IS so it is counted once.
           where ((m2 >= mifragthresh) .and. (m2 <= mifragupper) .and. &
                  (diam2 <= diam1))
              dintegral_collisional_breakup=0.0_wp
           end where
        endif

        dintegral_collisional_breakup=dintegral_collisional_breakup*nfrag

    end function dintegral_collisional_breakup
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!             ! type 1
!             elseif ((frimes>=0.9_wp).and.(frimel>=0.9_wp)) then
!
!                 ! collisions of hail and hail - no size constraint
!                     a0*oneovernine*abs(t-ttr-t0),0._wp)
!
!             ! types 2 or 3
!                 .and. (phis < 1._wp)) then
!                 ! collisions of ice /snow size 500 micron to 5 mm with any ice
!
!                 ! seems like columnar habits dont fragment?
!                     ! dendrites
!                         (1._wp+3.98e-5_wp/dsmall**1.5_wp)
!
!                 else
!                     ! spatial planar
!                         (1._wp+1.33e-4_wp/dsmall**1.5_wp)
!
!             else
!             ! CKE
!             ! finally apply equation 13
!             !---------------------------------------------------------------------------------
!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function limit1_coll(x)
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit1_coll
        limit1_coll=mithresh
    end function limit1_coll

    function limit1_coll_x(x)
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit1_coll_x
        limit1_coll_x=x
    end function limit1_coll_x

    function limit2_coll(x)
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit2_coll
        limit2_coll=miupper
    end function limit2_coll
!
    function limit1_mode1(x)
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit1_mode1
        limit1_mode1=milower
    end function limit1_mode1

    function limit2_mode1(x)
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit2_mode1
        limit2_mode1=min(x,miupper)
    end function limit2_mode1
!
    function limit1_mode2(x)
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit1_mode2
        limit1_mode2=x
    end function limit1_mode2

    function limit2_mode2(x)
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit2_mode2
        limit2_mode2=miupper
    end function limit2_mode2

  !-----------------------------------------------------------------------
  ! Quadrature for the Phillips et al. collision integrals
  !-----------------------------------------------------------------------
  !> Integral of f(x) from a to b with the 10-point Gauss-Legendre rule.
  function gl_quad_1d(f, a, b) result(total)

    implicit none

    procedure(integrand_1d) :: f
    real(wp), intent(in) :: a, b
    real(wp) :: total

    real(wp) :: centre, half_width
    real(wp) :: abscissa(2*n_half), weight(2*n_half)

    centre = 0.5_wp*(a + b)
    half_width = 0.5_wp*(b - a)

    abscissa(1:n_half) = centre - half_width*gl_node
    abscissa(n_half+1:2*n_half) = centre + half_width*gl_node
    weight(1:n_half) = gl_weight
    weight(n_half+1:2*n_half) = gl_weight

    total = half_width*sum(weight*f(abscissa))

  end function gl_quad_1d

  !> Integral over x from a to b, and over y from y_low(x) to y_high(x),
  !> of f(x,y), using the 10-point Gauss-Legendre rule in each direction.
  !> The outer coordinate is passed to the integrand explicitly, so the
  !> routine holds no module state and is safe to call from OpenMP threads.
  function gl_quad_2d(f, y_low, y_high, a, b) result(total)

    implicit none

    procedure(integrand_2d) :: f
    procedure(limit_function) :: y_low, y_high
    real(wp), intent(in) :: a, b
    real(wp) :: total

    real(wp) :: centre, half_width, x_outer, inner
    real(wp) :: y_centre, y_half_width
    real(wp) :: y_abscissa(2*n_half), weight(2*n_half)
    integer :: i

    weight(1:n_half) = gl_weight
    weight(n_half+1:2*n_half) = gl_weight

    centre = 0.5_wp*(a + b)
    half_width = 0.5_wp*(b - a)

    total = 0.0_wp
    do i = 1, 2*n_half
      if (i <= n_half) then
        x_outer = centre - half_width*gl_node(i)
      else
        x_outer = centre + half_width*gl_node(i-n_half)
      end if

      y_centre = 0.5_wp*(y_high(x_outer) + y_low(x_outer))
      y_half_width = 0.5_wp*(y_high(x_outer) - y_low(x_outer))
      y_abscissa(1:n_half) = y_centre - y_half_width*gl_node
      y_abscissa(n_half+1:2*n_half) = y_centre + y_half_width*gl_node

      inner = y_half_width*sum(weight*f(x_outer, y_abscissa))
      total = total + weight(i)*inner
    end do
    total = half_width*total

  end function gl_quad_2d

end module ice_multiplication
