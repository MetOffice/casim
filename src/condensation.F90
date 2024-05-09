module condensation
  use variable_precision, only: wp
  use passive_fields, only: rho, pressure, w, exner
  use mphys_switches, only: i_qv, i_ql, i_nl, i_th, i_qr, l_warm, &
       i_am4, i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am6, i_an6, i_am9, i_an11, i_an12,  &
       cloud_params, l_process, l_passivenumbers,l_passivenumbers_ice, aero_index, &
       l_cfrac_casim_diag_scheme
  use process_routines, only: process_rate, i_cond, i_aact
  use mphys_constants, only: Lv, cp
  use qsat_funs, only: qsaturation, dqwsatdt
  use thresholds, only: ql_small, ss_small, thresh_tidy
  use activation, only: activate
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use which_mode_to_use, only : which_mode
  use casim_runtime, only: casim_time, casim_smax, casim_smax_limit_time
  use casim_parent_mod, only: casim_parent, parent_um, parent_kid
  use cloud_frac_scheme, only: cloud_frac_casim_mphys

! #if DEF_MODEL==MODEL_KiD
!   use diagnostics, only: save_dg, i_dgtime, i_here, k_here
!   use runtime, only: time
!   Use namelists, only : smax, smax_limit_time
! #endif

  implicit none

  character(len=*), parameter, private :: ModuleName='CONDENSATION'

  private

      real(wp), allocatable :: dnccn_all(:),dmac_all(:)
      real(wp), allocatable :: dnccnd_all(:),dmad_all(:)

!$OMP THREADPRIVATE(dnccn_all, dmac_all, dnccnd_all, dmad_all)

  public condevp_initialise, condevp_finalise, condevp
!PRF
  public dnccn_all, dmac_all, dnccnd_all, dmad_all
!PRF

contains

  subroutine condevp_initialise()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Local variables
    character(len=*), parameter :: RoutineName='CONDEVP_INITIALISE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    allocate(dnccn_all(aero_index%nccn))
    allocate(dmac_all(aero_index%nccn))
    allocate(dnccnd_all(aero_index%nin))
    allocate(dmad_all(aero_index%nin))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine condevp_initialise  

  subroutine condevp_finalise()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Local variables
    character(len=*), parameter :: RoutineName='CONDEVP_FINALISE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    deallocate(dnccn_all)
    deallocate(dmac_all)
    deallocate(dnccnd_all)
    deallocate(dmad_all)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine condevp_finalise  

  subroutine condevp(ixy_inner, dt, nz, qfields, procs, aerophys, aerochem,   &
       aeroact, dustphys, dustchem, dustliq, aerosol_procs, rhcrit_lev)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! aerosol fields
    type(aerosol_phys), intent(in) :: aerophys(:)
    type(aerosol_chem), intent(in) :: aerochem(:)
    type(aerosol_active), intent(in) :: aeroact(:)
    type(aerosol_phys), intent(in) :: dustphys(:)
    type(aerosol_chem), intent(in) :: dustchem(:)
    type(aerosol_active), intent(in) :: dustliq(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), optional, target :: aerosol_procs(:,:)

    real(wp), intent(in) :: rhcrit_lev(:)

    ! Local variables
    real(wp) :: dmass, dnumber, dmac, dmad, dnumber_a, dnumber_d
    real(wp) :: dmac1, dmac2, dnac1, dnac2

    real(wp) :: th
    real(wp) :: qv
    real(wp) :: cloud_mass
    real(wp) :: cloud_number

    real(wp) :: qs, dqsdt, qsatfac

    real(wp) :: tau   ! timescale for adjustment of condensate
    real(wp) :: w_act ! vertical velocity to use for activation

    real(wp) :: smax,ait_ccn, acc_ccn, tot_ccn, activated_arg, &
         activated_cloud
    ! local variables for diagnostics cloud scheme (if needed)
    real(wp) :: cloud_mass_new, abs_liquid_t

    real(wp) :: cfrac, cfrac_old

    logical :: l_docloud  ! do we want to do the calculation of cond/evap
    integer :: k

    character(len=*), parameter :: RoutineName='CONDEVP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    tau=dt ! adjust instantaneously

    ! Initializations
    dmac=0.0
    dmad=0.0
    dnumber=0.0
    dnumber_a=0.0
    dnumber_d=0.0
    
    dnccn_all=0.0
    dmac_all=0.0

    do k = 1, nz

    cfrac=1.0
    cfrac_old=1.0

    ! Set pointers for convenience
    cloud_mass=qfields(k, i_ql)
    if (cloud_params%l_2m) cloud_number=qfields(k, i_nl)

    th=qfields(k, i_th)

    if (casim_parent == parent_um ) then

      if (l_cfrac_casim_diag_scheme ) then
        qv=qfields(k, i_qv)+cloud_mass
      else
        qv=qfields(k, i_qv)
      end if

      if (l_cfrac_casim_diag_scheme) then
        ! work out saturation vapour pressure/mixing ratio based on
        ! liquid water temperature
        abs_liquid_T=(th*exner(k,ixy_inner))-((lv * cloud_mass )/cp)
        qs=qsaturation(abs_liquid_T, pressure(k,ixy_inner)/100.0)
      else
        qs=qsaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)
      end if
    
    else ! casim_parent /= parent_um
      qv=qfields(k, i_qv)
      qs=qsaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)

    end if ! casim_parent == parent_um

    l_docloud=.true.
    if (qs==0.0) l_docloud=.false.

      
    if (casim_parent == parent_kid) then
      if ((qv/qs > 1.0 - ss_small .or. cloud_mass > 0.0) .and. l_docloud) then
!AH - following code limits the maximum supersaturation permitted. This 
!     is needed for the KiD-A 2d Sc case  - USE WITH CAUTION!!
        if ( (((qv/qs)-1.)*100.) > casim_smax .and. casim_time <= casim_smax_limit_time) then
          qs = qv/(1+(casim_smax/100.))
        endif
      endif
    endif ! casim_parent == parent_kid
! AH - set the cloud fraction for new calc of activation
    if (cloud_mass > epsilon(1.0_wp)) then 
       cfrac_old = 1.0_wp
    else 
       cfrac_old = 0.0_wp
    endif
       
    if ((qv/qs > 1.0 - ss_small .or. cloud_mass > 0.0 .or. l_cfrac_casim_diag_scheme) .and. l_docloud) then
! DPG - allow the cloud scheme to operate even if we are sub-saturated (since
! this is it's purpose!)
      if (l_cfrac_casim_diag_scheme .AND. casim_parent == parent_um ) then

        !Call Smith scheme before setting up microphysics vars, to work out
        ! cloud fraction, which is used to derive in-cloud mass and number
        !
        !IMPORTANT - qv is total water at this stage!
        call cloud_frac_casim_mphys(k, pressure(k,ixy_inner), th*exner(k,ixy_inner), abs_liquid_T, rhcrit_lev(k),  &
             qs, qv, cloud_mass, qfields(k,i_qr), cloud_mass_new )

        dmass=max(-cloud_mass, (cloud_mass_new-cloud_mass))/dt
      else
        dqsdt=dqwsatdt(qs, th*exner(k,ixy_inner))
        qsatfac=1.0/(1.0 + Lv/Cp*dqsdt)
        dmass=max(-cloud_mass, (qv-qs)*qsatfac )/dt
      end if ! l_cfrac_casim_diag_scheme

      if (dmass > 0.0_wp) then ! condensation
        if (dmass*dt + cloud_mass > ql_small) then ! is it worth bothering with?
           ! AH - if dmass > 0.0 there is a change in mass, so assume cloud fraction is 1.0
           ! this assumption is only valid with all-or-nothing scheme and no cloud fraction
           ! scheme
          cfrac = 1.0_wp
          if (cloud_params%l_2m) then
            ! If significant cloud formed then assume minimum velocity of 0.01m/s
            w_act=max(w(k,ixy_inner), 0.01_wp)

            call activate(tau, cloud_mass, cloud_number, w_act,         &
                 rho(k,ixy_inner), dnumber, dmac, th*exner(k,ixy_inner), pressure(k,ixy_inner),       &
                 cfrac, cfrac_old, aerophys(k), aerochem(k),            & 
                 aeroact(k), dustphys(k), dustchem(k), dustliq(k),      &
                 dnccn_all, dmac_all, dnumber_d, dmad,                  &
                 dnccnd_all, dmad_all, smax, ait_ccn, acc_ccn,          &
                 tot_ccn,activated_arg,activated_cloud)

            dnumber_a=dnumber
          end if
        else
          dmass=0.0 ! not worth doing anything
        end if
      else  ! evaporation
        if (cloud_mass > thresh_tidy(i_ql)) then ! anything significant to remove or just noise?
          if (dmass*dt + cloud_mass < ql_small) then  ! Remove all cloud
            ! Remove small quantities.
            dmass=-cloud_mass/dt
            ! liberate all number and aerosol
            if (cloud_params%l_2m) then
              dnumber=-cloud_number/dt

              !============================
              ! aerosol processing
              !============================
              if (l_process) then
                dmac=-aeroact(k)%mact1/dt
                dmad=-dustliq(k)%mact1/dt

                if (l_passivenumbers) then
                  dnumber_a=-aeroact(k)%nact1/dt
                else
                  dnumber_a=dnumber
                end if
                if (l_passivenumbers_ice) then
                  dnumber_d=-dustliq(k)%nact1/dt
                else
                  dnumber_d=dnumber
                end if

                if (aero_index%nin > 0) then 
                   dmad_all(aero_index%i_coarse_dust) = dmad
                   dnccnd_all(aero_index%i_coarse_dust) = dnumber_d
                endif

                if (aero_index%i_accum >0 .and. aero_index%i_coarse >0) then
                  ! We have both accumulation and coarse modes
                  if (dnumber_a*dmac<= 0) then
                    dnumber_a=dmac/1.0e-18/dt
                  end if
                  call which_mode(dmac, dnumber_a,                                 &
                       aerophys(k)%rd(aero_index%i_accum), aerophys(k)%rd(aero_index%i_coarse), &
                       aerochem(k)%density(aero_index%i_accum),     &
                       aerophys(k)%sigma(aero_index%i_accum),       &
                       dmac1, dmac2, dnac1, dnac2)

                  dmac_all(aero_index%i_accum)=dmac1  ! put it back into accumulation mode
                  dnccn_all(aero_index%i_accum)=dnac1
                  dmac_all(aero_index%i_coarse)=dmac2  ! put it back into coarse mode
                  dnccn_all(aero_index%i_coarse)=dnac2

                else

                  if (aero_index % i_accum > 0) then
                    dmac_all  ( aero_index % i_accum) = dmac
                    dnccn_all ( aero_index % i_accum) = dnumber_a
                  end if

                  if (aero_index % i_coarse > 0) then
                    dmac_all  (aero_index % i_coarse) = dmac
                    dnccn_all (aero_index % i_coarse) = dnumber_a
                  end if

                end if
              end if
            end if
          else ! Still some cloud will be left behind
            dnumber=0.0 ! we assume no change in number during evap
            dnccn_all=0.0 ! we assume no change in number during evap
            dmac=0.0 ! No aerosol processing required
            dmac_all=0.0 ! No aerosol processing required
            dnumber_a=0.0 ! No aerosol processing required
            dnumber_d=0.0 ! No aerosol processing required
            dnccnd_all = 0.0
            dmad_all = 0.0
          end if
        else  ! Nothing significant here to remove - the tidying routines will deal with this
          dmass=0.0 ! no need to do anything since this is now just numerical noise
          dnumber=0.0 ! we assume no change in number during evap
          dnccn_all=0.0 ! we assume no change in number during evap
          dmac=0.0 ! No aerosol processing required
          dmac_all=0.0 ! No aerosol processing required
          dnumber_a=0.0 ! No aerosol processing required
          dnumber_d=0.0 ! No aerosol processing required
          dnccnd_all = 0.0
          dmad_all = 0.0
        end if
      end if

      if (dmass /= 0.0_wp) then

        procs(i_qv, i_cond%id)%column_data(k)=-dmass
        procs(i_ql, i_cond%id)%column_data(k)=dmass

        if (cloud_params%l_2m) then
          procs(i_nl, i_cond%id)%column_data(k)=dnumber
        end if

        !============================
        ! aerosol processing
        !============================
        if (l_process) then

          aerosol_procs(i_am4, i_aact%id)%column_data(k)=dmac
          if (l_passivenumbers) aerosol_procs(i_an11, i_aact%id)%column_data(k)=dnumber_a
          if (l_passivenumbers_ice) aerosol_procs(i_an12, i_aact%id)%column_data(k)=dnumber_d

          if (aero_index%i_aitken > 0) then
            aerosol_procs(i_am1, i_aact%id)%column_data(k)=-dmac_all(aero_index%i_aitken)
            aerosol_procs(i_an1, i_aact%id)%column_data(k)=-dnccn_all(aero_index%i_aitken)
          end if
          if (aero_index%i_accum > 0) then
            aerosol_procs(i_am2, i_aact%id)%column_data(k)=-dmac_all(aero_index%i_accum)
            aerosol_procs(i_an2, i_aact%id)%column_data(k)=-dnccn_all(aero_index%i_accum)
          end if
          if (aero_index%i_coarse > 0) then
            aerosol_procs(i_am3, i_aact%id)%column_data(k)=-dmac_all(aero_index%i_coarse)
            aerosol_procs(i_an3, i_aact%id)%column_data(k)=-dnccn_all(aero_index%i_coarse)
          end if

          if (.not. l_warm .and. dmad /=0.0) then
            ! We may have some dust in the liquid...
             if (aero_index%nin > 0 ) then 
                aerosol_procs(i_am9, i_aact%id)%column_data(k)= dmad_all(aero_index%i_coarse_dust)
                aerosol_procs(i_am6, i_aact%id)%column_data(k)=-dmad_all(aero_index%i_coarse_dust)! < USING COARSE
                aerosol_procs(i_an6, i_aact%id)%column_data(k)=-dnccnd_all(aero_index%i_coarse_dust) ! < USING COARSE
             end if
          end if
        end if
       end if
     end if 
   enddo
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine condevp
end module condensation
