module condensation
  use variable_precision, only: wp
  use passive_fields, only: rho, pressure, w, exner, qws
  use mphys_switches, only: i_qv, i_ql, i_nl, i_th, i_qr, i_qi, i_qs, i_qg, hydro_complexity, l_warm, &
       i_am4, i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am6, i_an6, i_am9, i_an11, i_an12,  &
       cloud_params, l_process, l_passivenumbers,l_passivenumbers_ice, aero_index, &
       active_cloud, active_number, l_preventsmall, l_cfrac_casim_diag_scheme
  use process_routines, only: process_rate, i_cond, i_aact
  use mphys_constants, only: Lv, cp
  use qsat_funs, only: qsaturation, dqwsatdt
  use thresholds, only: ql_small, w_small, nl_small, ss_small, thresh_tidy
  use activation, only: activate
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use special, only: pi
  use which_mode_to_use, only : which_mode
  use casim_parent_mod, only: casim_parent, parent_um
  use cloud_frac_scheme, only: cloud_frac_casim_mphys

  implicit none

  character(len=*), parameter, private :: ModuleName='CONDENSATION'

  private

  logical :: l_notransfer=.true.  ! don't transfer aerosol from one mode to another.

      real(wp), allocatable :: dnccn_all(:),dmac_all(:)

  public condevp_initialise, condevp_finalise, condevp
contains

  subroutine condevp_initialise()
    implicit none

    character(len=*), parameter :: RoutineName='CONDEVP_INITIALISE'

    allocate(dnccn_all(aero_index%nccn))
    allocate(dmac_all(aero_index%nccn))
  end subroutine condevp_initialise  

  subroutine condevp_finalise()
    implicit none

    character(len=*), parameter :: RoutineName='CONDEVP_FINALISE'

    deallocate(dnccn_all)
    deallocate(dmac_all)
  end subroutine condevp_finalise  

  subroutine condevp(dt, k, qfields, aerofields, procs, aerophys, aerochem,   &
       aeroact, dustphys, dustchem, dustliq, aerosol_procs, rhcrit_lev)

    implicit none

    character(len=*), parameter :: RoutineName='CONDEVP'

    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    real(wp), intent(in), target :: qfields(:,:), aerofields(:,:)
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

    real(wp), intent(in) :: rhcrit_lev

    real(wp) :: dmass, dnumber, dmac, dmad, dnumber_a, dnumber_d
    real(wp) :: dmac1, dmac2, dnac1, dnac2

    real(wp) :: th
    real(wp) :: qv
    real(wp) :: cloud_mass
    real(wp) :: cloud_number
    type(process_rate), pointer :: this_proc
    type(process_rate), pointer :: aero_proc

    real(wp) :: qs, dqsdt, qsatfac

    real(wp) :: tau   ! timescale for adjustment of condensate
    real(wp) :: w_act ! vertical velocity to use for activation

    real(wp) :: rmean ! mean aerosol diameter

    ! local variables for diagnostics cloud scheme (if needed)
    real(wp) :: cloud_mass_new, abs_liquid_t

    logical :: l_docloud  ! do we want to do the calculation of cond/evap

    tau=dt ! adjust instantaneously

    ! Initializations
    dmac=0.0
    dmad=0.0
    dnumber=0.0
    dnumber_a=0.0
    dnumber_d=0.0
    
    dnccn_all=0.0
    dmac_all=0.0

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
        abs_liquid_T=(th*exner(k))-((lv * cloud_mass )/cp)
        qs=qsaturation(abs_liquid_T, pressure(k)/100.0)
      else
        qs=qsaturation(th*exner(k), pressure(k)/100.0)
      end if
    
    else ! casim_parent /= parent_um
      qv=qfields(k, i_qv)
      qs=qsaturation(th*exner(k), pressure(k)/100.0)

    end if ! casim_parent == parent_um

    l_docloud=.true.
    if (qs==0.0) l_docloud=.false.

    if ((qv/qs > 1.0 - ss_small .or. cloud_mass > 0.0) .and. l_docloud) then

      if (l_cfrac_casim_diag_scheme .AND. casim_parent == parent_um ) then

        !Call Smith scheme before setting up microphysics vars, to work out
        ! cloud fraction, which is used to derive in-cloud mass and number
        !
        !IMPORTANT - qv is total water at this stage!
        call cloud_frac_casim_mphys(k, pressure(k), abs_liquid_T, rhcrit_lev,         &
             qs, qv, cloud_mass, qfields(k,i_qr), cloud_mass_new )

        dmass=max(-cloud_mass, (cloud_mass_new-cloud_mass))/dt
      else
        dqsdt=dqwsatdt(qs, th*exner(k))
        qsatfac=1.0/(1.0 + Lv/Cp*dqsdt)
        dmass=max(-cloud_mass, (qv-qs)*qsatfac )/dt
      end if ! l_cfrac_casim_diag_scheme

      if (dmass > 0.0_wp) then ! condensation
        if (dmass*dt + cloud_mass > ql_small) then ! is it worth bothering with?
          if (cloud_params%l_2m) then
            ! If significant cloud formed then assume minimum velocity of 0.1m/s
            w_act=max(w(k), 0.1_wp)

            call activate(tau, cloud_mass, cloud_number, w_act, rho(k), dnumber, dmac, &
                 th*exner(k), pressure(k), aerophys(k), aerochem(k), aeroact(k),   &
                 dustphys(k), dustchem(k), dustliq(k),           &
                 dnccn_all, dmac_all, dnumber_d, dmad)

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

                if (aero_index%i_accum >0 .and. aero_index%i_coarse >0) then
                  ! We have both accumulation and coarse modes
                  if (dnumber_a*dmac<= 0) then
                    dnumber_a=dmac/1.0e-18/dt
                  end if
                  call which_mode(dmac, dnumber_a,                                 &
                       aerophys(k)%rd(aero_index%i_accum), aerophys(k)%rd(aero_index%i_coarse), &
                       aerochem(k)%density(aero_index%i_accum),     &
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
          end if
        else  ! Nothing significant here to remove - the tidying routines will deal with this
          dmass=0.0 ! no need to do anything since this is now just numerical noise
          dnumber=0.0 ! we assume no change in number during evap
          dnccn_all=0.0 ! we assume no change in number during evap
          dmac=0.0 ! No aerosol processing required
          dmac_all=0.0 ! No aerosol processing required
          dnumber_a=0.0 ! No aerosol processing required
          dnumber_d=0.0 ! No aerosol processing required
        end if
      end if

      if (dmass /= 0.0_wp) then
        this_proc=>procs(k, i_cond%id)

        this_proc%source(i_qv)=-dmass
        this_proc%source(i_ql)=dmass

        if (cloud_params%l_2m) then
          this_proc%source(i_nl)=dnumber
        end if
        nullify(this_proc)

        !============================
        ! aerosol processing
        !============================
        if (l_process) then
          aero_proc=>aerosol_procs(k, i_aact%id)

          aero_proc%source(i_am4)=dmac
          if (l_passivenumbers) aero_proc%source(i_an11)=dnumber_a
          if (l_passivenumbers_ice) aero_proc%source(i_an12)=dnumber_d

          if (aero_index%i_aitken > 0) then
            aero_proc%source(i_am1)=-dmac_all(aero_index%i_aitken)
            aero_proc%source(i_an1)=-dnccn_all(aero_index%i_aitken)
          end if
          if (aero_index%i_accum > 0) then
            aero_proc%source(i_am2)=-dmac_all(aero_index%i_accum)
            aero_proc%source(i_an2)=-dnccn_all(aero_index%i_accum)
          end if
          if (aero_index%i_coarse > 0) then
            aero_proc%source(i_am3)=-dmac_all(aero_index%i_coarse)
            aero_proc%source(i_an3)=-dnccn_all(aero_index%i_coarse)
          end if

          if (.not. l_warm .and. dmad /=0.0) then
            ! We may have some dust in the liquid...
            aero_proc%source(i_am9)=dmad
            aero_proc%source(i_am6)=-dmad ! < USING COARSE
            aero_proc%source(i_an6)=-dnumber_d ! < USING COARSE
          end if
          nullify(aero_proc)
        end if
      end if
    end if    
  end subroutine condevp
end module condensation
