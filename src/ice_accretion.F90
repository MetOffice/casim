module ice_accretion
  use variable_precision, only: wp, iwp
  use passive_fields, only: rho, pressure, w, exner, TdegC
  use mphys_switches, only: hydro_complexity,  i_ql, l_process,   &
       cloud_params, rain_params, i_am4, i_am7, i_am8, i_am9
  use type_process, only: process_name
  use process_routines, only: process_rate,   &
       i_iacw, i_sacw, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs, &
       i_diacw, i_dsacw, i_dgacw, i_dsacr, i_dgacr, i_draci
  use mphys_parameters, only: hydro_params
  use mphys_constants, only: Ls, cp
  use qsat_funs, only: qsaturation, qisaturation
  use thresholds, only: thresh_small
  use activation, only: activate
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active

  use m3_incs, only: m3_inc_type2, m3_inc_type3
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use sweepout_rate, only: sweepout, binary_collection

  use special, only: pi, Gammafunc

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='ICE_ACCRETION'

  public iacc
contains

  subroutine iacc(dt, k, params_X, params_Y, params_Z, qfields, procs,   &
       aeroact, dustliq, aerosol_procs)
    !
    !< CODE TIDYING: Move efficiencies into parameters

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    type(hydro_params), intent(in) :: params_X !< parameters for species which does the collecting
    type(hydro_params), intent(in) :: params_Y !< parameters for species which is collected
    type(hydro_params), intent(in) :: params_Z !< parameters for species to which resulting amalgamation is sent
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! aerosol fields
    type(aerosol_active), intent(in) :: aeroact(:)
    type(aerosol_active), intent(in) :: dustliq(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)

    ! Local variables
    type(process_name) :: iproc, iaproc  ! processes selected depending on
    ! which species we're depositing on.

    real(wp) :: dmass_Y, dnumber_Y, dmac, dmad
    real(wp) :: dmass_X, dmass_Z, dnumber_X, dnumber_Z

    real(wp) :: number_X, mass_X, m3_X, m1, m2, m3
    real(wp) :: mass_Y, number_Y, m3_Y
    type(process_rate), pointer :: this_proc, aero_proc

    real(wp) :: n0_X, lam_X, mu_X
    real(wp) :: n0_Y, lam_Y, mu_Y
    real(wp) :: dm1, dm2

    real(wp) :: Eff ! collection efficiencies need to re-evaluate these and put them in properly to mphys_parameters

    logical :: l_condition, l_alternate_Z, l_freezeall_X, l_freezeall_y
    logical :: l_slow ! fall speed is slow compared to interactive species
    logical :: l_aero ! If true then this process will modify aerosol

    character(len=*), parameter :: RoutineName='IACC'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    mass_Y=qfields(k, params_Y%i_1m)
    mass_X=qfields(k, params_X%i_1m)

    l_condition=mass_X > thresh_small(params_X%i_1m) .and. mass_Y > thresh_small(params_Y%i_1m)
    l_freezeall_X=.false.
    l_freezeall_Y=.false.

    if (l_condition) then
      ! initialize variables which may not have been set
      number_X=0.0
      m3_X=0.0
      dnumber_Y=0.0

      l_alternate_Z=params_X%id /= params_Z%id
      l_aero=.false.

      select case (params_Y%id)
      case(1_iwp) ! cloud liquid is collected
        l_slow=.true. ! Y catagory falls slowly compared to X
        select case (params_X%id)
        case (3_iwp) !ice collects
          iproc=i_iacw
          Eff=1.0
          iaproc=i_diacw
          l_aero=.true.
        case (4_iwp) !snow collects
          iproc=i_sacw
          Eff=1.0
          iaproc=i_dsacw
          l_aero=.true.
        case (5_iwp) !graupel collects
          iproc=i_gacw
          Eff=1.0
          iaproc=i_dgacw
          l_aero=.true.
        end select
      case(2_iwp) ! rain is collected
        l_slow=.false. ! Y catagory falls slowly compared to X
        select case (params_X%id)
        case (4_iwp) !snow collects
          iproc=i_sacr
          Eff=1.0
          iaproc=i_dsacr
          l_aero=.true.
        case (5_iwp) !graupel collects
          iproc=i_gacr
          Eff=1.0
          iaproc=i_dgacr
          l_aero=.true.
        end select
      case(3_iwp)! cloud ice is collected
        l_slow=.true. ! Y catagory falls slowly compared to X
        select case (params_X%id)
        case (2_iwp) !rain collects
          iproc=i_raci
          Eff=1.0
          iaproc=i_draci
          l_aero=.true.
        case (4_iwp) !snow collects
          iproc=i_saci
          !Eff=1.0
          Eff = min(1.0_wp, 0.2*exp(0.08*TdegC(k)))
        case (5_iwp) !graupel collects
          iproc=i_gaci
          !Eff=1.0
          Eff = min(1.0_wp, 0.2*exp(0.08*TdegC(k)))
        end select
      case(4_iwp) ! snow is collected
        l_slow=.false. ! Y catagory falls slowly compared to X
        select case (params_X%id)
        case (5_iwp) !graupel collects
          iproc=i_gacs
          !Eff=1.0
          Eff = min(1.0_wp, 0.2*exp(0.08*TdegC(k)))
        end select
      end select

      this_proc=>procs(k, iproc%id)

      if (params_X%l_2m) number_X=qfields(k, params_X%i_2m)
      if (params_X%l_3m) m3_X=qfields(k, params_X%i_3m)

      n0_X=dist_n0(k,params_X%id)
      mu_X=dist_mu(k,params_X%id)
      lam_X=dist_lambda(k,params_X%id)

      if (l_slow) then  ! collected species is approximated to have zero fallspeed
        dmass_Y=-Eff*sweepout(n0_X, lam_X, mu_X, params_X, rho(k))*mass_Y
        dmass_Y=max(dmass_Y, -mass_Y/dt)

        if (params_Y%l_2m) then
          number_Y=qfields(k, params_Y%i_2m)
          dnumber_Y=dmass_Y*number_Y/mass_Y
        end if

        if (l_alternate_Z) then ! We move resulting collision to another species.
          if (params_Y%l_2m) then
            number_Y=qfields(k, params_Y%i_2m)
          else
            ! we need an else in here
          end if
          dmass_X=-Eff*sweepout(n0_X, lam_X, mu_X, params_X, rho(k), mass_weight=.true.) * number_Y
          dmass_X=max(dmass_X, -mass_X/dt)
          dnumber_X=dnumber_Y !dmass_Y*number_X/mass_Y
          dmass_Z=-1.0*(dmass_X + dmass_Y)
          dnumber_Z=-1.0*dnumber_X
        end if
      else  ! both species have significant fall velocity
        if (params_Y%l_2m) then
          number_Y=qfields(k, params_Y%i_2m)
        end if
        if (params_Y%l_3m) m3_Y=qfields(k, params_Y%i_3m)

        n0_Y=dist_n0(k,params_Y%id)
        mu_Y=dist_mu(k,params_Y%id)
        lam_Y=dist_lambda(k,params_Y%id)

        dmass_Y=-Eff*binary_collection(n0_X, lam_X, mu_X, n0_Y, lam_Y, mu_Y, params_X, params_Y, rho(k), mass_weight=.true.)
        dmass_Y=max(dmass_Y, -mass_Y/dt)

        if (params_Y%l_2m) then
          dnumber_Y=-Eff*binary_collection(n0_X, lam_X, mu_X, n0_Y, lam_Y, mu_Y, params_X, params_Y, rho(k))
        end if

        if (l_alternate_Z) then ! We move resulting collision to another species.
          dmass_X=-Eff*binary_collection(n0_Y, lam_Y, mu_Y, n0_X, lam_X, mu_X, params_Y, params_X, rho(k), mass_weight=.true.)
          dmass_X=max(dmass_X, -mass_X/dt)
          dnumber_X=dnumber_Y
          dmass_Z=-1.0*(dmass_X + dmass_Y)
          dnumber_Z=-1.0*dnumber_X
        end if
      end if

      ! PRAGMATIC HACK - FIX ME (ALTHOUGH I QUITE LIKE IT)
      ! If most of the collected particles are to be removed then remove all of them
      if (-dmass_Y*dt >0.95*mass_Y .or. (params_Y%l_2m .and. -dnumber_Y*dt > 0.95*number_Y)) then
        dmass_Y=-mass_Y/dt
        dnumber_Y=-number_Y/dt
        l_freezeall_Y=.true.
      end if
      ! If most of the collecting particles are to be removed then remove all of them
      if (l_alternate_Z .and. (-dmass_X*dt >0.95*mass_X .or. (params_X%l_2m .and. -dnumber_X*dt > 0.95*number_X))) then
        dmass_X=-mass_X/dt
        dnumber_X=-number_X/dt
        l_freezeall_X=.true.
        dmass_Z = -1.0*( dmass_X + dmass_Y )
        dnumber_Z = -1.0*dnumber_X
      end if

      if (.not. l_alternate_Z) then
        dmass_X=-dmass_Y
        dnumber_X=0.0
        dnumber_Z=0.0
        dmass_Z=0.0
      end if

      this_proc%source(params_Y%i_1m)=dmass_Y
      if (params_Y%l_2m) then
        this_proc%source(params_Y%i_2m)=dnumber_Y
      end if

      this_proc%source(params_X%i_1m)=dmass_X
      if (params_X%l_2m) then
        this_proc%source(params_X%i_2m)=dnumber_X
      end if

      if (l_alternate_Z) then
        this_proc%source(params_Z%i_1m)=dmass_Z
        if (params_Z%l_2m) then
          this_proc%source(params_Z%i_2m)=dnumber_Z
        end if
      end if

      !----------------------
      ! Aerosol processing...
      !----------------------

      if (l_process .and. l_aero) then
        ! We note that all processes result in a source of frozen water
        ! i.e. params_Z is either ice, snow or graupel
        aero_proc=>aerosol_procs(k, iaproc%id)

        if (params_Y%id == cloud_params%id) then
          dmac=abs(dnumber_Y)*aeroact(k)%mact1_mean*aeroact(k)%nratio1
          dmad=abs(dnumber_Y)*dustliq(k)%mact1_mean*dustliq(k)%nratio1
        else if (params_Y%id == rain_params%id) then
          dmac=abs(dnumber_Y)*aeroact(k)%mact2_mean*aeroact(k)%nratio2
          dmad=abs(dnumber_Y)*dustliq(k)%mact2_mean*dustliq(k)%nratio2
        else if (params_X%id == cloud_params%id) then ! This is never the case !
          dmac=abs(dnumber_X)*aeroact(k)%mact1_mean*aeroact(k)%nratio1
          dmad=abs(dnumber_X)*dustliq(k)%mact1_mean*dustliq(k)%nratio1
        else if (params_X%id == rain_params%id) then
          dmac=abs(dnumber_X)*aeroact(k)%mact2_mean*aeroact(k)%nratio2
          dmad=abs(dnumber_X)*dustliq(k)%mact2_mean*dustliq(k)%nratio2
        end if

        aero_proc%source(i_am8)=dmac
        aero_proc%source(i_am4)=-dmac
        aero_proc%source(i_am7)=dmad
        aero_proc%source(i_am9)=-dmad
        nullify(aero_proc)
      end if
      nullify(this_proc)
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine iacc
end module ice_accretion
