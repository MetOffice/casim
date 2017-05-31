module graupel_embryo
  use variable_precision, only: wp
  use process_routines, only: process_rate, process_name, i_sacw
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use passive_fields, only: TdegC, TdegK, qws0, rho
  use mphys_parameters, only: ice_params, snow_params, graupel_params, rain_params, cloud_params
  use mphys_constants, only: rho0
  use thresholds, only: thresh_sig
  use special, only: pi, Gammafunc
  use m3_incs, only: m3_inc_type2, m3_inc_type3
  use ventfac, only: ventilation
  use distributions, only: dist_lambda, dist_mu, dist_n0

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='GRAUPEL_EMBRYO'

  public graupel_embryos
contains

  subroutine graupel_embryos(dt, k, qfields, procs, aerophys, aerochem, aeroact, aerosol_procs)

    !< Subroutine to convert some of the small rimed snow to graupel
    !< (Ikawa & Saito 1991)
    !<    !
    !< OPTIMISATION POSSIBILITIES: See gamma functions and distribution calculations
    !<
    !< AEROSOL: NOT DONE YET - internal category transfer

    implicit none

    character(len=*), parameter :: RoutineName='GRAUPEL_EMBRYOS'

    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! aerosol fields
    type(aerosol_phys), intent(in) :: aerophys(:)
    type(aerosol_chem), intent(in) :: aerochem(:)
    type(aerosol_active), intent(in) :: aeroact(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), optional :: aerosol_procs(:,:)

    real(wp) :: cloud_mass, snow_mass, snow_number
    real(wp) :: dmass, dnumber
    real(wp) :: m1, m2, m3, dm1, dm2, dm3
    real(wp) :: snow_n0, snow_lam, snow_mu ! distribution parameters

    real(wp) :: dnembryo ! rate of embryo creation
    real(wp) :: embryo_mass = 1.6e-10 ! mass(kg) of a new graupel embryo (should be in parameters)

    type(process_rate), pointer :: this_proc

    real(wp) :: pgsacw  ! Rate of mass transfer to graupel
    real(wp) :: Eff

    real(wp) :: alpha = 4.0 ! Tunable parameter see Reisner 1998, Ikawa+Saito 1991
    real(wp) :: rhogms      ! Difference between graupel density and snow density

    real(wp) :: Garg ! Argument for Gamma function
    integer  :: pid  ! Process id

    snow_mass=qfields(k, snow_params%i_1m)
    cloud_mass=qfields(k, cloud_params%i_1m)

    if (snow_mass > thresh_sig(snow_params%i_1m) .and. cloud_mass > thresh_sig(cloud_params%i_1m)) then

      if (snow_params%l_2m) snow_number=qfields(k, snow_params%i_2m)

      snow_n0=dist_n0(k,snow_params%id)
      snow_mu=dist_mu(k,snow_params%id)
      snow_lam=dist_lambda(k,snow_params%id)

      !< This efficiency should be the  same as used in sacw calculation, i.e.
      !< they should be defined consistently in the parameters
      Eff=1.0
      rhogms=graupel_params%density-snow_params%density

      Garg=2.0+2*snow_params%b_x + snow_mu
      pgsacw = (0.75*alpha*dt*pi/rhogms)*Eff*rho(k)*rho(k)*cloud_mass*cloud_mass   &
           *snow_params%a_x*snow_params%a_x*snow_n0*GammaFunc(Garg)*(2*snow_params%f_x + 2*snow_lam)**(-Garg) &
           *(rho(k)/rho0)**(2*snow_params%g_x)

      dnembryo=max(snow_params%density*pgsacw/rhogms/embryo_mass/rho(k), 0.0_wp)
      dnumber=min(dnembryo, 0.95*snow_number/dt)

      pid=i_sacw%id
      this_proc=>procs(k, pid)
      dmass=this_proc%source(snow_params%i_1m)-pgsacw
      this_proc%source(snow_params%i_1m)=dmass
      this_proc%source(graupel_params%i_1m)=pgsacw
      if (graupel_params%l_2m) then
        this_proc%source(graupel_params%i_2m)=dnumber
      end if
      if (snow_params%l_2m) then
        this_proc%source(snow_params%i_2m)=-dnumber
      end if

      !----------------------
      ! Aerosol processing...
      !----------------------
      nullify(this_proc)
    end if
  end subroutine graupel_embryos
end module graupel_embryo
