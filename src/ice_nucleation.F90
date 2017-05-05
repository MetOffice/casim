module ice_nucleation
  use mphys_die, only: throw_mphys_error
  use variable_precision, only: wp
  use passive_fields, only: rho, pressure, w, exner
  use mphys_switches, only: i_qv, i_ql, i_qi, i_ni, i_th , hydro_complexity, i_am4, i_am6, i_an2, l_2mi, l_2ms, l_2mg, &
       i_am8, i_am9, aerosol_option, i_nl, i_ns, i_ng, iopt_inuc, i_am7, i_an6, i_an12, l_process, l_passivenumbers, &
       l_passivenumbers_ice, active_number, active_ice, isol, iinsol, l_itotsg, contact_efficiency, immersion_efficiency
  use process_routines, only: process_rate, i_inuc, i_dnuc
  use mphys_parameters, only: nucleated_ice_mass, cloud_params, ice_params
  use mphys_constants, only: Ls, cp, pi, m3_to_cm3
  use qsat_funs, only: qsaturation, qisaturation
  use thresholds, only: ql_small, w_small, ni_tidy, nl_tidy
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active

  implicit none
contains

  !> Currently this routine considers heterogeneous nucleation
  !> notionally as a combination of deposition and/or condensation freezing.
  !> Immersion (i.e. freezing through preeixisting resident IN within a cloud drop) and
  !> Contact freezing (i.e. collision between cloud drop and IN) are not
  !> yet properly concidered.  Such freezing mechanisms should consider the
  !> processing of the aerosol in different ways.
  subroutine inuc(dt, k, qfields, procs, dustphys, dustchem, aeroact, dustliq, aerosol_procs)
    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    real(wp), intent(in), target :: qfields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! aerosol fields
    type(aerosol_phys), intent(in) :: dustphys(:)
    type(aerosol_chem), intent(in) :: dustchem(:)
    type(aerosol_active), intent(in) :: aeroact(:)
    type(aerosol_active), intent(in) :: dustliq(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), optional, target :: aerosol_procs(:,:)

    real(wp) :: dmass, dnumber, dmad, dmac, dmadl

    ! Liquid water and ice saturation for Meyers equation
    real(wp) :: lws_meyers, is_meyers

    !coefficients for Demott parametrization
    real(wp) :: a_demott, b_demott, c_demott, d_demott
    real(wp) :: Tp01 ! 273.16-Tk (or 0.01 - Tc)
    real(wp) :: Tc ! local temperature in C
    real(wp) :: Tk ! local temperature in K
    real(wp) :: th
    real(wp) :: qv
    real(wp) :: ice_number
    real(wp) :: cloud_number, cloud_mass
    real(wp) :: qs, qis, Si, Sw, limit, dN_imm, dN_contact, ql

    ! parameters for Meyers et al (1992)
    ! Meyers MP, DeMott PJ, Cotton WR (1992) New primary ice-nucleation
    ! parameterizations in an explicit cloud model. J Appl Meteorol 31:708–721
    real(wp), parameter :: meyers_a = -0.639 ! Meyers eq 2.4 coeff a
    real(wp), parameter :: meyers_b = 0.1296 ! Meyers eq 2.4 coeff b


    type(process_rate), pointer :: this_proc
    type(process_rate), pointer :: aero_proc

    integer :: kk
    logical :: l_condition

    qv=qfields(k, i_qv)
    th=qfields(k, i_th)

    ql=qfields(k, i_ql)
    Tk=th*exner(k)
    qs=qsaturation(Tk, pressure(k)/100.0)
    qis=qisaturation(Tk, pressure(k)/100.0)
    cloud_mass=qfields(k, i_ql)
    if (cloud_params%l_2m) then
      cloud_number=qfields(k, i_nl)
    else
      cloud_number=cloud_params%fix_N0
    end if

    if (qs==0.0 .or. qis==0.0) then
      call throw_mphys_error(2, 'ice_nucleation', 'Error in saturation calculation - qs or qis is zero')
    end if

    Sw=qv/qs - 1.0
    Si=qv/qis - 1.0
    Tc=Tk - 273.15

    ! What's the condition for ice nucleation...?
    ! This condition needs to be consistent with the mechanisms we're
    ! parametrized
    !l_condition=(( Sw >= -1.e-8 .and. TdegC(k) < -8)) .or. Si > 0.25

    select case(iopt_inuc)
    case default
      l_condition=(( Sw >= -0.001 .and. Tc < -8 .and. Tc > -38) .or. Si >= 0.08)
    case (2) 
      ! Meyers, same condition as DeMott
      l_condition=( cloud_number >= nl_tidy .and. Tc < 0)
    case (4)
      ! DeMott Depletion of dust (contact and immersion)
      l_condition=( cloud_number >= nl_tidy .and. Tc < 0)
    end select

    if (l_condition) then
      this_proc=>procs(k, i_inuc%id)

      if (ice_params%l_2m)then
        ice_number=qfields(k, i_ni)
      else
        ice_number=1.e3 ! PRAGMATIC SM HACK
      end if

      dN_contact=0.0
      dN_imm=0.0
      select case(iopt_inuc)
      case default
        ! Cooper
        ! Cooper WA (1986) Ice Initiation in Natural Clouds. Precipitation Enhancement -
        ! A Scientific Challenge. Meteor Monogr, (Am Meteor Soc, Boston, MA), 21, pp 29–32.
        dN_imm=5.0*exp(-0.304*Tc)/rho(k)
      case (2)
        ! Meyers
        ! Meyers MP, DeMott PJ, Cotton WR (1992) New primary ice-nucleation
        ! parameterizations in an explicit cloud model. J Appl Meteorol 31:708–721
        lws_meyers = 6.112 * exp(17.62*Tc/(243.12 + Tc))
        is_meyers  = 6.112 * exp(22.46*Tc/(272.62 + Tc))
        dN_imm     = 1.0e3 * exp(meyers_a + meyers_b *(100.0*(lws_meyers/is_meyers-1.0)))/rho(k)
        dN_imm     = MAX( dN_imm-ice_number, 0.0 )
        ! Applied just for water saturation, deposition freezing ignored

      case (3)
        ! Fletcher NH (1962) The Physics of Rain Clouds (Cambridge Univ Press, Cambridge, UK)
        dN_imm=0.01*exp(-0.6*Tc)/rho(k)
      case (4)
        ! DeMott Depletion of dust
        ! 'Predicting global atmospheric ice nuclei distributions and their impacts on climate',
        ! Proc. Natnl. Acad. Sci., 107 (25), 11217-11222, 2010, doi:10.1073/pnas.0910818107
        a_demott=5.94e-5
        b_demott=3.33
        c_demott=0.0264
        d_demott=0.0033
        Tp01=0.01-Tc

        if (dustphys(k)%N(1) > ni_tidy) then
          dN_contact=1.0e3/rho(k)*a_demott*(Tp01)**b_demott*                                &
                     (rho(k) * m3_to_cm3 * contact_efficiency*dustphys(k)%N(1))**(c_demott*Tp01+d_demott)
          dN_contact=min(.9*dustphys(k)%N(1), dN_contact)
        end if

        if (dustliq(k)%nact1 > ni_tidy) then
          dN_imm=1.0e3/rho(k)*a_demott*(Tp01)**b_demott*                                    &
                 (rho(k) * m3_to_cm3 * dustliq(k)%nact1)**(c_demott*Tp01+d_demott)
          dN_imm=immersion_efficiency*dN_imm
          dN_imm=min(dustliq(k)%nact, dN_imm)
        end if

      case (5)
        dN_imm=0.0
        dN_contact=max(0.0_wp, (dustphys(k)%N(1)-ice_number))
      end select

      if (cloud_params%l_2m) dN_imm=min(dN_imm, cloud_number)

      dN_imm=dN_imm/dt
      dN_contact=dN_contact/dt
      dnumber=dN_imm + dN_contact
      dnumber=min(dnumber, cloud_number/dt)

      if (dnumber > ni_tidy) then
        dmass=cloud_mass*dnumber/cloud_number
        this_proc%source(i_qi)=dmass

        if (l_2mi) then
          this_proc%source(i_ni)=dnumber
        end if
        if (cloud_params%l_2m) then
          this_proc%source(i_nl)=-dnumber
        end if
        this_proc%source(i_ql)=-dmass

        if (l_process) then
          aero_proc=>aerosol_procs(k, i_dnuc%id)

          ! New ice nuclei
          dmad=dN_contact*dustphys(k)%M(1)/dustphys(k)%N(1)

          ! Frozen soluble aerosol
          dmac=dnumber*aeroact(k)%mact1_mean*aeroact(k)%nratio1

          ! Dust already in the liquid phase
          dmadl=dN_imm*dustliq(k)%mact1_mean*dustliq(k)%nratio1

          aero_proc%source(i_am8)=dmac
          aero_proc%source(i_am4)=-dmac
          aero_proc%source(i_am9)=-dmadl

          aero_proc%source(i_am7)=dmad+dmadl
          aero_proc%source(i_am6)=-dmad    ! <WARNING: using coarse mode
          aero_proc%source(i_an6)=-dN_contact ! <WARNING: using coarse mode

          if (l_passivenumbers_ice) then
            ! we retain information on what'd been nucleated
            aero_proc%source(i_an12)=dN_contact
          end if
          nullify(aero_proc)
        end if
      end if
      nullify(this_proc)
    end if
  end subroutine inuc
end module ice_nucleation
