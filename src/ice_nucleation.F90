MODULE ice_nucleation

USE mphys_die, ONLY: throw_mphys_error
USE variable_precision, ONLY: wp
USE passive_fields, ONLY: rho, pressure, w, exner
USE mphys_switches, ONLY: i_qv, i_ql, i_qi, i_ni, i_th   &
     , hydro_complexity, i_am4, i_am6, i_an2, l_2mi, l_2ms, l_2mg, i_am8, i_am9&
     , aerosol_option, i_nl, i_ns, i_ng, iopt_inuc, i_am7, i_an6, i_an12 &
     , l_process, l_passivenumbers, l_passivenumbers_ice, active_number  &
     , active_ice, isol, iinsol, l_itotsg &
     , contact_efficiency, immersion_efficiency
USE process_routines, ONLY: process_rate, i_inuc,   &
     i_dnuc
USE mphys_parameters, ONLY: nucleated_ice_mass, cloud_params, ice_params
USE mphys_constants, ONLY: Ls, cp, pi
USE qsat_funs, ONLY: qsaturation, qisaturation
USE thresholds, ONLY: ql_small, w_small, ni_tidy, nl_tidy
USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, k_here
USE runtime, ONLY: time
#elif DEF_MODEL==MODEL_LEM
USE diaghelp_lem, ONLY: i_here, j_here
#elif DEF_MODEL==MODEL_UM
USE diaghelp_um, ONLY: i_here, j_here
USE UM_ParCore, ONLY: mype
USE timestep_mod,          ONLY : timestep_number
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here
#endif
IMPLICIT NONE

CONTAINS

SUBROUTINE inuc(dt, k, qfields, procs, dustphys, dustchem,   &
     aeroact, dustliq, aerosol_procs)

    ! Currently this routine considers heterogeneous nucleation
    ! notionally as a combination of deposition and/or condensation freezing.
    ! Immersion (i.e. freezing through preeixisting resident IN within a cloud drop) and
    ! Contact freezing (i.e. collision between cloud drop and IN) are not
    ! yet properly concidered.  Such freezing mechanisms should consider the
    ! processing of the aerosol in different ways.


REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: k
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! aerosol fields
TYPE(aerosol_phys), INTENT(IN) :: dustphys(:)
TYPE(aerosol_chem), INTENT(IN) :: dustchem(:)
TYPE(aerosol_active), INTENT(IN) :: aeroact(:)
TYPE(aerosol_active), INTENT(IN) :: dustliq(:)

    ! optional aerosol fields to be processed
TYPE(process_rate), INTENT(INOUT), OPTIONAL, TARGET :: aerosol_procs(:,:)

REAL(wp) :: dmass, dnumber, dmad, dmac, dmadl

    !coefficients for Demott parametrization
REAL(wp) :: a_demott, b_demott, c_demott, d_demott

REAL(wp) :: Tp01 ! 273.16-Tk (or 0.01 - Tc)
REAL(wp) :: Tc ! local temperature in C
REAL(wp) :: Tk ! local temperature in K

REAL(wp) :: th
REAL(wp) :: qv
REAL(wp) :: ice_number
REAL(wp) :: cloud_number, cloud_mass
TYPE(process_rate), POINTER :: this_proc
TYPE(process_rate), POINTER :: aero_proc

REAL(wp) :: qs, qis, Si, Sw, limit, dN_imm, dN_contact, ql

INTEGER :: kk

LOGICAL :: l_condition

qv = qfields(k, i_qv)
th = qfields(k, i_th)

ql = qfields(k, i_ql)
Tk = th*exner(k)
qs = qsaturation(Tk, pressure(k)/100.0)
qis = qisaturation(Tk, pressure(k)/100.0)
cloud_mass = qfields(k, i_ql)
IF (cloud_params%l_2m) THEN
  cloud_number = qfields(k, i_nl)
ELSE
  cloud_number=cloud_params%fix_N0
END IF

IF (qs==0.0 .OR. qis==0.0) THEN
#if DEF_MODEL==MODEL_UM
  DO kk=1,k
    PRINT*, 'DEBUG nuc', mype, i_here, j_here, kk, pressure(kk), w(kk), qfields(kk, i_th), qsaturation(qfields(kk, i_th)*exner(kk), pressure(kk)/100.0),qisaturation(qfields(kk, i_th)*exner(kk), pressure(kk)/100.0)
  END DO
#endif
  CALL throw_mphys_error(2, 'ice_nucleation', 'Error in saturation calculation')
END IF

Sw=qv/qs - 1.0
Si=qv/qis - 1.0
Tc=Tk - 273.15

    ! What's the condition for ice nucleation...?
    ! This condition needs to be consistent with the mechanisms we're
    ! parametrized
    !l_condition=(( Sw >= -1.e-8 .and. TdegC(k) < -8)) .or. Si > 0.25

SELECT CASE(iopt_inuc)
  CASE default
    l_condition=(( Sw >= -0.001 .AND. Tc < -8 .AND. Tc > -38) .OR. Si >= 0.08)
  CASE (4)
      ! DeMott Depletion of dust (contact and immersion)
    l_condition=( cloud_number >= nl_tidy .AND. Tc < 0)
END SELECT


IF (l_condition) THEN
  this_proc => procs(k, i_inuc%id)
  
  if (ice_params%l_2m)then
    ice_number = qfields(k, i_ni)
  else
    ice_number = 1.e3 ! PRAGMATIC SM HACK
  end if

  dN_contact=0.0
  dN_imm=0.0
  SELECT CASE(iopt_inuc)
    CASE default
        ! Cooper
        ! Cooper WA (1986) Ice Initiation in Natural Clouds. Precipitation Enhancement -
        ! A Scientific Challenge. Meteor Monogr, (Am Meteor Soc, Boston, MA), 21, pp 29–32.
      dN_imm = 5.0*EXP(-0.304*Tc)/rho(k)
    CASE (2)
        ! Meyers
        ! Meyers MP, DeMott PJ, Cotton WR (1992) New primary ice-nucleation
        ! parameterizations in an explicit cloud model. J Appl Meteorol 31:708–721
      dN_imm = EXP(4.108 - 0.262*Tc)/rho(k)
    CASE (3)
        ! Fletcher NH (1962) The Physics of Rain Clouds (Cambridge Univ Press, Cambridge, UK)
      dN_imm = 0.01*EXP(-0.6*Tc)/rho(k)
    CASE (4)
        ! DeMott Depletion of dust
        ! 'Predicting global atmospheric ice nuclei distributions and their impacts on climate',
        ! Proc. Natnl. Acad. Sci., 107 (25), 11217-11222, 2010, doi:10.1073/pnas.0910818107
      a_demott = 5.94e-5
      b_demott = 3.33
      c_demott = 0.0264
      d_demott = 0.0033
      Tp01 = 0.01 - Tc

      IF (dustphys(k)%N(1) > ni_tidy) THEN
        dN_contact=1.0e3/rho(k)*a_demott*(Tp01)**b_demott*(rho(k)*contact_efficiency*dustphys(k)%N(1))**(c_demott*Tp01+d_demott)
        dN_contact=MIN(.9*dustphys(k)%N(1), dN_contact)
      END IF

      IF (dustliq(k)%nact1 > ni_tidy) THEN
        dN_imm = 1.0e3/rho(k)*a_demott*(Tp01)**b_demott*(rho(k)*dustliq(k)%nact1)**(c_demott*Tp01+d_demott)
        dN_imm=immersion_efficiency*dN_imm
        dN_imm=MIN(dustliq(k)%nact, dN_imm)
      END IF

#if DEF_MODEL==MODEL_KiD
      CALL save_dg(k, dN_contact, 'dn_cnt', i_dgtime)
      CALL save_dg(k, dN_imm, 'dn_imm', i_dgtime)
#endif
    CASE (5)
      dN_imm = 0.0
      dN_contact = MAX(0.0_wp, (dustphys(k)%N(1) - ice_number))
  END SELECT

  IF (cloud_params%l_2m)dN_imm = MIN(dN_imm, cloud_number)

  dN_imm=dN_imm/dt
  dN_contact=dN_contact/dt
  dnumber = dN_imm + dN_contact
  dnumber = MIN(dnumber, cloud_number/dt)

  IF (dnumber > ni_tidy) THEN

    dmass=cloud_mass*dnumber/cloud_number

    this_proc%source(i_qi) = dmass

    IF (l_2mi) THEN
      this_proc%source(i_ni) = dnumber
    END IF
    IF (cloud_params%l_2m) THEN
      this_proc%source(i_nl) = -dnumber
    END IF
    this_proc%source(i_ql) = -dmass

    IF (l_process) THEN

      aero_proc => aerosol_procs(k, i_dnuc%id)

          ! New ice nuclei
      dmad = dN_contact*dustphys(k)%M(1)/dustphys(k)%N(1)

          ! Frozen soluble aerosol
      dmac  = dnumber*aeroact(k)%mact1_mean*aeroact(k)%nratio1

          ! Dust already in the liquid phase
      dmadl = dN_imm*dustliq(k)%mact1_mean*dustliq(k)%nratio1

      aero_proc%source(i_am8) = dmac
      aero_proc%source(i_am4) = -dmac
      aero_proc%source(i_am9) = -dmadl

      aero_proc%source(i_am7) = dmad + dmadl
      aero_proc%source(i_am6) = -dmad    ! <WARNING: using coarse mode
      aero_proc%source(i_an6) = -dN_contact ! <WARNING: using coarse mode

      IF (l_passivenumbers_ice) THEN
            ! we retain information on what'd been nucleated
        aero_proc%source(i_an12) = dN_contact
      END IF


      NULLIFY(aero_proc)

    END IF

  END IF
  NULLIFY(this_proc)
END IF


END SUBROUTINE inuc

END MODULE ice_nucleation
