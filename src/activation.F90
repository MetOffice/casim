MODULE activation

USE variable_precision, ONLY: wp
USE mphys_constants, ONLY: fixed_cloud_number
USE mphys_parameters, ONLY: C1, K1
USE mphys_switches, ONLY: iopt_act, iopt_rcrit, aero_index,   &
     i_am2, i_an2, i_am4, i_am1, i_an1, i_am3, i_an3, i_am5, i_am6, i_an6, i_am7, &
     l_warm
USE aerosol_routines, ONLY: aerosol_active, aerosol_phys, aerosol_chem   &
     , abdulRazzakGhan2000, invert_partial_moment, upperpartial_moment_logn &
     , invert_partial_moment_approx, invert_partial_moment_betterapprox
USE special, ONLY: erfinv, erf, pi
USE thresholds, ONLY: w_small, nl_tidy, ccn_tidy

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, k_here, i_here
USE runtime, ONLY: time
#elif DEF_MODEL==MODEL_LEM
USE diaghelp_lem, ONLY: k_here, i_here, j_here
#elif DEF_MODEL==MODEL_UM
USE UM_ParCore, ONLY: mype
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here
#endif
IMPLICIT NONE

CONTAINS

SUBROUTINE activate(dt, cloud_mass, cloud_number, w, rho, dnumber, dmac, T, p, &
     aerophys, aerochem, aeroact, dustphys, dustchem, dustliq, dnccn_all, dmac_all, &
     dnumber_d, dmass_d)

REAL(wp), INTENT(IN) :: dt
REAL(wp), INTENT(IN) :: cloud_mass, cloud_number, w, rho, T, p
REAL(wp), INTENT(OUT) ::  dnumber, dmac
TYPE(aerosol_phys), INTENT(IN) :: aerophys
TYPE(aerosol_chem), INTENT(IN) :: aerochem
TYPE(aerosol_active), INTENT(IN) :: aeroact
TYPE(aerosol_phys), INTENT(IN) :: dustphys
TYPE(aerosol_chem), INTENT(IN) :: dustchem
TYPE(aerosol_active), INTENT(IN) :: dustliq
REAL(wp), INTENT(OUT) :: dnccn_all(:),dmac_all(:)
REAL(wp), INTENT(OUT) :: dnumber_d, dmass_d ! activated dust number and mass
REAL(wp) :: active, nccn, Smax, ucrit, rcrit_ARG, rcrit, nccn_active
REAL(wp) :: Nd, rm, sigma, density

REAL :: work, dmac_i, diff
INTEGER :: imode
LOGICAL :: l_useactive

l_useactive=.FALSE.
dmac = 0.0
dnumber_d=0.0
dmass_d=0.0

SELECT CASE(iopt_act)
  CASE default
      ! fixed number
    active = fixed_cloud_number
  CASE(1)
      ! activate 100% aerosol
    active = SUM(aerophys%N(:))
  CASE(2)
      ! simple Twomey law Cs^k expressed as
      ! a function of w (Rogers and Yau 1989)
    active = 0.88*C1**(2.0/(K1+2.0))   &
         * (7.0E-2*(w*100.0)**1.5)**(K1/(K1+2.0)) &
         * 1.0e6/rho
  CASE(3)
      ! Use scheme of Abdul-Razzak and Ghan
    IF (w > w_small .AND. SUM(aerophys%N(:)) > ccn_tidy) THEN
!        call AbdulRazzakGhan(w, cloud_mass, p, T, aerophys, aerochem, nccn, Smax, rcrit_ARG)
      CALL AbdulRazzakGhan2000(w, p, T, aerophys, aerochem, dnccn_all, Smax, aeroact &
           , nccn_active, l_useactive)
      active   = SUM(dnccn_all(:))

      IF (Smax > 0.02 .AND. .NOT. l_warm) THEN
          ! need better model than this. Could do partitioning in the same way as CCN
        dnumber_d = .01*dustphys%N(1)/dt
        dmass_d = dnumber_d*dustphys%M(1)/dustphys%N(1)
      END IF
    ELSE
      active = 0.0
    END IF

#if DEF_MODEL==MODEL_KiD
    CALL save_dg(k_here, Smax, 'Smax', i_dgtime)
    CALL save_dg(k_here, dnumber_d, 'dnumber_d', i_dgtime)
    CALL save_dg(k_here, dmass_d, 'dmass_d', i_dgtime)
    CALL save_dg(k_here, dustphys%N(1), 'dustphys%N(1)', i_dgtime)
    CALL save_dg(k_here, dustphys%M(1), 'dustphys%M(1)', i_dgtime)
    CALL save_dg(k_here, i_here, SUM(aerophys%N(:)), 'sumN_ARG', i_dgtime)
    CALL save_dg(k_here, i_here, rcrit_ARG, 'rcrit_ARG', i_dgtime)
    CALL save_dg(k_here, i_here, active, 'nccn_ARG', i_dgtime)
#endif

END SELECT

IF (active < nl_tidy) THEN

  dnumber   = 0.0
  dmac      = 0.0
  dmac_all  = 0.0
  dnccn_all = 0.0

ELSE
#if DEF_MODEL==MODEL_KiD
  IF (cloud_Number > 0 .AND. iopt_act==3) THEN
    CALL save_dg(k_here, i_here, (cloud_number - nccn_active), 'nccn_active_extra', i_dgtime)
    CALL save_dg(k_here, i_here, (nccn_active/cloud_number), 'nccn_active_fraction', i_dgtime)
  END IF
#endif


      ! Need to make this consistent with all aerosol_options
  DO imode = 1,aero_index%nccn
    Nd      = aerophys%N(imode)

    IF (Nd > ccn_tidy) THEN
      rm      = aerophys%rd(imode)
      sigma   = aerophys%sigma(imode)
      density = aerochem%density(imode)

      SELECT CASE (iopt_rcrit)
        CASE default
          rcrit = invert_partial_moment_betterapprox(dnccn_all(imode), 0.0_wp, Nd, rm, sigma)
        CASE(2)
          rcrit = invert_partial_moment(Nd, dnccn_all(imode), 0.0_wp, rm, sigma)
        CASE(3)
          rcrit = invert_partial_moment_approx(dnccn_all(imode)/Nd, 0.0_wp, rm, sigma)
      END SELECT

      dmac_all(imode) = (4.0*pi*density/3.0)*(upperpartial_moment_logn(Nd, rm, sigma, 3.0_wp, rcrit))
      dmac_all(imode) = MIN(dmac_all(imode),0.999*aerophys%M(imode)) ! Don't remove more than 99.9%
      dmac = dmac + dmac_all(imode)
    END IF

  END DO
      ! Need to make this consistent with all aerosol_options
  dmac      = dmac/dt
  dmac_all  = dmac_all/dt
  dnccn_all = dnccn_all/dt
  IF (l_useactive) THEN
    dnumber   = active/dt
  ELSE
    dnumber   = MAX(0.0_wp,(active - cloud_number)/dt)
  END IF

END IF

END SUBROUTINE activate

END MODULE activation
