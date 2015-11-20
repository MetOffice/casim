MODULE sedimentation


USE mphys_die, ONLY: throw_mphys_error
USE variable_precision, ONLY: wp, iwp, defp
USE mphys_parameters, ONLY: nz, hydro_params, cloud_params, rain_params   &
     , ice_params, snow_params, graupel_params
USE passive_fields, ONLY: rho, rdz_on_rho, dz, z_half, z_centre
USE type_process, ONLY: process_name
USE mphys_switches, ONLY: hydro_complexity,   &
     l_abelshipway, l_sed_3mdiff, l_cons      &
     , i_an2, i_am2, i_am4                    &
     , l_ased, i_am5, i_am7, i_am8, i_am9, i_ql,            &
     l_process, l_passivenumbers, l_passivenumbers_ice,   &
     active_cloud, active_rain, isol, iinsol,  &
     active_ice, active_number, i_an11, i_an12, &
     l_separate_rain, l_warm, i_aerosed_method
USE mphys_constants, ONLY: rhow, rho0, cp
USE process_routines, ONLY: process_rate   &
     , i_psedr, i_asedr, i_asedl, i_psedl     &
     , i_pseds, i_psedi, i_psedg, i_dsedi, i_dseds, i_dsedg
USE thresholds, ONLY: ql_small, qr_small, nr_small, m3r_small,   &
     thresh_small
USE special, ONLY: Gammafunc

USE lookup, ONLY: moment
USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
USE aerosol_routines, ONLY: aerosol_active

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, k_here, i_here, nx
USE runtime, ONLY: l_dgstep, time
#elif DEF_MODEL==MODEL_LEM_DIAG
USE diaghelp_lem, ONLY: i_here, j_here, k_here
USE com_params, ONLY: time
USE extra_dgs
USE com_dgstore, ONLY: l_dodgs
#elif DEF_MODEL==MODEL_UM
USE diaghelp_um, ONLY: i_here, j_here, l_debug_um, debug_i, debug_j, debug_pe
USE UM_ParCore, ONLY: mype
USE timestep_mod, ONLY: timestep_number
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here
#endif
IMPLICIT NONE

CHARACTER*(2) :: qchar

CONTAINS

SUBROUTINE sedr(step_length, qfields, aerofields, aeroact, dustact,   &
     tend, params, procs, aerosol_procs, precip, l_doaerosol)

REAL(wp), INTENT(IN) :: step_length
REAL(wp), INTENT(IN), TARGET :: qfields(:,:), aerofields(:,:)
REAL(wp), INTENT(IN) :: tend(:,:)
TYPE(hydro_params), INTENT(IN) :: params
TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustact(:)
    ! NB for liquid phase (cloud/rain) dustact is dustliq
    ! for ice phase (ice/snow/graupel) aeroact is iceact
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)
REAL(wp), INTENT(OUT) :: precip
LOGICAL, OPTIONAL, INTENT(IN) :: l_doaerosol

CALL sedr_aero(step_length, qfields, aerofields, aeroact, dustact,       &
     tend, params, procs, aerosol_procs, precip, l_doaerosol)

END SUBROUTINE sedr

SUBROUTINE sedr_aero(step_length, qfields, aerofields, aeroact, dustact,   &
     tend, params, procs, aerosol_procs, precip, l_doaerosol)

REAL(wp), INTENT(IN) :: step_length
REAL(wp), INTENT(IN), TARGET :: qfields(:,:), aerofields(:,:)
REAL(wp), INTENT(IN) :: tend(:,:)
TYPE(hydro_params), INTENT(IN) :: params
TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustact(:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)
REAL(wp), INTENT(OUT) :: precip
LOGICAL, OPTIONAL, INTENT(IN) :: l_doaerosol

REAL(wp) :: dm1, dm2, dm3
REAL(wp) :: dn1, dn2, dn3, am2

REAL(wp) :: m1, m2, m3
REAL(wp) :: n1, n2, n3

REAL(wp) :: hydro_mass
TYPE(process_rate), POINTER :: this_proc
TYPE(process_rate), POINTER :: aero_proc
REAL(wp), ALLOCATABLE :: flux_n1(:)
REAL(wp), ALLOCATABLE :: flux_n2(:)
REAL(wp), ALLOCATABLE :: flux_n3(:)
REAL(wp), ALLOCATABLE :: Grho(:)
REAL(wp) :: n0, lam, mu, u1r, u2r, u3r, u1r2, u2r2, u3r2
REAL(wp) :: udp ! bulk fallspeed for mass (for precip diag)
INTEGER :: k

INTEGER :: i_1m, i_2m, i_3m
REAL(wp) :: p1, p2, p3
REAL(wp) :: sp1, sp2, sp3
REAL(wp) :: a_x, b_x, f_x, g_x, c_x, d_x
REAL(wp) :: a2_x, b2_x, f2_x
LOGICAL :: l_2m, l_3m

LOGICAL :: l_fluxin, l_fluxout

TYPE(process_name) :: iproc, iaproc  ! processes selected depending on
                                         ! which species we're depositing on.

REAL(wp) :: debug, mint, dmint, scal
REAL(wp) :: dmac, mac_mean_min=1.0e-20
REAL(wp) :: dmad

REAL(wp) :: dnumber_a, dnumber_d

LOGICAL :: l_da_local  ! local tranfer of l_doaerosol

REAL(wp) :: ratio ! used to scale back fluxes if the timestep is too small
                      ! If this is used, you can't trust the results.

l_da_local=.FALSE.
IF (PRESENT(l_doaerosol)) l_da_local=l_doaerosol

    ! precip diag
WRITE(qchar, '(i2.2)') params%id
precip=0.0

m1=0.0
m2=0.0
m3=0.0

ALLOCATE(flux_n1(nz))
flux_n1=0.0
i_1m = params%i_1m
i_2m = params%i_2m
i_3m = params%i_3m
p1 = params%p1
p2 = params%p2
p3 = params%p3

IF (l_sed_3mdiff) THEN
  sp1 = params%sp1
  sp2 = params%sp2
  sp3 = params%sp3
ELSE
  sp1 = params%p1
  sp2 = params%p2
  sp3 = params%p3
END IF

    ! we don't want flexible approach in this version....
IF (p1/=sp1 .OR. p2/=sp2 .OR. p3/=sp3) THEN
  PRINT*, 'Cannot have flexible sedimentation options with CASIM aerosol'
  CALL throw_mphys_error(1, 'sedr_aero', 'Cannot have flexible sedimentation options with CASIM aerosol')
END IF

a_x = params%a_x
b_x = params%b_x
f_x = params%f_x
g_x = params%g_x
c_x = params%c_x
d_x = params%d_x
a2_x = params%a2_x
b2_x = params%b2_x
f2_x = params%f2_x

IF (params%l_2m) THEN
  ALLOCATE(flux_n2(nz))
  flux_n2=0.0
END IF
IF (params%l_3m) THEN
  ALLOCATE(flux_n3(nz))
  flux_n3=0.0
END IF

ALLOCATE(Grho(nz)) ! may want to move this somewhere else
Grho = (rho0/rho)**g_x

dmint = 0.0
SELECT CASE (params%id)
  CASE (1_iwp) !cloud
    iproc=i_psedl
    iaproc=i_asedl
  CASE (2_iwp) !rain
    iproc=i_psedr
    iaproc=i_asedr
  CASE (3_iwp) !ice
    iproc=i_psedi
    iaproc=i_dsedi
  CASE (4_iwp) !snow
    iproc=i_pseds
    iaproc=i_dseds
  CASE (5_iwp) !graupel
    iproc=i_psedg
    iaproc=i_dsedg
END SELECT

DO k=nz-1,1,-1
#if DEF_MODEL==MODEL_KiD
  k_here=k
#endif
#if DEF_MODEL==MODEL_LEM_DIAG
  k_here=k
#endif
      ! initialize to zero
  dm1=0.0
  dm2=0.0
  dm3=0.0
  dn1=0.0
  dn2=0.0
  dn3=0.0
  n1=0.0
  n2=0.0
  n3=0.0
  m1=0.0
  m2=0.0
  m3=0.0

  this_proc => procs(k, iproc%id)
  IF (l_da_local)aero_proc => aerosol_procs(k, iaproc%id)

  hydro_mass = qfields(k, params%i_1m)
  IF (params%l_2m)m2 = qfields(k, params%i_2m)
  IF (params%l_3m)m3 = qfields(k, params%i_3m)

  l_fluxin=.FALSE.
  l_fluxout=.FALSE.
  IF (hydro_mass > thresh_small(params%i_1m))l_fluxout=.TRUE.
  IF (qfields(k+1, params%i_1m) > thresh_small(params%i_1m))l_fluxin=.TRUE.

  IF (l_fluxout) THEN
    m1 = (hydro_mass/c_x)

    n0 = dist_n0(k,params%id)
    mu = dist_mu(k,params%id)
    lam = dist_lambda(k,params%id)

    IF (l_sed_3mdiff) THEN
          !
          ! Moment transfer
          !
      n1 = moment(n0, lam, mu, sp1)
      n2 = moment(n0, lam, mu, sp2)
      n3 = moment(n0, lam, mu, sp3)
    ELSE
      n1=m1
      n2=m2
      n3=m3
    END IF

    u1r=a_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f_x)**(-(1.0+mu+sp1+b_x)))     &
           *(Gammafunc(1.0+mu+sp1+b_x)/Gammafunc(1.0+mu+sp1))

    IF (params%l_2m)     &
           u2r=a_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f_x)**(-(1.0+mu+sp2+b_x))) &
           *(Gammafunc(1.0+mu+sp2+b_x)/Gammafunc(1.0+mu+sp2))

    IF (params%l_3m)     &
           u3r=a_x*Grho(k)*(lam**(1.0+mu+sp3)*(lam+f_x)**(-(1.0+mu+sp3+b_x))) &
           *(Gammafunc(1.0+mu+sp3+b_x)/Gammafunc(1.0+mu+sp3))

    IF (l_abelshipway .AND. params%id==rain_params%id) THEN ! rain can use abel and shipway formulation

      u1r2=a2_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f2_x)**(-(1.0+mu+sp1+b2_x)))   &
             *(Gammafunc(1.0+mu+sp1+b2_x)/Gammafunc(1.0+mu+sp1))
      u1r=u1r+u1r2

      IF (params%l_2m) THEN
        u2r2=a2_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f2_x)**(-(1.0+mu+sp2+b2_x))) &
               *(Gammafunc(1.0+mu+sp2+b2_x)/Gammafunc(1.0+mu+sp2))
        u2r=u2r+u2r2
      END IF

      IF (params%l_3m) THEN
        u3r2=a2_x*Grho(k)*(lam**(1.0+mu+sp3)*(lam+f2_x)**(-(1.0+mu+sp3+b2_x))) &
             *(Gammafunc(1.0+mu+sp3+b2_x)/Gammafunc(1.0+mu+sp3))
        u3r=u3r+u3r2
      END IF
    END IF

        ! fall speeds shouldn't get too big...
    u1r=MIN(u1r,params%maxv)

    ! For clouds and ice, we only use a 1M representation of sedimentation
    ! so that spurious size sorting doesn't lead to overactive autoconversion
    IF (params%id==cloud_params%id .or. ice_params%id==cloud_params%id)then
      IF (params%l_2m)u2r=u1r
    ELSE
      IF (params%l_2m)u2r=MIN(u2r,params%maxv)
      IF (params%l_3m)u3r=MIN(u3r,params%maxv)
    END IF

        ! fall speeds shouldn't be negative (can happen with original AS formulation)
    u1r=MAX(u1r,0.0_wp)
    IF (params%l_2m)u2r=MAX(u2r,0.0_wp)
    IF (params%l_3m)u3r=MAX(u3r,0.0_wp)

#if DEF_MODEL==MODEL_KiD
    IF (nx==1) THEN
      CALL save_dg(k, u1r, 'u1r_'//qchar, i_dgtime)
      IF (params%l_2m)CALL save_dg(k, u2r, 'u2r_'//qchar, i_dgtime)
      IF (params%l_3m)CALL save_dg(k, u3r, 'u3r_'//qchar, i_dgtime)
    ELSE
      CALL save_dg(k, i_here, u1r, 'u1r_'//qchar, i_dgtime)
      IF (params%l_2m)CALL save_dg(k, i_here, u2r, 'u2r_'//qchar, i_dgtime)
      IF (params%l_3m)CALL save_dg(k, i_here, u3r, 'u3r_'//qchar, i_dgtime)
    END IF
#elif DEF_MODEL==MODEL_LEM_DIAG
    IF (l_dodgs .AND. j_here>=jstartdg .AND.     &
           j_here<=jenddg) THEN
      IF (params%id == rain_params%id) THEN
        idgproc = req_dgproc('VR')
      ELSE IF (params%id == ice_params%id) THEN
        idgproc = req_dgproc('VI')
      ELSE IF (params%id == snow_params%id) THEN
        idgproc = req_dgproc('VS')
      ELSE IF (params%id == graupel_params%id) THEN
        idgproc = req_dgproc('VG')
      END IF
      IF (idgproc > 0) THEN
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = u1r
      END IF
      IF (params%id == rain_params%id) THEN
        idgproc = req_dgproc('VRN')
      ELSE IF (params%id == ice_params%id) THEN
        idgproc = req_dgproc('VIN')
      ELSE IF (params%id == snow_params%id) THEN
        idgproc = req_dgproc('VSN')
      ELSE IF (params%id == graupel_params%id) THEN
        idgproc = req_dgproc('VGN')
      END IF
      IF (idgproc > 0) THEN
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = u2r
      END IF
      IF (params%id == rain_params%id) THEN
        idgproc = req_dgproc('LAMR')
      ELSE IF (params%id == ice_params%id) THEN
        idgproc = req_dgproc('LAMI')
      ELSE IF (params%id == snow_params%id) THEN
        idgproc = req_dgproc('LAMS')
      ELSE IF (params%id == graupel_params%id) THEN
        idgproc = req_dgproc('LAMG')
      END IF
      IF (idgproc > 0) THEN
        dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = lam
      END IF
    END IF
#endif

    flux_n1(k)=n1*u1r

    IF (params%l_2m)     &
           flux_n2(k)=n2*u2r

    IF (params%l_3m)     &
           flux_n3(k)=n3*u3r

#if DEF_MODEL==MODEL_UM
        !===========================================================
        ! Cludge to slow down the rain if CFL condition is broken
        ! Should be replaced with a semi-lagrangian scheme
        !===========================================================
    IF (k>1) THEN
      IF (c_x*flux_n1(k)*rdz_on_rho(k)*step_length > 0.99*hydro_mass .OR.      &
                flux_n2(k)*rdz_on_rho(k)*step_length > 0.99*m2) THEN

        ratio=MIN(0.95*hydro_mass/c_x/flux_n1(k)/rdz_on_rho(k)/step_length, 0.95*m2/flux_n2(k)/rdz_on_rho(k)/step_length)

        PRINT*, 'WARNING: slowing down precip - your timestep is too long'
        PRINT*, 'INFO:', timestep_number, i_here, j_here, k, step_length,      &
                   ratio, qfields(k, params%i_1m), c_x*flux_n1(k)*rdz_on_rho(k)&
                   , qfields(k, params%i_2m), flux_n2(k)*rdz_on_rho(k), u1r, u2r

        flux_n1(k)=flux_n1(k)*ratio

        IF (params%l_2m) flux_n2(k)=flux_n2(k)*ratio

        IF (params%l_3m) flux_n3(k)=flux_n3(k)*ratio
      END IF
    END IF
#endif

        ! diagnostic for precip
    IF (k==1) THEN
      udp=a_x*Grho(k)*(lam**d_x*(lam+f_x)**(-(d_x+b_x)))     &
             *(Gammafunc(1.0+mu+d_x+b_x)/Gammafunc(1.0+mu+d_x))
      precip = (rho(k)*hydro_mass)*udp
    END IF

  END IF

  dmac=0.0
  dmad=0.0
  dnumber_a=0.0
  dnumber_d=0.0

  IF (l_fluxout) THEN !flux out (flux(k+1) will be zero if no flux in)
    dn1 = (flux_n1(k+1) - flux_n1(k))* rdz_on_rho(k)

    IF (params%l_2m)     &
           dn2 = (flux_n2(k+1) - flux_n2(k)) * rdz_on_rho(k)

    IF (params%l_3m)     &
           dn3 = (flux_n3(k+1) - flux_n3(k)) * rdz_on_rho(k)

        !============================
        ! aerosol processing
        !============================
    IF (l_ased .AND. l_da_local) THEN
      IF (params%id == cloud_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean -    &
                 flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)* rdz_on_rho(k)
        IF (l_passivenumbers) THEN
          dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio1 - flux_n2(k)*aeroact(k)%nratio1)* rdz_on_rho(k)
        END IF
        IF (.NOT. l_warm) THEN
          dmad = (flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean -  &
                   flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)* rdz_on_rho(k)
          IF (l_passivenumbers_ice .AND. dustact(k)%mact_mean > 0.0) THEN
            dnumber_d = (flux_n2(k+1)*dustact(k+1)%nratio1 -      &
                      flux_n2(k)*dustact(k)%nratio1)* rdz_on_rho(k)
          END IF
        END IF
      ELSE IF (params%id == rain_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean -    &
                 flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)* rdz_on_rho(k)
        IF (l_passivenumbers .AND. aeroact(k)%mact_mean > 0.0) THEN
          dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio2 - flux_n2(k)*aeroact(k)%nratio2)* rdz_on_rho(k)
        END IF
        IF (.NOT. l_warm) THEN
          dmad = (flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean -  &
                    flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)* rdz_on_rho(k)
          IF (l_passivenumbers_ice) THEN
            dnumber_d = (flux_n2(k+1)*dustact(k+1)%nratio2 -    &
                    flux_n2(k)*dustact(k)%nratio2)* rdz_on_rho(k)
          END IF
        END IF
      END IF

      IF (params%id == ice_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean -    &
                 flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)* rdz_on_rho(k)
        dmad = (flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean -    &
                 flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)* rdz_on_rho(k)
        IF (l_passivenumbers) THEN
          dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio1 -      &
                    flux_n2(k)*aeroact(k)%nratio1)* rdz_on_rho(k)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d = (flux_n2(k+1)*dustact(k+1)%nratio1 -     &
                   flux_n2(k)*dustact(k)%nratio1)* rdz_on_rho(k)
        END IF
      ELSE IF (params%id == snow_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean -    &
                 flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)* rdz_on_rho(k)
        dmad = (flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean -    &
                 flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)* rdz_on_rho(k)
        IF (l_passivenumbers) THEN
          dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio2 - flux_n2(k)*aeroact(k)%nratio2)* rdz_on_rho(k)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d = (flux_n2(k+1)*dustact(k+1)%nratio2 - flux_n2(k)*dustact(k)%nratio2)* rdz_on_rho(k)
        END IF
      ELSE IF (params%id == graupel_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean -    &
                  flux_n2(k)*aeroact(k)%nratio3*aeroact(k)%mact3_mean)* rdz_on_rho(k)
        dmad = (flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean -    &
                  flux_n2(k)*dustact(k)%nratio3*dustact(k)%mact3_mean)* rdz_on_rho(k)
#if DEF_MODEL==MODEL_KiD
        CALL save_dg(k, dmad, 'dmac_gs', i_dgtime)
        CALL save_dg(k, flux_n2(k+1)*dustact(k+1)%mact3_mean, 'flux_in', i_dgtime)
        CALL save_dg(k, flux_n2(k)*dustact(k)%mact3_mean, 'flux_out', i_dgtime)
        CALL save_dg(k, n2*dustact(k)%mact3_mean, 'n2Xmean3', i_dgtime)
        CALL save_dg(k, dustact(k)%mact3, 'mact', i_dgtime)
#endif
        IF (l_passivenumbers) THEN
          dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio3 - flux_n2(k)*aeroact(k)%nratio3)* rdz_on_rho(k)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d = (flux_n2(k+1)*dustact(k+1)%nratio3 - flux_n2(k)*dustact(k)%nratio3)* rdz_on_rho(k)
        END IF

      END IF

    END IF


  ELSE IF (l_fluxin) THEN !flux in, but not out

    dn1 = flux_n1(k+1)* rdz_on_rho(k)

    IF (params%l_2m)     &
           dn2 = flux_n2(k+1)* rdz_on_rho(k)

    IF (params%l_3m)     &
           dn3 = flux_n3(k+1)* rdz_on_rho(k)


        !============================
        ! aerosol processing
        !============================
    IF (l_ased .AND. l_da_local) THEN
      IF (params%id == cloud_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)* rdz_on_rho(k)
        IF (l_passivenumbers) THEN
          dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio1)* rdz_on_rho(k)
        END IF
        IF (.NOT. l_warm) THEN
          dmad = (flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)* rdz_on_rho(k)
          IF (l_passivenumbers_ice) THEN
            dnumber_d = (flux_n2(k+1)*dustact(k+1)%nratio1)* rdz_on_rho(k)
          END IF
        END IF
      ELSE IF (params%id == rain_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)* rdz_on_rho(k)
        IF (l_passivenumbers) THEN
          dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio2)* rdz_on_rho(k)
        END IF
        IF (.NOT. l_warm) THEN
          dmad = (flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)* rdz_on_rho(k)
          IF (l_passivenumbers_ice) THEN
            dnumber_d = flux_n2(k+1)*dustact(k+1)%nratio2* rdz_on_rho(k)
          END IF
        END IF
      END IF

      IF (params%id == ice_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)* rdz_on_rho(k)
        dmad = (flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)* rdz_on_rho(k)
        IF (l_passivenumbers) THEN
          dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio1)* rdz_on_rho(k)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d = (flux_n2(k+1)*dustact(k+1)%nratio1)* rdz_on_rho(k)
        END IF
      ELSE IF (params%id == snow_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)* rdz_on_rho(k)
        dmad = (flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)* rdz_on_rho(k)
        IF (l_passivenumbers) THEN
          dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio2)* rdz_on_rho(k)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d = (flux_n2(k+1)*dustact(k+1)%nratio2)* rdz_on_rho(k)
        END IF
      ELSE IF (params%id == graupel_params%id) THEN
        IF (i_aerosed_method==1) THEN
          dmac = (flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean)* rdz_on_rho(k)
          dmad = (flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean)* rdz_on_rho(k)
          IF (l_passivenumbers) THEN
            dnumber_a = (flux_n2(k+1)*aeroact(k+1)%nratio3)* rdz_on_rho(k)
        END IF
        IF (l_passivenumbers_ice) THEN
            dnumber_d = (flux_n2(k+1)*dustact(k+1)%nratio3)* rdz_on_rho(k)
          END IF
        ELSE
          PRINT*, 'ERROR: GET RID OF i_aerosed_method variable!'
        END IF

      END IF

    END IF

  END IF

      ! Store the aerosol process terms...
  IF (l_da_local) THEN
    IF (params%id == cloud_params%id .OR. params%id == rain_params%id) THEN
          !liquid phase
      IF (l_separate_rain .AND. params%id == rain_params%id) THEN
        aero_proc%source(i_am5) = dmac
      ELSE
        aero_proc%source(i_am4) = dmac
      END IF
      IF (.NOT. l_warm)aero_proc%source(i_am9) = dmad
      IF (l_passivenumbers) THEN
        aero_proc%source(i_an11) = dnumber_a
      END IF
      IF (l_passivenumbers_ice) THEN
        aero_proc%source(i_an12) = dnumber_d
      END IF

    ELSE
          !ice phase
      aero_proc%source(i_am7) = dmad
      aero_proc%source(i_am8) = dmac
      IF (l_passivenumbers) THEN
        aero_proc%source(i_an11) = dnumber_a
      END IF
      IF (l_passivenumbers_ice) THEN
        aero_proc%source(i_an12) = dnumber_d
      END IF
    END IF
  END IF

#if DEF_MODEL==MODEL_KiD

  IF (nx==1) THEN
    CALL save_dg(k, dmac, 'dmac', i_dgtime)
    CALL save_dg(k, (c_x*flux_n1(k+1)*9.8/cp), 'frictional_heating'//qchar, i_dgtime)
  ELSE
    CALL save_dg(k, i_here, dmac, 'dmac', i_dgtime)
    CALL save_dg(k, i_here, (c_x*flux_n1(k+1)*9.8/cp), 'frictional_heating'//qchar, i_dgtime)
  END IF
#endif

  dm1=dn1
  dm2=dn2
  dm3=dn3

  this_proc%source(params%i_1m) = c_x*dm1

  IF (params%l_2m) THEN
    this_proc%source(params%i_2m) = dm2
  END IF

  IF (params%l_3m) THEN
    this_proc%source(params%i_3m) = dm3
  END IF


  NULLIFY(this_proc)
  IF (l_da_local)NULLIFY(aero_proc)

#if DEF_MODEL==MODEL_UM
  IF (l_debug_um .AND. mype==debug_pe .AND. i_here==debug_i .AND. j_here==debug_j .AND. (l_fluxin .OR. l_fluxout)) THEN
    PRINT*, 'DEBUG sed' , timestep_number, k, l_fluxin, l_fluxout, qfields(k, params%i_1m), qfields(k, params%i_2m), c_x*dm1, dm2, n0, lam, mu, flux_n1(k), rdz_on_rho(k)
  END IF
#endif

END DO


DEALLOCATE(Grho)

DEALLOCATE(flux_n1)
IF (params%l_2m)DEALLOCATE(flux_n2)
IF (params%l_3m)DEALLOCATE(flux_n3)

END SUBROUTINE sedr_aero

END MODULE sedimentation
