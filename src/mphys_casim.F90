MODULE mphys_casim
#if DEF_MODEL==MODEL_KiD

USE variable_precision, ONLY: wp
USE micro_main, ONLY: shipway_microphysics
USE mphys_switches, ONLY: set_mphys_switches, option, aerosol_option,       &
     nq_l, nq_r, nq_i, nq_s, nq_g                                           &
     , i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am4, i_am5               &
     , i_am6, i_an6, i_am7, i_am8 , i_am9, i_am10, i_an10, i_an11, i_an12
USE initialize, ONLY: mphys_init

  !KiD variables
USE physconst, ONLY : p0, r_on_cp
USE column_variables, ONLY:                              &
     hydrometeors, dhydrometeors_div, dhydrometeors_adv, &
     dtheta_adv, dtheta_div, dqv_adv, dqv_div,           &
     aerosol, daerosol_div, daerosol_adv,                &
     dhydrometeors_mphys, daerosol_mphys,                &
     dtheta_mphys, dqv_mphys,                            &
     exner_in=>exner, w_in=>w,                           &
     theta_in=>theta, qv_in=>qv, dz_in=>dz, rho_in=>rho, &
     z_centre_in=>z

USE parameters, ONLY : nx, nz, dt

IMPLICIT NONE

  !Logical switches
LOGICAL :: micro_unset=.TRUE.

CONTAINS

SUBROUTINE mphys_casim_interface

    !< OPTIMIZATIONS: A bit unwieldy, but should probably put the loops inside the
    ! if tests for each microphysics variable

REAL(wp) :: theta(nz,nx,1), pressure(nz,nx,1),      &
       z_half(0:nz,nx,1), z_centre(nz,nx,1), dz(nz,nx,1), qv(nz,nx,1),qc(nz,nx,1) &
       , nc(nz,nx,1), qr(nz,nx,1), nr(nz,nx,1), m3r(nz,nx,1),rho(nz,nx,1) &
       , exner(nz,nx,1), w(nz,nx,1), tke(nz,nx,1)                              &
       , qi(nz,nx,1), ni(nz,nx,1), qs(nz,nx,1), ns(nz,nx,1), m3s(nz,nx,1) &
       , qg(nz,nx,1), ng(nz,nx,1), m3g(nz,nx,1)

REAL(wp) :: AccumSolMass(nz,nx,1), AccumSolNumber(nz,nx,1) ! Accumulation mode aerosol
REAL(wp) :: ActiveSolLiquid(nz,nx,1)                      ! Activated aerosol
REAL(wp) :: AitkenSolMass(nz,nx,1), AitkenSolNumber(nz,nx,1) ! Aitken mode aerosol
REAL(wp) :: CoarseSolMass(nz,nx,1), CoarseSolNumber(nz,nx,1) ! Course mode aerosol
REAL(wp) :: ActiveSolRain(nz,nx,1)                      ! Activeated aerosol in rain
REAL(wp) :: CoarseDustMass(nz,nx,1), CoarseDustNumber(nz,nx,1) ! Coarse Dust
REAL(wp) :: ActiveInsolIce(nz,nx,1)                      ! Activeated dust
REAL(wp) :: ActiveSolIce(nz,nx,1)                      ! Activeated aerosol in ice
REAL(wp) :: ActiveInsolLiquid(nz,nx,1)                      ! Activeated dust in cloud
REAL(wp) :: AccumInsolMass(nz,nx,1)                      ! Accum mode dust mass
REAL(wp) :: AccumInsolNumber(nz,nx,1)                      ! Accum mode dust number
REAL(wp) :: ActiveSolNumber(nz,nx,1)                      ! Activated soluble number (if we need a tracer)
REAL(wp) :: ActiveInsolNumber(nz,nx,1)                      ! Activated insoluble number (if we need a tracer)


    ! Tendency from other physics/advection/forcing
    ! NB This is tendency (/s) not an increment over the timestep
REAL(wp) :: dqv(nz,nx,1), dth(nz,nx,1), dqc(nz,nx,1), dnc(nz,nx,1)             &
       , dqr(nz,nx,1), dnr(nz,nx,1), dm3r(nz,nx,1)                             &
       , dqi(nz,nx,1), dni(nz,nx,1), dqs(nz,nx,1), dns(nz,nx,1), dm3s(nz,nx,1) &
       , dqg(nz,nx,1), dng(nz,nx,1), dm3g(nz,nx,1)

REAL(wp) :: dAccumSolMass(nz,nx,1), dAccumSolNumber(nz,nx,1) ! Accumulation mode aerosol
REAL(wp) :: dActiveSolLiquid(nz,nx,1)                       ! Activated aerosol
REAL(wp) :: dAitkenSolMass(nz,nx,1), dAitkenSolNumber(nz,nx,1) ! Aitken mode aerosol
REAL(wp) :: dCoarseSolMass(nz,nx,1), dCoarseSolNumber(nz,nx,1) ! Course mode aerosol
REAL(wp) :: dActiveSolRain(nz,nx,1)                       ! Activeated aerosol in rain
REAL(wp) :: dCoarseDustMass(nz,nx,1), dCoarseDustNumber(nz,nx,1) ! Dust
REAL(wp) :: dActiveInsolIce(nz,nx,1)                       ! Activeated dust
REAL(wp) :: dActiveSolIce(nz,nx,1)                       ! Activeated aerosol in ice
REAL(wp) :: dActiveInsolLiquid(nz,nx,1)                       ! Activeated dust in cloud
REAL(wp) :: dAccumInsolMass(nz,nx,1)                      ! Accum mode dust mass
REAL(wp) :: dAccumInsolNumber(nz,nx,1)                      ! Accum mode dust number
REAL(wp) :: dActiveSolNumber(nz,nx,1)                      ! Activated soluble number (if we need a tracer)
REAL(wp) :: dActiveInsolNumber(nz,nx,1)                      ! Activated insoluble number (if we need a tracer)

INTEGER ::   ids=1,ide=nx, jds=1,jde=1, kds=1,kde=nz ,     &
       ims=1,ime=nx, jms=1,jme=1, kms=1,kme=nz , &
       its=1,ite=nx, jts=1,jte=1, kts=1,kte=nz

INTEGER :: i,k

REAL(wp) :: dt_wp

dt_wp=dt

    ! Initialize aerosol fields in case...
AitkenSolMass = 0.0
dAitkenSolMass = 0.0
AitkenSolNumber = 0.0
dAitkenSolNumber = 0.0
AccumSolMass = 0.0
dAccumSolMass = 0.0
AccumSolNumber = 0.0
dAccumSolNumber = 0.0
CoarseSolMass = 0.0
dCoarseSolMass = 0.0
CoarseSolNumber = 0.0
dCoarseSolNumber = 0.0
ActiveSolLiquid = 0.0
dActiveSolLiquid = 0.0
CoarseDustMass = 0.0
dCoarseDustMass = 0.0
CoarseDustNumber = 0.0
dCoarseDustNumber = 0.0
ActiveInsolIce = 0.0
dActiveInsolIce = 0.0
ActiveSolIce = 0.0
dActiveSolIce = 0.0
ActiveInsolLiquid = 0.0
dActiveInsolLiquid = 0.0
AccumInsolMass = 0.0
dAccumInsolMass = 0.0
AccumInsolNumber = 0.0
dAccumInsolNumber = 0.0
ActiveSolNumber = 0.0
dActiveSolNumber = 0.0
ActiveInsolNumber = 0.0
dActiveInsolNumber = 0.0


    ! Initialise microphysics
IF (micro_unset) THEN
  CALL set_mphys_switches(option, aerosol_option)
  CALL mphys_init
  micro_unset=.FALSE.
END IF

DO i=1,nx
  DO k=1,nz
    theta(k,i,1) = theta_in(k,i)
    dth(k,i,1) = (dtheta_adv(k,i)+dtheta_div(k,i))
    exner(k,i,1) = exner_in(k,i)
    pressure(k,i,1) = p0*exner_in(k,i)**(1.0/r_on_cp)
    z_centre(k,i,1) = z_centre_in(k)
    dz(k,i,1) = dz_in(k)
        !z_half(k,i,1) = z_half_in(k)
    z_half(k,i,1) = z_centre_in(k)+0.5*dz(k,i,1)
    rho(k,i,1) = rho_in(k)
    w(k,i,1) = w_in(k,i)
    tke(k,i,1) = .1 ! Test value

    qv(k,i,1) = qv_in(k,i)
    dqv(k,i,1) = (dqv_adv(k,i)+dqv_div(k,i))

    IF (nq_l > 0)qc(k,i,1) = hydrometeors(k,i,1)%moments(1,1)
    IF (nq_l > 0)dqc(k,i,1) =     &
           + (dhydrometeors_adv(k,i,1)%moments(1,1) &
           + dhydrometeors_div(k,i,1)%moments(1,1))

    IF (nq_r > 0)qr(k,i,1) = hydrometeors(k,i,2)%moments(1,1)
    IF (nq_r > 0)dqr(k,i,1) =     &
           + (dhydrometeors_adv(k,i,2)%moments(1,1) &
           + dhydrometeors_div(k,i,2)%moments(1,1))

    IF (nq_l > 1)nc(k,i,1) = hydrometeors(k,i,1)%moments(1,2)
    IF (nq_l > 1)dnc(k,i,1) =     &
           + (dhydrometeors_adv(k,i,1)%moments(1,2) &
           + dhydrometeors_div(k,i,1)%moments(1,2))

    IF (nq_r > 1)nr(k,i,1) = hydrometeors(k,i,2)%moments(1,2)
    IF (nq_r > 1)dnr(k,i,1)=     &
           + (dhydrometeors_adv(k,i,2)%moments(1,2) &
           + dhydrometeors_div(k,i,2)%moments(1,2))

    IF (nq_r > 2)m3r(k,i,1) = hydrometeors(k,i,2)%moments(1,3)
    IF (nq_r > 2)dm3r(k,i,1) =     &
           + (dhydrometeors_adv(k,i,2)%moments(1,3) &
           + dhydrometeors_div(k,i,2)%moments(1,3))

        ! Ice microphysical fields
    IF (nq_i > 0)  qi(k,i,1) = hydrometeors(k,i,3)%moments(1,1)
    IF (nq_s > 0)  qs(k,i,1) = hydrometeors(k,i,4)%moments(1,1)
    IF (nq_g > 0)  qg(k,i,1) = hydrometeors(k,i,5)%moments(1,1)
    IF (nq_i > 1)  ni(k,i,1) = hydrometeors(k,i,3)%moments(1,2)
    IF (nq_s > 1)  ns(k,i,1) = hydrometeors(k,i,4)%moments(1,2)
    IF (nq_g > 1)  ng(k,i,1) = hydrometeors(k,i,5)%moments(1,2)
    IF (nq_s > 2)  m3s(k,i,1) = hydrometeors(k,i,4)%moments(1,3)
    IF (nq_g > 2)  m3g(k,i,1) = hydrometeors(k,i,5)%moments(1,3)
    IF (nq_i > 0)  dqi(k,i,1) =     &
           + (dhydrometeors_adv(k,i,3)%moments(1,1) &
           + dhydrometeors_div(k,i,3)%moments(1,1))
    IF (nq_s > 0)  dqs(k,i,1) =     &
           + (dhydrometeors_adv(k,i,4)%moments(1,1) &
           + dhydrometeors_div(k,i,4)%moments(1,1))
    IF (nq_g > 0)  dqg(k,i,1) =     &
           + (dhydrometeors_adv(k,i,5)%moments(1,1) &
           + dhydrometeors_div(k,i,5)%moments(1,1))
    IF (nq_i > 1)  dni(k,i,1) =     &
           + (dhydrometeors_adv(k,i,3)%moments(1,2) &
           + dhydrometeors_div(k,i,3)%moments(1,2))
    IF (nq_s > 1)  dns(k,i,1) =     &
           + (dhydrometeors_adv(k,i,4)%moments(1,2) &
           + dhydrometeors_div(k,i,4)%moments(1,2))
    IF (nq_g > 1)  dng(k,i,1) =     &
           + (dhydrometeors_adv(k,i,5)%moments(1,2) &
           + dhydrometeors_div(k,i,5)%moments(1,2))
    IF (nq_s > 2)  dm3s(k,i,1) =     &
           + (dhydrometeors_adv(k,i,4)%moments(1,3) &
           + dhydrometeors_div(k,i,4)%moments(1,3))
    IF (nq_g > 2)  dm3g(k,i,1) =     &
           + (dhydrometeors_adv(k,i,5)%moments(1,3) &
           + dhydrometeors_div(k,i,5)%moments(1,3))

        ! Aerosol fields
    IF (i_am1 > 0)AitkenSolMass(k,i,1) = aerosol(k,i,1)%moments(1,2)
    IF (i_an1 > 0)AitkenSolNumber(k,i,1) = aerosol(k,i,1)%moments(1,1)
    IF (i_am2 > 0)AccumSolMass(k,i,1) = aerosol(k,i,2)%moments(1,2)
    IF (i_an2 > 0)AccumSolNumber(k,i,1) = aerosol(k,i,2)%moments(1,1)
    IF (i_am3 > 0)CoarseSolMass(k,i,1) = aerosol(k,i,3)%moments(1,2)
    IF (i_an3 > 0)CoarseSolNumber(k,i,1) = aerosol(k,i,3)%moments(1,1)
    IF (i_am4 > 0)ActiveSolLiquid(k,i,1) = aerosol(k,i,4)%moments(1,1)
    IF (i_am5 > 0)ActiveSolRain(k,i,1) = aerosol(k,i,5)%moments(1,1)
    IF (i_am6 > 0)CoarseDustMass(k,i,1) = aerosol(k,i,6)%moments(1,2)
    IF (i_an6 > 0)CoarseDustNumber(k,i,1) = aerosol(k,i,6)%moments(1,1)
    IF (i_am7 > 0)ActiveInsolIce(k,i,1) = aerosol(k,i,7)%moments(1,1)
    IF (i_am8 > 0)ActiveSolIce(k,i,1) = aerosol(k,i,8)%moments(1,1)
    IF (i_am9 > 0)ActiveInsolLiquid(k,i,1) = aerosol(k,i,9)%moments(1,1)
    IF (i_am10 > 0)AccumInsolMass(k,i,1) = aerosol(k,i,10)%moments(1,2)
    IF (i_an10 > 0)AccumInsolNumber(k,i,1) = aerosol(k,i,10)%moments(1,1)
    IF (i_an11 > 0)ActiveSolNumber(k,i,1) = aerosol(k,i,11)%moments(1,1)
    IF (i_an12 > 0)ActiveInsolNumber(k,i,1) = aerosol(k,i,12)%moments(1,1)
    IF (i_am1 > 0)dAitkenSolMass(k,i,1)=    &
           + (daerosol_adv(k,i,1)%moments(1,2) &
           + daerosol_div(k,i,1)%moments(1,2))
    IF (i_an1 > 0)dAitkenSolNumber(k,i,1)=    &
           + (daerosol_adv(k,i,1)%moments(1,1) &
           + daerosol_div(k,i,1)%moments(1,1))
    IF (i_am2 > 0)dAccumSolMass(k,i,1)=    &
           + (daerosol_adv(k,i,2)%moments(1,2) &
           + daerosol_div(k,i,2)%moments(1,2))
    IF (i_an2 > 0)dAccumSolNumber(k,i,1)=    &
           + (daerosol_adv(k,i,2)%moments(1,1) &
           + daerosol_div(k,i,2)%moments(1,1))
    IF (i_am3 > 0)dCoarseSolMass(k,i,1)=    &
           + (daerosol_adv(k,i,3)%moments(1,2) &
           + daerosol_div(k,i,3)%moments(1,2))
    IF (i_an3 > 0)dCoarseSolNumber(k,i,1)=    &
           + (daerosol_adv(k,i,3)%moments(1,1) &
           + daerosol_div(k,i,3)%moments(1,1))
    IF (i_am4 > 0)dActiveSolLiquid(k,i,1)=    &
           + (daerosol_adv(k,i,4)%moments(1,1) &
           + daerosol_div(k,i,4)%moments(1,1))
    IF (i_am5 > 0)dActiveSolRain(k,i,1)=    &
           + (daerosol_adv(k,i,4)%moments(1,1) &
           + daerosol_div(k,i,4)%moments(1,1))
    IF (i_am6 > 0)dCoarseDustMass(k,i,1)=    &
           + (daerosol_adv(k,i,6)%moments(1,2) &
           + daerosol_div(k,i,6)%moments(1,2))
    IF (i_an6 > 0)dCoarseDustNumber(k,i,1)=    &
           + (daerosol_adv(k,i,6)%moments(1,1) &
           + daerosol_div(k,i,6)%moments(1,1))
    IF (i_am7 > 0)dActiveInsolIce(k,i,1)=    &
           + (daerosol_adv(k,i,7)%moments(1,1) &
           + daerosol_div(k,i,7)%moments(1,1))
    IF (i_am8 > 0)dActiveSolIce(k,i,1)=    &
           + (daerosol_adv(k,i,8)%moments(1,1) &
           + daerosol_div(k,i,8)%moments(1,1))
    IF (i_am9 > 0)dActiveInsolLiquid(k,i,1)=    &
           + (daerosol_adv(k,i,9)%moments(1,1) &
           + daerosol_div(k,i,9)%moments(1,1))
    IF (i_am10 > 0)dAccumInsolMass(k,i,1)=    &
           + (daerosol_adv(k,i,10)%moments(1,2) &
           + daerosol_div(k,i,10)%moments(1,2))
    IF (i_an10 > 0)dAccumInsolNumber(k,i,1)=    &
           + (daerosol_adv(k,i,10)%moments(1,1) &
           + daerosol_div(k,i,10)%moments(1,1))
    IF (i_an11 > 0)dActiveSolNumber(k,i,1)=    &
           + (daerosol_adv(k,i,11)%moments(1,1) &
           + daerosol_div(k,i,11)%moments(1,1))
    IF (i_an12 > 0)dActiveInsolNumber(k,i,1)=    &
           + (daerosol_adv(k,i,12)%moments(1,1) &
           + daerosol_div(k,i,12)%moments(1,1))

  END DO
      !z_half(0,i,1) = z_half_in(0)
  z_half(0,i,1) = z_centre_in(1)-0.5*dz(1,i,1)
END DO

CALL shipway_microphysics(                     &
! in
     its, ite,                                  &
     jts, jte,                                  &
     kts, kte,                                  &
     dt_wp,                                              &
     qv, qc, qr,    &
     nc, nr, m3r,     &
     qi, qs, qg,    &
     ni, ns, ng,    &
     m3s, m3g,     &
     theta,                         &
     AitkenSolMass, AitkenSolNumber,            &
     AccumSolMass, AccumSolNumber,            &
     CoarseSolMass, CoarseSolNumber,            &
     ActiveSolLiquid,                                      &
     ActiveSolRain,                                      &
     CoarseDustMass, CoarseDustNumber,            &
     ActiveInsolIce,                                      &
     ActiveSolIce,                                      &
     ActiveInsolLiquid,                                      &
     AccumInsolMass, AccumInsolNumber,            &
     ActiveSolNumber, ActiveInsolNumber,            &
     exner,                                           &
     pressure, rho,     &
     w, tke,       &
     z_half, z_centre, dz,                        &
! in/out
     dqv, dqc, dqr, dnc, dnr, dm3r,              &
     dqi, dqs, dqg, dni, dns, dng, dm3s, dm3g,     &
     dth,                                  &
     dAitkenSolMass, dAitkenSolNumber,            &
     dAccumSolMass, dAccumSolNumber,            &
     dCoarseSolMass, dCoarseSolNumber,            &
     dActiveSolLiquid,                        &
     dActiveSolRain,                        &
     dCoarseDustMass, dCoarseDustNumber,           &
     dActiveInsolIce,                        &
     dActiveSolIce,                        &
     dActiveInsolLiquid,                         &
     dAccumInsolMass, dAccumInsolNumber,            &
     dActiveSolNumber, dActiveInsolNumber            &
     )

DO i=1,nx
  DO k=1,nz
    dtheta_mphys(k,i)=dth(k,i,1)

    dqv_mphys(k,i)=dqv(k,i,1)

        ! Retreive warm fields
    IF (nq_l > 0)dhydrometeors_mphys(k,i,1)%moments(1,1)=     &
           dqc(k,i,1)

    IF (nq_r > 0)dhydrometeors_mphys(k,i,2)%moments(1,1)=     &
           dqr(k,i,1)

    IF (nq_l > 1)dhydrometeors_mphys(k,i,1)%moments(1,2)=     &
           dnc(k,i,1)

    IF (nq_r > 1)dhydrometeors_mphys(k,i,2)%moments(1,2)=     &
           dnr(k,i,1)

    IF (nq_r > 2)dhydrometeors_mphys(k,i,2)%moments(1,3)=     &
           dm3r(k,i,1)
        ! Retreive ice fields
    IF (nq_i > 0)  dhydrometeors_mphys(k,i,3)%moments(1,1)=     &
             dqi(k,i,1)

    IF (nq_i > 1)  dhydrometeors_mphys(k,i,3)%moments(1,2)=     &
             dni(k,i,1)

    IF (nq_s > 0)  dhydrometeors_mphys(k,i,4)%moments(1,1)=     &
             dqs(k,i,1)

    IF (nq_s > 1)  dhydrometeors_mphys(k,i,4)%moments(1,2)=     &
             dns(k,i,1)

    IF (nq_s > 2)  dhydrometeors_mphys(k,i,4)%moments(1,3)=     &
             dm3s(k,i,1)

    IF (nq_g > 0)  dhydrometeors_mphys(k,i,5)%moments(1,1)=     &
             dqg(k,i,1)

    IF (nq_g > 1)  dhydrometeors_mphys(k,i,5)%moments(1,2)=     &
             dng(k,i,1)

    IF (nq_g > 2)  dhydrometeors_mphys(k,i,5)%moments(1,3)=     &
             dm3g(k,i,1)

        ! Retreive aerosol fields...
    IF (i_am1 > 0)daerosol_mphys(k,i,1)%moments(1,2)= dAitkenSolMass(k,i,1)
    IF (i_an1 > 0)daerosol_mphys(k,i,1)%moments(1,1)= dAitkenSolNumber(k,i,1)
    IF (i_am2 > 0)daerosol_mphys(k,i,2)%moments(1,2)= dAccumSolMass(k,i,1)
    IF (i_an2 > 0)daerosol_mphys(k,i,2)%moments(1,1)= dAccumSolNumber(k,i,1)
    IF (i_am3 > 0)daerosol_mphys(k,i,3)%moments(1,2)= dCoarseSolMass(k,i,1)
    IF (i_an3 > 0)daerosol_mphys(k,i,3)%moments(1,1)= dCoarseSolNumber(k,i,1)
    IF (i_am4 > 0)daerosol_mphys(k,i,4)%moments(1,1)= dActiveSolLiquid(k,i,1)
    IF (i_am5 > 0)daerosol_mphys(k,i,5)%moments(1,1)= dActiveSolRain(k,i,1)
    IF (i_am6 > 0)daerosol_mphys(k,i,6)%moments(1,2)= dCoarseDustMass(k,i,1)
    IF (i_an6 > 0)daerosol_mphys(k,i,6)%moments(1,1)= dCoarseDustNumber(k,i,1)
    IF (i_am7 > 0)daerosol_mphys(k,i,7)%moments(1,1)= dActiveInsolIce(k,i,1)
    IF (i_am8 > 0)daerosol_mphys(k,i,8)%moments(1,1)= dActiveSolIce(k,i,1)
    IF (i_am9 > 0)daerosol_mphys(k,i,9)%moments(1,1)= dActiveInsolLiquid(k,i,1)
    IF (i_am10 > 0)daerosol_mphys(k,i,10)%moments(1,2)= dAccumInsolMass(k,i,1)
    IF (i_an10 > 0)daerosol_mphys(k,i,10)%moments(1,1)= dAccumInsolNumber(k,i,1)
    IF (i_an11 > 0)daerosol_mphys(k,i,11)%moments(1,1)= dActiveSolNumber(k,i,1)
    IF (i_an12 > 0)daerosol_mphys(k,i,12)%moments(1,1)= dActiveInsolNumber(k,i,1)

  END DO
END DO

END SUBROUTINE mphys_casim_interface

#endif
END MODULE mphys_casim

