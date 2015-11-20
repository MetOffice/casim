MODULE mphys_CASIM_lem
#if DEF_MODEL==MODEL_LEM

USE variable_precision, ONLY: wp
USE micro_main, ONLY: shipway_microphysics
USE mphys_switches, ONLY: set_mphys_switches, aerosol_option, option,   &
     nq_l, nq_r, nq_i, nq_s, nq_g                                       &
     , i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am4, i_am5           &
     , i_am6, i_an6, i_am7, i_am8 , i_am9, i_am10, i_an10, i_an11, i_an12

USE initialize, ONLY: mphys_init


  ! LEM
USE com_params, ONLY: p0=>PSFR
USE com_prametr, ONLY : ilatentdgp, isurfdg, iclassp, kkp, jjp, nqp
USE com_qcnt
USE com_grid1, ONLY : zn, z, dzn, prefrcp, prefn, thref, rprefrcp, rhon
USE com_physcon, ONLY : r_on_cp=>rlvap_on_cp

USE diaghelp_lem, ONLY: koff, i_here

IMPLICIT NONE

  !Logical switches
LOGICAL :: micro_unset=.TRUE.

CONTAINS

SUBROUTINE mphys_CASIM_interface(I, RDT, Q, TH, SQ, STH, W_IN)

    ! Input

INTEGER, INTENT(IN) ::   I    ! X loop counter
REAL, INTENT(IN) ::   RDT  ! 1/timestep called RDTM higher in the model
REAL, INTENT(IN), DIMENSION(0:JJP+1,KKP,NQP) ::     Q
    ! All q variables only changed through SQ
REAL, INTENT(IN), DIMENSION(0:JJP+1,KKP) ::   TH       &
    ! Liquid water static temperature perturbation
    ! only changed through STH
         , W_IN    ! vertical velocity
REAL, INTENT(INOUT), DIMENSION(0:JJP+1,KKP,NQP) ::   SQ
    ! Source term for Q fields
REAL, INTENT(INOUT), DIMENSION(0:JJP+1,KKP) ::   STH
    ! Source term for TH variable

INTEGER :: its, ite

INTEGER, PARAMETER :: jts=0, jte=JJP+1
INTEGER, PARAMETER :: kts=1, kte=KKP
INTEGER, PARAMETER :: ils=1, ile=1
INTEGER, PARAMETER :: jls=1, jle=JJP
INTEGER, PARAMETER :: kls=2, kle=KKP-5

REAL :: theta(kts:kte,1,jts:jte), pressure(kts:kte,1,jts:jte), dz(kts:kte,1,jts:jte), qv(kts:kte,1,jts:jte),qc(kts:kte,1,jts:jte) &
       ,nc(kts:kte,1,jts:jte), qr(kts:kte,1,jts:jte), nr(kts:kte,1,jts:jte), m3r(kts:kte,1,jts:jte),rho(kts:kte,1,jts:jte) &
       , pii(kts:kte,1,jts:jte), w(kts:kte,1,jts:jte), tke(kts:kte,1,jts:jte) &
       , qi(kts:kte,1,jts:jte), qs(kts:kte,1,jts:jte), qg(kts:kte,1,jts:jte)  &
       , ni(kts:kte,1,jts:jte), ns(kts:kte,1,jts:jte), ng(kts:kte,1,jts:jte)  &
       , m3s(kts:kte,1,jts:jte), m3g(kts:kte,1,jts:jte)

REAL :: z_half(kts-1:kte,1,jts:jte), z_centre(kts:kte,1,jts:jte)

REAL :: AccumSolMass(kts:kte,1,jts:jte), AccumSolNumber(kts:kte,1,jts:jte) ! Accumulation mode aerosol
REAL :: ActiveSolLiquid(kts:kte,1,jts:jte)                                ! Activated aerosol
REAL :: AitkenSolMass(kts:kte,1,jts:jte), AitkenSolNumber(kts:kte,1,jts:jte) ! Aitken mode aerosol
REAL :: CoarseSolMass(kts:kte,1,jts:jte), CoarseSolNumber(kts:kte,1,jts:jte) ! Course mode aerosol
REAL :: ActiveSolRain(kts:kte,1,jts:jte)                                ! Activeated aerosol in rain
REAL :: CoarseDustMass(kts:kte,1,jts:jte), CoarseDustNumber(kts:kte,1,jts:jte) ! Dust
REAL :: ActiveInsolIce(kts:kte,1,jts:jte)                                ! Activeated dust
REAL :: ActiveSolIce(kts:kte,1,jts:jte)                                ! Activeated aerosol in ice
REAL :: ActiveInsolLiquid(kts:kte,1,jts:jte)                                ! Activeated dust in cloud
REAL :: AccumInsolMass(kts:kte,1,jts:jte)                                ! Accumulation mode dust mass
REAL :: AccumInsolNumber(kts:kte,1,jts:jte)                                ! Accumulation mode dust number
REAL :: ActiveSolNumber(kts:kte,1,jts:jte)                                ! Active soluble number (if we need to keep as a tracer)
REAL :: ActiveInsolNumber(kts:kte,1,jts:jte)                                ! Active insoluble number (if we need to keep as a tracer)


REAL(wp) :: dqv(kts:kte,1,jts:jte), dth(kts:kte,1,jts:jte), dqc(kts:kte,1,jts:jte), dnc(kts:kte,1,jts:jte) &
       , dqr(kts:kte,1,jts:jte), dnr(kts:kte,1,jts:jte), dm3r(kts:kte,1,jts:jte) &
       , dqi(kts:kte,1,jts:jte), dqs(kts:kte,1,jts:jte), dqg(kts:kte,1,jts:jte)&
       , dni(kts:kte,1,jts:jte), dns(kts:kte,1,jts:jte), dng(kts:kte,1,jts:jte)&
       , dm3s(kts:kte,1,jts:jte), dm3g(kts:kte,1,jts:jte)

REAL(wp) :: dAccumSolMass(kts:kte,1,jts:jte), dAccumSolNumber(kts:kte,1,jts:jte) ! Accumulation mode aerosol
REAL(wp) :: dActiveSolLiquid(kts:kte,1,jts:jte)                                 ! Activated aerosol
REAL(wp) :: dAitkenSolMass(kts:kte,1,jts:jte), dAitkenSolNumber(kts:kte,1,jts:jte) ! Aitken mode aerosol
REAL(wp) :: dCoarseSolMass(kts:kte,1,jts:jte), dCoarseSolNumber(kts:kte,1,jts:jte) ! Course mode aerosol
REAL(wp) :: dActiveSolRain(kts:kte,1,jts:jte)                                 ! Activeated aerosol in rain
REAL(wp) :: dCoarseDustMass(kts:kte,1,jts:jte), dCoarseDustNumber(kts:kte,1,jts:jte) ! Dust
REAL(wp) :: dActiveInsolIce(kts:kte,1,jts:jte)                                 ! Activeated dust
REAL(wp) :: dActiveSolIce(kts:kte,1,jts:jte)                                 ! Activeated aerosol in ice
REAL(wp) :: dActiveInsolLiquid(kts:kte,1,jts:jte)                                 ! Activeated dust in cloud
REAL :: dAccumInsolMass(kts:kte,1,jts:jte)                                ! Accumulation mode dust mass
REAL :: dAccumInsolNumber(kts:kte,1,jts:jte)                                ! Accumulation mode dust number
REAL :: dActiveSolNumber(kts:kte,1,jts:jte)                                ! Active soluble number (if we need to keep as a tracer)
REAL :: dActiveInsolNumber(kts:kte,1,jts:jte)                                ! Active insoluble number (if we need to keep as a tracer)

INTEGER :: j, k, iqx

REAL(wp) :: dt

REAL :: ratio ! multiplicative factor

    ! Initialise microphysics
IF (micro_unset) THEN
  CALL set_mphys_switches(option, aerosol_option)
  CALL mphys_init
  micro_unset=.FALSE.
END IF

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
ActiveSolRain = 0.0
dActiveSolRain = 0.0
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

qi = 0.0
qs = 0.0
qg = 0.0
ni = 0.0
ns = 0.0
ng = 0.0
m3s = 0.0
m3g = 0.0
dqi = 0.0
dqs = 0.0
dqg = 0.0
dni = 0.0
dns = 0.0
dng = 0.0
dm3s = 0.0
dm3g = 0.0

its=1
ite=1
dt=1.0/rdt

i_here=i
koff=kts-1

DO k=kts,kte
  DO j=jts, jte

    theta(k,its,j)  = (thref(k) + th(j,k))
    dth(k,its,j)   = sth(j,k)
    pii(k,its,j) = rprefrcp(k)
    pressure(k,its,j)   = prefn(k)
    dz(k,its,j)  = dzn(k)
    z_half(k,its,j)  = z(k)
    z_centre(k,its,j)  = zn(k)

    rho(k,its,j) = rhon(k)
    w(k,its,j)   = w_in(j,k)
    tke(k,its,j) = .1 ! Test value

    iqx=iqv
    qv(k,its,j)  =  MAX(0.0, q(j,k,iqx))
    dqv(k,its,j)   = sq(j,k,iqx)
    iqx=iql
    IF (nq_l > 0)qc(k,its,j)  = MAX(0.0, q(j,k,iqx))
    IF (nq_l > 0)dqc(k,its,j)   = sq(j,k,iqx)
    iqx=iqr
    IF (nq_r > 0)qr(k,its,j)  = MAX(0.0, q(j,k,iqx))
    IF (nq_r > 0)dqr(k,its,j)   = sq(j,k,iqx)

    iqx=iqnr
    IF (nq_r > 1)nr(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_r > 1)dnr(k,its,j) = sq(j,k,iqx)
    iqx=iqnl
    IF (nq_l > 1)nc(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_l > 1)dnc(k,its,j) =  sq(j,k,iqx)
    iqx=iqm3r
    IF (nq_r > 2)m3r(k,its,j) = MAX(0.0, q(j,k,iqx) )
    IF (nq_r > 2)dm3r(k,its,j) = sq(j,k,iqx)


          ! Set ice microphysical fields
    iqx=iqi
    IF (nq_i > 0)qi(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_i > 0)dqi(k,its,j) = sq(j,k,iqx)
    iqx=iqs
    IF (nq_s > 0)qs(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_s > 0)dqs(k,its,j) = sq(j,k,iqx)
    iqx=iqg
    IF (nq_g > 0)qg(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_g > 0)dqg(k,its,j) = sq(j,k,iqx)
    iqx=iqn
    IF (nq_i > 1)ni(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_i > 1)dni(k,its,j) = sq(j,k,iqx)
    iqx=iqns
    IF (nq_s > 1)ns(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_s > 1)dns(k,its,j) = sq(j,k,iqx)
    iqx=iqng
    IF (nq_g > 1)ng(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_g > 1)dng(k,its,j) = sq(j,k,iqx)
          ! Not coupled in 3rd moment yet
    iqx=iqm3s
    IF (nq_s > 2)m3s(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_s > 2)dm3s(k,its,j) = sq(j,k,iqx)
    iqx=iqm3g
    IF (nq_g > 2)m3g(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (nq_g > 2)dm3g(k,its,j) = sq(j,k,iqx)

    iqx=iam1
    IF (i_am1 >0)AitkenSolMass(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am1 >0)dAitkenSolMass(k,its,j)= sq(j,k,iqx)
    iqx=ian1
    IF (i_an1 >0)AitkenSolNumber(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_an1 >0)dAitkenSolNumber(k,its,j)= sq(j,k,iqx)
    iqx=iam2
    IF (i_am2 >0)AccumSolMass(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am2 >0)dAccumSolMass(k,its,j)= sq(j,k,iqx)
    iqx=ian2
    IF (i_an2 >0)AccumSolNumber(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_an2 >0)dAccumSolNumber(k,its,j)= sq(j,k,iqx)
    iqx=iam3
    IF (i_am3 >0)CoarseSolMass(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am3 >0)dCoarseSolMass(k,its,j)= sq(j,k,iqx)
    iqx=ian3
    IF (i_an3 >0)CoarseSolNumber(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_an3 >0)dCoarseSolNumber(k,its,j)= sq(j,k,iqx)
    iqx=iam4
    IF (i_am4 >0)ActiveSolLiquid(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am4 >0)dActiveSolLiquid(k,its,j)= sq(j,k,iqx)
    iqx=iam5
    IF (i_am5 >0)ActiveSolRain(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am5 >0)dActiveSolRain(k,its,j)= sq(j,k,iqx)
    iqx=iam6
    IF (i_am6 >0)CoarseDustMass(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am6 >0)dCoarseDustMass(k,its,j)= sq(j,k,iqx)
    iqx=ian6
    IF (i_an6 >0)CoarseDustNumber(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_an6 >0)dCoarseDustNumber(k,its,j)= sq(j,k,iqx)
    iqx=iam7
    IF (i_am7 >0)ActiveInsolIce(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am7 >0)dActiveInsolIce(k,its,j)= sq(j,k,iqx)
    iqx=iam8
    IF (i_am8 >0)ActiveSolIce(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am8 >0)dActiveSolIce(k,its,j)= sq(j,k,iqx)
    iqx=iam9
    IF (i_am9 >0)ActiveInsolLiquid(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am9 >0)dActiveInsolLiquid(k,its,j)= sq(j,k,iqx)
    iqx=iam10
    IF (i_am10 >0)AccumInsolMass(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_am10 >0)dAccumInsolMass(k,its,j)= sq(j,k,iqx)
    iqx=ian10
    IF (i_an10 >0)AccumInsolNumber(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_an10 >0)dAccumInsolNumber(k,its,j)= sq(j,k,iqx)
    iqx=ian11
    IF (i_an11 >0)ActiveSolNumber(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_an11 >0)dActiveSolNumber(k,its,j)= sq(j,k,iqx)
    iqx=ian12
    IF (i_an12 >0)ActiveInsolNumber(k,its,j) = MAX(0.0, q(j,k,iqx))
    IF (i_an12 >0)dActiveInsolNumber(k,its,j)= sq(j,k,iqx)

  END DO
END DO


CALL shipway_microphysics(                             &
! in
     its, ite,                                  &
     jts, jte,                                  &
     kts, kte,                                  &
     dt,                                              &
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
     AccumInsolMass,                                 &
     AccumInsolNumber,                               &
     ActiveSolNumber,                                &
     ActiveInsolNumber,                              &
     pii,                                           &
     pressure, rho,     &
     w, tke,       &
     z_half, z_centre,                             &
     dz,                        &
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
     dActiveInsolLiquid,                        &
     dAccumInsolMass,                                 &
     dAccumInsolNumber,                               &
     dActiveSolNumber,                                &
     dActiveInsolNumber,                              &
     ils, ile,                                  &
     jls, jle,                                  &
     kls, kle,                                  &
     l_tendency=.TRUE.                          &
     )

DO k=kls,kle
       ! save tendencies
  sth(jls:jle,k) = sth(jls:jle,k) + dth(k,1,jls:jle)
  iqx=iqv
  sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dqv(k,1,jls:jle)
  iqx=iql
  IF (nq_l > 0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dqc(k,1,jls:jle)
  iqx=iqr
  IF (nq_r > 0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dqr(k,1,jls:jle)

  iqx=iqnr
  IF (nq_r > 1)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dnr(k,1,jls:jle)
  iqx=iqnl
  IF (nq_l > 1)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dnc(k,1,jls:jle)
  iqx=iqm3r
  IF (nq_r > 2)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dm3r(k,1,jls:jle)

  iqx=iqi
  IF (nq_i > 0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dqi(k,1,jls:jle)
  iqx=iqs
  IF (nq_s > 0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dqs(k,1,jls:jle)
  iqx=iqg
  IF (nq_g > 0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dqg(k,1,jls:jle)
  iqx=iqn
  IF (nq_i > 1)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dni(k,1,jls:jle)
  iqx=iqns
  IF (nq_s > 1)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dns(k,1,jls:jle)
  iqx=iqng
  IF (nq_g > 1)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dng(k,1,jls:jle)
  iqx=iqm3s
  IF (nq_s > 2)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dm3s(k,1,jls:jle)
  iqx=iqm3g
  IF (nq_g > 2)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dm3g(k,1,jls:jle)

  iqx=iam1
  IF (i_am1>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dAitkenSolMass(k,1,jls:jle)
  iqx=ian1
  IF (i_an1>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dAitkenSolNumber(k,1,jls:jle)
  iqx=iam2
  IF (i_am2>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dAccumSolMass(k,1,jls:jle)
  iqx=ian2
  IF (i_an2>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dAccumSolNumber(k,1,jls:jle)
  iqx=iam3
  IF (i_am3>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dCoarseSolMass(k,1,jls:jle)
  iqx=ian3
  IF (i_an3>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dCoarseSolNumber(k,1,jls:jle)
  iqx=iam4
  IF (i_am4>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dActiveSolLiquid(k,1,jls:jle)
  iqx=iam5
  IF (i_am5>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dActiveSolLiquid(k,1,jls:jle)
  iqx=iam6
  IF (i_am6>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dCoarseDustMass(k,1,jls:jle)
  iqx=ian6
  IF (i_an6>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dCoarseDustNumber(k,1,jls:jle)
  iqx=iam7
  IF (i_am7>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dActiveInsolIce(k,1,jls:jle)
  iqx=iam8
  IF (i_am8>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dActiveSolIce(k,1,jls:jle)
  iqx=iam9
  IF (i_am9>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dActiveInsolLiquid(k,1,jls:jle)
  iqx=iam10
  IF (i_am10>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dAccumInsolMass(k,1,jls:jle)
  iqx=ian10
  IF (i_an10>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dAccumInsolNumber(k,1,jls:jle)
  iqx=ian11
  IF (i_an11>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dActiveSolNumber(k,1,jls:jle)
  iqx=ian12
  IF (i_an12>0)sq(jls:jle,k,iqx) = sq(jls:jle,k,iqx) + dActiveInsolNumber(k,1,jls:jle)

END DO

        ! Prevent small negative numbers being generated through
        ! rounding error
ratio=(1.0-SPACING(dt))*rdt
DO k=kts,kte
  DO j=jts,jte
    DO iqx=1,nqp
      IF (q(j,k,iqx)+dt*sq(j,k,iqx) < 0.0 ) THEN
        sq(j,k,iqx) = - ratio*q(j,k,iqx)
      END IF
    END DO
  END DO
END DO

        ! bottom boundary
DO iqx=1,nqp
  DO j=0,jjp
    sq(j,1,iqx)=sq(j, 2, iqx)
  END DO
END DO

        ! Do we need to wrap things??
DO iqx=1,nqp
  DO k=1,kkp
    sq(0,k,iqx)=sq(jjp, k, iqx)
    sq(jjp+1,k,iqx)=sq(1, k, iqx)
  END DO
END DO

CALL ClassSet(i, q)

END SUBROUTINE mphys_CASIM_interface
#endif

END MODULE mphys_CASIM_lem


