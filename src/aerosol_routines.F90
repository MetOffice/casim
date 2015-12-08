MODULE aerosol_routines

USE mphys_die, ONLY: throw_mphys_error
USE type_aerosol, ONLY: aerosol_type, aerosol_phys, aerosol_chem, aerosol_active
USE variable_precision, ONLY: wp
USE mphys_constants, ONLY: Mw, zetasa, Ru, rhow   &
     , g, Lv, eps, cp, Rd, Rv, ka, Dv, pi
USE special, ONLY: erfc, erfinv
USE mphys_switches, ONLY: i_am2, i_an2, i_am4, i_nl, i_nr     &
     , i_am1, i_an1, i_am3, i_an3, i_am5, i_am6, i_an6, i_am7 &
     , i_am8, i_am9, i_an11, i_an12                           &
     , i_ni, i_ns, i_ng, i_ql, i_qr               &
     , i_qi, i_qs, i_qg                                       &
     , l_active_inarg2000, aero_index, l_warm                 &
     , active_cloud, active_rain, active_ice, isol, iinsol    &
     , l_process, l_passivenumbers, l_passivenumbers_ice, l_separate_rain
USE thresholds, ONLY: nr_tidy, nl_tidy, ni_tidy, ccn_tidy, qr_tidy     &
     , aeromass_small, aeronumber_small
USE mphys_parameters, ONLY: sigma_arc, nz

USE preconditioning, ONLY: precondition

  ! Really want to be allocating this elsewhere

#if DEF_MODEL==MODEL_KiD
USE diagnostics, ONLY: save_dg, i_dgtime, k_here, i_here, j_here, nx
USE runtime, ONLY: time
#elif  DEF_MODEL==MODEL_LEM
USE diaghelp_lem, ONLY: k_here, i_here, j_here
USE com_params, ONLY: time
#elif  DEF_MODEL==MODEL_UM
USE timestep_mod, ONLY: timestep_number
USE diaghelp_um, ONLY: k_here, i_here, j_here
USE UM_ParCore, ONLY: mype
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here
#endif

IMPLICIT NONE

CONTAINS

SUBROUTINE allocate_aerosol(aerophys, aerochem, nmodes, initphys, initchem)
    !
    ! Allocate space for aerosol_chem and aerosol_phys
    !
TYPE(aerosol_phys), INTENT(INOUT) :: aerophys(:)
TYPE(aerosol_chem), INTENT(INOUT) :: aerochem(:)
INTEGER, INTENT(IN) :: nmodes    ! number of modes
TYPE(aerosol_phys), INTENT(IN), OPTIONAL :: initphys
TYPE(aerosol_chem), INTENT(IN), OPTIONAL :: initchem

INTEGER :: k

    !------------------------
    ! First allocate aerophys
    !------------------------
DO k=1,SIZE(aerophys)
  aerophys(k)%nmodes=nmodes
  ALLOCATE(aerophys(k)%N(MAX(1,nmodes)))
  ALLOCATE(aerophys(k)%M(MAX(1,nmodes)))
  ALLOCATE(aerophys(k)%rd(MAX(1,nmodes)))
  ALLOCATE(aerophys(k)%sigma(MAX(1,nmodes)))
  ALLOCATE(aerophys(k)%rpart(MAX(1,nmodes)))

  IF (PRESENT(initphys)) aerophys(k)=initphys
END DO

    !------------------------
    ! Then allocate aerochem
    !------------------------
DO k=1,SIZE(aerochem)
  aerochem(k)%nmodes=nmodes
  ALLOCATE(aerochem(k)%vantHoff(MAX(1,nmodes)))
  ALLOCATE(aerochem(k)%massMole(MAX(1,nmodes)))
  ALLOCATE(aerochem(k)%density(MAX(1,nmodes)))
  ALLOCATE(aerochem(k)%epsv(MAX(1,nmodes)))
  ALLOCATE(aerochem(k)%beta(MAX(1,nmodes)))

  IF (PRESENT(initphys)) aerochem(k)=initchem
END DO

END SUBROUTINE allocate_aerosol

SUBROUTINE deallocate_aerosol(aerophys, aerochem)
    !
    ! Deallocate space for aerosol_chem and aerosol_phys
    !
TYPE(aerosol_phys), INTENT(INOUT) :: aerophys(:)
TYPE(aerosol_chem), INTENT(INOUT) :: aerochem(:)

INTEGER :: k

    !------------------------
    ! First deallocate aerophys
    !------------------------
DO k=1,SIZE(aerophys)
  DEALLOCATE(aerophys(k)%N)
  DEALLOCATE(aerophys(k)%M)
  DEALLOCATE(aerophys(k)%rd)
  DEALLOCATE(aerophys(k)%sigma)
  DEALLOCATE(aerophys(k)%rpart)
END DO

    !------------------------
    ! Then deallocate aerochem
    !------------------------
DO k=1,SIZE(aerochem)
  DEALLOCATE(aerochem(k)%vantHoff)
  DEALLOCATE(aerochem(k)%massMole)
  DEALLOCATE(aerochem(k)%density)
  DEALLOCATE(aerochem(k)%epsv)
  DEALLOCATE(aerochem(k)%beta)
END DO

END SUBROUTINE deallocate_aerosol

SUBROUTINE examine_aerosol(aerofields, qfields,    &
     aerophys, aerochem, aeroact,                  &
     dustphys, dustchem, dustact,                  &
     aeroice, dustliq, icall)
    !
    !< This routine examines the aerosol input and derives physical and
    !< chemical properties in
    !<     - aerophys: information about physical composition of aerosol
    !<     - aerochem: information about chemical composition of aerosol (currently not derived)
    !<     - aeroact:  information about physical composition activated aerosol
    !<     - dustphys: information about physical composition of dust
    !<     - dustchem: information about chemical composition of dust (currently not derived)
    !<     - dustact:  information about physical composition activated dust
    !<     - aeroice:  information about physical composition activated aerosol in ice
    !<     - dustliq:  information about physical composition activated dust in liquid water
    !
REAL(wp), INTENT(IN), TARGET :: aerofields(:,:)
TYPE(aerosol_phys), INTENT(INOUT), TARGET :: aerophys(:)
TYPE(aerosol_chem), INTENT(IN), TARGET :: aerochem(:)
TYPE(aerosol_active), INTENT(INOUT) :: aeroact(:)
TYPE(aerosol_phys), INTENT(INOUT), TARGET, OPTIONAL :: dustphys(:)
TYPE(aerosol_chem), INTENT(IN), TARGET, OPTIONAL :: dustchem(:)
TYPE(aerosol_active), INTENT(INOUT), OPTIONAL :: dustact(:)
TYPE(aerosol_active), INTENT(INOUT), OPTIONAL :: aeroice(:)
TYPE(aerosol_active), INTENT(INOUT), OPTIONAL :: dustliq(:)
REAL(wp), INTENT(INOUT), TARGET :: qfields(:,:)
INTEGER, INTENT(IN), OPTIONAL :: icall

REAL(wp) :: n, m, mac, mar, mad, maai, madl, cloud_number, rain_number
REAL(wp) :: cloud_mass, rain_mass
REAL(wp) :: ice_number, snow_number, graupel_number
REAL(wp) :: ice_mass, snow_mass, graupel_mass
REAL(wp) :: mact2, rcrit
REAL(wp) :: Nd, rm, sigma, density
INTEGER :: k

REAL(wp) :: mtot, tmp, m_ratio, n_ratio, dcloud_number, rm_arc
CHARACTER(2) :: chcall, chmode

REAL(wp) :: ratio_l, ratio_r, nratio_l, nratio_r, mratio_l, mratio_r
REAL(wp) :: ratio_i, ratio_s, ratio_g, nratio_i, nratio_s, nratio_g, mratio_i, mratio_s, mratio_g
REAL(wp) :: ratio_dil, ratio_dli, ratio_ali, ratio_ail

INTEGER :: imode
REAL(wp) :: mode_N, mode_M
REAL(wp) :: nitot, ntot, nhtot, mitot

REAL(wp) :: mact, rcrit2

REAL(wp) :: foo

LOGICAL :: l_condition, l_condition_r

    ! Some diagnostic strings
IF (PRESENT(icall)) THEN
  WRITE(chcall,'(A1,i1)') '_', icall
ELSE
  chcall=''
END IF

DO k=1,nz

  IF (precondition(k)) THEN


    cloud_number=qfields(k,i_nl)
    rain_number=qfields(k,i_nr)
    cloud_mass=qfields(k,i_ql)
    rain_mass=qfields(k,i_qr)


    density=aerochem(k)%density(1)

    DO imode=1,aero_index%nccn
      mode_N=aerofields(k,aero_index%ccn_n(imode))
      mode_M=aerofields(k,aero_index%ccn_m(imode))
      IF (mode_N > ccn_tidy .AND. mode_M > ccn_tidy*1e-18) THEN ! FiX This
        aerophys(k)%N(imode) = mode_N
        aerophys(k)%M(imode) = mode_M
        aerophys(k)%rd(imode)= MNtoRm(mode_M,mode_N,density,aerophys(k)%sigma(imode))
      ELSE
        aerophys(k)%N(imode) = 0.0
        aerophys(k)%M(imode) = 0.0
        aerophys(k)%rd(imode)= 0.0
      END IF
    END DO

    aeroact(k)%nact            = 0.0
    aeroact(k)%mact            = 0.0
    aeroact(k)%rcrit           = 999.0
    aeroact(k)%mact_mean       = 0.0
    aeroact(k)%nact2           = 0.0
    aeroact(k)%nact1           = 0.0
    aeroact(k)%mact1           = 0.0
    aeroact(k)%rcrit1          = 999.0
    aeroact(k)%nact2           = 0.0
    aeroact(k)%rcrit2          = 999.0
    aeroact(k)%mact2           = 0.0
    aeroact(k)%mact2_mean      = 0.0
    aeroact(k)%mact1_mean      = 0.0

    aeroact(k)%nratio1         = 0.0
    aeroact(k)%nratio2         = 0.0

    IF (l_process) THEN
          ! Examine activated aerosol
      mac=aerofields(k,i_am4)

      IF (l_separate_rain) THEN
        mar=aerofields(k,i_am5)
        mact=mac+mar
        l_condition = mac > 0.0 .AND. cloud_number > nl_tidy
        l_condition_r = mar > 0.0 .AND. rain_number > nr_tidy .AND. rain_mass > qr_tidy
      ELSE
        mact=mac
        l_condition = mac > 0.0 .AND. cloud_number + rain_number > nl_tidy
        l_condition_r = mac > 0.0 .AND. rain_number > nr_tidy .AND. rain_mass > qr_tidy
      END IF

      if (l_warm)then
         maai = 0.0
      else
         maai=aerofields(k,i_am8)
      end if
      
      ratio_ali=mact/(mact+maai+1e-38)

      nhtot=cloud_number+rain_number

      IF (l_passivenumbers) THEN
        ntot = aerofields(k,i_an11)*ratio_ali
      ELSE
        ntot = nhtot
      END IF

      l_condition=l_condition .and. ntot > 0.0

      nhtot=nhtot+1e-38 ! prevent possible later division by zero
      nratio_l=cloud_number/nhtot
      nratio_r=rain_number/nhtot

      IF (l_condition) THEN

        aeroact(k)%nact               = ntot
        aeroact(k)%mact               = mac
        aeroact(k)%rcrit              = 0.0
        aeroact(k)%mact_mean          = aeroact(k)%mact /(aeroact(k)%nact + TINY(ntot))

            ! Get mean radius of distribution
        rm_arc = MNtoRm(aeroact(k)%mact,aeroact(k)%nact,density,sigma_arc)

        aeroact(k)%rd    = rm_arc
        aeroact(k)%sigma = sigma_arc

      END IF

      IF (l_condition_r) THEN
        IF (l_separate_rain) THEN! Separate rain category
          rcrit2 = 0
          mact2  = mar
        ELSE IF (.NOT. l_condition) THEN ! No cloud here
          rcrit2 = 0
          mact2  = mac
        ELSE ! Diagnostic partitioning of aerosol between cloud and rain.
          rcrit2= 0.0
          mact2 = aeroact(k)%mact*rain_mass/(cloud_mass + rain_mass)
        END IF

        aeroact(k)%nact2            = ntot*nratio_r
        aeroact(k)%rcrit2           = rcrit2
        aeroact(k)%mact2            = mact2
        aeroact(k)%mact2_mean       = aeroact(k)%mact2/(aeroact(k)%nact2 + TINY(ntot))
      END IF


      aeroact(k)%nact1            = MAX(0.0_wp, aeroact(k)%nact - aeroact(k)%nact2)
      aeroact(k)%mact1            = MAX(0.0_wp, aeroact(k)%mact - aeroact(k)%mact2)
      aeroact(k)%rcrit1           = 0.0
      aeroact(k)%mact1_mean       = aeroact(k)%mact1 / (aeroact(k)%nact1 + TINY(mar))

      IF (cloud_number > 0.0)       aeroact(k)%nratio1 = aeroact(k)%nact1/cloud_number
      IF (rain_number  > 0.0)       aeroact(k)%nratio2 = aeroact(k)%nact2/rain_number

    END IF

!        if (time > 114 .and. k==20 .and. i_here==64)then
!           print*, 'DEBUG CLA1 ', cloud_mass, cloud_number, rain_mass, rain_number, aeroact(k)%mact1, aeroact(k)%mact2
!        end if

  END IF
END DO

IF (.NOT. l_warm) THEN
      ! Activated dust
  DO k=1,nz

    density=dustchem(k)%density(1)

    ice_number=qfields(k,i_ni)
    snow_number=qfields(k,i_ns)
    nhtot=ice_number+snow_number
    IF (i_ng > 0) THEN
      graupel_number=qfields(k,i_ng)
      nhtot=nhtot+graupel_number
    END IF

    IF (l_process) THEN
      mad=aerofields(k,i_am7)
      madl=aerofields(k,i_am9)
      ratio_dil = mad/(madl+mad+TINY(mad))
    END IF

    IF (l_passivenumbers_ice) THEN
      nitot=aerofields(k,i_an12)*ratio_dil
    ELSE
      nitot=nhtot
    END IF

    nhtot=nhtot+TINY(nhtot) ! prevent possible later division by zero
    nratio_i=ice_number/nhtot
    nratio_s=snow_number/nhtot
    nratio_g=graupel_number/nhtot

!        mratio_i=ice_mass/mitot
!        mratio_s=snow_mass/mitot
!        mratio_g=graupel_mass/mitot

        ! Use nratios...
    ratio_i=nratio_i
    ratio_s=nratio_s
    ratio_g=nratio_g

        ! Examine interstitial dust
    DO imode=1,aero_index%nin
      mode_N=aerofields(k,aero_index%in_n(imode))
      mode_M=aerofields(k,aero_index%in_m(imode))
      IF (mode_m < 0 .OR. mode_n < 0 ) THEN
        PRINT*, aero_index%in_m(imode), aerofields(k,:)
        PRINT*, 'Problem! Using pragmatic hack to continue.  If you see this message, report it to Ben!!!'
        PRINT*, k, chcall, imode, mode_m, mode_n, nitot, ice_number, snow_number, graupel_number
        IF (mode_m < 0)mode_m=aeromass_small
        IF (mode_n < 0)mode_n=aeronumber_small
      END IF

      dustphys(k)%N(imode) = mode_N
      dustphys(k)%M(imode) = mode_M

      dustphys(k)%rd(imode)= MNtoRm(mode_M,mode_N,density,dustphys(k)%sigma(imode))
    END DO


        ! Examine activated dust
        ! Initialize to zero/defaults
    dustact(k)%nact             = 0.0
    dustact(k)%mact             = 0.0
    dustact(k)%rcrit            = 999.0
    dustact(k)%mact_mean        = 0.0
    dustact(k)%nact1            = 0.0
    dustact(k)%nratio1          = 0.0
    dustact(k)%mact1            = 0.0
    dustact(k)%rcrit1           = 999.0
    dustact(k)%mact1_mean       = 0.0
    dustact(k)%mact2            = 0.0
    dustact(k)%nact2            = 0.0
    dustact(k)%nratio2          = 0.0
    dustact(k)%rcrit2           = 999.0
    dustact(k)%mact2_mean       = 0.0
    dustact(k)%mact3            = 0.0
    dustact(k)%nact3            = 0.0
    dustact(k)%nratio3          = 0.0
    dustact(k)%rcrit3           = 999.0
    dustact(k)%mact3_mean       = 0.0


    IF (l_process) THEN

      IF (mad > 0.0 .AND. nitot > ni_tidy) THEN

            ! Equal partitioning of dust distribution across all ice species
            ! This will cause problems, e.g. when ice splinters, there will be a
            ! spurious increase in the distributed dust mass in the ice category
            ! (This isn't a problem if there is no snow or graupel)

        dustact(k)%nact            = nitot
        dustact(k)%mact            = mad
        dustact(k)%rcrit           = 0.0
        dustact(k)%mact_mean       = dustact(k)%mact /(dustact(k)%nact+TINY(mad))

        IF (ice_number > 0.0) THEN
          dustact(k)%nact1           = nitot * ratio_i
          dustact(k)%mact1           = mad * dustact(k)%nact1/dustact(k)%nact
          dustact(k)%rcrit1          = 0.0
          dustact(k)%mact1_mean      = dustact(k)%mact1/(dustact(k)%nact1+TINY(nhtot))
          dustact(k)%nratio1         = dustact(k)%nact1/(ice_number+TINY(nhtot))
        END IF

        IF (snow_number > 0.0) THEN
          dustact(k)%nact2           = nitot * ratio_s
          dustact(k)%rcrit2          = 0.0
          dustact(k)%mact2           = mad * dustact(k)%nact2/(dustact(k)%nact+TINY(nhtot))
          dustact(k)%mact2_mean      = dustact(k)%mact2/(dustact(k)%nact2+TINY(nhtot))
          dustact(k)%nratio2         = dustact(k)%nact2/(snow_number+TINY(nhtot))
        END IF

        IF (i_ng > 0 .AND. graupel_number > ni_tidy) THEN
          dustact(k)%nact3           = nitot * ratio_g
          dustact(k)%rcrit3          = 0.0
          dustact(k)%mact3           = mad * dustact(k)%nact3/(dustact(k)%nact+TINY(nhtot))
          dustact(k)%mact3_mean      = dustact(k)%mact3/(dustact(k)%nact3+TINY(nhtot))
          dustact(k)%nratio3         = dustact(k)%nact3/(graupel_number+TINY(nhtot))
        END IF

      END IF

    END IF

  END DO

  IF (PRESENT(aeroice)) THEN
    IF (l_process) THEN
          ! Activated aerosol in ice
      DO k=1,nz

        density=aerochem(k)%density(1)

        maai=aerofields(k,i_am8)
        mac=aerofields(k,i_am4)
        ratio_ail = maai/(maai+mac+1e-38)


        ice_number=qfields(k,i_ni)
        snow_number=qfields(k,i_ns)
        ice_mass=qfields(k,i_qi)
        snow_mass=qfields(k,i_qs)
        graupel_mass=qfields(k,i_qg)

        nhtot=ice_number+snow_number
        mitot=ice_mass+snow_mass
        IF (i_ng > 0) THEN
          graupel_number=qfields(k,i_ng)
          nhtot=nhtot+graupel_number
          mitot=mitot+graupel_mass
        END IF

        IF (l_passivenumbers) THEN
          nitot=aerofields(k,i_an11)*ratio_ail
        ELSE
          nitot=nhtot
        END IF


        nhtot=nhtot+1e-38 ! prevent possible later division by zero
        nratio_i=ice_number/nhtot
        nratio_s=snow_number/nhtot
        nratio_g=graupel_number/nhtot
!            mratio_i=ice_mass/mitot
!            mratio_s=snow_mass/mitot
!            mratio_g=graupel_mass/mitot

            ! Use nratios...
        ratio_i=nratio_i
        ratio_s=nratio_s
        ratio_g=nratio_g

            ! Initialize to zero/defaults
        aeroice(k)%nact             = 0.0
        aeroice(k)%mact             = 0.0
        aeroice(k)%rcrit            = 999.0
        aeroice(k)%mact_mean        = 0.0
        aeroice(k)%nact1            = 0.0
        aeroice(k)%nratio1          = 0.0
        aeroice(k)%mact1            = 0.0
        aeroice(k)%rcrit1           = 999.0
        aeroice(k)%mact1_mean       = 0.0
        aeroice(k)%mact2            = 0.0
        aeroice(k)%nact2            = 0.0
        aeroice(k)%nratio2          = 0.0
        aeroice(k)%rcrit2           = 999.0
        aeroice(k)%mact2_mean       = 0.0
        aeroice(k)%mact3            = 0.0
        aeroice(k)%nact3            = 0.0
        aeroice(k)%nratio3          = 0.0
        aeroice(k)%rcrit3           = 999.0
        aeroice(k)%mact3_mean       = 0.0

        IF (maai > 0.0 .AND. nitot > ni_tidy) THEN

              ! Equal partitioning of dust distribution across all ice species

          aeroice(k)%nact            = nitot
          aeroice(k)%mact            = maai
          aeroice(k)%rcrit           = 0.0
          aeroice(k)%mact_mean       = aeroice(k)%mact /(aeroice(k)%nact+TINY(nhtot))

          IF (ratio_i > 0.0) THEN
            aeroice(k)%nact1           = nitot * ratio_i
            aeroice(k)%mact1           = maai * ratio_i
            aeroice(k)%rcrit1          = 0.0
            aeroice(k)%mact1_mean      = aeroice(k)%mact1/(aeroice(k)%nact1+TINY(nhtot))
            aeroice(k)%nratio1         = aeroice(k)%nact1/(ice_number+TINY(nhtot))
          END IF

          IF (ratio_s > 0.0) THEN
            aeroice(k)%nact2           = nitot * ratio_s
            aeroice(k)%rcrit2          = 0.0
            aeroice(k)%mact2           = maai * ratio_s
            aeroice(k)%mact2_mean      = aeroice(k)%mact2/(aeroice(k)%nact2+TINY(nhtot))
            aeroice(k)%nratio2         = aeroice(k)%nact2/(snow_number+TINY(nhtot))
          END IF

          IF (ratio_g > 0.0) THEN
            aeroice(k)%nact3           = nitot * ratio_g
            aeroice(k)%rcrit3          = 0.0
            aeroice(k)%mact3           = maai * ratio_g
            aeroice(k)%mact3_mean      = aeroice(k)%mact3/(aeroice(k)%nact3+TINY(nhtot))
            aeroice(k)%nratio3         = aeroice(k)%nact3/(graupel_number+TINY(nhtot))
          END IF

        END IF

      END DO

    END IF
  END IF

  IF (PRESENT(dustliq)) THEN
    IF (l_process) THEN
          ! Activated aerosol in ice
      DO k=1,nz

        density=dustchem(k)%density(1)

        madl=aerofields(k,i_am9)
        mad=aerofields(k,i_am7)
        ratio_dli = madl/(mad+madl+1e-38)

        cloud_number=qfields(k,i_nl)
        rain_number=qfields(k,i_nr)

        nratio_l=cloud_number/(cloud_number+rain_number+TINY(cloud_number))
        nratio_r=rain_number/(cloud_number+rain_number+TINY(cloud_number))
        mratio_l=cloud_mass/(cloud_mass+rain_mass+TINY(cloud_number))
        mratio_r=rain_mass/(cloud_mass+rain_mass+TINY(cloud_number))

            ! Use nratios...
        ratio_l=nratio_l
        ratio_r=nratio_r

        IF (l_passivenumbers_ice) THEN
          ntot=aerofields(k,i_an12)*ratio_dli
        ELSE
          ntot=cloud_number + rain_number
        END IF

            ! Initialize to zero/defaults
        dustliq(k)%nact             = 0.0
        dustliq(k)%mact             = 0.0
        dustliq(k)%rcrit            = 999.0
        dustliq(k)%mact_mean        = 0.0
        dustliq(k)%nact1            = 0.0
        dustliq(k)%nratio1          = 0.0
        dustliq(k)%mact1            = 0.0
        dustliq(k)%rcrit1           = 999.0
        dustliq(k)%mact1_mean       = 0.0
        dustliq(k)%mact2            = 0.0
        dustliq(k)%nact2            = 0.0
        dustliq(k)%nratio2          = 0.0
        dustliq(k)%rcrit2           = 999.0
        dustliq(k)%mact2_mean       = 0.0
        dustliq(k)%mact3            = 0.0
        dustliq(k)%nact3            = 0.0
        dustliq(k)%nratio3          = 0.0
        dustliq(k)%rcrit3           = 999.0
        dustliq(k)%mact3_mean       = 0.0

        IF (madl > 0.0 .AND. ntot > nr_tidy) THEN

              ! Equal partitioning of aerosol distribution across all liquid species

          dustliq(k)%nact            = ntot
          dustliq(k)%mact            = madl
          dustliq(k)%rcrit           = 0.0
          dustliq(k)%mact_mean       = dustliq(k)%mact/(dustliq(k)%nact+TINY(nhtot))

          IF (ratio_l > 0.0) THEN
            dustliq(k)%nact1           = ntot * ratio_l
            dustliq(k)%mact1           = madl * ratio_l
            dustliq(k)%rcrit1          = 0.0
            dustliq(k)%mact1_mean      = dustliq(k)%mact1/(dustliq(k)%nact1+TINY(nhtot))
            dustliq(k)%nratio1         = dustliq(k)%nact1/(cloud_number+TINY(nhtot))
          END IF

          IF (ratio_r > 0.0) THEN
            dustliq(k)%nact2           = ntot * ratio_r
            dustliq(k)%rcrit2          = 0.0
            dustliq(k)%mact2           = madl * ratio_r
            dustliq(k)%mact2_mean      = dustliq(k)%mact2/(dustliq(k)%nact2+TINY(nhtot))
            dustliq(k)%nratio2         = dustliq(k)%nact2/(rain_number+TINY(nhtot))
          END IF
        END IF

      END DO

    END IF
  END IF

END IF

#if DEF_MODEL==MODEL_KiD
DO k=1,nz
  IF (nx==1) THEN

    CALL save_dg(k, dustact(k)%mact3_mean, 'dustact%mact3_mean'//chcall, i_dgtime )
    CALL save_dg(k, dustact(k)%mact_mean, 'dustact%mact_mean'//chcall, i_dgtime )
    CALL save_dg(k, dustact(k)%mact, 'dustact%mact'//chcall, i_dgtime)
    CALL save_dg(k, dustliq(k)%mact3_mean, 'dustliq%mact3_mean'//chcall, i_dgtime )
    CALL save_dg(k, dustliq(k)%mact_mean, 'dustliq%mact_mean'//chcall, i_dgtime )
    CALL save_dg(k, dustliq(k)%mact, 'dustliq%mact'//chcall, i_dgtime)
    DO imode=1,aerophys(k)%nmodes
      WRITE(chmode, '(i2)') imode
      CALL save_dg(k, aerophys(k)%N(imode), 'aerophys%N'//chmode//chcall, i_dgtime )
      CALL save_dg(k, aerophys(k)%M(imode), 'aerophys%M'//chmode//chcall, i_dgtime )
      CALL save_dg(k, aerophys(k)%rd(imode), 'aerophys%rd'//chmode//chcall, i_dgtime )
    END DO
    CALL save_dg(k, dustphys(k)%N(1), 'dustphys%N(1)'//chcall, i_dgtime )
    CALL save_dg(k, aeroact(k)%nact1, 'aeroact%nact1'//chcall, i_dgtime )
    CALL save_dg(k, aeroact(k)%mact1, 'aeroact%mact1'//chcall, i_dgtime )
    CALL save_dg(k, aeroact(k)%rcrit1, 'aeroact%rcrit1'//chcall, i_dgtime )
    CALL save_dg(k, aeroact(k)%mact1_mean, 'aeroact%mact1_mean'//chcall, i_dgtime )
    CALL save_dg(k, aeroact(k)%nact2, 'aeroact%nact2'//chcall, i_dgtime )
    CALL save_dg(k, aeroact(k)%rcrit2, 'aeroact%rcrit2'//chcall, i_dgtime )
    CALL save_dg(k, aeroact(k)%mact2, 'aeroact%mact2'//chcall, i_dgtime )
    CALL save_dg(k, aeroact(k)%mact2_mean, 'aeroact%mact2_mean'//chcall, i_dgtime )
  ELSE
    DO imode=1,aerophys(k)%nmodes
      WRITE(chmode, '(i2)') imode
      CALL save_dg(k, i_here, aerophys(k)%N(imode), 'aerophys%N'//chmode//chcall, i_dgtime )
      CALL save_dg(k, i_here, aerophys(k)%M(imode), 'aerophys%M'//chmode//chcall, i_dgtime )
      CALL save_dg(k, i_here, aerophys(k)%rd(imode), 'aerophys%rd'//chmode//chcall, i_dgtime )
    END DO
    CALL save_dg(k, i_here, dustphys(k)%N(1), 'dustphys%N(1)'//chcall, i_dgtime )
    CALL save_dg(k, i_here, aeroact(k)%nact1, 'aeroact%nact1'//chcall, i_dgtime )
    CALL save_dg(k, i_here, aeroact(k)%mact1, 'aeroact%mact1'//chcall, i_dgtime )
    CALL save_dg(k, i_here, aeroact(k)%rcrit1, 'aeroact%rcrit1'//chcall, i_dgtime )
    CALL save_dg(k, i_here, aeroact(k)%mact1_mean, 'aeroact%mact1_mean'//chcall, i_dgtime )
    CALL save_dg(k, i_here, aeroact(k)%nact2, 'aeroact%nact2'//chcall, i_dgtime )
    CALL save_dg(k, i_here, aeroact(k)%rcrit2, 'aeroact%rcrit2'//chcall, i_dgtime )
    CALL save_dg(k, i_here, aeroact(k)%mact2, 'aeroact%mact2'//chcall, i_dgtime )
    CALL save_dg(k, i_here, aeroact(k)%mact2_mean, 'aeroact%mact2_mean'//chcall, i_dgtime )
  END IF
END DO
#endif


END SUBROUTINE examine_aerosol


SUBROUTINE AbdulRazzakGhan(w, qc, p, T, phys, chem, nccn, Smax, rcrit, tau)

    ! Note: despite appearences this is only coded up
    ! for a single mode so far...

REAL(wp), INTENT(IN) :: w ! vertical velocity (ms-1)
REAL(wp), INTENT(IN) :: qc ! Pre-existing cloud mass
REAL(wp), INTENT(IN) :: p ! pressure (Pa)
REAL(wp), INTENT(IN) :: T ! temperature (K)
TYPE(aerosol_phys), INTENT(IN) :: phys
TYPE(aerosol_chem), INTENT(IN) :: chem
REAL(wp), INTENT(OUT) :: Nccn ! number of activated aerosol

REAL(wp), INTENT(OUT) :: smax ! peak supersaturation
REAL(wp), INTENT(OUT) :: rcrit ! radius of smallest activated aerosol
REAL(wp), INTENT(OUT) :: tau   ! equilibrium adjustment timescale


REAL(wp), ALLOCATABLE :: s_cr(:)
REAL(wp) :: Tc, Ak, Bk, eta, alpha, gamma, es, bigG, Galt
REAL(wp) :: f1, f2, zeta, error_func
INTEGER :: i

ALLOCATE(s_cr(phys%nmodes))

Tc = T-273.15 ! Temperature (C)

    ! surface tension of water, not solute (units are N/m)
    !zetasa = (76.1-(0.155*Tc)) * 1e-3

Ak = 2.0*Mw*zetasa/(Ru*T*rhow)

alpha = g*(Lv/(eps*cp*T)-1)/(T*Rd)

es = (100.0*6.1121)*EXP((18.678-Tc/(234.5))*Tc/(257.14+Tc))

gamma= eps*p/es+Lv**2/(Rv*T**2*cp)

    !Dv = (0.211*((T/273.15)**(1.94))*((100000.)/(p))) / 1.e4

bigG = 1.0/(rhow*(Rv*T/(es*Dv)+Lv*(Lv/(Rv*T)-1)/(ka*T)))

DO i=1,1!phys%nmodes

  Bk = chem%vantHoff(i)*Mw*chem%density(i)/      &
            (chem%massMole(i)*rhow)

  s_cr(i) = (2.0/SQRT(Bk))*(Ak/(3.0*phys%rd(i)))**1.5

  eta = (w*alpha/BigG)**1.5/      &
            (2.0*pi*rhow*gamma*phys%N(i))

  zeta = (2.0/3.0)*Ak*(w*alpha/BigG)**0.5

  f1 = 1.5*EXP(2.25*(LOG(phys%sigma(i)))**2)

  f2 = 1.0 + 0.25*LOG(phys%sigma(i))

  error_func = 1.0 - erf(LOG(f1*(zeta/eta)**1.5 +                  &
            f2*((s_cr(i)*s_cr(i))/(eta+3.0*zeta))**.75)&
            /(3.0*SQRT(2.0)*LOG(phys%sigma(i))))


  nccn = 0.5*phys%N(i)*error_func

  smax = s_cr(i)/(f1*(zeta/eta)**1.5      &
            + f2*(s_cr(i)*s_cr(i)/eta)**.75)**.5

  rcrit = phys%rd(i) * (s_cr(i)/smax)**(.66667)

END DO

tau = 1.0/(alpha*w + gamma*4*pi*rhow*bigG*qc)

DEALLOCATE(s_cr)

END SUBROUTINE AbdulRazzakGhan


SUBROUTINE AbdulRazzakGhan2000(w, p, T, phys, chem, nccn, Smax, active_phys   &
     , nccn_active, l_useactive )

    ! Note: despite appearences this is only coded up
    ! for a single mode so far...

REAL(wp), INTENT(IN) :: w ! vertical velocity (ms-1)
REAL(wp), INTENT(IN) :: p ! pressure (Pa)
REAL(wp), INTENT(IN) :: T ! temperature (K)
TYPE(aerosol_phys), INTENT(IN) :: phys
TYPE(aerosol_chem), INTENT(IN) :: chem
TYPE(aerosol_active), INTENT(IN), OPTIONAL :: active_phys
REAL(wp), INTENT(OUT) :: nccn_active ! notional nccn that would be generated by activated aerosol
REAL(wp), INTENT(OUT) :: Nccn(:) ! number of activated aerosol
                                     ! for each mode

REAL(wp), INTENT(OUT) :: smax ! peak supersaturation
LOGICAL, INTENT(OUT) :: l_useactive ! Do we use consider already activated aerosol

REAL(wp), ALLOCATABLE :: s_cr(:)
REAL(wp) :: s_cr_active
REAL(wp) :: Tc, Ak, Bk, eta, alpha, gamma, es, bigG, Galt
REAL(wp) :: f1, f2, zeta, error_func
REAL(wp) :: rsmax2  ! 1/(smax*smax)
REAL(wp) :: diff
INTEGER :: i
REAL(wp) :: sum

l_useactive=.FALSE.
IF (PRESENT(active_phys)) l_useactive = active_phys%nact > ccn_tidy .AND. l_active_inarg2000

nccn_active=0.0

ALLOCATE(s_cr(phys%nmodes))

Tc = T-273.15 ! Temperature (C)

    ! surface tension of water, not solute (units are N/m)
    !zetasa = (76.1-(0.155*TC)) * 1e-3

Ak = 2.0*Mw*zetasa/(Ru*T*rhow)

alpha = g*(Lv/(eps*cp*T)-1)/(T*Rd)

es = (100.0*6.1121)*EXP((18.678-Tc/(234.5))*Tc/(257.14+Tc))

gamma= eps*p/es+Lv**2/(Rv*T**2*cp)

    !Dv = (0.211*((T/273.15)**(1.94))*((100000.)/(p))) / 1.e4

bigG = 1.0/(rhow*(Rv*T/(es*Dv)+Lv*(Lv/(Rv*T)-1)/(ka*T)))

zeta = (2.0/3.0)*Ak*(w*alpha/BigG)**0.5
rsmax2 = 0.0

DO i=1,phys%nmodes

  IF (phys%N(i) > ccn_tidy) THEN

    Bk = chem%vantHoff(i)*Mw*chem%density(i)/     &
           (chem%massMole(i)*rhow)

    s_cr(i) = (2.0/SQRT(Bk))*(Ak/(3.0*phys%rd(i)))**1.5

    eta = (w*alpha/BigG)**1.5/     &
           (2.0*pi*rhow*gamma*phys%N(i))

    f1 = 0.5*EXP(2.5*(LOG(phys%sigma(i)))**2)

    f2 = 1.0 + 0.25*LOG(phys%sigma(i))

    rsmax2 = rsmax2 + (f1*(zeta/eta)**1.5     &
           + f2*(s_cr(i)*s_cr(i)/(eta+3.0*zeta))**.75)/(s_cr(i)*s_cr(i))

  END IF

END DO

IF (l_useactive) THEN  ! We should include all the activated aerosol

  IF (active_phys%nact > ccn_tidy) THEN

    i=aero_index%i_accum ! use accumulation mode chem (should do this properly)

    Bk = chem%vantHoff(i)*Mw*chem%density(i)/     &
           (chem%massMole(i)*rhow)

    s_cr_active = (2.0/SQRT(Bk))*(Ak/(3.0*active_phys%rd))**1.5

    eta = (w*alpha/BigG)**1.5/     &
           (2.0*pi*rhow*gamma*active_phys%nact)

    f1 = 0.5*EXP(2.5*(LOG(active_phys%sigma))**2)

    f2 = 1.0 + 0.25*LOG(active_phys%sigma)


    rsmax2 = rsmax2 + (f1*(zeta/eta)**1.5     &
           + f2*(s_cr_active*s_cr_active/(eta+3.0*zeta))**.75)/(s_cr_active*s_cr_active)

  END IF

END IF
smax=SQRT(1.0/rsmax2)

nccn = 0.0

IF (l_useactive) THEN  ! We should include all the activated aerosol

  IF (active_phys%nact > ccn_tidy) THEN

    error_func = 1.0 - erf(2.0*LOG(s_cr_active/smax)/(3.0*SQRT(2.0)*LOG(active_phys%sigma)))

    nccn_active = 0.5*active_phys%nact*error_func

  END IF

END IF

DO i=1,phys%nmodes

  IF (phys%N(i) > ccn_tidy) THEN

    error_func = 1.0 - erf(2.0*LOG(s_cr(i)/smax)/(3.0*SQRT(2.0)*LOG(phys%sigma(i))))

    nccn(i) = 0.5*phys%N(i)*error_func

        ! Make sure we don't activate too many...
    nccn(i) = MIN(nccn(i), .999*phys%N(i))
        ! And don't bother if we don't do much
    IF (nccn(i) < ccn_tidy) nccn(i)=0.0

  END IF

END DO

IF (l_useactive) THEN  ! Remove already activated aerosol
                          ! from newly activated
                          ! Start from smallest mode
  diff  = active_phys%nact - nccn_active

  sum=0.0
  DO i=1,SIZE(nccn)
    sum = nccn(i) + sum
    nccn(i) = MIN(nccn(i), MAX(0.0_wp, sum - diff))
  END DO

      !nccn(1) = max(0.0, nccn(1) - diff)
      !nccn(2) = min(nccn(2), max(0.0, nccn(2) + nccn(1) - diff))
      !nccn(3) = min(nccn(3), max(0.0, nccn(3) + nccn(2) + nccn(1) - diff))

END IF

nccn=.99*nccn ! Don't allow all aerosol to be removed.


DEALLOCATE(s_cr)

END SUBROUTINE AbdulRazzakGhan2000

FUNCTION moment_logn(N, rm, sigma, p)
    !
    ! Calculate the moments of a lognormal distribution
    !
    ! int_0^infinty r^p*n(r)dr
    !
    ! where
    !
    !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
    !            *exp(-.5*(log(r/rm)/log(sigma))^2)
    !
    ! rm = median radius
    ! log(rm) = mean of log(n(r)/N)
    ! log(sigma) = s.d. of log(n(r)/N)
    !

REAL(wp) :: N, rm, sigma
REAL(wp) :: p ! calculate pth moment
REAL(wp) :: moment_logn

moment_logn=N*rm**p*EXP(.5*p*p*LOG(sigma)**2)

END FUNCTION moment_logn

FUNCTION upperpartial_moment_logn(N, rm, sigma, p, rcrit)
    !
    ! Calculate the upper partial moment of a
    ! lognormal distribution
    !
    ! int_rcrit^infinty r^p*n(r)dr
    !
    ! where
    !
    !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
    !            *exp(-.5*(log(r/rm)/log(sigma))^2)
    !
    ! rm = median radius
    ! log(rm) = mean of log(n(r)/N)
    ! log(sigma) = s.d. of log(n(r)/N)
    !

REAL(wp) :: N, rm, sigma
REAL(wp) :: p ! calculate pth moment
REAL(wp) :: rcrit ! lower threshold for partial moment
REAL(wp) :: upperpartial_moment_logn

IF (rcrit==0.0) THEN
  upperpartial_moment_logn=moment_logn(N, rm, sigma, p)
ELSE
  upperpartial_moment_logn=N*rm**p*EXP(.5*p*p*LOG(sigma)**2)     &
         * .5*erfc((LOG(rcrit/rm)/LOG(sigma) - p*LOG(sigma))/SQRT(2.0))
END IF

END FUNCTION upperpartial_moment_logn


FUNCTION invert_partial_moment(m, mup, p, rm, sigma)
    !
    ! Calculate the lower bound of integral given
    ! a complete and upper partial moment
    !
    ! m = int_0^infinty r^p*n(r)dr
    ! mup = int_x^infinty r^p*n(r)dr
    !
    ! where
    !
    !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
    !            *exp(-.5*(log(r/rm)/log(sigma))^2)
    !
    ! rm = median radius
    ! log(rm) = mean of log(n(r)/N)
    ! log(sigma) = s.d. of log(n(r)/N)
    !

REAL(wp) :: m, mup
REAL(wp) :: p ! pth moment
REAL(wp) :: rm, sigma
REAL(wp) :: x, invert_partial_moment

x=rm * EXP(p*LOG(sigma)**2                 &
               + SQRT(2.0)*LOG(sigma)*erfinv(1.0-2.0*mup/m))

invert_partial_moment = x

END FUNCTION invert_partial_moment

FUNCTION invert_partial_moment_approx(mup, p, rm, sigma)
    !
    ! Calculate the lower bound of integral given
    ! a complete and upper partial moment
    !
    ! mup = int_x^infinty r^p*n(r)dr
    !
    ! where
    !
    !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
    !            *exp(-.5*(log(r/rm)/log(sigma))^2)
    !
    ! rm = median radius
    ! log(rm) = mean of log(n(r)/N)
    ! log(sigma) = s.d. of log(n(r)/N)
    !
    ! This uses an approximation for erfc(log(x)) to achieve this

REAL(wp) :: mup
REAL(wp) :: p ! pth moment
REAL(wp) :: rm, sigma
REAL(wp) :: beta, c, lsig, mbeta
REAL(wp) :: x, invert_partial_moment_approx

c=-0.1021 ! constant in approximation
lsig=SQRT(2.0)*LOG(sigma)

beta=mup*EXP(-0.25*p*p*lsig*lsig)/rm**p
mbeta=1.0-beta

IF (beta < 1.0e-4) THEN
  x=rm *(-c)**(-lsig)
ELSE IF (beta >.9999) THEN
  x=0.0
ELSE
  x=rm * (sigma**(-0.5*p)*(     &
         ((mbeta)*c+SQRT((mbeta)**2*c*c + 4.0*beta*(mbeta)))/(2.0*beta) &
         ))**(lsig)
END IF

invert_partial_moment_approx = x

END FUNCTION invert_partial_moment_approx


FUNCTION invert_partial_moment_betterapprox(mup, p, N, rm, sigma)
    !
    ! Calculate the lower bound of integral given
    ! a complete and upper partial moment
    !
    ! mup = int_x^infinty r^p*n(r)dr
    !
    ! where
    !
    !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
    !            *exp(-.5*(log(r/rm)/log(sigma))^2)
    !
    ! rm = median radius
    ! log(rm) = mean of log(n(r)/N)
    ! log(sigma) = s.d. of log(n(r)/N)
    !
    ! This uses an approximation 26.2.23 from A&S

REAL(wp) :: mup
REAL(wp) :: p ! pth moment
REAL(wp) :: rm, sigma, N
REAL(wp) :: frac, x, invert_partial_moment_betterapprox
REAL(wp) :: small_frac=1e-6

IF (mup==0.0) THEN
  frac=EPSILON(x)
ELSE


!              if (mype==240 .and. timestep_number>193 .and. timestep_number<195)write(mype+3000,*)' DEBUG mup', N,rm,p

  frac = mup/(N*rm**p)*EXP(-.5*p*p*LOG(sigma)*LOG(sigma))
END IF

IF (frac>=1.0 - EPSILON(x))frac=.999

IF (frac > small_frac) THEN

  x=normal_quantile(frac)

  invert_partial_moment_betterapprox=rm*EXP(x*LOG(sigma)+p*LOG(sigma)*LOG(sigma))

ELSE

  invert_partial_moment_betterapprox=0.0

END IF

END FUNCTION invert_partial_moment_betterapprox

FUNCTION normal_quantile(p)
    ! This uses an approximation 26.2.23 from A&S
REAL(wp) :: p, pp

REAL(wp) :: normal_quantile
REAL(wp) :: t,x
REAL(wp) :: c0,c1,c2,d1,d2,d3

c0=2.515517
c1=0.802853
c2=0.010328
d1=1.432788
d2=0.189269
d3=0.001308

pp=p

IF (p<=EPSILON(x)) THEN

  x=HUGE(x)

ELSE IF (p >=1.0 - EPSILON(x)) THEN

  x=-HUGE(x)

ELSE

  IF (p>0.5)pp=1.0-p

  t=SQRT(LOG(1.0/(pp*pp)))

  x = t - (c0+c1*t+c2*t*t)/(1.0+d1*t+d2*t*t+d3*t*t*t)

  IF (p>0.5)x=-x

END IF

normal_quantile = x

END FUNCTION normal_quantile

FUNCTION MNtoRm(M, N, density, sigma)
    !
    ! Calculate mean radius of lognormal distribution
    ! given Mass and number
    !

REAL(wp) :: M, N, density, sigma
REAL(wp) :: MNtoRm
IF (N==0 .OR. M==0) THEN ! shouldn't really be here
  MNtoRm = 0.0
ELSE
  MNtoRm = (3.0*M*EXP(-4.5*LOG(sigma)**2)/(4.0*N*pi*density))**(1.0/3.0)
END IF

END FUNCTION MNtoRm

END MODULE aerosol_routines

