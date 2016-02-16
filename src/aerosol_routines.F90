module aerosol_routines
  use mphys_die, only: throw_mphys_error
  use type_aerosol, only: aerosol_type, aerosol_phys, aerosol_chem, aerosol_active
  use variable_precision, only: wp
  use mphys_constants, only: Mw, zetasa, Ru, rhow, g, Lv, eps, cp, Rd, Rv, ka, Dv, pi
  use special, only: erfc, erfinv
  use mphys_switches, only: i_am2, i_an2, i_am4, i_nl, i_nr, i_am1, i_an1, i_am3, i_an3, &
       i_am5, i_am6, i_an6, i_am7, i_am8, i_am9, i_an11, i_an12, i_ni, i_ns, i_ng, i_ql, &
       i_qr, i_qi, i_qs, i_qg, l_active_inarg2000, aero_index, l_warm, active_cloud, active_rain, &
       active_ice, isol, iinsol, l_process, l_passivenumbers, l_passivenumbers_ice, l_separate_rain
  use thresholds, only: nr_tidy, nl_tidy, ni_tidy, ccn_tidy, qr_tidy, aeromass_small, aeronumber_small
  use mphys_parameters, only: sigma_arc, nz
  use preconditioning, only: precondition

  ! Really want to be allocating this elsewhere

#if DEF_MODEL==MODEL_KiD
  use diagnostics, only: save_dg, i_dgtime, k_here, i_here, j_here, nx
  use runtime, only: time
#elif  DEF_MODEL==MODEL_LEM
  use diaghelp_lem, only: k_here, i_here, j_here
  use com_params, only: time
#elif  DEF_MODEL==MODEL_UM
  use timestep_mod, only: timestep_number
  use diaghelp_um, only: k_here, i_here, j_here
  use UM_ParCore, only: mype
#elif  DEF_MODEL==MODEL_MONC
  use diaghelp_monc, only: i_here, j_here
#endif

  implicit none
  private

  public aerosol_active, aerosol_phys, aerosol_chem, abdulRazzakGhan2000, invert_partial_moment, upperpartial_moment_logn, &
       invert_partial_moment_approx, invert_partial_moment_betterapprox, examine_aerosol, allocate_aerosol, &
       deallocate_aerosol, MNtoRm
contains
  !
  ! Allocate space for aerosol_chem and aerosol_phys
  !
  subroutine allocate_aerosol(aerophys, aerochem, nmodes, initphys, initchem)
    type(aerosol_phys), intent(inout) :: aerophys(:)
    type(aerosol_chem), intent(inout) :: aerochem(:)
    integer, intent(in) :: nmodes    ! number of modes
    type(aerosol_phys), intent(in), optional :: initphys
    type(aerosol_chem), intent(in), optional :: initchem

    integer :: k

    !------------------------
    ! First allocate aerophys
    !------------------------
    do k=1, size(aerophys)
      aerophys(k)%nmodes=nmodes
      allocate(aerophys(k)%N(max(1,nmodes)))
      allocate(aerophys(k)%M(max(1,nmodes)))
      allocate(aerophys(k)%rd(max(1,nmodes)))
      allocate(aerophys(k)%sigma(max(1,nmodes)))
      allocate(aerophys(k)%rpart(max(1,nmodes)))

      if (present(initphys)) aerophys(k)=initphys
    end do

    !------------------------
    ! Then allocate aerochem
    !------------------------
    do k=1, size(aerochem)
      aerochem(k)%nmodes=nmodes
      allocate(aerochem(k)%vantHoff(max(1,nmodes)))
      allocate(aerochem(k)%massMole(max(1,nmodes)))
      allocate(aerochem(k)%density(max(1,nmodes)))
      allocate(aerochem(k)%epsv(max(1,nmodes)))
      allocate(aerochem(k)%beta(max(1,nmodes)))

      if (present(initphys)) aerochem(k)=initchem
    end do
  end subroutine allocate_aerosol

  !
  ! Deallocate space for aerosol_chem and aerosol_phys
  !
  subroutine deallocate_aerosol(aerophys, aerochem)  
    type(aerosol_phys), intent(inout) :: aerophys(:)
    type(aerosol_chem), intent(inout) :: aerochem(:)

    integer :: k

    !------------------------
    ! First deallocate aerophys
    !------------------------
    do k=1, size(aerophys)
      deallocate(aerophys(k)%N)
      deallocate(aerophys(k)%M)
      deallocate(aerophys(k)%rd)
      deallocate(aerophys(k)%sigma)
      deallocate(aerophys(k)%rpart)
    end do

    !------------------------
    ! Then deallocate aerochem
    !------------------------
    do k=1, size(aerochem)
      deallocate(aerochem(k)%vantHoff)
      deallocate(aerochem(k)%massMole)
      deallocate(aerochem(k)%density)
      deallocate(aerochem(k)%epsv)
      deallocate(aerochem(k)%beta)
    end do

  end subroutine deallocate_aerosol

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
  subroutine examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact, &
       dustphys, dustchem, dustact, aeroice, dustliq, icall)
    real(wp), intent(in), target :: aerofields(:,:)
    type(aerosol_phys), intent(inout), target :: aerophys(:)
    type(aerosol_chem), intent(in), target :: aerochem(:)
    type(aerosol_active), intent(inout) :: aeroact(:)
    type(aerosol_phys), intent(inout), target, optional :: dustphys(:)
    type(aerosol_chem), intent(in), target, optional :: dustchem(:)
    type(aerosol_active), intent(inout), optional :: dustact(:)
    type(aerosol_active), intent(inout), optional :: aeroice(:)
    type(aerosol_active), intent(inout), optional :: dustliq(:)
    real(wp), intent(inout), target :: qfields(:,:)
    integer, intent(in), optional :: icall

    real(wp) :: n, m, mac, mar, mad, maai, madl, cloud_number, rain_number
    real(wp) :: cloud_mass, rain_mass
    real(wp) :: ice_number, snow_number, graupel_number
    real(wp) :: ice_mass, snow_mass, graupel_mass
    real(wp) :: mact2, rcrit
    real(wp) :: Nd, rm, sigma, density
    integer :: k

    real(wp) :: mtot, tmp, m_ratio, n_ratio, dcloud_number, rm_arc
    character(2) :: chcall, chmode
    real(wp) :: ratio_l, ratio_r, nratio_l, nratio_r, mratio_l, mratio_r
    real(wp) :: ratio_i, ratio_s, ratio_g, nratio_i, nratio_s, nratio_g, mratio_i, mratio_s, mratio_g
    real(wp) :: ratio_dil, ratio_dli, ratio_ali, ratio_ail

    integer :: imode
    real(wp) :: mode_N, mode_M
    real(wp) :: nitot, ntot, nhtot, mitot
    real(wp) :: mact, rcrit2
    real(wp) :: foo

    logical :: l_condition, l_condition_r

    ! Some diagnostic strings
    if (present(icall)) then
      write(chcall,'(A1,i1)') '_', icall
    else
      chcall=''
    end if

    do k=1, nz
      if (precondition(k)) then
        cloud_number=qfields(k,i_nl)
        rain_number=qfields(k,i_nr)
        cloud_mass=qfields(k,i_ql)
        rain_mass=qfields(k,i_qr)
        density=aerochem(k)%density(1)

        do imode=1, aero_index%nccn
          mode_N=aerofields(k,aero_index%ccn_n(imode))
          mode_M=aerofields(k,aero_index%ccn_m(imode))
          if (mode_N > ccn_tidy .and. mode_M > ccn_tidy*1e-18) then ! FiX This
            aerophys(k)%N(imode)=mode_N
            aerophys(k)%M(imode)=mode_M
            aerophys(k)%rd(imode)=MNtoRm(mode_M,mode_N,density,aerophys(k)%sigma(imode))
          else
            aerophys(k)%N(imode)=0.0
            aerophys(k)%M(imode)=0.0
            aerophys(k)%rd(imode)=0.0
          end if
        end do

        aeroact(k)%nact=0.0
        aeroact(k)%mact=0.0
        aeroact(k)%rcrit=999.0
        aeroact(k)%mact_mean=0.0
        aeroact(k)%nact2=0.0
        aeroact(k)%nact1=0.0
        aeroact(k)%mact1=0.0
        aeroact(k)%rcrit1=999.0
        aeroact(k)%nact2=0.0
        aeroact(k)%rcrit2=999.0
        aeroact(k)%mact2=0.0
        aeroact(k)%mact2_mean=0.0
        aeroact(k)%mact1_mean=0.0

        aeroact(k)%nratio1=0.0
        aeroact(k)%nratio2=0.0

        if (l_process) then
          ! Examine activated aerosol
          mac=aerofields(k,i_am4)

          if (l_separate_rain) then
            mar=aerofields(k,i_am5)
            mact=mac+mar
            l_condition=mac > 0.0 .and. cloud_number > nl_tidy
            l_condition_r=mar > 0.0 .and. rain_number > nr_tidy .and. rain_mass > qr_tidy
          else
            mact=mac
            l_condition=mac > 0.0 .and. cloud_number + rain_number > nl_tidy
            l_condition_r=mac > 0.0 .and. rain_number > nr_tidy .and. rain_mass > qr_tidy
          end if

          if (l_warm)then
            maai=0.0
          else
            maai=aerofields(k,i_am8)
          end if

          ratio_ali=mact/(mact+maai+1e-38)

          nhtot=cloud_number+rain_number

          if (l_passivenumbers) then
            ntot=aerofields(k,i_an11)*ratio_ali
          else
            ntot=nhtot
          end if

          l_condition=l_condition .and. ntot > 0.0

          nhtot=nhtot+1e-38 ! prevent possible later division by zero
          nratio_l=cloud_number/nhtot
          nratio_r=rain_number/nhtot

          if (l_condition) then
            aeroact(k)%nact=ntot
            aeroact(k)%mact=mac
            aeroact(k)%rcrit=0.0
            aeroact(k)%mact_mean=aeroact(k)%mact /(aeroact(k)%nact + tiny(ntot))

            ! Get mean radius of distribution
            rm_arc=MNtoRm(aeroact(k)%mact,aeroact(k)%nact,density,sigma_arc)

            aeroact(k)%rd=rm_arc
            aeroact(k)%sigma=sigma_arc

          end if

          if (l_condition_r) then
            if (l_separate_rain) then! Separate rain category
              rcrit2=0
              mact2=mar
            else if (.not. l_condition) then ! No cloud here
              rcrit2=0
              mact2=mac
            else ! Diagnostic partitioning of aerosol between cloud and rain.
              rcrit2=0.0
              mact2=aeroact(k)%mact*rain_mass/(cloud_mass + rain_mass)
            end if

            aeroact(k)%nact2=ntot*nratio_r
            aeroact(k)%rcrit2=rcrit2
            aeroact(k)%mact2=mact2
            aeroact(k)%mact2_mean=aeroact(k)%mact2/(aeroact(k)%nact2 + tiny(ntot))
          end if

          aeroact(k)%nact1=max(0.0_wp, aeroact(k)%nact-aeroact(k)%nact2)
          aeroact(k)%mact1=max(0.0_wp, aeroact(k)%mact-aeroact(k)%mact2)
          aeroact(k)%rcrit1=0.0
          aeroact(k)%mact1_mean=aeroact(k)%mact1/(aeroact(k)%nact1+tiny(mar))

          if (cloud_number > 0.0) aeroact(k)%nratio1=aeroact(k)%nact1/cloud_number
          if (rain_number > 0.0) aeroact(k)%nratio2=aeroact(k)%nact2/rain_number
        end if
      end if
    end do

    if (.not. l_warm) then
      ! Activated dust
      do k=1,nz
        density=dustchem(k)%density(1)
        ice_number=qfields(k,i_ni)
        snow_number=qfields(k,i_ns)
        nhtot=ice_number+snow_number
        if (i_ng > 0) then
          graupel_number=qfields(k,i_ng)
          nhtot=nhtot+graupel_number
        end if

        if (l_process) then
          mad=aerofields(k,i_am7)
          madl=aerofields(k,i_am9)
          ratio_dil=mad/(madl+mad+tiny(mad))
        end if

        if (l_passivenumbers_ice) then
          nitot=aerofields(k,i_an12)*ratio_dil
        else
          nitot=nhtot
        end if

        nhtot=nhtot+tiny(nhtot) ! prevent possible later division by zero
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
        do imode=1, aero_index%nin
          mode_N=aerofields(k,aero_index%in_n(imode))
          mode_M=aerofields(k,aero_index%in_m(imode))
          if (mode_m < 0 .or. mode_n < 0 ) then
            print*, aero_index%in_m(imode), aerofields(k,:)
            print*, 'Problem! Using pragmatic hack to continue.  If you see this message, report it to Ben!!!'
            print*, k, chcall, imode, mode_m, mode_n, nitot, ice_number, snow_number, graupel_number
            if (mode_m < 0) mode_m=aeromass_small
            if (mode_n < 0) mode_n=aeronumber_small
          end if

          dustphys(k)%N(imode)=mode_N
          dustphys(k)%M(imode)=mode_M

          dustphys(k)%rd(imode)=MNtoRm(mode_M,mode_N,density,dustphys(k)%sigma(imode))
        end do

        ! Examine activated dust
        ! Initialize to zero/defaults
        dustact(k)%nact=0.0
        dustact(k)%mact=0.0
        dustact(k)%rcrit=999.0
        dustact(k)%mact_mean=0.0
        dustact(k)%nact1=0.0
        dustact(k)%nratio1=0.0
        dustact(k)%mact1=0.0
        dustact(k)%rcrit1=999.0
        dustact(k)%mact1_mean=0.0
        dustact(k)%mact2=0.0
        dustact(k)%nact2=0.0
        dustact(k)%nratio2=0.0
        dustact(k)%rcrit2=999.0
        dustact(k)%mact2_mean=0.0
        dustact(k)%mact3=0.0
        dustact(k)%nact3=0.0
        dustact(k)%nratio3=0.0
        dustact(k)%rcrit3=999.0
        dustact(k)%mact3_mean=0.0

        if (l_process) then
          if (mad > 0.0 .and. nitot > ni_tidy) then

            ! Equal partitioning of dust distribution across all ice species
            ! This will cause problems, e.g. when ice splinters, there will be a
            ! spurious increase in the distributed dust mass in the ice category
            ! (This isn't a problem if there is no snow or graupel)

            dustact(k)%nact=nitot
            dustact(k)%mact=mad
            dustact(k)%rcrit=0.0
            dustact(k)%mact_mean=dustact(k)%mact/(dustact(k)%nact+tiny(mad))

            if (ice_number > 0.0) then
              dustact(k)%nact1=nitot * ratio_i
              dustact(k)%mact1=mad * dustact(k)%nact1/dustact(k)%nact
              dustact(k)%rcrit1=0.0
              dustact(k)%mact1_mean=dustact(k)%mact1/(dustact(k)%nact1+tiny(nhtot))
              dustact(k)%nratio1=dustact(k)%nact1/(ice_number+tiny(nhtot))
            end if

            if (snow_number > 0.0) then
              dustact(k)%nact2=nitot*ratio_s
              dustact(k)%rcrit2=0.0
              dustact(k)%mact2=mad*dustact(k)%nact2/(dustact(k)%nact+tiny(nhtot))
              dustact(k)%mact2_mean=dustact(k)%mact2/(dustact(k)%nact2+tiny(nhtot))
              dustact(k)%nratio2=dustact(k)%nact2/(snow_number+tiny(nhtot))
            end if

            if (i_ng > 0 .and. graupel_number > ni_tidy) then
              dustact(k)%nact3=nitot*ratio_g
              dustact(k)%rcrit3=0.0
              dustact(k)%mact3=mad*dustact(k)%nact3/(dustact(k)%nact+tiny(nhtot))
              dustact(k)%mact3_mean=dustact(k)%mact3/(dustact(k)%nact3+tiny(nhtot))
              dustact(k)%nratio3=dustact(k)%nact3/(graupel_number+tiny(nhtot))
            end if
          end if
        end if
      end do

      if (present(aeroice)) then
        if (l_process) then
          ! Activated aerosol in ice
          do k=1, nz
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
            if (i_ng > 0) then
              graupel_number=qfields(k,i_ng)
              nhtot=nhtot+graupel_number
              mitot=mitot+graupel_mass
            end if

            if (l_passivenumbers) then
              nitot=aerofields(k,i_an11)*ratio_ail
            else
              nitot=nhtot
            end if

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
            aeroice(k)%nact=0.0
            aeroice(k)%mact=0.0
            aeroice(k)%rcrit=999.0
            aeroice(k)%mact_mean=0.0
            aeroice(k)%nact1=0.0
            aeroice(k)%nratio1=0.0
            aeroice(k)%mact1=0.0
            aeroice(k)%rcrit1=999.0
            aeroice(k)%mact1_mean=0.0
            aeroice(k)%mact2=0.0
            aeroice(k)%nact2=0.0
            aeroice(k)%nratio2=0.0
            aeroice(k)%rcrit2=999.0
            aeroice(k)%mact2_mean=0.0
            aeroice(k)%mact3=0.0
            aeroice(k)%nact3=0.0
            aeroice(k)%nratio3=0.0
            aeroice(k)%rcrit3=999.0
            aeroice(k)%mact3_mean=0.0

            if (maai > 0.0 .and. nitot > ni_tidy) then
              ! Equal partitioning of dust distribution across all ice species
              aeroice(k)%nact=nitot
              aeroice(k)%mact=maai
              aeroice(k)%rcrit=0.0
              aeroice(k)%mact_mean=aeroice(k)%mact/(aeroice(k)%nact+tiny(nhtot))

              if (ratio_i > 0.0) then
                aeroice(k)%nact1=nitot*ratio_i
                aeroice(k)%mact1=maai*ratio_i
                aeroice(k)%rcrit1=0.0
                aeroice(k)%mact1_mean=aeroice(k)%mact1/(aeroice(k)%nact1+tiny(nhtot))
                aeroice(k)%nratio1=aeroice(k)%nact1/(ice_number+tiny(nhtot))
              end if

              if (ratio_s > 0.0) then
                aeroice(k)%nact2=nitot*ratio_s
                aeroice(k)%rcrit2=0.0
                aeroice(k)%mact2=maai*ratio_s
                aeroice(k)%mact2_mean=aeroice(k)%mact2/(aeroice(k)%nact2+tiny(nhtot))
                aeroice(k)%nratio2=aeroice(k)%nact2/(snow_number+tiny(nhtot))
              end if

              if (ratio_g > 0.0) then
                aeroice(k)%nact3=nitot*ratio_g
                aeroice(k)%rcrit3=0.0
                aeroice(k)%mact3=maai*ratio_g
                aeroice(k)%mact3_mean=aeroice(k)%mact3/(aeroice(k)%nact3+tiny(nhtot))
                aeroice(k)%nratio3=aeroice(k)%nact3/(graupel_number+tiny(nhtot))
              end if
            end if
          end do
        end if
      end if

      if (present(dustliq)) then
        if (l_process) then
          ! Activated aerosol in ice
          do k=1, nz
            density=dustchem(k)%density(1)
            madl=aerofields(k,i_am9)
            mad=aerofields(k,i_am7)
            ratio_dli = madl/(mad+madl+1e-38)

            cloud_number=qfields(k,i_nl)
            rain_number=qfields(k,i_nr)

            nratio_l=cloud_number/(cloud_number+rain_number+tiny(cloud_number))
            nratio_r=rain_number/(cloud_number+rain_number+tiny(cloud_number))
            mratio_l=cloud_mass/(cloud_mass+rain_mass+tiny(cloud_number))
            mratio_r=rain_mass/(cloud_mass+rain_mass+tiny(cloud_number))

            ! Use nratios...
            ratio_l=nratio_l
            ratio_r=nratio_r

            if (l_passivenumbers_ice) then
              ntot=aerofields(k,i_an12)*ratio_dli
            else
              ntot=cloud_number + rain_number
            end if

            ! Initialize to zero/defaults
            dustliq(k)%nact=0.0
            dustliq(k)%mact=0.0
            dustliq(k)%rcrit=999.0
            dustliq(k)%mact_mean=0.0
            dustliq(k)%nact1=0.0
            dustliq(k)%nratio1=0.0
            dustliq(k)%mact1=0.0
            dustliq(k)%rcrit1=999.0
            dustliq(k)%mact1_mean=0.0
            dustliq(k)%mact2=0.0
            dustliq(k)%nact2=0.0
            dustliq(k)%nratio2=0.0
            dustliq(k)%rcrit2=999.0
            dustliq(k)%mact2_mean=0.0
            dustliq(k)%mact3=0.0
            dustliq(k)%nact3=0.0
            dustliq(k)%nratio3=0.0
            dustliq(k)%rcrit3=999.0
            dustliq(k)%mact3_mean=0.0

            if (madl > 0.0 .and. ntot > nr_tidy) then
              ! Equal partitioning of aerosol distribution across all liquid species

              dustliq(k)%nact=ntot
              dustliq(k)%mact=madl
              dustliq(k)%rcrit=0.0
              dustliq(k)%mact_mean=dustliq(k)%mact/(dustliq(k)%nact+tiny(nhtot))

              if (ratio_l > 0.0) then
                dustliq(k)%nact1=ntot*ratio_l
                dustliq(k)%mact1=madl*ratio_l
                dustliq(k)%rcrit1=0.0
                dustliq(k)%mact1_mean=dustliq(k)%mact1/(dustliq(k)%nact1+tiny(nhtot))
                dustliq(k)%nratio1=dustliq(k)%nact1/(cloud_number+tiny(nhtot))
              end if

              if (ratio_r > 0.0) then
                dustliq(k)%nact2=ntot*ratio_r
                dustliq(k)%rcrit2=0.0
                dustliq(k)%mact2=madl*ratio_r
                dustliq(k)%mact2_mean=dustliq(k)%mact2/(dustliq(k)%nact2+tiny(nhtot))
                dustliq(k)%nratio2=dustliq(k)%nact2/(rain_number+tiny(nhtot))
              end if
            end if
          end do
        end if
      end if
    end if

#if DEF_MODEL==MODEL_KiD
    do k=1, nz
      if (nx==1) then
        call save_dg(k, dustact(k)%mact3_mean, 'dustact%mact3_mean'//chcall, i_dgtime )
        call save_dg(k, dustact(k)%mact_mean, 'dustact%mact_mean'//chcall, i_dgtime )
        call save_dg(k, dustact(k)%mact, 'dustact%mact'//chcall, i_dgtime)
        call save_dg(k, dustliq(k)%mact3_mean, 'dustliq%mact3_mean'//chcall, i_dgtime )
        call save_dg(k, dustliq(k)%mact_mean, 'dustliq%mact_mean'//chcall, i_dgtime )
        call save_dg(k, dustliq(k)%mact, 'dustliq%mact'//chcall, i_dgtime)
        do imode=1, aerophys(k)%nmodes
          write(chmode, '(i2)') imode
          call save_dg(k, aerophys(k)%N(imode), 'aerophys%N'//chmode//chcall, i_dgtime )
          call save_dg(k, aerophys(k)%M(imode), 'aerophys%M'//chmode//chcall, i_dgtime )
          call save_dg(k, aerophys(k)%rd(imode), 'aerophys%rd'//chmode//chcall, i_dgtime )
        end do
        call save_dg(k, dustphys(k)%N(1), 'dustphys%N(1)'//chcall, i_dgtime )
        call save_dg(k, aeroact(k)%nact1, 'aeroact%nact1'//chcall, i_dgtime )
        call save_dg(k, aeroact(k)%mact1, 'aeroact%mact1'//chcall, i_dgtime )
        call save_dg(k, aeroact(k)%rcrit1, 'aeroact%rcrit1'//chcall, i_dgtime )
        call save_dg(k, aeroact(k)%mact1_mean, 'aeroact%mact1_mean'//chcall, i_dgtime )
        call save_dg(k, aeroact(k)%nact2, 'aeroact%nact2'//chcall, i_dgtime )
        call save_dg(k, aeroact(k)%rcrit2, 'aeroact%rcrit2'//chcall, i_dgtime )
        call save_dg(k, aeroact(k)%mact2, 'aeroact%mact2'//chcall, i_dgtime )
        call save_dg(k, aeroact(k)%mact2_mean, 'aeroact%mact2_mean'//chcall, i_dgtime )
      else
        do imode=1, aerophys(k)%nmodes
          write(chmode, '(i2)') imode
          call save_dg(k, i_here, aerophys(k)%N(imode), 'aerophys%N'//chmode//chcall, i_dgtime )
          call save_dg(k, i_here, aerophys(k)%M(imode), 'aerophys%M'//chmode//chcall, i_dgtime )
          call save_dg(k, i_here, aerophys(k)%rd(imode), 'aerophys%rd'//chmode//chcall, i_dgtime )
        end do
        call save_dg(k, i_here, dustphys(k)%N(1), 'dustphys%N(1)'//chcall, i_dgtime )
        call save_dg(k, i_here, aeroact(k)%nact1, 'aeroact%nact1'//chcall, i_dgtime )
        call save_dg(k, i_here, aeroact(k)%mact1, 'aeroact%mact1'//chcall, i_dgtime )
        call save_dg(k, i_here, aeroact(k)%rcrit1, 'aeroact%rcrit1'//chcall, i_dgtime )
        call save_dg(k, i_here, aeroact(k)%mact1_mean, 'aeroact%mact1_mean'//chcall, i_dgtime )
        call save_dg(k, i_here, aeroact(k)%nact2, 'aeroact%nact2'//chcall, i_dgtime )
        call save_dg(k, i_here, aeroact(k)%rcrit2, 'aeroact%rcrit2'//chcall, i_dgtime )
        call save_dg(k, i_here, aeroact(k)%mact2, 'aeroact%mact2'//chcall, i_dgtime )
        call save_dg(k, i_here, aeroact(k)%mact2_mean, 'aeroact%mact2_mean'//chcall, i_dgtime )
      end if
    end do
#endif
  end subroutine examine_aerosol

  ! Note: despite appearences this is only coded up
  ! for a single mode so far...
  subroutine AbdulRazzakGhan(w, qc, p, T, phys, chem, nccn, Smax, rcrit, tau)
    real(wp), intent(in) :: w ! vertical velocity (ms-1)
    real(wp), intent(in) :: qc ! Pre-existing cloud mass
    real(wp), intent(in) :: p ! pressure (Pa)
    real(wp), intent(in) :: T ! temperature (K)
    type(aerosol_phys), intent(in) :: phys
    type(aerosol_chem), intent(in) :: chem
    real(wp), intent(out) :: Nccn ! number of activated aerosol
    real(wp), intent(out) :: smax ! peak supersaturation
    real(wp), intent(out) :: rcrit ! radius of smallest activated aerosol
    real(wp), intent(out) :: tau   ! equilibrium adjustment timescale

    real(wp), allocatable :: s_cr(:)
    real(wp) :: Tc, Ak, Bk, eta, alpha, gamma, es, bigG, Galt
    real(wp) :: f1, f2, zeta, error_func
    integer :: i

    allocate(s_cr(phys%nmodes))

    Tc=T-273.15 ! Temperature (C)

    ! surface tension of water, not solute (units are N/m)
    !zetasa = (76.1-(0.155*Tc)) * 1e-3
    Ak=2.0*Mw*zetasa/(Ru*T*rhow)

    alpha=g*(Lv/(eps*cp*T)-1)/(T*Rd)

    es=(100.0*6.1121)*exp((18.678-Tc/(234.5))*Tc/(257.14+Tc))

    gamma=eps*p/es+Lv**2/(Rv*T**2*cp)

    !Dv = (0.211*((T/273.15)**(1.94))*((100000.)/(p))) / 1.e4

    bigG=1.0/(rhow*(Rv*T/(es*Dv)+Lv*(Lv/(Rv*T)-1)/(ka*T)))

    do i=1, 1!phys%nmodes
      Bk=chem%vantHoff(i)*Mw*chem%density(i)/(chem%massMole(i)*rhow)
      s_cr(i)=(2.0/sqrt(Bk))*(Ak/(3.0*phys%rd(i)))**1.5
      eta=(w*alpha/BigG)**1.5/(2.0*pi*rhow*gamma*phys%N(i))
      zeta=(2.0/3.0)*Ak*(w*alpha/BigG)**0.5
      f1=1.5*exp(2.25*(log(phys%sigma(i)))**2)
      f2=1.0 + 0.25*log(phys%sigma(i))

      error_func=1.0-erf(log(f1*(zeta/eta)**1.5+&
           f2*((s_cr(i)*s_cr(i))/(eta+3.0*zeta))**.75)/(3.0*sqrt(2.0)*log(phys%sigma(i))))

      nccn=0.5*phys%N(i)*error_func
      smax=s_cr(i)/(f1*(zeta/eta)**1.5+f2*(s_cr(i)*s_cr(i)/eta)**.75)**.5
      rcrit=phys%rd(i)*(s_cr(i)/smax)**(.66667)
    end do
    tau=1.0/(alpha*w+gamma*4*pi*rhow*bigG*qc)
    deallocate(s_cr)
  end subroutine AbdulRazzakGhan

  ! Note: despite appearences this is only coded up
  ! for a single mode so far...
  subroutine AbdulRazzakGhan2000(w, p, T, phys, chem, nccn, Smax, active_phys, nccn_active, l_useactive )
    real(wp), intent(in) :: w ! vertical velocity (ms-1)
    real(wp), intent(in) :: p ! pressure (Pa)
    real(wp), intent(in) :: T ! temperature (K)
    type(aerosol_phys), intent(in) :: phys
    type(aerosol_chem), intent(in) :: chem
    type(aerosol_active), intent(in), optional :: active_phys
    real(wp), intent(out) :: nccn_active ! notional nccn that would be generated by activated aerosol
    real(wp), intent(out) :: Nccn(:) ! number of activated aerosol for each mode
    real(wp), intent(out) :: smax ! peak supersaturation
    logical, intent(out) :: l_useactive ! Do we use consider already activated aerosol

    real(wp), allocatable :: s_cr(:)
    real(wp) :: s_cr_active
    real(wp) :: Tc, Ak, Bk, eta, alpha, gamma, es, bigG, Galt
    real(wp) :: f1, f2, zeta, error_func
    real(wp) :: rsmax2  ! 1/(smax*smax)
    real(wp) :: diff
    integer :: i
    real(wp) :: sum

    l_useactive=.false.
    if (present(active_phys)) l_useactive=active_phys%nact > ccn_tidy .and. l_active_inarg2000

    nccn_active=0.0

    allocate(s_cr(phys%nmodes))

    Tc=T-273.15 ! Temperature (C)

    ! surface tension of water, not solute (units are N/m)
    !zetasa = (76.1-(0.155*TC)) * 1e-3

    Ak=2.0*Mw*zetasa/(Ru*T*rhow)
    alpha=g*(Lv/(eps*cp*T)-1)/(T*Rd)
    es=(100.0*6.1121)*exp((18.678-Tc/(234.5))*Tc/(257.14+Tc))
    gamma=eps*p/es+Lv**2/(Rv*T**2*cp)
    !Dv = (0.211*((T/273.15)**(1.94))*((100000.)/(p))) / 1.e4
    bigG=1.0/(rhow*(Rv*T/(es*Dv)+Lv*(Lv/(Rv*T)-1)/(ka*T)))
    zeta=(2.0/3.0)*Ak*(w*alpha/BigG)**0.5
    rsmax2=0.0

    do i=1, phys%nmodes
      if (phys%N(i) > ccn_tidy) then
        Bk=chem%vantHoff(i)*Mw*chem%density(i)/(chem%massMole(i)*rhow)
        s_cr(i)=(2.0/sqrt(Bk))*(Ak/(3.0*phys%rd(i)))**1.5
        eta=(w*alpha/BigG)**1.5/(2.0*pi*rhow*gamma*phys%N(i))
        f1=0.5*exp(2.5*(log(phys%sigma(i)))**2)
        f2=1.0+0.25*log(phys%sigma(i))
        rsmax2=rsmax2+(f1*(zeta/eta)**1.5+f2*(s_cr(i)*s_cr(i)/(eta+3.0*zeta))**.75)/(s_cr(i)*s_cr(i))
      end if
    end do

    if (l_useactive) then  ! We should include all the activated aerosol
      if (active_phys%nact > ccn_tidy) then
        i=aero_index%i_accum ! use accumulation mode chem (should do this properly)
        Bk=chem%vantHoff(i)*Mw*chem%density(i)/(chem%massMole(i)*rhow)
        s_cr_active=(2.0/sqrt(Bk))*(Ak/(3.0*active_phys%rd))**1.5
        eta=(w*alpha/BigG)**1.5/(2.0*pi*rhow*gamma*active_phys%nact)
        f1=0.5*exp(2.5*(log(active_phys%sigma))**2)
        f2=1.0+0.25*log(active_phys%sigma)
        rsmax2=rsmax2+(f1*(zeta/eta)**1.5+f2*(s_cr_active*s_cr_active/(eta+3.0*zeta))**.75)/(s_cr_active*s_cr_active)
      end if
    end if
    smax=sqrt(1.0/rsmax2)
    nccn=0.0

    if (l_useactive) then  ! We should include all the activated aerosol
      if (active_phys%nact > ccn_tidy) then
        error_func=1.0-erf(2.0*log(s_cr_active/smax)/(3.0*sqrt(2.0)*log(active_phys%sigma)))
        nccn_active=0.5*active_phys%nact*error_func
      end if
    end if

    do i=1, phys%nmodes
      if (phys%N(i) > ccn_tidy) then
        error_func=1.0-erf(2.0*log(s_cr(i)/smax)/(3.0*sqrt(2.0)*log(phys%sigma(i))))
        nccn(i)=0.5*phys%N(i)*error_func

        ! Make sure we don't activate too many...
        nccn(i)=min(nccn(i), .999*phys%N(i))
        ! And don't bother if we don't do much
        if (nccn(i) < ccn_tidy) nccn(i)=0.0
      end if
    end do

    if (l_useactive) then  ! Remove already activated aerosol
      ! from newly activated
      ! Start from smallest mode
      diff =active_phys%nact-nccn_active
      sum=0.0
      do i=1, size(nccn)
        sum=nccn(i)+sum
        nccn(i)=min(nccn(i), max(0.0_wp, sum - diff))
      end do
      !nccn(1) = max(0.0, nccn(1) - diff)
      !nccn(2) = min(nccn(2), max(0.0, nccn(2) + nccn(1) - diff))
      !nccn(3) = min(nccn(3), max(0.0, nccn(3) + nccn(2) + nccn(1) - diff))
    end if
    nccn=.99*nccn ! Don't allow all aerosol to be removed.
    deallocate(s_cr)
  end subroutine AbdulRazzakGhan2000

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
  function moment_logn(N, rm, sigma, p)
    real(wp), intent(in) :: N, rm, sigma
    real(wp), intent(in) :: p ! calculate pth moment
    real(wp) :: moment_logn

    moment_logn=N*rm**p*exp(.5*p*p*log(sigma)**2)
  end function moment_logn

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
  function upperpartial_moment_logn(N, rm, sigma, p, rcrit)
    real(wp), intent(in) :: N, rm, sigma
    real(wp), intent(in) :: p ! calculate pth moment
    real(wp), intent(in) :: rcrit ! lower threshold for partial moment
    real(wp) :: upperpartial_moment_logn

    if (rcrit==0.0) then
      upperpartial_moment_logn=moment_logn(N, rm, sigma, p)
    else
      upperpartial_moment_logn=N*rm**p*exp(.5*p*p*log(sigma)**2)     &
           * .5*erfc((log(rcrit/rm)/log(sigma) - p*log(sigma))/sqrt(2.0))
    end if
  end function upperpartial_moment_logn

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
  function invert_partial_moment(m, mup, p, rm, sigma)
    real(wp), intent(in) :: m, mup
    real(wp), intent(in) :: p ! pth moment
    real(wp), intent(in) :: rm, sigma
    real(wp) :: x, invert_partial_moment

    x=rm*exp(p*log(sigma)**2+sqrt(2.0)*log(sigma)*erfinv(1.0-2.0*mup/m))
    invert_partial_moment=x
  end function invert_partial_moment

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
  function invert_partial_moment_approx(mup, p, rm, sigma)
    real(wp), intent(in) :: mup, rm, sigma, p ! pth moment
    real(wp) :: invert_partial_moment_approx

    real(wp) :: beta, c, lsig, mbeta
    real(wp) :: x

    c=-0.1021 ! constant in approximation
    lsig=sqrt(2.0)*log(sigma)

    beta=mup*exp(-0.25*p*p*lsig*lsig)/rm**p
    mbeta=1.0-beta

    if (beta < 1.0e-4) then
      x=rm *(-c)**(-lsig)
    else if (beta >.9999) then
      x=0.0
    else
      x=rm * (sigma**(-0.5*p)*(((mbeta)*c+sqrt((mbeta)**2*c*c + 4.0*beta*(mbeta)))/(2.0*beta)))**(lsig)
    end if
    invert_partial_moment_approx=x
  end function invert_partial_moment_approx

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
  function invert_partial_moment_betterapprox(mup, p, N, rm, sigma)
    real(wp), intent(in) :: mup, rm, sigma, N, p ! pth moment
    real(wp) :: invert_partial_moment_betterapprox

    real(wp) :: frac, x
    real(wp) :: small_frac=1e-6

    if (mup==0.0) then
      frac=epsilon(x)
    else
      frac = mup/(N*rm**p)*exp(-.5*p*p*log(sigma)*log(sigma))
    end if

    if (frac>=1.0 - epsilon(x)) frac=.999
    if (frac > small_frac) then
      x=normal_quantile(frac)
      invert_partial_moment_betterapprox=rm*exp(x*log(sigma)+p*log(sigma)*log(sigma))
    else
      invert_partial_moment_betterapprox=0.0
    end if
  end function invert_partial_moment_betterapprox

  ! This uses an approximation 26.2.23 from A&S
  function normal_quantile(p)
    real(wp), intent(in) :: p
    real(wp) :: normal_quantile

    real(wp) :: pp, t,x
    real(wp) :: c0,c1,c2,d1,d2,d3

    c0=2.515517
    c1=0.802853
    c2=0.010328
    d1=1.432788
    d2=0.189269
    d3=0.001308

    pp=p

    if (p<=epsilon(x)) then
      x=huge(x)
    else if (p >=1.0 - epsilon(x)) then
      x=-huge(x)
    else
      if (p>0.5) pp=1.0-p
      t=sqrt(log(1.0/(pp*pp)))
      x=t-(c0+c1*t+c2*t*t)/(1.0+d1*t+d2*t*t+d3*t*t*t)
      if (p>0.5)x=-x
    end if
    normal_quantile=x
  end function normal_quantile

  !
  ! Calculate mean radius of lognormal distribution
  ! given Mass and number
  !
  function MNtoRm(M, N, density, sigma)    
    real(wp), intent(in) :: M, N, density, sigma
    real(wp) :: MNtoRm

    if (N==0 .or. M==0) then ! shouldn't really be here
      MNtoRm=0.0
    else
      MNtoRm=(3.0*M*exp(-4.5*log(sigma)**2)/(4.0*N*pi*density))**(1.0/3.0)
    end if
  end function MNtoRm
end module aerosol_routines
