module sedimentation
  use mphys_die, only: throw_mphys_error
  use variable_precision, only: wp, iwp, defp
  use mphys_parameters, only: nz, hydro_params, cloud_params, rain_params   &
       , ice_params, snow_params, graupel_params
  use passive_fields, only: rho, rdz_on_rho, dz, z_half, z_centre
  use type_process, only: process_name
  use mphys_switches, only: hydro_complexity, l_abelshipway, l_sed_3mdiff, l_cons, &
       i_an2, i_am2, i_am4, l_ased, i_am5, i_am7, i_am8, i_am9, i_ql, l_process, l_passivenumbers, l_passivenumbers_ice,   &
       active_cloud, active_rain, isol, iinsol, active_ice, active_number, i_an11, i_an12, &
       l_separate_rain, l_warm, i_aerosed_method
  use mphys_constants, only: rhow, rho0, cp
  use process_routines, only: process_rate, i_psedr, i_asedr, i_asedl, i_psedl, &
       i_pseds, i_psedi, i_psedg, i_dsedi, i_dseds, i_dsedg
  use thresholds, only: ql_small, qr_small, nr_small, m3r_small, thresh_small
  use special, only: Gammafunc

  use lookup, only: moment
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use aerosol_routines, only: aerosol_active

#if DEF_MODEL==MODEL_KiD
  use diagnostics, only: save_dg, i_dgtime, k_here, i_here, nx
  use runtime, only: l_dgstep, time
#elif DEF_MODEL==MODEL_LEM_DIAG
  use diaghelp_lem, only: i_here, j_here, k_here
  use com_params, only: time
  use extra_dgs
  use com_dgstore, only: l_dodgs
#elif DEF_MODEL==MODEL_UM
  use diaghelp_um, only: i_here, j_here, l_debug_um, debug_i, debug_j, debug_pe
  use UM_ParCore, only: mype
  use timestep_mod, only: timestep_number
#elif  DEF_MODEL==MODEL_MONC
  use diaghelp_monc, only: i_here, j_here
#endif
  implicit none
  private

  character*(2) :: qchar

  public sedr
contains

  subroutine sedr(step_length, qfields, aerofields, aeroact, dustact,   &
       tend, params, procs, aerosol_procs, precip, l_doaerosol)
    real(wp), intent(in) :: step_length
    real(wp), intent(in), target :: qfields(:,:), aerofields(:,:)
    real(wp), intent(in) :: tend(:,:)
    type(hydro_params), intent(in) :: params
    type(aerosol_active), intent(in) :: aeroact(:), dustact(:)
    ! NB for liquid phase (cloud/rain) dustact is dustliq
    ! for ice phase (ice/snow/graupel) aeroact is iceact
    type(process_rate), intent(inout), target :: procs(:,:)
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)
    real(wp), intent(out) :: precip
    logical, optional, intent(in) :: l_doaerosol

    call sedr_aero(step_length, qfields, aerofields, aeroact, dustact,       &
         tend, params, procs, aerosol_procs, precip, l_doaerosol)
  end subroutine sedr

  subroutine sedr_aero(step_length, qfields, aerofields, aeroact, dustact,   &
       tend, params, procs, aerosol_procs, precip, l_doaerosol)
    real(wp), intent(in) :: step_length
    real(wp), intent(in), target :: qfields(:,:), aerofields(:,:)
    real(wp), intent(in) :: tend(:,:)
    type(hydro_params), intent(in) :: params
    type(aerosol_active), intent(in) :: aeroact(:), dustact(:)
    type(process_rate), intent(inout), target :: procs(:,:)
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)
    real(wp), intent(out) :: precip
    logical, optional, intent(in) :: l_doaerosol

    real(wp) :: dm1, dm2, dm3
    real(wp) :: dn1, dn2, dn3, am2
    real(wp) :: m1, m2, m3
    real(wp) :: n1, n2, n3
    real(wp) :: hydro_mass
    type(process_rate), pointer :: this_proc
    type(process_rate), pointer :: aero_proc
    real(wp), allocatable :: flux_n1(:)
    real(wp), allocatable :: flux_n2(:)
    real(wp), allocatable :: flux_n3(:)
    real(wp), allocatable :: Grho(:)
    real(wp) :: n0, lam, mu, u1r, u2r, u3r, u1r2, u2r2, u3r2
    real(wp) :: udp ! bulk fallspeed for mass (for precip diag)
    integer :: k

    integer :: i_1m, i_2m, i_3m
    real(wp) :: p1, p2, p3
    real(wp) :: sp1, sp2, sp3
    real(wp) :: a_x, b_x, f_x, g_x, c_x, d_x
    real(wp) :: a2_x, b2_x, f2_x
    logical :: l_2m, l_3m
    logical :: l_fluxin, l_fluxout

    type(process_name) :: iproc, iaproc  ! processes selected depending on
    ! which species we're depositing on.

    real(wp) :: debug, mint, dmint, scal
    real(wp) :: dmac, mac_mean_min=1.0e-20
    real(wp) :: dmad
    real(wp) :: dnumber_a, dnumber_d
    logical :: l_da_local  ! local tranfer of l_doaerosol
    real(wp) :: ratio ! used to scale back fluxes if the timestep is too small
    ! If this is used, you can't trust the results.

    l_da_local=.false.
    if (present(l_doaerosol)) l_da_local=l_doaerosol

    ! precip diag
    write(qchar, '(i2.2)') params%id
    precip=0.0

    m1=0.0
    m2=0.0
    m3=0.0

    allocate(flux_n1(nz))
    flux_n1=0.0
    i_1m=params%i_1m
    i_2m=params%i_2m
    i_3m=params%i_3m
    p1=params%p1
    p2=params%p2
    p3=params%p3

    if (l_sed_3mdiff) then
      sp1=params%sp1
      sp2=params%sp2
      sp3=params%sp3
    else
      sp1=params%p1
      sp2=params%p2
      sp3=params%p3
    end if

    ! we don't want flexible approach in this version....
    if (p1/=sp1 .or. p2/=sp2 .or. p3/=sp3) then
      print*, 'Cannot have flexible sedimentation options with CASIM aerosol'
      call throw_mphys_error(1, 'sedr_aero', 'Cannot have flexible sedimentation options with CASIM aerosol')
    end if

    a_x=params%a_x
    b_x=params%b_x
    f_x=params%f_x
    g_x=params%g_x
    c_x=params%c_x
    d_x=params%d_x
    a2_x=params%a2_x
    b2_x=params%b2_x
    f2_x=params%f2_x

    if (params%l_2m) then
      allocate(flux_n2(nz))
      flux_n2=0.0
    end if
    if (params%l_3m) then
      allocate(flux_n3(nz))
      flux_n3=0.0
    end if

    allocate(Grho(nz)) ! may want to move this somewhere else
    Grho=(rho0/rho)**g_x

    dmint=0.0
    select case (params%id)
    case (1_iwp) !cloud
      iproc=i_psedl
      iaproc=i_asedl
    case (2_iwp) !rain
      iproc=i_psedr
      iaproc=i_asedr
    case (3_iwp) !ice
      iproc=i_psedi
      iaproc=i_dsedi
    case (4_iwp) !snow
      iproc=i_pseds
      iaproc=i_dseds
    case (5_iwp) !graupel
      iproc=i_psedg
      iaproc=i_dsedg
    end select

    do k=nz-1, 1, -1
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

      this_proc=>procs(k, iproc%id)
      if (l_da_local) aero_proc=>aerosol_procs(k, iaproc%id)

      hydro_mass=qfields(k, params%i_1m)
      if (params%l_2m) m2=qfields(k, params%i_2m)
      if (params%l_3m) m3=qfields(k, params%i_3m)

      l_fluxin=.false.
      l_fluxout=.false.
      if (hydro_mass > thresh_small(params%i_1m)) l_fluxout=.true.
      if (qfields(k+1, params%i_1m) > thresh_small(params%i_1m)) l_fluxin=.true.

      if (l_fluxout) then
        m1=(hydro_mass/c_x)

        n0=dist_n0(k,params%id)
        mu=dist_mu(k,params%id)
        lam=dist_lambda(k,params%id)

        if (l_sed_3mdiff) then
          ! Moment transfer
          n1=moment(n0, lam, mu, sp1)
          n2=moment(n0, lam, mu, sp2)
          n3=moment(n0, lam, mu, sp3)
        else
          n1=m1
          n2=m2
          n3=m3
        end if

        u1r=a_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f_x)**(-(1.0+mu+sp1+b_x)))     &
             *(Gammafunc(1.0+mu+sp1+b_x)/Gammafunc(1.0+mu+sp1))

        if (params%l_2m)     &
             u2r=a_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f_x)**(-(1.0+mu+sp2+b_x))) &
             *(Gammafunc(1.0+mu+sp2+b_x)/Gammafunc(1.0+mu+sp2))

        if (params%l_3m)     &
             u3r=a_x*Grho(k)*(lam**(1.0+mu+sp3)*(lam+f_x)**(-(1.0+mu+sp3+b_x))) &
             *(Gammafunc(1.0+mu+sp3+b_x)/Gammafunc(1.0+mu+sp3))

        if (l_abelshipway .and. params%id==rain_params%id) then ! rain can use abel and shipway formulation
          u1r2=a2_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f2_x)**(-(1.0+mu+sp1+b2_x)))   &
               *(Gammafunc(1.0+mu+sp1+b2_x)/Gammafunc(1.0+mu+sp1))
          u1r=u1r+u1r2

          if (params%l_2m) then
            u2r2=a2_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f2_x)**(-(1.0+mu+sp2+b2_x))) &
                 *(Gammafunc(1.0+mu+sp2+b2_x)/Gammafunc(1.0+mu+sp2))
            u2r=u2r+u2r2
          end if

          if (params%l_3m) then
            u3r2=a2_x*Grho(k)*(lam**(1.0+mu+sp3)*(lam+f2_x)**(-(1.0+mu+sp3+b2_x))) &
                 *(Gammafunc(1.0+mu+sp3+b2_x)/Gammafunc(1.0+mu+sp3))
            u3r=u3r+u3r2
          end if
        end if

        ! fall speeds shouldn't get too big...
        u1r=min(u1r,params%maxv)

        ! For clouds and ice, we only use a 1M representation of sedimentation
        ! so that spurious size sorting doesn't lead to overactive autoconversion
        if (params%id==cloud_params%id .or. ice_params%id==cloud_params%id) then
          if (params%l_2m)u2r=u1r
        else
          if (params%l_2m)u2r=min(u2r,params%maxv)
          if (params%l_3m)u3r=min(u3r,params%maxv)
        end if

        ! fall speeds shouldn't be negative (can happen with original AS formulation)
        u1r=max(u1r,0.0_wp)
        if (params%l_2m)u2r=max(u2r,0.0_wp)
        if (params%l_3m)u3r=max(u3r,0.0_wp)

#if DEF_MODEL==MODEL_KiD
        if (nx==1) then
          call save_dg(k, u1r, 'u1r_'//qchar, i_dgtime)
          if (params%l_2m)call save_dg(k, u2r, 'u2r_'//qchar, i_dgtime)
          if (params%l_3m)call save_dg(k, u3r, 'u3r_'//qchar, i_dgtime)
        else
          call save_dg(k, i_here, u1r, 'u1r_'//qchar, i_dgtime)
          if (params%l_2m)call save_dg(k, i_here, u2r, 'u2r_'//qchar, i_dgtime)
          if (params%l_3m)call save_dg(k, i_here, u3r, 'u3r_'//qchar, i_dgtime)
        end if
#elif DEF_MODEL==MODEL_LEM_DIAG
        if (l_dodgs .and. j_here>=jstartdg .and.     &
             j_here<=jenddg) then
          if (params%id == rain_params%id) then
            idgproc = req_dgproc('VR')
          else if (params%id == ice_params%id) then
            idgproc = req_dgproc('VI')
          else if (params%id == snow_params%id) then
            idgproc = req_dgproc('VS')
          else if (params%id == graupel_params%id) then
            idgproc = req_dgproc('VG')
          end if
          if (idgproc > 0) then
            dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = u1r
          end if
          if (params%id == rain_params%id) then
            idgproc = req_dgproc('VRN')
          else if (params%id == ice_params%id) then
            idgproc = req_dgproc('VIN')
          else if (params%id == snow_params%id) then
            idgproc = req_dgproc('VSN')
          else if (params%id == graupel_params%id) then
            idgproc = req_dgproc('VGN')
          end if
          if (idgproc > 0) then
            dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = u2r
          end if
          if (params%id == rain_params%id) then
            idgproc = req_dgproc('LAMR')
          else if (params%id == ice_params%id) then
            idgproc = req_dgproc('LAMI')
          else if (params%id == snow_params%id) then
            idgproc = req_dgproc('LAMS')
          else if (params%id == graupel_params%id) then
            idgproc = req_dgproc('LAMG')
          end if
          if (idgproc > 0) then
            dgprocs(j_here-jstartdg+1,k,i_here,idgproc) = lam
          end if
        end if
#endif

        flux_n1(k)=n1*u1r

        if (params%l_2m)     &
             flux_n2(k)=n2*u2r

        if (params%l_3m)     &
             flux_n3(k)=n3*u3r

#if DEF_MODEL==MODEL_UM
        !===========================================================
        ! Cludge to slow down the rain if CFL condition is broken
        ! Should be replaced with a semi-lagrangian scheme
        !===========================================================
        if (k>1) then
          if (c_x*flux_n1(k)*rdz_on_rho(k)*step_length > 0.99*hydro_mass .or.      &
               flux_n2(k)*rdz_on_rho(k)*step_length > 0.99*m2) then

            ratio=min(0.95*hydro_mass/c_x/flux_n1(k)/rdz_on_rho(k)/step_length, 0.95*m2/flux_n2(k)/rdz_on_rho(k)/step_length)

            print*, 'WARNING: slowing down precip - your timestep is too long'
            print*, 'INFO:', timestep_number, i_here, j_here, k, step_length,      &
                 ratio, qfields(k, params%i_1m), c_x*flux_n1(k)*rdz_on_rho(k)&
                 , qfields(k, params%i_2m), flux_n2(k)*rdz_on_rho(k), u1r, u2r

            flux_n1(k)=flux_n1(k)*ratio

            if (params%l_2m) flux_n2(k)=flux_n2(k)*ratio

            if (params%l_3m) flux_n3(k)=flux_n3(k)*ratio
          end if
        end if
#endif

        ! diagnostic for precip
        if (k==1) then
          udp=a_x*Grho(k)*(lam**d_x*(lam+f_x)**(-(d_x+b_x)))*(Gammafunc(1.0+mu+d_x+b_x)/Gammafunc(1.0+mu+d_x))
          precip=(rho(k)*hydro_mass)*udp
        end if
      end if

      dmac=0.0
      dmad=0.0
      dnumber_a=0.0
      dnumber_d=0.0

      if (l_fluxout) then !flux out (flux(k+1) will be zero if no flux in)
        dn1=(flux_n1(k+1)-flux_n1(k))*rdz_on_rho(k)

        if (params%l_2m) dn2=(flux_n2(k+1)-flux_n2(k))*rdz_on_rho(k)

        if (params%l_3m) dn3=(flux_n3(k+1)-flux_n3(k))*rdz_on_rho(k)

        !============================
        ! aerosol processing
        !============================
        if (l_ased .and. l_da_local) then
          if (params%id == cloud_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean -    &
                 flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)* rdz_on_rho(k)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1 - flux_n2(k)*aeroact(k)%nratio1)* rdz_on_rho(k)
            end if
            if (.not. l_warm) then
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-  &
                   flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k)
              if (l_passivenumbers_ice .and. dustact(k)%mact_mean > 0.0) then
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k)
              end if
            end if
          else if (params%id == rain_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-    &
                 flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k)
            if (l_passivenumbers .and. aeroact(k)%mact_mean > 0.0) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k)
            end if
            if (.not. l_warm) then
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-  &
                   flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k)
              if (l_passivenumbers_ice) then
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k)
              end if
            end if
          end if

          if (params%id == ice_params%id) then
            dmac = (flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean-    &
                 flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)*rdz_on_rho(k)
            dmad = (flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-    &
                 flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1-flux_n2(k)*aeroact(k)%nratio1)*rdz_on_rho(k)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k)
            end if
          else if (params%id == snow_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-    &
                 flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-    &
                 flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k)
            end if
          else if (params%id == graupel_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean-    &
                 flux_n2(k)*aeroact(k)%nratio3*aeroact(k)%mact3_mean)*rdz_on_rho(k)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean-    &
                 flux_n2(k)*dustact(k)%nratio3*dustact(k)%mact3_mean)*rdz_on_rho(k)
#if DEF_MODEL==MODEL_KiD
            call save_dg(k, dmad, 'dmac_gs', i_dgtime)
            call save_dg(k, flux_n2(k+1)*dustact(k+1)%mact3_mean, 'flux_in', i_dgtime)
            call save_dg(k, flux_n2(k)*dustact(k)%mact3_mean, 'flux_out', i_dgtime)
            call save_dg(k, n2*dustact(k)%mact3_mean, 'n2Xmean3', i_dgtime)
            call save_dg(k, dustact(k)%mact3, 'mact', i_dgtime)
#endif
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3-flux_n2(k)*aeroact(k)%nratio3)*rdz_on_rho(k)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3-flux_n2(k)*dustact(k)%nratio3)*rdz_on_rho(k)
            end if
          end if
        end if
      else if (l_fluxin) then !flux in, but not out
        dn1=flux_n1(k+1)*rdz_on_rho(k)
        if (params%l_2m) dn2=flux_n2(k+1)*rdz_on_rho(k)
        if (params%l_3m) dn3=flux_n3(k+1)*rdz_on_rho(k)

        !============================
        ! aerosol processing
        !============================
        if (l_ased .and. l_da_local) then
          if (params%id == cloud_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k)
            end if
            if (.not. l_warm) then
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k)
              if (l_passivenumbers_ice) then
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k)
              end if
            end if
          else if (params%id == rain_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k)
            end if
            if (.not. l_warm) then
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k)
              if (l_passivenumbers_ice) then
                dnumber_d=flux_n2(k+1)*dustact(k+1)%nratio2*rdz_on_rho(k)
              end if
            end if
          end if

          if (params%id == ice_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k)
            end if
          else if (params%id == snow_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2)*rdz_on_rho(k)
            end if
          else if (params%id == graupel_params%id) then
            if (i_aerosed_method==1) then
              dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean)*rdz_on_rho(k)
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean)*rdz_on_rho(k)
              if (l_passivenumbers) then
                dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3)*rdz_on_rho(k)
              end if
              if (l_passivenumbers_ice) then
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3)*rdz_on_rho(k)
              end if
            else
              print*, 'ERROR: GET RID OF i_aerosed_method variable!'
            end if
          end if
        end if
      end if

      ! Store the aerosol process terms...
      if (l_da_local) then
        if (params%id == cloud_params%id .or. params%id == rain_params%id) then
          !liquid phase
          if (l_separate_rain .and. params%id == rain_params%id) then
            aero_proc%source(i_am5)=dmac
          else
            aero_proc%source(i_am4)=dmac
          end if
          if (.not. l_warm) aero_proc%source(i_am9) =dmad
          if (l_passivenumbers) then
            aero_proc%source(i_an11)=dnumber_a
          end if
          if (l_passivenumbers_ice) then
            aero_proc%source(i_an12)=dnumber_d
          end if

        else
          !ice phase
          aero_proc%source(i_am7)=dmad
          aero_proc%source(i_am8)=dmac
          if (l_passivenumbers) then
            aero_proc%source(i_an11)=dnumber_a
          end if
          if (l_passivenumbers_ice) then
            aero_proc%source(i_an12)=dnumber_d
          end if
        end if
      end if

#if DEF_MODEL==MODEL_KiD

      if (nx==1) then
        call save_dg(k, dmac, 'dmac', i_dgtime)
        call save_dg(k, (c_x*flux_n1(k+1)*9.8/cp), 'frictional_heating'//qchar, i_dgtime)
      else
        call save_dg(k, i_here, dmac, 'dmac', i_dgtime)
        call save_dg(k, i_here, (c_x*flux_n1(k+1)*9.8/cp), 'frictional_heating'//qchar, i_dgtime)
      end if
#endif

      dm1=dn1
      dm2=dn2
      dm3=dn3

      this_proc%source(params%i_1m)=c_x*dm1

      if (params%l_2m) this_proc%source(params%i_2m)=dm2
      if (params%l_3m) this_proc%source(params%i_3m) = dm3

      nullify(this_proc)
      if (l_da_local) nullify(aero_proc)
    end do

    deallocate(Grho)
    deallocate(flux_n1)
    if (params%l_2m)deallocate(flux_n2)
    if (params%l_3m)deallocate(flux_n3)
  end subroutine sedr_aero
end module sedimentation
