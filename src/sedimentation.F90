module sedimentation
  use mphys_die, only: throw_mphys_error, incorrect_opt, std_msg
  use variable_precision, only: wp, iwp, defp
  use mphys_parameters, only: nz, hydro_params, cloud_params, rain_params   &
       , ice_params, snow_params, graupel_params
  use passive_fields, only: rho, rdz_on_rho, dz, z_half, z_centre
  use type_process, only: process_name
  use mphys_switches, only: hydro_complexity, l_abelshipway, l_sed_3mdiff, l_cons, &
       i_an2, i_am2, i_am4, l_ased, i_am5, i_am7, i_am8, i_am9, i_ql, l_process, l_passivenumbers, l_passivenumbers_ice,   &
       active_cloud, active_rain, isol, iinsol, active_ice, active_number, i_an11, i_an12, &
       l_separate_rain, l_warm, i_aerosed_method, l_sed_icecloud_as_1m
  use mphys_constants, only: rhow, rho0, cp
  use process_routines, only: process_rate, i_psedr, i_asedr, i_asedl, i_psedl, &
       i_pseds, i_psedi, i_psedg, i_dsedi, i_dseds, i_dsedg
  use thresholds, only: ql_small, qr_small, nr_small, m3r_small, thresh_small
  use special, only: Gammafunc

  use lookup, only: moment
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use aerosol_routines, only: aerosol_active

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='SEDIMENTATION'

  character*(2) :: qchar
  real(wp), allocatable :: flux_n1(:)
  real(wp), allocatable :: flux_n2(:)
  real(wp), allocatable :: flux_n3(:)
  real(wp), allocatable :: Grho(:)

  public sedr, initialise_sedr, finalise_sedr
contains

  subroutine initialise_sedr()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='INITIALISE_SEDR'


    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    allocate(flux_n1(nz), flux_n2(nz), flux_n3(nz), Grho(nz))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine initialise_sedr

  subroutine finalise_sedr()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='FINALISE_SEDR'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    deallocate(flux_n1, flux_n2, flux_n3, Grho)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine finalise_sedr  

  subroutine sedr(step_length, qfields, aerofields, aeroact, dustact,   &
       tend, params, procs, aerosol_procs, precip, l_doaerosol)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SEDR'

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

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    call sedr_aero(step_length, qfields, aerofields, aeroact, dustact,       &
         tend, params, procs, aerosol_procs, precip, l_doaerosol)


    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine sedr

  subroutine sedr_aero(step_length, qfields, aerofields, aeroact, dustact,   &
       tend, params, procs, aerosol_procs, precip, l_doaerosol)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SEDR_AERO'

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

    real(wp) :: mint, dmint, scal
    real(wp) :: dmac, mac_mean_min=1.0e-20
    real(wp) :: dmad
    real(wp) :: dnumber_a, dnumber_d
    logical :: l_da_local  ! local tranfer of l_doaerosol
    real(wp) :: ratio ! used to scale back fluxes if the timestep is too small
    ! If this is used, you can't trust the results.

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_da_local=.false.
    if (present(l_doaerosol)) l_da_local=l_doaerosol

    ! precip diag
    write(qchar, '(i2.2)') params%id
    precip=0.0

    m1=0.0
    m2=0.0
    m3=0.0
    
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

    Grho=(rho0/rho)**params%g_x

    ! we don't want flexible approach in this version....
    if (p1/=sp1 .or. p2/=sp2 .or. p3/=sp3) then
      write(std_msg, '(A)') 'Cannot have flexible sedimentation options '//&
                            'with CASIM aerosol'
      call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                             std_msg)
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

    if (params%l_2m) flux_n2=0.0     
    if (params%l_3m) flux_n3=0.0      

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
          n1=m1*rho(k)
          n2=m2*rho(k)
          n3=m3*rho(k)
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
          if (l_sed_icecloud_as_1m)then
            if (params%l_2m)u2r=u1r
          else
            if (params%l_2m)u2r=min(u2r,params%maxv)
          end if
        else
          if (params%l_2m)u2r=min(u2r,params%maxv)
          if (params%l_3m)u3r=min(u3r,params%maxv)
        end if

        ! fall speeds shouldn't be negative (can happen with original AS formulation)
        u1r=max(u1r,0.0_wp)
        if (params%l_2m)u2r=max(u2r,0.0_wp)
        if (params%l_3m)u3r=max(u3r,0.0_wp)

        flux_n1(k)=n1*u1r

        if (params%l_2m) flux_n2(k)=n2*u2r

        if (params%l_3m) flux_n3(k)=n3*u3r

        ! diagnostic for precip
        if (k==1) then
          precip=flux_n1(k)*c_x
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
              write(std_msg, '(A)') 'ERROR: GET RID OF i_aerosed_method variable!'
              call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                                     std_msg )
            end if
          end if
        end if
      end if

      ! Store the aerosol process terms...
      if (l_da_local) then
        if (params%id == cloud_params%id .or. params%id == rain_params%id) then
          !liquid phase
          if (l_separate_rain .and. params%id == rain_params%id) then
            aerosol_procs(k, iaproc%id)%source(i_am5)=dmac
          else
            aerosol_procs(k, iaproc%id)%source(i_am4)=dmac
          end if
          if (.not. l_warm) aerosol_procs(k, iaproc%id)%source(i_am9) =dmad
          if (l_passivenumbers) then
            aerosol_procs(k, iaproc%id)%source(i_an11)=dnumber_a
          end if
          if (l_passivenumbers_ice) then
            aerosol_procs(k, iaproc%id)%source(i_an12)=dnumber_d
          end if

        else
          !ice phase
          aerosol_procs(k, iaproc%id)%source(i_am7)=dmad
          aerosol_procs(k, iaproc%id)%source(i_am8)=dmac
          if (l_passivenumbers) then
            aerosol_procs(k, iaproc%id)%source(i_an11)=dnumber_a
          end if
          if (l_passivenumbers_ice) then
            aerosol_procs(k, iaproc%id)%source(i_an12)=dnumber_d
          end if
        end if
      end if

      dm1=dn1
      dm2=dn2
      dm3=dn3

      procs(k, iproc%id)%source(params%i_1m)=c_x*dm1

      if (params%l_2m) procs(k, iproc%id)%source(params%i_2m)=dm2
      if (params%l_3m) procs(k, iproc%id)%source(params%i_3m)=dm3
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine sedr_aero
end module sedimentation
