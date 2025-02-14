module sedimentation
  use mphys_die, only: throw_mphys_error, incorrect_opt, std_msg
  use variable_precision, only: wp, iwp
  use mphys_parameters, only: nz, hydro_params, cloud_params, rain_params   &
       , ice_params, snow_params, graupel_params, a_s, b_s
  use passive_fields, only: rho, rdz_on_rho, dz
  use type_process, only: process_name
  use mphys_switches, only: l_abelshipway, l_sed_3mdiff, &
       i_am4, l_ased, i_am5, i_am7, i_am8, i_am9, l_passivenumbers, l_passivenumbers_ice,   &
       i_an11, i_an12, &
       l_separate_rain, l_warm, i_aerosed_method, l_sed_icecloud_as_1m, l_sed_rain_1m, l_sed_snow_1m, &
       l_sed_graupel_1m, l_kfsm, & 
       cfl_vt_max, l_sed_eulexp
  use mphys_constants, only: rho0
  use process_routines, only: process_rate, i_psedr, i_asedr, i_asedl, i_psedl, &
       i_pseds, i_psedi, i_psedg, i_dsedi, i_dseds, i_dsedg
  use thresholds, only: thresh_small
  use special, only: Gammafunc

  use lookup, only: moment
  use distributions, only: dist_lambda, dist_mu, dist_n0, dist_lams
  use aerosol_routines, only: aerosol_active

  use lsp_sedim_eulexp_mod, only: lsp_sedim_eulexp
  ! lsp_sedim_eulexp_mod is a UM module, so can not be stored in CASIM or MONC repo. 
  ! if required please contact Adrian Hill. With MONC this uses a dummy routine in the CASIM
  ! component

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='SEDIMENTATION'

  real(wp), allocatable :: flux_n1(:)
  real(wp), allocatable :: flux_n2(:)
  real(wp), allocatable :: flux_n3(:)
  real(wp), allocatable :: Grho(:)

!$OMP THREADPRIVATE(flux_n1, flux_n2, flux_n3, Grho)

  public sedr, initialise_sedr, finalise_sedr, sedr_1M_2M, terminal_velocity_CFL
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

  subroutine sedr(ixy_inner, qfields, aeroact, dustact,   &
       params, procs, aerosol_procs, precip1d, l_doaerosol)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SEDR'

    integer, intent(in) :: ixy_inner

    real(wp), intent(in), target :: qfields(:,:)
    type(hydro_params), intent(in) :: params
    type(aerosol_active), intent(in) :: aeroact(:), dustact(:)
    type(process_rate), intent(inout), target :: procs(:,:)
    type(process_rate), intent(inout), target :: aerosol_procs(:,:)
    real(wp), intent(out) :: precip1d(nz)
    logical, optional, intent(in) :: l_doaerosol

    real(wp) :: dm1, dm2, dm3
    real(wp) :: dn1, dn2, dn3
    real(wp) :: m1, m2, m3
    real(wp) :: n1, n2, n3
    real(wp) :: hydro_mass
    real(wp) :: n0, lam, mu, u1r, u2r, u1r2, u2r2
!   real(wp) :: u3r, u3r2
    integer :: k

    real(wp) :: p1, p2, p3
    real(wp) :: sp1, sp2, sp3
    real(wp) :: a_x, b_x, f_x, c_x
    real(wp) :: a2_x, b2_x, f2_x
    logical :: l_fluxin, l_fluxout

    type(process_name) :: iproc, iaproc  ! processes selected depending on
    ! which species we're depositing on.

    real(wp) :: dmac
    real(wp) :: dmad
    real(wp) :: dnumber_a, dnumber_d
    logical :: l_sedim_generic
    logical :: l_da_local  ! local tranfer of l_doaerosol
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

    if (l_kfsm) then
      l_sedim_generic = (params%id==ice_params%id)
    else
      l_sedim_generic = .false.
    end if

    ! precip diag
    do k = 1, nz
      precip1d(k) = 0.0
      flux_n1(k)  = 0.0
      Grho(k)=(rho0/rho(k,ixy_inner))**params%g_x
    end do

    m1=0.0
    m2=0.0
    m3=0.0

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
      write(std_msg, '(A)') 'Cannot have flexible sedimentation options '//&
                            'with CASIM aerosol'
      call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                             std_msg)
    end if

    if (.not. l_kfsm) then
      a_x=params%a_x
      b_x=params%b_x
    end if

    f_x=params%f_x
    c_x=params%c_x
    a2_x=params%a2_x
    b2_x=params%b2_x
    f2_x=params%f2_x

    if (params%l_2m) flux_n2(:)=0.0     
!!$    if (params%l_3m) flux_n3=0.0      

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
!!$      if (params%l_3m) m3=qfields(k, params%i_3m)

      l_fluxin=.false.
      l_fluxout=.false.
      if (hydro_mass > thresh_small(params%i_1m)) l_fluxout=.true.
      if (qfields(k+1, params%i_1m) > thresh_small(params%i_1m)) l_fluxin=.true.

      if (l_sedim_generic) then
        a_x=a_s(k)
        b_x=b_s(k)
      else
        a_x=params%a_x
        b_x=params%b_x
      end if

      if (l_fluxout) then
        m1=(hydro_mass/c_x)

        n0=dist_n0(k,params%id)
        mu=dist_mu(k,params%id)

        if (l_sedim_generic) then
          lam=dist_lams(k,params%id,1)
        else
          lam=dist_lambda(k,params%id)
        end if

        if (l_sed_3mdiff) then
          ! Moment transfer
          n1=moment(n0, lam, mu, sp1)
          n2=moment(n0, lam, mu, sp2)
          n3=moment(n0, lam, mu, sp3)
        else
          n1=m1*rho(k,ixy_inner)
          n2=m2*rho(k,ixy_inner)
          n3=m3*rho(k,ixy_inner)
        end if

        u1r=a_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f_x)**(-(1.0+mu+sp1+b_x)))     &
             *(Gammafunc(1.0+mu+sp1+b_x)/Gammafunc(1.0+mu+sp1))

        if (params%l_2m)     &
             u2r=a_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f_x)**(-(1.0+mu+sp2+b_x))) &
             *(Gammafunc(1.0+mu+sp2+b_x)/Gammafunc(1.0+mu+sp2))

!!$        if (params%l_3m)     &
!!$             u3r=a_x*Grho(k)*(lam**(1.0+mu+sp3)*(lam+f_x)**(-(1.0+mu+sp3+b_x))) &
!!$             *(Gammafunc(1.0+mu+sp3+b_x)/Gammafunc(1.0+mu+sp3))

        if (l_abelshipway .and. params%id==rain_params%id) then ! rain can use abel and shipway formulation
          u1r2=a2_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f2_x)**(-(1.0+mu+sp1+b2_x)))   &
               *(Gammafunc(1.0+mu+sp1+b2_x)/Gammafunc(1.0+mu+sp1))
          u1r=u1r+u1r2

          if (params%l_2m) then
            u2r2=a2_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f2_x)**(-(1.0+mu+sp2+b2_x))) &
                 *(Gammafunc(1.0+mu+sp2+b2_x)/Gammafunc(1.0+mu+sp2))
            u2r=u2r+u2r2
          end if

!!$          if (params%l_3m) then
!!$            u3r2=a2_x*Grho(k)*(lam**(1.0+mu+sp3)*(lam+f2_x)**(-(1.0+mu+sp3+b2_x))) &
!!$                 *(Gammafunc(1.0+mu+sp3+b2_x)/Gammafunc(1.0+mu+sp3))
!!$            u3r=u3r+u3r2
!!$          end if
        end if

        ! fall speeds shouldn't get too big...
        u1r=min(u1r,params%maxv)

        ! For clouds and ice, we only use a 1M representation of sedimentation
        ! so that spurious size sorting doesn't lead to overactive autoconversion
        if (params%id==cloud_params%id .or. params%id==ice_params%id) then
          if (l_sed_icecloud_as_1m)then
            if (params%l_2m)u2r=u1r
          else
            if (params%l_2m)u2r=min(u2r,params%maxv)
          end if
        else
          if (params%l_2m)u2r=min(u2r,params%maxv)
!!$          if (params%l_3m)u3r=min(u3r,params%maxv)
        end if

        ! fall speeds shouldn't be negative (can happen with original AS formulation)
        u1r=max(u1r,0.0_wp)
        if (params%l_2m)u2r=max(u2r,0.0_wp)
!!$        if (params%l_3m)u3r=max(u3r,0.0_wp)

         flux_n1(k)=n1*u1r
           
         if (params%l_2m) flux_n2(k)=n2*u2r
        
        !if (params%l_3m) flux_n3(k)=n3*u3r

        precip1d(k) = flux_n1(k)*c_x
        
     end if

      dmac=0.0
      dmad=0.0
      dnumber_a=0.0
      dnumber_d=0.0

      if (l_fluxout) then !flux out (flux(k+1) will be zero if no flux in)
        dn1=(flux_n1(k+1)-flux_n1(k))*rdz_on_rho(k,ixy_inner)
        if (params%l_2m) dn2=(flux_n2(k+1)-flux_n2(k))*rdz_on_rho(k,ixy_inner)
        !if (params%l_3m) dn3=(flux_n3(k+1)-flux_n3(k))*rdz_on_rho(k)

        !============================
        ! aerosol processing
        !============================
        if (l_ased .and. l_da_local) then
          if (params%id == cloud_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean -    &
                 flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)* rdz_on_rho(k,ixy_inner)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1 - flux_n2(k)*aeroact(k)%nratio1)* rdz_on_rho(k,ixy_inner)
            end if
            if (.not. l_warm) then
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-  &
                   flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
              if (l_passivenumbers_ice .and. dustact(k)%mact_mean > 0.0) then
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
              end if
            end if
          else if (params%id == rain_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-    &
                 flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            if (l_passivenumbers .and. aeroact(k)%mact_mean > 0.0) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
            end if
            if (.not. l_warm) then
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-  &
                   flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
              if (l_passivenumbers_ice) then
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
              end if
            end if
          end if

          if (params%id == ice_params%id) then
            dmac = (flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean-    &
                 flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            dmad = (flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-    &
                 flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1-flux_n2(k)*aeroact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
            end if
          else if (params%id == snow_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-    &
                 flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-    &
                 flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
            end if
          else if (params%id == graupel_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean-    &
                 flux_n2(k)*aeroact(k)%nratio3*aeroact(k)%mact3_mean)*rdz_on_rho(k,ixy_inner)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean-    &
                 flux_n2(k)*dustact(k)%nratio3*dustact(k)%mact3_mean)*rdz_on_rho(k,ixy_inner)

            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3-flux_n2(k)*aeroact(k)%nratio3)*rdz_on_rho(k,ixy_inner)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3-flux_n2(k)*dustact(k)%nratio3)*rdz_on_rho(k,ixy_inner)
            end if
          end if
        end if
      else if (l_fluxin) then !flux in, but not out
        dn1=flux_n1(k+1)*rdz_on_rho(k,ixy_inner)
        if (params%l_2m) dn2=flux_n2(k+1)*rdz_on_rho(k,ixy_inner)
        !if (params%l_3m) dn3=flux_n3(k+1)*rdz_on_rho(k)

        !============================
        ! aerosol processing
        !============================
        if (l_ased .and. l_da_local) then
          if (params%id == cloud_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
            end if
            if (.not. l_warm) then
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
              if (l_passivenumbers_ice) then
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
              end if
            end if
          else if (params%id == rain_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
            end if
            if (.not. l_warm) then
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
              if (l_passivenumbers_ice) then
                dnumber_d=flux_n2(k+1)*dustact(k+1)%nratio2*rdz_on_rho(k,ixy_inner)
              end if
            end if
          end if

          if (params%id == ice_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
            end if
          else if (params%id == snow_params%id) then
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            if (l_passivenumbers) then
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
            end if
            if (l_passivenumbers_ice) then
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
            end if
          else if (params%id == graupel_params%id) then
            if (i_aerosed_method==1) then
              dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean)*rdz_on_rho(k,ixy_inner)
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean)*rdz_on_rho(k,ixy_inner)
              if (l_passivenumbers) then
                dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3)*rdz_on_rho(k,ixy_inner)
              end if
              if (l_passivenumbers_ice) then
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3)*rdz_on_rho(k,ixy_inner)
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
            aerosol_procs(i_am5, iaproc%id)%column_data(k)=dmac
          else
            aerosol_procs(i_am4, iaproc%id)%column_data(k)=dmac
          end if
          if (.not. l_warm) aerosol_procs(i_am9, iaproc%id)%column_data(k) =dmad
          if (l_passivenumbers) then
            aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
          end if
          if (l_passivenumbers_ice) then
            aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d
          end if

        else
          !ice phase
          aerosol_procs(i_am7, iaproc%id)%column_data(k)=dmad
          aerosol_procs(i_am8, iaproc%id)%column_data(k)=dmac
          if (l_passivenumbers) then
            aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
          end if
          if (l_passivenumbers_ice) then
            aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d
          end if
        end if
      end if

      dm1=dn1
      dm2=dn2
      !dm3=dn3

      procs(params%i_1m, iproc%id)%column_data(k)=c_x*dm1

      if (params%l_2m) procs(params%i_2m, iproc%id)%column_data(k)=dm2

      !if (params%l_3m) procs(params%i_3m, iproc%id)%column_data(k)=dm3
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine sedr

subroutine sedr_1M_2M(ixy_inner, step_length, qfields, aeroact, dustact,   &
       params, procs, aerosol_procs, precip1d, l_doaerosol)

!!! This routine has the same functionality as sedr, except it will only work with 
!!! single moment or double moment settings, since the fallspeeds are calced using 
!!! pre-calculated gamma functions with fixed shape. This routine is default when 
!!! the third prognostic moment is switched off. 

! UM dummy routines
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE um_types, ONLY: real_lsprec

use casim_parent_mod, only: casim_parent, parent_um

implicit none

character(len=*), parameter :: RoutineName='SEDR_1M_2M'

integer, intent(in) :: ixy_inner

real(wp), intent(in) :: step_length
real(wp), intent(in), target :: qfields(:,:)
type(hydro_params), intent(in) :: params
type(aerosol_active), intent(in) :: aeroact(:), dustact(:)
type(process_rate), intent(inout), target :: procs(:,:)
type(process_rate), intent(inout), target :: aerosol_procs(:,:)
real(wp), intent(out) :: precip1d(nz)
logical, optional, intent(in) :: l_doaerosol

real(wp) :: dn1, dn2
real(wp) :: m1, m2
real(wp) :: n1, n2
real(wp) :: hydro_mass
real(wp) :: n0, lam, mu, u1r, u2r, u1r2, u2r2
integer :: k

integer :: i_1m, i_2m
real(wp) :: p1, p2 
real(wp) :: sp1, sp2 
real :: a_x, b_x, f_x, g_x, c_x, d_x
real(wp) :: a2_x, b2_x, f2_x
logical :: l_fluxin, l_fluxout

! AH: variables needed for lsp_sedim_eulexp. All declared using real_lsprec from the UM or 
!     the dummy real_lsprec from MONC casim component, which is wp
real(real_lsprec) :: lsp_sedim_c_x
real(real_lsprec) :: u1r_above, u2r_above
!   real(real_lsprec) :: u1w, u2w
!   real(real_lsprec) :: mixingratio_fromabove, mixingratio
!   real(real_lsprec) :: numberconc, numberconc_fromabove

real(real_lsprec) :: m0 
integer :: points
real(real_lsprec) :: dhi(1), dhir(1),rhor(1),rhoin(1),mixratio_thislayer(1), &
     flux_fromabove(1),fallspeed_fromabove(1), fallspeed_thislayer(1),total_flux_out(1)
! End declaration of variables for lsp_sedim_eulexp

type(process_name) :: iproc, iaproc  ! processes selected depending on
! which species we're depositing on.

real(wp) :: dmac
real(wp) :: dmad
real(wp) :: dnumber_a, dnumber_d
logical :: l_da_local  ! local tranfer of l_doaerosol
logical :: l_sedim_generic
! If this is used, you can't trust the results.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (present(l_doaerosol)) then
  l_da_local=l_doaerosol
else
  l_da_local=.false.
end if

if (l_kfsm) then
  l_sedim_generic = (params%id==ice_params%id)
else
  l_sedim_generic = .false.
  a_x=params%a_x
  b_x=params%b_x
end if
    
! precip diag
do k = 1, nz
  precip1d(k) = 0.0
  flux_n1(k)  = 0.0
  Grho(k)=(rho0/rho(k,ixy_inner))**params%g_x
end do

i_1m=params%i_1m
i_2m=params%i_2m
p1=params%p1
p2=params%p2

sp1=params%p1
sp2=params%p2

! we don't want flexible approach in this version....
if (p1/=sp1 .or. p2/=sp2) then
  write(std_msg, '(A)') 'Cannot have flexible sedimentation options '//&
                        'with CASIM aerosol'
  call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                         std_msg)
end if

f_x=params%f_x
g_x=params%g_x
c_x=params%c_x
lsp_sedim_c_x = c_x
d_x=params%d_x
a2_x=params%a2_x
b2_x=params%b2_x
f2_x=params%f2_x

if (params%l_2m) flux_n2(:)=0.0          

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

! AH : initialise terminal velocity for level+1 before loop starts
!     
u1r_above=0.0
u2r_above=0.0

do k=nz-1, 1, -1

  ! initialize to zero
  dn1=0.0
  dn2=0.0
  n1=0.0
  n2=0.0
  m2=0.0

  hydro_mass=qfields(k, params%i_1m)
  if (params%l_2m) m2=qfields(k, params%i_2m)

  if (casim_parent == parent_um .and. l_sed_eulexp ) then
    ! AH - the code below is required to run in the UM (GA) and 
    !      it over-rides the above conditions. This code should be wrapped in 
    !      in a condition that relates to l_sed_eulexp
    !      NOTE: at this stage it is not clear to me why this has to 
    !            be set to true when sed_eulexp is true!
    l_fluxin=.true. 
    l_fluxout=.true.
  else
    l_fluxin=.false.
    l_fluxout=.false.
    if (hydro_mass > thresh_small(params%i_1m)) l_fluxout=.true.
    if (qfields(k+1, params%i_1m) > thresh_small(params%i_1m)) l_fluxin=.true.
  endif

  if (l_sedim_generic) then
    a_x=a_s(k)
    b_x=b_s(k)
  else
    a_x=params%a_x
    b_x=params%b_x
  end if

  if (l_fluxout) then
    m1=(hydro_mass/c_x)

    n0=dist_n0(k,params%id)
    mu=dist_mu(k,params%id)
    lam=dist_lambda(k,params%id)

    n1=m1*rho(k,ixy_inner)
    n2=m2*rho(k,ixy_inner)

    if (lam > 0.0_wp) then 
      if (l_kfsm) then 
        u1r=a_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f_x)**(-(1.0+mu+sp1+b_x)))     &
            *(Gammafunc(1.0+mu+sp1+b_x)/Gammafunc(1.0+mu+sp1))
      else 
        u1r=a_x*Grho(k)*(lam**(params%exp_1_mu_sp1)*(lam+f_x)**(-(params%exp_1_mu_sp1_bx)))     &
            *(params%gam_1_mu_sp1_bx/params%gam_1_mu_sp1)
        if (params%l_2m)     &
            u2r=a_x*Grho(k)*(lam**(params%exp_1_mu_sp2)*(lam+f_x)**(-(params%exp_1_mu_sp2_bx))) &
                *(params%gam_1_mu_sp2_bx/params%gam_1_mu_sp2)
      endif
        
      if (l_abelshipway .and. params%id==rain_params%id) then ! rain can use abel and shipway formulation
        u1r2=a2_x*Grho(k)*(lam**(params%exp_1_mu_sp1)*(lam+f2_x)**(-(params%exp_1_mu_sp1_b2x)))   &
             *(params%gam_1_mu_sp1_b2x/params%gam_1_mu_sp1)
        u1r=u1r+u1r2
           
        if (params%l_2m) then
          !u2r2=a2_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f2_x)**(-(params%exp_1_mu_sp1_b2x))) &
          !     *(Gammafunc(1.0+mu+sp2+b2_x)/Gammafunc(1.0+mu+sp2))
          u2r2=a2_x*Grho(k)*(lam**(params%exp_1_mu_sp2)*(lam+f2_x)**(-(params%exp_1_mu_sp2_b2x))) &
               *(params%gam_1_mu_sp2_b2x/params%gam_1_mu_sp2)
          u2r=u2r+u2r2
        end if
      endif

      ! fall speeds shouldn't get too big...
      u1r=min(u1r,params%maxv)

      ! For clouds and ice, we only use a 1M representation of sedimentation
      ! so that spurious size sorting doesn't lead to overactive autoconversion
      if (params%id==cloud_params%id .or. params%id==ice_params%id) then
        if (l_sed_icecloud_as_1m)then
          if (params%l_2m)u2r=u1r
        else
          if (params%l_2m)u2r=min(u2r,params%maxv)
        end if
      else
        if (params%l_2m)u2r=min(u2r,params%maxv)
      end if

      if ( (params%id==rain_params%id .and. l_sed_rain_1m) .or. &
           (params%id==snow_params%id .and. l_sed_snow_1m) .or. &
           (params%id==graupel_params%id .and. l_sed_graupel_1m) )  then
        if (params%l_2m)u2r=u1r
      end if

      ! fall speeds shouldn't be negative (can happen with original AS formulation)
      u1r=max(u1r,0.0_wp)
      if (params%l_2m)u2r=max(u2r,0.0_wp)
    else 
      u1r=0.0_wp
      if (params%l_2m)u2r=0.0_wp
    endif
    flux_n1(k)=n1*u1r
        
    if (params%l_2m) flux_n2(k)=n2*u2r
           

    ! AH : derive the sedimentation flux
    if (casim_parent == parent_um .and. l_sed_eulexp) then 
            ! use the method, which is based on UM and is most appropriate 
            ! (and stable) for  long timesteps
      !PRF use WB method -overwrite everything for now - put logicals in
      !print *, 'sed_eulexp called', k, hydro_mass
      if (hydro_mass < thresh_small(params%i_1m)) then
        u1r=0.0
        n1=0.0
        if (params%l_2m) then
          u2r=0.0
          n2=0.0
        endif
      endif
           
      points=1
      m0=1e-10/lsp_sedim_c_x  !for the mass - does not apply to number - 
                              !fix zero in one but not the other below)
      dhi(1)=step_length/dz(k,ixy_inner)
      dhir(1)=dz(k,ixy_inner)/step_length
      rhoin(1)=rho(k,ixy_inner)
      rhor(1)=1.0/rho(k,ixy_inner)
      mixratio_thislayer(1)=lsp_sedim_c_x*n1   /rho(k,ixy_inner)
      flux_fromabove(1)=lsp_sedim_c_x*flux_n1(k+1)
      fallspeed_fromabove(1)=u1r_above
      fallspeed_thislayer(1)=u1r
      total_flux_out(1)=0.0

      call lsp_sedim_eulexp(points,m0,dhi,dhir,rhoin,rhor,                     &
                            flux_fromabove, fallspeed_thislayer,               &
                            mixratio_thislayer, fallspeed_fromabove,           &
                            total_flux_out)
      flux_n1(k)=total_flux_out(1)/lsp_sedim_c_x
      u1r_above=fallspeed_thislayer(1)
           
      if (params%l_2m) then
        points=1
        m0=1e-10/lsp_sedim_c_x  !for the mass - does not apply to number - 
                                !fix zero in one but not the other below)
        dhi(1)=step_length/dz(k,ixy_inner)
        dhir(1)=dz(k,ixy_inner)/step_length
        rhoin(1)=rho(k,ixy_inner)
        rhor(1)=1.0/rho(k,ixy_inner)
        mixratio_thislayer(1)=n2/rho(k,ixy_inner)
        flux_fromabove(1)=flux_n2(k+1)
        fallspeed_fromabove(1)=u2r_above
        fallspeed_thislayer(1)=u2r
        total_flux_out(1)=0.0
        
        call lsp_sedim_eulexp(points,m0,dhi,dhir,rhoin,rhor,                   &
                              flux_fromabove, fallspeed_thislayer,             &
                              mixratio_thislayer, fallspeed_fromabove,         &
                              total_flux_out)
        flux_n2(k)=total_flux_out(1)
        u2r_above=fallspeed_thislayer(1)
      endif
           
      if (flux_n1(k) .lt. epsilon(1.0_wp) ) then
        flux_n1(k)=0.0
        u1r_above=0.0
      endif

      ! if the fallspeed is within a typical number value of the precision
      ! used in lsp_sedim_eulexp, there is the chance a floating point
      ! error will occur on the next level down, hence stop sedimentation     
      if (params%l_2m .and. (flux_n2(k) .lt. epsilon(1.0_wp) &
           .or. u2r_above .lt. 1.0e5*tiny(1.0_real_lsprec))) then
        flux_n2(k)=0.0
        u2r_above=0.0
      endif    
        
!!PRF
    endif
!         else if (parent_model
!                       ! ! Ported code that produces very different results to standard sed
!             flux_fromabove(1) = c_x * flux_n1(k+1)
!             mixingratio_fromabove = (flux_fromabove(1)*step_length/dz(k)*(1/rho(k)))
!             mixingratio = c_x * n1 * (1.0/rho(k))
           
!             if (mixingratio_fromabove + mixingratio > epsilon(1.0_wp)) then 
!                u1w=(u1r*mixingratio + &
!                     u1r_above*mixingratio_fromabove)/ & 
!                     (mixingratio_fromabove + mixingratio) 
!             else

!               u1w = 0.0_wp
              
!            endif

!            if (u1w > 0.0_wp) then 

!               flux_n1(k)= flux_fromabove(1) + (dz(k)/step_length) * &
!                    (rho(k)*mixingratio - flux_fromabove(1)/u1w) * &
!                    (1.0_wp - exp(-u1w*(step_length/dz(k))))
              
!               flux_n1(k) = flux_n1(k)/c_x
              
!            else
              
!               flux_n1(k) = 0.0
             
!            endif
              
!            u1r_above = u1w

!            if (flux_n1(k) < epsilon(1.0_wp)) then
!               flux_n1(k)=0.0
!               u1r_above=0.0
!            endif
           
!            if  (params%l_2m) Then
              
!               numberconc = n2 * (1.0/rho(k))
!               numberconc_fromabove = (flux_n2(k+1)*step_length/dz(k)*(1/rho(k)))

!               if (numberconc_fromabove + numberconc > epsilon(1.0_wp)) then 
!                  u2w=(u2r*numberconc + & 
!                       u2r_above*numberconc_fromabove)/ &
!                       (numberconc+numberconc_fromabove) ! PF:number also mass wtd 
!                                                         ! AH: no it is not, it is number weighted
!               else
                 
!                  u2w = 0.0_wp
                 
!               endif
              
!               if ( u2w > 0.0_wp ) then 
                 
!                  flux_n2(k)=flux_n2(k+1) + (dz(k)/step_length) * &
!                       (rho(k)*numberconc - flux_n2(k+1)/u2w) * &
!                       (1.0-exp(-u2w*(step_length/dz(k))))
!               else
                 
!                  flux_n2(k) = 0.0
                 
!               endif
                 
!               u2r_above = u2w

!               if (flux_n2(k) < epsilon(1.0_wp)) then
!                  flux_n2(k)=0.0
!                  u2r_above=0.0
!               endif
!            endif
! #endif
!        endif ! L_sed_eulexp

    precip1d(k) = flux_n1(k)*c_x
        ! diagnostic for precip

  end if

  dmac=0.0
  dmad=0.0
  dnumber_a=0.0
  dnumber_d=0.0

  if (l_fluxout) then !flux out (flux(k+1) will be zero if no flux in)
    dn1=(flux_n1(k+1)-flux_n1(k))*rdz_on_rho(k,ixy_inner)
    if (params%l_2m) dn2=(flux_n2(k+1)-flux_n2(k))*rdz_on_rho(k,ixy_inner)

    !============================
    ! aerosol processing
    !============================
    if (l_ased .and. l_da_local) then
      if (params%id == cloud_params%id) then
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean -      &
             flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)* rdz_on_rho(k,ixy_inner)
        if (l_passivenumbers) then
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1 -                       &
                    flux_n2(k)*aeroact(k)%nratio1)* rdz_on_rho(k,ixy_inner)
        end if
        if (.not. l_warm) then
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-     &
               flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
          if (l_passivenumbers_ice .and. dustact(k)%mact_mean > 0.0) then
            dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-                      &
                      flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
          end if
        end if
      else if (params%id == rain_params%id) then
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-       &
             flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        if (l_passivenumbers .and. aeroact(k)%mact_mean > 0.0) then
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-                        &
                    flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
        end if
        if (.not. l_warm) then
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-     &
               flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
          if (l_passivenumbers_ice) then
            dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-                      &
                      flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
          end if
        end if
      end if

      if (params%id == ice_params%id) then
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean-     &
             flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        dmad = (flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-     &
             flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        if (l_passivenumbers) then
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1-                        &
                    flux_n2(k)*aeroact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
        end if
        if (l_passivenumbers_ice) then
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-                        &
                    flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
        end if
      else if (params%id == snow_params%id) then
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-       &
             flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-       &
             flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        if (l_passivenumbers) then
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-                        &
                    flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
        end if
        if (l_passivenumbers_ice) then
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-                        &
                    flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
        end if
      else if (params%id == graupel_params%id) then
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean-       &
             flux_n2(k)*aeroact(k)%nratio3*aeroact(k)%mact3_mean)*rdz_on_rho(k,ixy_inner)
        dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean-       &
             flux_n2(k)*dustact(k)%nratio3*dustact(k)%mact3_mean)*rdz_on_rho(k,ixy_inner)

        if (l_passivenumbers) then
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3-                        &
                    flux_n2(k)*aeroact(k)%nratio3)*rdz_on_rho(k,ixy_inner)
        end if
        if (l_passivenumbers_ice) then
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3-                        &
                    flux_n2(k)*dustact(k)%nratio3)*rdz_on_rho(k,ixy_inner)
        end if
      end if
    end if
  else if (l_fluxin) then !flux in, but not out
    dn1=flux_n1(k+1)*rdz_on_rho(k,ixy_inner)
    if (params%l_2m) dn2=flux_n2(k+1)*rdz_on_rho(k,ixy_inner)

    !============================
    ! aerosol processing
    !============================
    if (l_ased .and. l_da_local) then
      if (params%id == cloud_params%id) then
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        if (l_passivenumbers) then
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
        end if
        if (.not. l_warm) then
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
          if (l_passivenumbers_ice) then
            dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
          end if
        end if
      else if (params%id == rain_params%id) then
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        if (l_passivenumbers) then
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
        end if
        if (.not. l_warm) then
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
          if (l_passivenumbers_ice) then
            dnumber_d=flux_n2(k+1)*dustact(k+1)%nratio2*rdz_on_rho(k,ixy_inner)
          end if
        end if
      end if

      if (params%id == ice_params%id) then
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        if (l_passivenumbers) then
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
        end if
        if (l_passivenumbers_ice) then
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
        end if
      else if (params%id == snow_params%id) then
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        if (l_passivenumbers) then
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
        end if
        if (l_passivenumbers_ice) then
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
        end if
      else if (params%id == graupel_params%id) then
        if (i_aerosed_method==1) then
          dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean)*rdz_on_rho(k,ixy_inner)
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean)*rdz_on_rho(k,ixy_inner)
          if (l_passivenumbers) then
            dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3)*rdz_on_rho(k,ixy_inner)
          end if
          if (l_passivenumbers_ice) then
            dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3)*rdz_on_rho(k,ixy_inner)
          end if
        else
          write(std_msg, '(A)') 'ERROR: GET RID OF i_aerosed_method variable!'
          call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, std_msg )
        end if
      end if
    end if
  end if

  ! Store the aerosol process terms...
  if (l_da_local) then
    if (params%id == cloud_params%id .or. params%id == rain_params%id) then
      !liquid phase
      if (l_separate_rain .and. params%id == rain_params%id) then
        aerosol_procs(i_am5, iaproc%id)%column_data(k)=dmac
      else
        aerosol_procs(i_am4, iaproc%id)%column_data(k)=dmac
      end if
      if (.not. l_warm) aerosol_procs(i_am9, iaproc%id)%column_data(k) =dmad
      if (l_passivenumbers) then
        aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
      end if
      if (l_passivenumbers_ice) then
        aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d
      end if

    else
      !ice phase
      aerosol_procs(i_am7, iaproc%id)%column_data(k)=dmad
      aerosol_procs(i_am8, iaproc%id)%column_data(k)=dmac
      if (l_passivenumbers) then
        aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
      end if
      if (l_passivenumbers_ice) then
        aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d
      end if
    end if
  end if

  procs(params%i_1m, iproc%id)%column_data(k)=c_x*dn1

  if (params%l_2m) procs(params%i_2m, iproc%id)%column_data(k)=dn2

end do

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine sedr_1M_2M


  subroutine terminal_velocity_CFL(dt, vt, nsubseds_hydro, sed_length_hydro, nsubseds_max, & 
    sed_length_max, mindz)

    ! work out the number of substeps for each hydrometeor
    ! using the minimum dz and maxv to determine a CFL limit
    ! NOTE: if the parent_model has a fixed timestep this calc
    !       should be done in initialise_sedr. If there is a dynamic
    !       timestep, this calc has to be done on every timestep
    
    implicit none

    real(wp), intent(in) :: dt  ! timestep from parent model
    real(wp), intent(in) :: vt  ! max hydrometeor terminal velocity for a hydrometeor
                                ! (comes from mphys_params)
    ! Both of the following varaibles are derived from mphys_switches
    integer, intent(in)  ::  nsubseds_max  ! number of sedimentation substeps for sedimentation
    real(wp), intent(in) :: sed_length_max ! sedimentation substep length sedimentation
    ! cfl_vt_max is declared and set in mphys_switches
    real(wp), intent(in) :: mindz  ! pass it in 
    

    integer, intent(out)  ::  nsubseds_hydro  ! number of sedimentation substeps for a hydrometeor
    real(wp), intent(out) :: sed_length_hydro ! sedimentation substep length for a hydormeteor

    ! local variables
    real(wp) :: cfl_vt ! CFL value for the terminal velocity argument
    

    
    cfl_vt =( dt * vt ) / mindz

    !!write(6,*) 'PRF1s ',cfl_vt, cfl_vt_max

    if ( cfl_vt > cfl_vt_max ) then
       nsubseds_hydro = max(1, ceiling(cfl_vt / cfl_vt_max))
       sed_length_hydro = dt / (real(nsubseds_hydro, kind=wp))
       
       if (sed_length_hydro.lt. sed_length_max) then
          nsubseds_hydro = nsubseds_max
          sed_length_hydro = sed_length_max
       endif

       
    else    
       nsubseds_hydro = 1
       sed_length_hydro = dt
    endif

   !! write(6,*) 'PRF2s ',dt,nsubseds_hydro,sed_length_hydro,nsubseds_max,sed_length_max
  
  end subroutine terminal_velocity_CFL

end module sedimentation
