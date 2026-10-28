!===============================================================================
! MODULE: hail_diagnostic_fast_mod

!===============================================================================
module hail_diagnostic_fast_mod

! Some of the content of this change has been produced with the assistance
! of Claude Met Office Github Copilot Enterprise

  USE precision, ONLY: wp

  USE umPrintMgr, ONLY: umPrint, umMessage

  USE mphys_constants, ONLY: pi, Rd, Rv, Lf, Lv, ka, Dv, rho0

  implicit none
  private
  public :: diagnose_hail_fast

  ! ---- physical constants ----

  real(wp), parameter :: Kw0    = 0.56_wp         ! W/m/K, thermal conductivity of water (~0C)
  real(wp), parameter :: T0     = 273.15_wp      ! K, melting point
  real(wp), parameter :: es0    = 611.2_wp       ! Pa, saturation vapor pressure at 0C
  real(wp), parameter :: Pr_no  = 0.71_wp        ! Prandtl number of air
  real(wp), parameter :: Sc_no  = 0.60_wp        ! Schmidt number of air
  real(wp), parameter :: visc   = 1.7e-5_wp      ! visc kg/m/s

  ! ---- terminal velocity power law: V = a_v * (D_cm)^b_v * sqrt(rho0_ref/rho_air) ----
  real(wp), parameter :: a_v = 14.0_wp  ! for cm; MASON 1971
  real(wp), parameter :: b_v = 0.5_wp

  ! ---- minimum thresholds ----
  real(wp), parameter :: qg_min = 1.0e-6_wp      ! kg/kg, minimum graupel to bother

contains

  subroutine diagnose_hail_fast( ixy_inner, nz, nq, qfields,                   &
                                  D_max_sfc, precip_rate,                      &
                                  hail_flag, ierr )

    USE distributions,        ONLY: query_distributions, dist_lambda,          &
                                    dist_mu, dist_n0
    USE mphys_parameters,     ONLY: graupel_params
    USE passive_fields,       ONLY: Tdegk, dz, pressure, rho
    use mphys_switches, only: i_qv, i_qg
    
    real(wp), parameter :: Nthresh=1e-6

    integer,  intent(in)  :: nz, nq
    integer,  intent(in)  :: ixy_inner

    real(wp), intent(in) :: qfields(nz,nq)
!    real(wp), intent(in) :: cffields(nz,nq) !not used at the moment. 

    real(wp), intent(out) :: D_max_sfc   ! largest hailstone diameter at surface [m]
    real(wp), intent(out) :: precip_rate ! surface hail mass flux [kg/m2/s]
                                          ! (x3600 for mm/hr liquid-equivalent)
    real(wp), intent(out) :: hail_flag   ! .true. if any hail reaches the surface
    integer,  intent(out) :: ierr        ! 0=ok, 1=no graupel, 2=no melting level

    real(wp) :: rhoa(nz), qv(nz), t(nz), qg(nz), dz_in(nz)
    real(wp) :: z_ml
    real(wp) :: N_t, lam
    real(wp) :: rho_g, mu_g
    integer  :: k_ml, nn, k
    real(wp) :: a, a0, v, da

    integer, parameter :: nh = 7 !ensure that this matches the haildedge array length
    real(wp), parameter :: rhoi = 900.0
    real(wp) :: haild(nh)
    real(wp) :: haildf(nh)
    real(wp) :: haildedge(nh+1)
    real(wp) :: hailprecipf(nh)
    real(wp) :: hailbin(nh)
    real(wp) :: hailconc(nh)

    data haildedge /1e-3, 3e-3,6e-3,1e-2,3e-2,6e-2,1e-1, 3e-1/  ! initial hail size edges at melting layer diameter in m
        
    rho_g=graupel_params%density
    qv=qfields(:, i_qv)
    qg=qfields(:, i_qg)
    t(:)=TdegK(:,ixy_inner)
    rhoa(:)=rho(:,ixy_inner)
    dz_in(:)=dz(:,ixy_inner)

    D_max_sfc = 0.0_wp
    precip_rate = 0.0_wp
    hail_flag = 0.0_wp
    ierr = 0

    if (maxval(qg) < qg_min) then
      ierr = 1
      return
    end if


    call find_melting_level(nz, dz_in, t, k_ml, z_ml, ierr)
    if (ierr /= 0) return
    mu_g=dist_mu(k_ml,graupel_params%id)
    lam=dist_lambda(k_ml, graupel_params%id)
    N_t=dist_n0(k_ml, graupel_params%id)

    do nn=1,nh
      haild(nn)=(haildedge(nn+1)+haildedge(nn))/2.0  !initial diameters
      haildf(nn)=0.0_wp      !final diameter
      hailprecipf(nn)=0.0_wp   ! precip in bins
      hailbin(nn)=(haildedge(nn+1)-haildedge(nn))  ! bin width
      hailconc(nn)=N_t*lam**(mu_g+1)/gamma(mu_g+1)*haild(nn)**mu_g*exp(-lam*haild(nn))*hailbin(nn)  !conc from graupel dist
    end do

    if ( maxval(hailconc) < Nthresh)  then
      ierr = 2
      return
    end if
    
    do nn=1,nh
      a0=haild(nn)/2.0 ! convert to radius at melting level
      a=a0
      do k = k_ml,1,-1
        v=a_v*(2.0*a*100.0)**b_v * sqrt(rho0/rhoa(k)) !diam in cm
        da=mason_melt(a,a0,dz_in(k),t(k)-T0,v,rhoi,rhoa(k))
        a=a-da
      end do
      if ( a .gt. 0.0 .and. hailconc(nn) .gt. Nthresh) then
        haildf(nn)=2.0*a
        hailprecipf(nn)=hailconc(nn)*v*pi/6.0*haildf(nn)**3*rhoi*rhoa(1)
        D_max_sfc=2.0*a
        hail_flag=1.0
      else
        haildf(nn)=0.0
        hailprecipf(nn)=0.0
      end if
    end do
    precip_rate=sum(hailprecipf(:)) !kg m-2 s-1

    


  end subroutine diagnose_hail_fast

  subroutine find_melting_level(nz, dz, t, k_ml, z_ml, ierr)

    integer,  intent(in)  :: nz
    real(wp), intent(in)  :: dz(nz), t(nz)
    integer,  intent(out) :: k_ml
    real(wp), intent(out) :: z_ml
    integer,  intent(out) :: ierr
    integer :: k
    real(wp) :: frac

    if (t(1) >= T0) then
      do k = nz-1, 1, -1

        if (t(k+1) < T0 .and. t(k) >= T0) then
          frac = (T0 - t(k)) / (t(k+1) - t(k))
          z_ml = SUM(dz(1:k-1)) + frac*dz(k)
          k_ml = k
          ierr = 0

          return
        end if
      end do
    end if

    if (t(1) < T0) then
      ! surface already below freezing -- no melting layer above the surface
      z_ml = 0.0_wp
      k_ml = 1
      ierr = 2
    end if

  end subroutine find_melting_level

  !==================
  ! Mason melting
  !==================
  function mason_melt(a,a0,dz,t,v,rhoi,rhoa) result(da)
    real(wp), intent(in) :: a,a0,dz,t,rhoi,rhoa,v
    real(wp) :: da, C, Re, beta
    beta=0.0 ! ignore condensation/evap
    
    Re=v*(2*a)/visc*rhoa
    C=1.6+0.3*Re**0.5
    da=(Kw0*t*dz/Lf/rhoi/v)/ &
      ((a0-a)*a/a0+(Kw0/(C*(Ka+Lv*Dv*beta)))*a**2/a0)


  end function mason_melt



end module hail_diagnostic_fast_mod



