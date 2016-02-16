!Fields which aren't updated during the microphysics calculation
module passive_fields
  use mphys_die, only: throw_mphys_error
  use variable_precision, only: wp
  use qsat_funs, only: qsaturation
  use mphys_switches, only: i_th

#if DEF_MODEL==MODEL_KiD
  use diagnostics, only: i_here, j_here
#elif  DEF_MODEL==MODEL_UM
  use diaghelp_um, only: i_here, j_here
#elif  DEF_MODEL==MODEL_LEM
  use diaghelp_lem, only: i_here, j_here
#elif  DEF_MODEL==MODEL_MONC
  use diaghelp_monc, only: i_here, j_here
#endif

  implicit none
  private

  real(wp), allocatable :: rho(:), pressure(:), z(:), exner(:), rexner(:)
  real(wp), allocatable :: z_half(:), z_centre(:), dz(:)
  real(wp), allocatable :: qws(:), qws0(:), TdegC(:), TdegK(:), w(:), tke(:)

  real(wp), allocatable :: rdz_on_rho(:)
  real(wp) :: dt
  integer :: kl, ku, nz

  public set_passive_fields, rho, pressure, initialise_passive_fields, z, qws, exner, rexner, z_half, z_centre, dz, qws0, &
       TdegC, TdegK, w, tke, rdz_on_rho
contains

  subroutine initialise_passive_fields(kl_arg, ku_arg)
    integer, intent(in) :: kl_arg, ku_arg

    kl=kl_arg
    ku=ku_arg
    nz=ku-kl+1
    allocate(rho(nz))
    allocate(pressure(nz))
    allocate(exner(nz))
    allocate(rexner(nz))
    allocate(dz(nz))
    allocate(z_half(0:nz))
    allocate(z_centre(nz))
    allocate(w(nz))
    allocate(tke(nz))
    allocate(rdz_on_rho(nz))
    allocate(qws(nz))
    allocate(qws0(nz))
    allocate(TdegK(nz))
    allocate(TdegC(nz))
  end subroutine initialise_passive_fields

  subroutine set_passive_fields(dt_in, rho_in, p_in, exner_in,   &
       z_half_in, z_centre_in, dz_in, w_in, tke_in, qfields)
    real(wp), intent(in) :: dt_in
    real(wp), intent(in) :: rho_in(kl:ku), p_in(kl:ku), exner_in(kl:ku)
    real(wp), intent(in) :: z_half_in(kl-1:ku),z_centre_in(kl:ku),dz_in(kl:ku)
    real(wp), intent(in) :: w_in(kl:ku), tke_in(kl:ku)
    real(wp), intent(in), target :: qfields(:,:)
    integer :: k

    real(wp), pointer :: theta(:)

    theta=>qfields(:,i_th)

    dt=dt_in

    rho(:)=rho_in(kl:ku)
    pressure(:)=p_in(kl:ku)
    exner(:)=exner_in(kl:ku)
    rexner(:)=1.0/exner(:)
    w(:)=w_in(kl:ku)
    tke(:)=tke_in(kl:ku)
    dz(:)=dz_in(kl:ku)
    z_half(0:)=z_half_in(kl-1:ku)
    z_centre(:)=z_centre_in(kl:ku)
    rdz_on_rho(:)=1.0/(dz_in(kl:ku)*rho_in(kl:ku))
    do k=1,nz
      TdegK(k)=theta(k)*exner(k)
      TdegC(k)=TdegK(k)-273.15
      qws(k)=qsaturation(TdegK(k), pressure(k)/100.0)
      qws0(k)=qsaturation(273.15_wp, pressure(k)/100.0)
    end do

    qws(nz)=1.0e-8

    if (any(qws==0.0)) then
      do k=1, nz
        print *, 'DEBUG qws', i_here, j_here, k, theta(k), TdegK(k), qws(k),rdz_on_rho(k)
      end do
      call throw_mphys_error(2, 'passive_fields', 'Error in saturation calculation')
    end if
    nullify(theta)
  end subroutine set_passive_fields
end module passive_fields
