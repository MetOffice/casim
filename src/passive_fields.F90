!Fields which aren't updated during the microphysics calculation
module passive_fields
  use mphys_die, only: throw_mphys_error, bad_values, std_msg
  use variable_precision, only: wp
  use qsat_funs, only: qsaturation
  use mphys_switches, only: i_th

  implicit none
  private

  real(wp), allocatable :: rho(:), pressure(:), z(:), exner(:), rexner(:)
  real(wp), allocatable :: dz(:)
  real(wp), allocatable :: qws(:), qws0(:), TdegC(:), TdegK(:), w(:), tke(:)

!$OMP THREADPRIVATE(rho, pressure, z, exner, rexner, dz, &
!$OMP               qws, qws0, TdegC, TdegK, w, tke)

  real(wp), allocatable :: rhcrit_1d(:)
!$OMP THREADPRIVATE(rhcrit_1d)

  real(wp), allocatable :: rdz_on_rho(:)
!$OMP THREADPRIVATE(rdz_on_rho)

  real(wp) :: dt
  real(wp) :: min_dz ! minimum vertical resolution
  integer :: kl, ku, nz

!$OMP THREADPRIVATE(dt,min_dz,nz,kl,ku)

  character(len=*), parameter, private :: ModuleName='PASSIVE_FIELDS'

  public set_passive_fields, rho, pressure, initialise_passive_fields, z, qws, &
       exner, rexner, dz, qws0, &
       TdegC, TdegK, w, tke, rhcrit_1d, rdz_on_rho, min_dz
contains

  subroutine initialise_passive_fields(kl_arg, ku_arg)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='INITIALISE_PASSIVE_FIELDS'

    integer, intent(in) :: kl_arg, ku_arg

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    kl=kl_arg
    ku=ku_arg
    nz=ku-kl+1
    allocate(rho(nz))
    allocate(pressure(nz))
    allocate(exner(nz))
    allocate(rexner(nz))
    allocate(dz(nz))
    allocate(w(nz))
    allocate(tke(nz))
    allocate(rdz_on_rho(nz))
    allocate(qws(nz))
    allocate(qws0(nz))
    allocate(TdegK(nz))
    allocate(TdegC(nz))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine initialise_passive_fields

  subroutine set_passive_fields(dt_in, rho_in, p_in, exner_in,   &
                                dz_in, w_in, tke_in, qfields)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SET_PASSIVE_FIELDS'

    real(wp), intent(in) :: dt_in
    real(wp), intent(in) :: rho_in(kl:ku), p_in(kl:ku), exner_in(kl:ku)
    real(wp), intent(in) :: dz_in(kl:ku)
    real(wp), intent(in) :: w_in(kl:ku), tke_in(kl:ku)
    real(wp), intent(in), target :: qfields(:,:)
    integer :: k

    real(wp), pointer :: theta(:)

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    theta=>qfields(:,i_th)

    dt=dt_in

    rho(:)=rho_in(kl:ku)
    pressure(:)=p_in(kl:ku)
    exner(:)=exner_in(kl:ku)
    rexner(:)=1.0/exner(:)
    w(:)=w_in(kl:ku)
    tke(:)=tke_in(kl:ku)
    dz(:)=dz_in(kl:ku)
    rdz_on_rho(:)=1.0/(dz_in(kl:ku)*rho_in(kl:ku))
    do k=1,nz
      TdegK(k)=theta(k)*exner(k)
      TdegC(k)=TdegK(k)-273.15
      qws(k)=qsaturation(TdegK(k), pressure(k)/100.0)
      qws0(k)=qsaturation(273.15_wp, pressure(k)/100.0)
    end do

    qws(nz)=1.0e-8

    min_dz = minval(dz)

    if (any(qws==0.0)) then
      write(std_msg, '(A)') 'Error in saturation calculation - qws is zero'
      call throw_mphys_error( bad_values, ModuleName//':'//RoutineName,     &
                              std_msg )
    end if
    nullify(theta)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine set_passive_fields
end module passive_fields
