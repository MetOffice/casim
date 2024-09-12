!Fields which aren't updated during the microphysics calculation
module passive_fields
  use mphys_die, only: throw_mphys_error, bad_values, std_msg
  use variable_precision, only: wp
  use qsat_funs, only: qsaturation
  use mphys_switches, only: i_th
  use mphys_parameters, only: nxy_inner

  implicit none
  private

  real(wp), allocatable :: rho(:,:), pressure(:,:), z(:), exner(:,:), rexner(:,:)
  real(wp), allocatable :: dz(:,:)
  real(wp), allocatable :: qws(:,:), qws0(:,:), TdegC(:,:), TdegK(:,:), w(:,:), tke(:)

!$OMP THREADPRIVATE(rho, pressure, z, exner, rexner, dz, &
!$OMP               qws, qws0, TdegC, TdegK, w, tke)

  real(wp), allocatable :: rhcrit_1d(:)
!$OMP THREADPRIVATE(rhcrit_1d)

  real(wp), allocatable :: rdz_on_rho(:,:)
!$OMP THREADPRIVATE(rdz_on_rho)

  real(wp) :: dt
  real(wp), allocatable :: min_dz(:) ! minimum vertical resolution
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
    allocate(rho(nz, nxy_inner))
    allocate(pressure(nz, nxy_inner))
    allocate(exner(nz, nxy_inner))
    allocate(rexner(nz, nxy_inner))
    allocate(dz(nz, nxy_inner))
    allocate(w(nz, nxy_inner))
    allocate(tke(nz)) !! tke(:) is not used in other subroutines. So no need to have inner dimension
    allocate(rdz_on_rho(nz, nxy_inner))
    allocate(qws(nz, nxy_inner))
    allocate(qws0(nz, nxy_inner))
    allocate(TdegK(nz, nxy_inner))
    allocate(TdegC(nz, nxy_inner))
    allocate(min_dz(nxy_inner))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine initialise_passive_fields

  subroutine set_passive_fields(nxy_inner_loop, ixy_outer, is_in, js_in, je_in, &
                                dt_in, rho_in, p_in, exner_in,   &
                                dz_in, w_in, tke_in, qfields)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='SET_PASSIVE_FIELDS'

    integer, intent(in) :: nxy_inner_loop
    integer, intent(in) :: ixy_outer
    integer, intent(in) :: is_in, js_in, je_in

    real(wp), intent(in) :: dt_in
    real(wp), intent(in) :: rho_in(:,:,:), p_in(:,:,:), exner_in(:,:,:)
    real(wp), intent(in) :: dz_in(:,:,:)
    real(wp), intent(in) :: w_in(:,:,:), tke_in(:,:,:)
    real(wp), intent(in), target :: qfields(:,:,:)
    integer :: k, ixy_inner, ixy, jy, ix

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dt=dt_in

    do ixy_inner=1, nxy_inner_loop

      ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
      jy = modulo(ixy-1,(je_in-js_in+1))+js_in
      ix = (ixy-1)/(je_in-js_in+1)+is_in


      do k=1,nz
        rho(k,ixy_inner)=rho_in(k,ix,jy)
        pressure(k, ixy_inner)=p_in(k,ix,jy)
        exner(k, ixy_inner)=exner_in(k,ix,jy)
        rexner(k, ixy_inner)=1.0/exner(k, ixy_inner)
        w(k, ixy_inner)=w_in(k,ix,jy)
        !tke(k)=tke_in(kl:ku)
        dz(k, ixy_inner)=dz_in(k,ix,jy)
        rdz_on_rho(k, ixy_inner)=1.0/(dz_in(k,ix,jy)*rho_in(k,ix,jy))
        !  do k=1,nz
        TdegK(k, ixy_inner)=qfields(k,i_th,ixy_inner)*exner(k, ixy_inner)
        TdegC(k, ixy_inner)=TdegK(k, ixy_inner)-273.15
        qws(k, ixy_inner)=qsaturation(TdegK(k, ixy_inner), pressure(k, ixy_inner)/100.0)
        qws0(k, ixy_inner)=qsaturation(273.15_wp, pressure(k, ixy_inner)/100.0)
      end do ! k

      qws(nz, ixy_inner)=1.0e-8

      min_dz(ixy_inner) = minval(dz(:,ixy_inner))

      if (any(qws(:,ixy_inner)==0.0)) then
        write(std_msg, '(A)') 'Error in saturation calculation - qws is zero'
        call throw_mphys_error( bad_values, ModuleName//':'//RoutineName,     &
                                std_msg )

      end if
    end do ! ixy_inner
    !nullify(theta)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine set_passive_fields
end module passive_fields
