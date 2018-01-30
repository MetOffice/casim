module m3_incs
  use variable_precision, only: wp
  use lookup, only: Gfunc
  use mphys_die, only: throw_mphys_error, warn, std_msg
  implicit none
  private

  character(len=*), parameter, private :: ModuleName='M3_INCS'

  public m3_inc_type4, m3_inc_type3, m3_inc_type2, m3_inc_type1
contains

  ! for changes of phase from category y to category x
  ! dM3_x = -(c_y/c_x)**(p3/3) * dM3_y
  ! c_x, c_y are densities of x and y respectively
  ! p3 is the third moment (must be the same for both x and y)
  ! This assumes both categories are spherical.
  subroutine m3_inc_type4(dm3_y, c_x, c_y, p3, dm3_x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='M3_INC_TYPE4'

    real(wp), intent(in) :: dm3_y
    real(wp), intent(in) :: c_x,c_y
    real(wp), intent(in) :: p3
    real(wp), intent(out) :: dm3_x

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dm3_x = -(c_y/c_x)**(p3/3.0) * dm3_y

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine m3_inc_type4

  ! tendency of m3 as a function of dm1, dm2
  ! assuming initial value for mu
  subroutine m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='M3_INC_TYPE3'

    real(wp), intent(in) :: dm1, dm2
    real(wp), intent(in) :: p1, p2, p3
    real(wp), intent(in) :: mu
    real(wp), intent(out) :: dm3
    
    real(wp) :: k1, k2, k3

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    k3 = p2-p1
    k1 = (p3-p2)/k3
    k2 = (p1-p3)/k3

    dm3 = (Gfunc(mu, p1, p2, p3)*dm1**(-k1)*dm2**(-k2))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine m3_inc_type3

  ! tendency of m3 as a function of m1,m2,dm1,dm2
  ! assuming shape parameter does not vary
  subroutine m3_inc_type2(m1, m2, m3, p1, p2, p3, dm1, dm2, dm3, mu_init)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='M3_INC_TYPE2'
   
    real(wp), intent(in) :: m1, m2, m3, dm1, dm2
    real(wp), intent(in) :: p1, p2, p3
    real(wp), intent(out) :: dm3
    real(wp), intent(in), optional :: mu_init ! initial mu for type3
    
    real(wp) :: k1, k2, k3
    real(wp) :: fac

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (m1 > 0.0) then
      k1 = p3-p2
      k2 = p1-p3
      k3 = p2-p1

      if (abs(p3-6.0)< spacing(p3)) then ! p3=6
        ! since we've hardwired p1=3 and p2=0, we can reduce the calculation to...
        fac=((m1+dm1)/m1)*((m1+dm1)/m1)*(m2/(m2+dm2))-1.0
      else if (abs(p3-4.5)< spacing(p3)) then ! p3=6
        ! since we've hardwired p1=3 and p2=0, we can reduce the calculation to.
        fac=sqrt(((m1+dm1)/m1)*((m1+dm1)/m1)*((m1+dm1)/m1)*(m2/(m2+dm2)))-1.0
      else if (abs(p3-1.5)< spacing(p3)) then ! p3=6
        ! since we've hardwired p1=3 and p2=0, we can reduce the calculation to...
        fac=sqrt(((m1+dm1)/m1)*((m2+dm2)/m2))-1.0
      else
        !      fac=((m1/(m1+dm1))**(k1/k3)*(m2/(m2+dm2))**(k2/k3)-1.)
        fac=exp((k1*log(m1/(m1+dm1)) + k2*log(m2/(m2+dm2)))/k3)-1.0
      end if
      if (fac < -1.0) then

        write(std_msg, *) 'm3inc_2 ERROR:', fac, m1,m2,m3,dm1,dm2,k1,k2,k3, &
                            sqrt(fac)

        call throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)

      end if
      dm3 = m3*fac
    else ! If there is no pre-existing mass, then use type 3
      call m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu_init)
    end if

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine m3_inc_type2

  ! tendency of m2 as a function of m1,dm1
  ! assuming shape parameter and slope do not vary
  subroutine m3_inc_type1(m1, m2, dm1, dm2)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='M3_INC_TYPE1'

    real(wp), intent(in) :: m1, m2, dm1
    real(wp), intent(out) :: dm2

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dm2 = dm1*m1/m2

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine m3_inc_type1
end module m3_incs
