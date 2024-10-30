module shipway_erf

  use variable_precision, only: wp
  Use shipway_constants, only: pi

  Implicit None

  character(len=*), parameter, private :: ModuleName='SHIPWAY_ERF'

contains

  function erfg(x,c)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none
    
    real(wp), intent(in) :: x
    integer, intent(in) :: c ! 0 gives erf(x)
                             ! 1 gives erfc(x)
    real(wp) :: erfg, f, z
    integer :: j, cc

    real(wp) :: t, a1, a2, a3, a4, a5, p
    integer :: sign_x

    character(len=*), parameter :: RoutineName='ERFG'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    z=x
    cc=c

    if (abs(z) < 1.5)then 
       j=3 + int(9*abs(z))
       f=1
       do while(j /= 0)
          f=1.+f*z**2*(.5-j)/(j*(.5+j))
          j=j-1
       end do
       f=cc+f*z*(2.-4.*cc)/sqrt(pi)
    else
      ! Something wrong with this for x<~-2
      ! so use A&S instead
      ! cc=cc*sign_x(1.,z)
      ! j=3 + int(32/abs(z))
      ! f=0
      ! do while(j /= 0)
      !    f=1./(f*j+sqrt(2*z*z))
      !    j=j-1
      ! end do
      ! f=f*(cc*cc+cc-1.)*sqrt(2./pi)*exp(-z*z) + (1-cc)
      !---------------------
      ! Following from A&S 
      ! save the sign_x of x
      sign_x = 1
      if (x <= 0) sign_x = -1
      z = abs(x)

      !# constants
      a1 =  0.254829592
      a2 = -0.284496736
      a3 =  1.421413741
      a4 = -1.453152027
      a5 =  1.061405429
      p  =  0.3275911

      ! A&S formula 7.1.26
      t = 1.0/(1.0 + p*x)
      f = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
      if (c==1)f=1-f
    end if
    
    erfg=f

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function erfg

end module shipway_erf
