module ventfac
  use variable_precision, only: wp
  use passive_fields, only: rho
  use special, only: pi, GammaFunc
  use mphys_parameters, only: hydro_params, vent_1, vent_2
  use mphys_constants, only: visair, rho0, Dv

  implicit none
  private

  public ventilation
contains

  subroutine ventilation(k, V, n0, lam, mu, params)
    integer, intent(in) :: k
    real(wp), intent(in) :: n0, lam, mu
    real(wp), intent(out) :: V  ! bulk ventilation factor
    type(hydro_params), intent(in) :: params

    real(wp) :: T1, T2
    real(wp) :: Sc ! Schmidt number
    real(wp) :: a_x, b_x, f_x

    a_x=params%a_x
    b_x=params%b_x
    f_x=params%f_x

    Sc=visair/Dv

    T1=vent_1*(mu+1.0)/lam
    T2=vent_2*Sc**(1.0/3.0)*(a_x*(rho(k)/rho0)**(0.5)*rho(k)/visair)**(0.5)

    V=2.0*pi*n0/rho(k)*(T1+T2*GammaFunc(.5*b_x+mu+2.5)/GammaFunc(1.0+mu) &
         *(1.0 + .5*f_x/lam)**(-(.5*b_x + mu + 2.5))*lam**(-.5*b_x-1.5))
    ! for computational efficiency/accuracy changed from...
    ![' ']    !*(lam + .5*f_x)**(-(.5*b_x + mu + 2.5))*lam**(1.+mu)) &
  end subroutine ventilation
end module ventfac
