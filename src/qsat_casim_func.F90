module qsat_funs
  use variable_precision, only: wp
  implicit none
  private

  public Qsaturation, Qisaturation, dqwsatdt

contains
  ! Function to return the saturation mr over water
  ! Based on tetans formular
  ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
  function Qsaturation (T, p)

    implicit none

    real(wp), intent(IN) :: T, p
    real(wp) :: Qsaturation
    ! Temperature in Kelvin
    ! Pressure in mb
    real(wp), parameter ::tk0c = 273.15, qsa1 = 3.8, qsa2 = - 17.2693882, qsa3 = 35.86, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------

    if (T > qsa3 .and. p * exp (qsa2 * (t - tk0c) / (T - qsa3)) > qsa4) then
      Qsaturation=qsa1/(p*exp(qsa2*(t-tk0c)/(T-qsa3))-qsa4)
    else
      qsaturation=999.0
    end if

  end function Qsaturation

  ! Function to return the saturation mr over ice
  ! Based on tetans formular
  ! QS=3.8/(P*EXP(-21.8745584*(T-273.15)/(T-7.66))-6.109)
  function Qisaturation(T, p)

    implicit none

    real(wp), intent(IN) ::  T, p
    real(wp) :: qisaturation
    ! Temperature in Kelvin
    ! Pressure in mb
    real(wp), parameter :: tk0c = 273.15, qis1 = 3.8, qis2 = -21.8745584 , qis3 = 7.66  , qis4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qisat equation
    ! Constant in qisat equation
    ! Constant in qisat equation

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------

    if (T > qis3 .and. p * exp (qis2 * (t - tk0c) / (T - qis3)) > qis4) then
       qisaturation = qis1/(p*exp(qis2*(T - tk0c)/(T - qis3)) - qis4)
    else
       qisaturation = 999.0
    end if

  end function Qisaturation

  ! Function to return the rate of change with temperature
  ! of saturation mixing ratio over liquid water.
  !
  ! Based on tetans formular
  ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
  function dqwsatdt (qsat, T)

    implicit none

    real(wp) , intent(IN) ::qsat, T
    real(wp) ::dqwsatdt
    ! Saturatio mixing ratio
    ! Temperature in Kelvin

    real(wp), parameter ::tk0c = 273.15, qsa1 = 3.8, qsa2 = - 17.2693882, qsa3 = 35.86, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------

    dqwsatdt=-qsa2*(TK0C-qsa3)*(1.0+qsa4*qsat/qsa1)*qsat*(T-qsa3)**(- 2.0)

  end function dqwsatdt
end module qsat_funs
