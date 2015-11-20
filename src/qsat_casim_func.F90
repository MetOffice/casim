MODULE qsat_funs

USE variable_precision, ONLY: wp

IMPLICIT NONE

CONTAINS

FUNCTION Qsaturation (T, p)

    !
    ! Function to return the saturation mr over water
    ! Based on tetans formular
    ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
    !
REAL(wp) :: Qsaturation
REAL(wp), INTENT(IN) :: T, p
    ! Temperature in Kelvin
    ! Pressure in mb
REAL(wp), PARAMETER ::tk0c = 273.15, qsa1 = 3.8, qsa2 = - 17.2693882,    &
         qsa3 = 35.86, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !
IF (T > qsa3 .AND. p * EXP (qsa2 * (t - tk0c) / (T - qsa3) )     &
            > qsa4) THEN
  Qsaturation = qsa1 / (p * EXP (qsa2 * (t - tk0c) / (T - qsa3) )      &
            - qsa4)
ELSE
  qsaturation=999.0
END IF

END FUNCTION Qsaturation

FUNCTION Qisaturation(T, p)
    !
    ! Function to return the saturation mr over ice
    ! Based on tetans formular
    ! QS=3.8/(P*EXP(-21.8745584*(T-273.15)/(T-7.66))-6.109)
    !
REAL(wp) :: qisaturation
REAL(wp), INTENT(IN) ::  T, p
    ! Temperature in Kelvin
    ! Pressure in mb
REAL(wp), PARAMETER :: tk0c = 273.15, qis1 = 3.8, qis2 = -21.8745584     &
       , qis3 = 7.66  , qis4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qisat equation
    ! Constant in qisat equation
    ! Constant in qisat equation

qisaturation = qis1/(p*EXP(qis2*(T - tk0c)/(T - qis3)) - qis4)

END FUNCTION Qisaturation

FUNCTION dqwsatdt (qsat, T)
IMPLICIT NONE
    !
    ! Function to return the rate of change with temperature
    ! of saturation mixing ratio over liquid water.
    !
    ! Based on tetans formular
    ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
    !
REAL(wp) ::dqwsatdt

REAL(wp) , INTENT(IN) ::qsat, T
    ! Saturatio mixing ratio
    ! Temperature in Kelvin

REAL(wp), PARAMETER ::tk0c = 273.15, qsa1 = 3.8, qsa2 = - 17.2693882,    &
         qsa3 = 35.86, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !
dqwsatdt = - qsa2 * (TK0C - qsa3) * (1.0 + qsa4 * qsat / qsa1)     &
         * qsat * (T - qsa3) ** ( - 2.0)
    !
    !

END FUNCTION dqwsatdt

END MODULE qsat_funs
