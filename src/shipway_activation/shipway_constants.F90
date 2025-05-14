module shipway_constants
  use variable_precision, only: wp
  use mphys_constants, only: pi

  implicit none

  character(len=*), parameter, private :: ModuleName='SHIPWAY_CONSTANTS'

  real(wp), parameter :: &
       Ru = 8.314472    & ! Universal gas constant
       ,Mw = 0.18015e-1 & ! Molecular weight of water 
       ,zetasa = 0.8e-1 & ! Surface tension at solution-air interface
       ,rhow = 1000.    & ! Density of water
       ,rho = 1           ! Density of  air 
  
  real(wp), parameter :: &
         eps = 1.608   &  !(Rv/Rd)
         ,Rd = 287.05   &
         ,Rv = 461.5    &
         ,Dv = 0.226e-4 &
         ,Lv = 0.2501e7 &
         ,cp = 1005.    &
         ,ka = 0.243e-1 

  real(wp) :: alpha_c = 0.05 ! set do default value which we may be changed elsewhere
                         ! kinetic parameter, expressing probability
                         ! of water vapour molecules being incorporated
                         ! into droplet upon collision + dissolution 
                         ! kinetics (see fountoukis & nenes 2007 and 
                         ! Asa-Awuku and Nenes,2007)
  
  real(wp) :: &
       Dp_big  & ! droplet size upper bound (m)
       ,Dp_low & ! droplet size lower bound (m)
       ,delDp    ! difference between upper and lower bound
  

contains

  real(wp) function Dv_mean(T, alpha_c)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(in):: T ! temperature (Kelvin)
    real(wp), intent(in):: alpha_c ! condensation coefficient
                               ! (Sassen & Dodd 1988 use 0.05)

    ! Local Variables

    real(wp) :: B ! a combination of variables

    character(len=*), parameter :: RoutineName='DV_MEAN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    Dp_big=5.e-6
    Dp_low=min(0.207683*alpha_c**(-0.33048)*1e-6,Dp_big-1e-15)
    delDp=Dp_big-Dp_low

    B=(2.*Dv/alpha_c)*sqrt(2*pi*Mw/(Ru*T))
    Dv_mean=Dv*(delDp-B*log((Dp_big+B)/(Dp_low+B)))/delDp

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function Dv_mean
       
end module shipway_constants
