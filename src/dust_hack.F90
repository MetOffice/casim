module dust_hack

! No ModuleName required as no subroutines are present

  implicit none
  ! Hacks to use 3rd moments as aerosol fields
  ! These are plumbed through in mphys_casim_um
  !                      0-2km        2-4km          4-6km
  ! dust /cc               0.0466311    0.0506125  0.000694496
  ! soluble /cc              543.755      176.647      185.030
  ! soluble mass kg/m3   1.27056e-09  5.96604e-10  1.01946e-10

  real :: n_dust_1=0.0466311*1.0e6
  real :: n_dust_2=0.0506125*1.0e6
  real :: n_dust_3=0.000694496*1.0e6

  real :: n_sol_1=543.755*1.0e6
  real :: n_sol_2=176.647*1.0e6
  real :: n_sol_3=185.030*1.0e6

  real :: m_sol_1=1.27056e-09
  real :: m_sol_2=5.96604e-10
  real :: m_sol_3=1.01946e-10

  real :: dust_factor = 1.0
  real :: sol_factor = 1.0

end module dust_hack
