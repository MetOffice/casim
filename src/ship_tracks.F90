#if DEF_MODEL==MODEL_UM
MODULE ship_tracksCasim

USE missing_data_mod, ONLY: RMDI
IMPLICIT NONE

  ! Some variables used for ship tracks
LOGICAL :: l_shiptrack=.FALSE.  ! Logical to switch on the shiptrack
LOGICAL :: l_shipfromfile=.FALSE. ! read in ship location from a binary file
INTEGER :: nships = 1 ! how many ships
INTEGER, PARAMETER :: maxships = 10 ! maximum number of ships
REAL :: ship_z=20.0    ! exhaust height
REAL :: ship_lat(maxships) = RMDI  ! latitude
REAL :: ship_lon(maxships) = RMDI ! longitude

REAL :: ship_lat_start(maxships) = RMDI  ! start latitude
REAL :: ship_lon_start(maxships) = RMDI ! stop longitude
REAL :: ship_lat_end(maxships) = RMDI  ! start latitude
REAL :: ship_lon_end(maxships) = RMDI ! stop longitude
REAL :: ship_speed(maxships) = RMDI ! ship speed (knots)

REAL :: ship_dndt = RMDI ! rate of production number
REAL :: ship_dmdt = RMDI ! rate of production mass

REAL :: default_mean_mass = 3.0e-18 ! default mean aerosol mass

INTEGER :: iship ! ship counter

END MODULE ship_tracksCasim
#endif
