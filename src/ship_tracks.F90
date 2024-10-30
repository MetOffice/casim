MODULE ship_tracks_casim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
ModuleName = 'SHIP_TRACKS_CASIM'

REAL, PARAMETER    :: RMDI     = -32768.0*32768.0 
! Borrowed from UM mod for now - consider whether we should 
! include this routine or not for the future.

  ! Some variables used for ship tracks
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

CONTAINS

SUBROUTINE include_ship_source(i_start, i_end, j_start, j_end, k_start, k_end, timestep, &
                               timestep_number, height,          &
                               r_theta_levels,                   &
                               true_latitude, true_longitude,    &
                               dAccumSolMass, dAccumSolNumber )

USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
USE mphys_constants,  ONLY: pi

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_start, i_end, j_start, j_end, k_start, k_end

REAL, INTENT(IN) :: timestep
INTEGER, INTENT(IN) :: timestep_number

REAL, INTENT(IN) :: height(k_start:k_end, i_start:i_end, j_start:j_end)
REAL, INTENT(IN) :: r_theta_levels(i_start:i_end, j_start:j_end, 0:k_end)
REAL, INTENT(IN) :: true_latitude(i_start:i_end, j_start:j_end)
REAL, INTENT(IN) :: true_longitude(i_start:i_end,j_start:j_end)

REAL, INTENT(INOUT) :: dAccumSolMass(k_start:k_end, i_start:i_end, j_start:j_end)
REAL, INTENT(INOUT) :: dAccumSolNumber(k_start:k_end, i_start:i_end, j_start:j_end)

! Local Variables

CHARACTER(LEN=*), PARAMETER :: RoutineName='INCLUDE_SHIP_SOURCE'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


REAL :: dttime, rspeed, theta, true_lat, true_lon

REAL :: xarray1(k_start:k_end)
REAL :: xarray2(i_start:i_end,j_start:j_end)

INTEGER :: ijval(2)
INTEGER :: i, j, kval(1), k

REAL :: dlat, dlon, dr3, boxsize, rrho

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

dttime = timestep * timestep_number

DO iship = 1, nships

  IF (ship_lat_start(iship) > rmdi) THEN ! check data has been 
                                         ! entered for this ship

    rspeed = 2.0*pi*0.51444*ship_speed(iship)/40075000.0 

    theta  = ASIN((ship_lat_end(iship)-ship_lat_start(iship)) / &
                  (ship_lon_end(iship)-ship_lon_start(iship)))

    true_lat = ship_lat_start(iship)/360.0 * 2 * pi + &
               dttime * rspeed * SIN(theta)

    true_lon = ship_lon_start(iship)/360.0 * 2 * pi + &
               dttime*rspeed*COS(theta)

    true_lon = MOD(true_lon+2*pi, 2*pi)

    xarray2 = (true_latitude(:,:)  - true_lat) &
            * (true_latitude(:,:)  - true_lat) &
            + (true_longitude(:,:) - true_lon) &
            * (true_longitude(:,:) - true_lon)

    ijval = MINLOC(xarray2)
    i     = ijval(1)
    j     = ijval(2)

! The if test determines if the nearest point is within 2
! grid cells of the edge
! AND to prevent duplication across pes requires i>1 and j>1

    IF ( MINVAL(xarray2) < 0.25 * (true_latitude(i,j) - &
         true_latitude(i,j-1)) * (true_longitude(i,j) - &
         true_longitude(i-1,j)) .AND. i>1 .AND. j> 1 ) THEN

      xarray1(:) = (height(:,i,j) - ship_z) ** 2

      kval = MINLOC(xarray1)
      ! Ensure not model level 1 or top of model
      k    = MAX(2, MIN(kval(1), k_end-1) )

      IF ( ship_dmdt < 0.0 ) ship_dmdt = ship_dndt * default_mean_mass
      IF ( ship_dndt < 0.0 ) ship_dndt = ship_dmdt / default_mean_mass

      dlat = true_latitude(i,j)  - true_latitude(i,j-1)
      dlon = true_longitude(i,j) - true_longitude(i-1,j)
      dr3  = r_theta_levels(i,j,k+1)**3 - r_theta_levels(i,j,k)**3

! Assume rotated pole puts use near the equator,
! i.e. sin(theta)=1
! r*r*sin(theta)*dlat*pi/180*dlon*pi/180

      boxsize = r_theta_levels(i,j,k)**2 * dlat * dlon * pi * pi/32400.0

      dAccumSolMass(k,i,j)   = dAccumSolMass(k,i,j) + ship_dmdt * rrho * &
                               timestep/boxsize

      dAccumSolNumber(k,i,j) = dAccumSolNumber(k,i,j) + ship_dndt * rrho * &
                               timestep/boxsize

    END IF ! Minval(xarray)

  END IF ! ship_lat_start

END DO ! iship


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE include_ship_source

END MODULE ship_tracks_casim
