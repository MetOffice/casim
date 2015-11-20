MODULE mphys_end

USE lookup, ONLY: mu_i, mu_g, mu_i_sed, mu_g_sed
USE passive_fields, ONLY:   &
    rho, pressure, exner, dz, w, &
    tke, rdz_on_rho, qws, TdegK

IMPLICIT NONE

CONTAINS

SUBROUTINE cleanup

!     ! Allocated in mphys_init
!     deallocate(mu_i)
!     deallocate(mu_g)
!     deallocate(mu_i_sed)
!     deallocate(mu_g_sed)

!    ! Allocated in set_passive_fields
!    deallocate(rho)
!    deallocate(pressure)
!    deallocate(exner)
!    deallocate(dz)
!    deallocate(w)
!    deallocate(tke)
!    deallocate(rdz_on_rho)
!    deallocate(qws)
!    deallocate(TdegK)
!    deallocate(qws0)
!    deallocate(TdegC)

END SUBROUTINE cleanup

END MODULE mphys_end
