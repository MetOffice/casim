#if  DEF_MODEL==MODEL_UM
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_um_store_progs

CONTAINS

SUBROUTINE eg_store_casim_progs

USE atm_fields_bounds_mod, ONLY : tdims,pdims

USE mphys_casim_switches, ONLY:                                          &
         l_mp_CloudNumber, l_mp_RainNumber,       &
         l_mp_Rain3mom, l_mp_IceNumber,           &
         l_mp_SnowNumber, l_mp_Snow3mom,          &
         l_mp_GraupNumber, l_mp_Graup3mom,        &
       ! Switches for the activated aerosols
         l_mp_ActiveSolLiquid, l_mp_ActiveSolRain,&
         l_mp_ActiveInsolIce, l_mp_ActiveSolIce,  &
         l_mp_ActiveInsolLiquid,                  &
         l_mp_ActiveSolNumber,                    &
         l_mp_ActiveInSolNumber
USE mphys_casim_prognostics, ONLY:                                       &
         ! prognostics for CASIM cloud and ice microphysics
         CloudNumber, RainNumber, Rain3mom,    &
         IceNumber, SnowNumber,                &
         Snow3mom, GraupNumber, Graup3mom,     &
         ! prognostic increments from CASIM microphysics
         CloudNumber_inc, RainNumber_inc,      &
         Rain3mom_inc, IceNumber_inc,          &
         SnowNumber_inc, Snow3mom_inc,         &
         GraupNumber_inc, Graup3mom_inc,  &
         ! ENDGAME prognostic store values from casim
         CloudNumber_eg_store, RainNumber_eg_store,      &
         Rain3mom_eg_store, IceNumber_eg_store,          &
         SnowNumber_eg_store, Snow3mom_eg_store,         &
         GraupNumber_eg_store, Graup3mom_eg_store,  &
         ! prognostics for activated aerosol (if used)
         ActiveSolLiquid, ActiveSolRain,       &
         ActiveInsolIce, ActiveSolIce,         &
         ActiveInsolLiquid, ActiveSolNumber,   &
         ActiveInsolNumber,                    &
         ! prognostic increments for activated aerosol
         dActiveSolLiquid, dActiveSolRain,          &
         dActiveInsolIce, dActiveSolIce,            &
         dActiveInsolLiquid, dActiveSolNumber,      &
         dActiveInsolNumber,                        &
         ! ENDGAME prognostic store values from casim (in-cloud aersol)
         ActiveSolLiquid_eg_store, ActiveSolRain_eg_store,       &
         ActiveInsolIce_eg_store, ActiveSolIce_eg_store,         &
         ActiveInsolLiquid_eg_store, ActiveSolNumber_eg_store,   &
         ActiveInsolNumber_eg_store

IMPLICIT NONE

INTEGER :: i, j, k

    ! Store the end of atmosphys1 value for iteration in ENDGAME.
    ! Zero the increments to prevent double counting.
    ! NOTE: loops are not required, I think an array operation will work
    !
IF (l_mp_CloudNumber) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        CloudNumber_eg_store(i,j,k) = CloudNumber(i,j,k) +           &
                 CloudNumber_inc(i,j,k)
        CloudNumber_inc(i,j,k)=0.0
      END DO
    END DO
  END DO
END IF !l_mp_CloudNumber

IF (l_mp_RainNumber) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        RainNumber_eg_store(i,j,k) = RainNumber(i,j,k) +        &
                 RainNumber_inc(i,j,k)
        RainNumber_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_RainNumber

IF (l_mp_Rain3mom) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        Rain3mom_eg_store(i,j,k) = Rain3mom(i,j,k) +          &
                    Rain3mom_inc(i,j,k)
        Rain3mom_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_Rain3mom

IF (l_mp_IceNumber) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IceNumber_eg_store(i,j,k) = IceNumber(i,j,k) +         &
                    IceNumber_inc(i,j,k)
        IceNumber_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_IceNumber

IF (l_mp_SnowNumber) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        SnowNumber_eg_store(i,j,k) = SnowNumber(i,j,k) +        &
                    SnowNumber_inc(i,j,k)
        SnowNumber_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_SnowNumber

IF (l_mp_Snow3mom) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        Snow3mom_eg_store(i,j,k) = Snow3mom(i,j,k) +          &
                    Snow3mom_inc(i,j,k)
        Snow3mom_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_Snow3mom

IF (l_mp_GraupNumber) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        GraupNumber_eg_store(i,j,k) = GraupNumber(i,j,k) +       &
                    GraupNumber_inc(i,j,k)
        GraupNumber_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_GraupNumber

IF (l_mp_Graup3mom) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        Graup3mom_eg_store(i,j,k) = Graup3mom(i,j,k) +         &
                    Graup3mom_inc(i,j,k)
        Graup3mom_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_Graup3mom
!
! End of the CASIM extra cloud and ice moment prognostics
! Begin adding the activated aerosol to the advection super array
!
IF (l_mp_ActiveSolLiquid) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ActiveSolLiquid_eg_store(i,j,k) = ActiveSolLiquid(i,j,k) +        &
                    dActiveSolLiquid(i,j,k)
        dActiveSolLiquid(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_ActiveSolLiquid

IF (l_mp_ActiveSolRain) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ActiveSolRain_eg_store(i,j,k) = ActiveSolRain(i,j,k) +          &
                    dActiveSolRain(i,j,k)
        dActiveSolRain(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_ActiveSolRain

IF (l_mp_ActiveInSolIce) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ActiveInSolIce_eg_store(i,j,k) = ActiveInSolIce(i,j,k) +         &
                    dActiveInSolIce(i,j,k)
        dActiveInSolIce(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_ActiveInSolIce

IF (l_mp_ActiveSolIce) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ActiveSolIce(i,j,k) = ActiveSolIce(i,j,k) +           &
                    dActiveSolIce(i,j,k)
        dActiveSolIce(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_ActiveSolIce

IF (l_mp_ActiveInSolLiquid) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ActiveInSolLiquid_eg_store(i,j,k) = ActiveInSolLiquid(i,j,k) +         &
                    dActiveInSolLiquid(i,j,k)
        dActiveInSolLiquid(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_ActiveInSolLiquid

IF (l_mp_ActiveSolNumber) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ActiveSolNumber_eg_store(i,j,k) = ActiveSolNumber(i,j,k) +         &
                    dActiveSolNumber(i,j,k)
        dActiveSolNumber(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_ActiveSolNumber

IF (l_mp_ActiveInSolNumber) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ActiveInSolNumber_eg_store(i,j,k) = ActiveInSolNumber(i,j,k) +         &
                 dActiveInSolNumber(i,j,k)
        dActiveInSolNumber(i,j,k) = 0.0
      END DO
    END DO
  END DO
END IF !l_mp_ActiveInSolNumber

END SUBROUTINE eg_store_casim_progs


END MODULE EG_UM_STORE_PROGS

#endif
