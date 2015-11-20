#if DEF_MODEL==MODEL_UM
MODULE mphys_um_update_fields
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

CONTAINS

SUBROUTINE update_mphys_fields

USE atm_fields_bounds_mod, ONLY: tdims
USE mphys_casim_switches, ONLY:                                        &
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
USE mphys_casim_prognostics, ONLY:                                     &
! prognostics for CASIM cloud and ice microphysics
       CloudNumber, RainNumber, Rain3mom,    &
       IceNumber, SnowNumber,                &
       Snow3mom, GraupNumber, Graup3mom,     &
       ! prognostics for activated aerosol (if used)
       ActiveSolLiquid, ActiveSolRain,       &
       ActiveInsolIce, ActiveSolIce,         &
       ActiveInsolLiquid, ActiveSolNumber,   &
       ActiveInsolNumber,                    &
! updated CASIM fields following advection and microphysics
       CloudNumber_star, RainNumber_star,    &
       Rain3mom_star, IceNumber_star,        &
       SnowNumber_star, Snow3mom_star,       &
       GraupNumber_star, Graup3mom_star,     &
       ActiveSolLiquid_star,                 &
       ActiveSolRain_star,                   &
       ActiveInsolIce_star,                  &
       ActiveSolIce_star,                    &
       ActiveInsolLiquid_star,               &
       ActiveSolNumber_star,                 &
       ActiveInsolNumber_star

IMPLICIT NONE

IF (l_mp_CloudNumber)                                                     &
       CloudNumber(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) =  &
       CloudNumber_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_RainNumber)                                                   &
       RainNumber(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) =   &
       RainNumber_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_Rain3mom)                                                  &
       Rain3mom(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) =     &
       Rain3mom_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_IceNumber)                                                    &
       IceNumber(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) =    &
       IceNumber_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_SnowNumber)                                                   &
       SnowNumber(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) =   &
       SnowNumber_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_Snow3mom)                                                  &
       Snow3mom(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) =     &
       Snow3mom_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_GraupNumber)                                                  &
       GraupNumber(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) =  &
       GraupNumber_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_Graup3mom)                                                 &
       Graup3mom(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) =    &
       Graup3mom_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
  ! Update the activate aerosol fields from CASIM
IF (l_mp_ActiveSolLiquid)                                         &
       ActiveSolLiquid(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) = &
       ActiveSolLiquid_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_ActiveSolRain)                                           &
       ActiveSolRain(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) =&
       ActiveSolRain_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_ActiveInSolIce)                                          &
       ActiveInSolIce(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) = &
       ActiveInSolIce_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_ActiveSolIce)                                            &
       ActiveSolIce(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) = &
       ActiveSolIce_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_ActiveInSolLiquid)                                       &
       ActiveInSolLiquid(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) = &
       ActiveInSolLiquid_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_ActiveInSolNumber)                                       &
       ActiveInSolNumber(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) = &
       ActiveInSolNUmber_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)
IF (l_mp_ActiveSolNumber)                                       &
       ActiveSolNumber(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :) = &
       ActiveSolNUmber_star(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, :)


END SUBROUTINE update_mphys_fields

END MODULE mphys_um_update_fields
#endif
