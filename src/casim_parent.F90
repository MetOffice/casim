! Module to contain CASIM parent information

MODULE casim_parent_mod

IMPLICIT NONE

! No RoutineName as no subroutines or functions in this module

! Parameters for each model CASIM is driven from
INTEGER, PARAMETER :: parent_unset = 0
INTEGER, PARAMETER :: parent_kid   = 1
INTEGER, PARAMETER :: parent_monc  = 2
INTEGER, PARAMETER :: parent_um    = 3

! Variable containing casim parent - initialised
! as unset for now, but should be set by each
! parent model
INTEGER :: casim_parent = parent_unset

END MODULE casim_parent_mod
