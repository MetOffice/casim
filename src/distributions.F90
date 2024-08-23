module distributions
  use variable_precision, only: wp
  use mphys_parameters, only: hydro_params, nz, a_s, b_s
! use mphys_parameters, only: ice_params
  use mphys_switches, only: i_th, l_limit_psd, &
                            max_mu, max_mu_frac, l_kfsm, l_prf_cfrac,    &
                            i_cfl, i_cfr, i_cfi, i_cfs, i_cfg,           &
                            l_adjust_D0
  use lookup, only: get_slope_generic, get_slope_generic_kf
! use lookup, only: moment, get_lam_n0
  use thresholds, only: cfliq_small
! use thresholds, only: thresh_large
  use mphys_die, only: throw_mphys_error, bad_values, std_msg
! use mphys_die, only: warn
! use m3_incs, only: m3_inc_type3
  use passive_fields, only: exner

  ! VERBOSE = 1 for verbose print statements, otherwise dont display
#define VERBOSE 0

  implicit none
  private

  character(len=*), parameter, private :: ModuleName='DISTRIBUTIONS'

  real(wp), allocatable :: dist_lambda(:,:), dist_mu(:,:), dist_n0(:,:)
  real(wp), allocatable :: dist_lams(:,:,:)
  real(wp), dimension(:), allocatable :: m1,m2, m3, m3_old

!$OMP THREADPRIVATE(dist_lambda, dist_mu, dist_n0, dist_lams, m1,m2,m3,m3_old)

  public query_distributions, initialise_distributions, dist_lambda, dist_mu, dist_n0, dist_lams
contains

  subroutine initialise_distributions(nz, nspecies)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: nz, nspecies

    ! Local variables
    character(len=*), parameter :: RoutineName='INITIALISE_DISTRIBUTIONS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    allocate(dist_lambda(nz,nspecies), dist_mu(nz,nspecies), dist_n0(nz,nspecies), &
             dist_lams(nz, nspecies, 2), m1(nz), m2(nz), m3(nz), m3_old(nz),       &
             a_s(nz), b_s(nz))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine initialise_distributions

  ! Any changes in number should be applied to the prognostic variable
  ! rather than just these parameters.  Currently this is not done.
  subroutine query_distributions(ixy_inner, params, qfields, cffields)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    
    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    type(hydro_params), intent(inout) :: params !< species parameters
    real(wp), intent(inout) :: qfields(:,:)
    real(wp), intent(in) :: cffields(:,:)

    ! Local variables
    integer :: k    
    integer(wp) :: i1,i2,i3,ispec 
    real(wp) :: alpha, D0, mu_maxes_calc, Tk
!   real(wp) :: mu_pass=1.0
#if VERBOSE==1
    real(wp) :: n0_old, mu_old, lam_old
#endif

    real(wp) :: cf , q1_in, m2_in
    integer :: i_cf


    character(len=*), parameter :: RoutineName='QUERY_DISTRIBUTIONS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    mu_maxes_calc=max_mu_frac*max_mu

    ispec=params%id
    i1=params%i_1m
    i2=params%i_2m
    i3=params%i_3m

    dist_lambda(:,ispec)=0.0
    dist_mu(:,ispec)=0.0
    dist_n0(:,ispec)=0.0
    dist_lams(:,ispec,1) = 0.0
    dist_lams(:,ispec,2) = 0.0

    m1(:)=qfields(:, i1)/params%c_x                
    if (params%l_2m) m2(:)=qfields(:, i2)          
    ! if (params%l_3m) then
    !   m3(:)=qfields(:, i3)
    ! else
    !   m3=0.0
    ! end if
    !   m3_old=m3

    ! If 3rd moment has become too small then recalculate using max_mu
    ! This shouldn't happen except in small number situations - be careful
!     if (params%l_3m) then          
!       do k=1, nz
!         if (qfields(k, i1) .gt. 0.0) then
!           if (m3(k) < spacing(m3(k))) then
!             call m3_inc_type3(params%p1, params%p2, params%p3, m1(k), m2(k), m3(k), max_mu)
! #if VERBOSE==1                
!             write(std_msg,*) 'WARNING: resetting negative third moment',  &
!                                params%id, m3(k), m3_old(k), 'm1 and m2 are: ',  &
!                                qfields(k, i1)/params%c_x, m2(k)

!             call throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)
! #endif                
!           else if (m3(k) > thresh_large(params%i_3m)) then
!             call m3_inc_type3(params%p1, params%p2, params%p3, m1(k), m2(k), m3(k), 0.0_wp)
! #if VERBOSE==1                 
!             write(std_msg,*) 'WARNING: resetting large third moment',    &
!                                params%id, m3(k), m3_old(k), 'm1 and m2 are: ', &
!                                qfields(k, i1)/params%c_x, m2(k)

!             call throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)
! #endif                
!           end if            
!           qfields(k,i3)=m3(k)
!         end if
!       end do
!    end if

    if (l_prf_cfrac) then 
      select case (ispec)
      case (1) !cloud
        i_cf=i_cfl
      case (2) !rain
        i_cf=i_cfr
      case (3) !ice
        i_cf=i_cfi
      case (4) !snow
        i_cf=i_cfs
      case (5) !graupel
        i_cf=i_cfg
      case default
        STOP
      end select
    end if

    do k=1, nz
      if (l_prf_cfrac) then
        if (cffields(k,i_cf) .gt. cfliq_small) then
          cf=cffields(k,i_cf)
        else
          cf=cfliq_small !nonzero value - maybe move cf test higher up
        endif
      else
          cf=1.0
      endif

      q1_in=qfields(k, i1) / cf  !incloud mass                 
      m2_in=m2(k)          / cf  !incloud number
    
    
      if (qfields(k, i1) > 0.0) then
        if (l_kfsm) then
          Tk = qfields(k, i_th) * exner(k,ixy_inner)
          call get_slope_generic_kf(ixy_inner, k, params, dist_n0(k,ispec),     &
                                    dist_lambda(k,ispec), dist_mu(k,ispec),     &
                                    dist_lams(k,ispec,:), q1_in, Tk,   &
                                    m2_in, m3(k))

        else
          call get_slope_generic(params, dist_n0(k,ispec),                      &
                                 dist_lambda(k,ispec), dist_mu(k,ispec),        &
                                 q1_in, m2_in, m3(k))

        end if ! l_kfsm

      end if ! qfields > 0

      !KF.< moved back to lookup.F90
      !!if ( l_kfsm .and. ispec == ice_params%id) then
      !!  a_s(k) = params%a_x
      !!  b_s(k) = params%b_x
      !!end if
      !KF.>

    end do

!    ! If we diagnose a mu out of bounds, then reset m3
!     if (params%l_3m) then
!       do k=1, nz
!         if (qfields(k, i1) .gt. 0.0) then
!           if (dist_mu(k,ispec)<0.0 .or. mu_pass < 0.0) then
!             mu_old=dist_mu(k,ispec)
!             dist_mu(k,ispec)=0.0
!             call get_lam_n0(m1(k), m2(k), m3(k), params, dist_mu(k,ispec), dist_lambda(k,ispec), &
!                  dist_n0(k,ispec))
!             m3(k)=moment(dist_n0(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
!             qfields(k,i3)=m3(k)
! #if VERBOSE==1
!             write(std_msg, *) 'WARNING: resetting negative mu',  mu_old, m1(k), m2(k), m3(k), m3_old(k)
!             call throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)
! #endif
!           else if (.not. l_limit_psd .and. &
!                (dist_mu(k,ispec) + epsilon(1.0) > max_mu .or. mu_pass +epsilon(1.0) > max_mu)) then
!             mu_old=dist_mu(k,ispec)
!             dist_mu(k,ispec)=max_mu
!             call get_lam_n0(m1(k), m2(k), m3(k), params, dist_mu(k,ispec), dist_lambda(k,ispec), &
!                  dist_n0(k,ispec))
!             m3(k)=moment(dist_mu(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
!             qfields(k,i3)=m3(k)
! #if VERBOSE==1
!             write(std_msg, *) 'WARNING: resetting large mu',  mu_old, m1(k), m2(k), m3(k), m3_old(k)
!             call throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)
! #endif
!           end if
!         end if
!       end do
!     end if

    if (l_limit_psd .and. params%l_2m) then

 !      if (params%l_3m) then
!         do k=1, nz
!           if (qfields(k, i1) .gt. 0.0) then
!             if (dist_mu(k,ispec) > mu_maxes_calc) then
!               !-----------------------
!               ! Adjust mu/m3 necessary
!               !-----------------------              
!               dist_mu(k,ispec)=(dist_mu(k,ispec) + mu_maxes_calc)*0.5
!               call get_lam_n0(m1(k), m2(k), m3(k), params, dist_mu(k,ispec), dist_lambda(k,ispec), &
!                    dist_n0(k,ispec))
!               m3_old(k)=m3(k)
!               m3(k)=moment(dist_n0(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
!               qfields(k,i3) = m3(k)
! #if VERBOSE==1
!                write(std_msg, *) 'WARNING: adjusting m3 with large mu',  params%id, dist_mu(k,ispec), mu_old, &
!                                    m1(k), m2(k), m3(k), m3_old(k)
!                call throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)
! #endif
!             end if
!           end if
!         end do
!       end if

      !-----------------------
      ! Adjust D0 if necessary
      ! only conserves mass - not number
      !-----------------------
      if (l_adjust_D0) then  
      ! Adjust the PSD to limit to the Dmax hydrometeor, default is true
        !D0 = (m1/m2)**(1./(params%p1-params%p2))
         do k=1, nz
            if (qfields(k, i1) .gt. 0.0) then
               D0=(1+dist_mu(k,ispec))/dist_lambda(k,ispec)
               if (D0 > params%Dmax) then
#if VERBOSE==1
                  mu_old=dist_mu(k,ispec)
                  lam_old=dist_lambda(k,ispec)
                  n0_old=dist_n0(k,ispec)
#endif
                  alpha=D0/params%Dmax
                  dist_lambda(k,ispec)=alpha*dist_lambda(k,ispec)
                  dist_n0(k,ispec)=alpha**(params%p1)*dist_n0(k,ispec)
                  
                  qfields(k,i2)=dist_n0(k,ispec)*cf  !grid mean
                  !             if (params%l_3m) then
                  !               m3(k)=moment(dist_n0(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
                  !               qfields(k,i3)=m3(k)
                  !             end if
                  ! #if VERBOSE==1
                  !             write(std_msg,*) 'WARNING: adjusting number and m3',  params%id, n0_old,     &
                  !                               dist_n0(k,ispec), m3_old(k), m3(k), 'new m1, m2, m3 are: ', &
                  !                               m1(k), m2(k), m3(k)
                  !             call throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)
                  ! #endif
               end if
               if (D0 < params%Dmin) then
#if VERBOSE==1
                  mu_old=dist_mu(k,ispec)
                  lam_old=dist_lambda(k,ispec)
                  n0_old=dist_n0(k,ispec)
#endif
                  alpha=D0/params%Dmin
                  dist_lambda(k,ispec)=alpha*dist_lambda(k,ispec)
                  dist_n0(k,ispec)=alpha**(params%p1)*dist_n0(k,ispec)
                  
                  qfields(k,i2)=dist_n0(k,ispec)  *cf  !grid mean
                  !             if (params%l_3m) then
                  !                m3(k)=moment(dist_n0(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
                  !                qfields(k,i3)=m3(k)
                  !             end if
                  ! #if VERBOSE==1
                  !             write(std_msg,*) 'WARNING: adjusting number and m3',  params%id, n0_old,      &
                  !                               dist_n0(k,ispec), m3_old(k), m3(k), 'new m1, m2, m3 are: ',  &
                  !                               m1(k), m2(k), m3(k)
                  !             call throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)
                  ! #endif
                  
               end if
            end if
         end do
      endif ! end adjust PSD
      
   end if
      




    ! Final check that distributions do actually make sense.
   do k = 1, nz

      if ( m1(k) > 0.0 ) then
         if (params % l_1m .and. dist_lambda(k, ispec) <= 0.0) then
            write(std_msg, '(A,F7.2)')'Unexpected zero or negative lambda: lambda =', dist_lambda(k, ispec)
            call throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)
         end if
         
         if (params % l_2m .and. dist_n0(k, ispec) <= 0.0) then
            write(std_msg, '(A,F7.2)')'Unexpected zero or negative n0: n0 =', dist_n0(k, ispec)
            call throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)
         end if
         
         ! if (params % l_3m .and. dist_mu(k, ispec) < 0.0) then
         !   write(std_msg, '(A,F7.2)')'Unexpected zero or negative mu: mu =', dist_mu(k, ispec)
         !   call throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)
         ! end if
      end if ! m1(k) > 0
      
   end do ! k

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

 end subroutine query_distributions

end module distributions
