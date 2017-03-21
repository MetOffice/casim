module distributions
  use variable_precision, only: wp
  use mphys_parameters, only: hydro_params, nz, rain_params
  use mphys_switches, only: hydro_names, l_limit_psd, l_passive3m, max_mu, max_mu_frac
  use lookup, only: get_slope_generic, get_n0, moment, get_mu, get_lam_n0
  use special, only: GammaFunc
  use thresholds, only: thresh_tidy, thresh_sig, thresh_large

  use m3_incs, only: m3_inc_type3

  ! VERBOSE = 1 for verbose print statements, otherwise dont display
#define VERBOSE 0

#if DEF_MODEL==MODEL_KiD
  ! Kid modules
  use diagnostics, only: save_dg, i_dgtime, i_here, j_here, k_here
  use runtime, only: time
  use parameters, only: nx
#elif  DEF_MODEL==MODEL_MONC
  use diaghelp_monc, only: i_here, j_here, time
#endif
  implicit none
  private

  real(wp), allocatable :: dist_lambda(:,:), dist_mu(:,:), dist_n0(:,:)
  real(wp), dimension(:), allocatable :: m1,m2, m3, m3_old

  public query_distributions, initialise_distributions, dist_lambda, dist_mu, dist_n0
contains

  subroutine initialise_distributions(nz, nspecies)
    integer, intent(in) :: nz, nspecies
    
    allocate(dist_lambda(nz,nspecies), dist_mu(nz,nspecies), dist_n0(nz,nspecies), m1(nz), m2(nz), m3(nz), m3_old(nz))
  end subroutine initialise_distributions  

  ! Any changes in number should be applied to the prognostic variable
  ! rather than just these parameters.  Currently this is not done.
  subroutine query_distributions(params, qfields, icall)
    type(hydro_params), intent(in) :: params !< species parameters
    real(wp), intent(inout) :: qfields(:,:)
    integer, intent(in), optional :: icall

    integer :: k    
    integer(wp) :: i1,i2,i3,ispec 
    real(wp) :: n0_old, mu_old, lam_old, alpha, D0, mu_pass=1.0, k1,k2,k3, mu_maxes_calc

    character(2) :: chcall

    ! Some diagnostic strings
    if (present(icall)) then
      write(chcall,'(A1,i1)') '_', icall
    else
      chcall=''
    end if

    mu_maxes_calc=max_mu_frac*max_mu

    ispec=params%id
    i1=params%i_1m
    i2=params%i_2m
    i3=params%i_3m

    !do k=1, nz 
     
      dist_lambda(:,ispec)=0.0
      dist_mu(:,ispec)=0.0
      dist_n0(:,ispec)=0.0

      !if (qfields(k, i1) .gt. 0.0) then
        m1(:)=qfields(:, i1)/params%c_x
        if (params%l_2m) m2(:)=qfields(:, i2)
        if (params%l_3m) then
          m3(:)=qfields(:, i3)
        else
          m3=0.0
        end if
        m3_old=m3

        ! If 3rd moment has become too small then recalculate using max_mu
        ! This shouldn't happen except in small number situations - be careful
        if (params%l_3m) then          
          do k=1, nz
            if (qfields(k, i1) .gt. 0.0) then
              if (m3(k) < spacing(m3(k))) then
                call m3_inc_type3(params%p1, params%p2, params%p3, m1(k), m2(k), m3(k), max_mu)
#if VERBOSE==1                
                print*, 'WARNING: resetting negative third moment', time, params%id, m3(k), m3_old(k)
                print*, 'm1 and m2 are: ', qfields(k, i1)/params%c_x, m2(k)
#endif                
              else if (m3(k) > thresh_large(params%i_3m)) then
                call m3_inc_type3(params%p1, params%p2, params%p3, m1(k), m2(k), m3(k), 0.0_wp)
#if VERBOSE==1                 
                print*, 'WARNING: resetting large third moment', time, params%id, m3(k), m3_old(k)
                print*, 'm1 and m2 are: ', qfields(k, i1)/params%c_x, m2(k)
#endif                
              end if            
              qfields(k,i3)=m3(k)
          end if
        end do
      end if

#if DEF_MODEL==MODEL_KiD
        ! as a test calculate what mu should be
        if (l_passive3m) then
          do k=1, nz
            if (qfields(k, i1) .gt. 0.0) then
              call get_mu(m1(k), m2(k), m3(k), params%p1, params%p2, params%p3, mu_pass)
              if (nx>1) then
                call save_dg(k, i_here, mu_pass, 'save_mu_'//trim(hydro_names(params%i_1m)), i_dgtime)
              else
                call save_dg(k, mu_pass, 'save_mu_'//trim(hydro_names(params%i_1m)), i_dgtime)
              end if
            end if
          end do
        end if
#endif
        do k=1, nz
          if (qfields(k, i1) .gt. 0.0) then
            call get_slope_generic(k, params, dist_n0(k,ispec), dist_lambda(k,ispec), dist_mu(k,ispec), &
                 qfields(k, i1), m2(k), m3(k))
          end if
        end do

        ! If we diagnose a mu out of bounds, then reset m3
        if (params%l_3m) then
          do k=1, nz
            if (qfields(k, i1) .gt. 0.0) then
              if (dist_mu(k,ispec)<0.0 .or. mu_pass < 0.0) then
                mu_old=dist_mu(k,ispec)
                dist_mu(k,ispec)=0.0
                call get_lam_n0(m1(k), m2(k), params%p1, params%p2, dist_mu(k,ispec), dist_lambda(k,ispec), dist_n0(k,ispec))
                m3(k)=moment(dist_n0(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
                qfields(k,i3)=m3(k)
#if VERBOSE==1
                print*, 'WARNING: resetting negative mu', time, mu_old, m1(k), m2(k), m3(k), m3_old(k)
#endif
              else if (.not. l_limit_psd .and. (dist_mu(k,ispec) + epsilon(1.0) > max_mu .or. mu_pass +epsilon(1.0) > max_mu)) then
                mu_old=dist_mu(k,ispec)
                dist_mu(k,ispec)=max_mu
                call get_lam_n0(m1(k), m2(k), params%p1, params%p2, dist_mu(k,ispec), dist_lambda(k,ispec), dist_n0(k,ispec))
                m3(k)=moment(dist_mu(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
                qfields(k,i3)=m3(k)
#if VERBOSE==1
                print*, 'WARNING: resetting large mu', time, mu_old, m1(k), m2(k), m3(k), m3_old(k)
#endif
              end if
            end if
          end do
        end if

        if (l_limit_psd .and. params%l_2m) then
          if (params%l_3m) then
            do k=1, nz
              if (qfields(k, i1) .gt. 0.0) then
                if (dist_mu(k,ispec) > mu_maxes_calc) then
                  !-----------------------
                  ! Adjust mu/m3 necessary
                  !-----------------------              
                  dist_mu(k,ispec)=(dist_mu(k,ispec) + mu_maxes_calc)*0.5
                  call get_lam_n0(m1(k), m2(k), params%p1, params%p2, dist_mu(k,ispec), dist_lambda(k,ispec), dist_n0(k,ispec))
                  m3_old(k)=m3(k)
                  m3(k)=moment(dist_n0(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
                  qfields(k,i3) = m3(k)
#if VERBOSE==1
                  print*, 'WARNING: adjusting m3 with large mu', time, params%id, dist_mu(k,ispec), mu_old, &
                       m1(k), m2(k), m3(k), m3_old(k)
#endif
                end if
              end if
            end do
          end if

          !-----------------------
          ! Adjust D0 if necessary
          !-----------------------
          !D0 = (m1/m2)**(1./(params%p1-params%p2))
          do k=1, nz
            if (qfields(k, i1) .gt. 0.0) then
              D0=(1+dist_mu(k,ispec))/dist_lambda(k,ispec)
              if (D0 > params%Dmax) then
                mu_old=dist_mu(k,ispec)
                lam_old=dist_lambda(k,ispec)
                n0_old=dist_n0(k,ispec)

                alpha=D0/params%Dmax
                dist_lambda(k,ispec)=alpha*dist_lambda(k,ispec)
                dist_n0(k,ispec)=alpha**(params%p1)*dist_n0(k,ispec)

                qfields(k,i2)=dist_n0(k,ispec)
                if (params%l_3m) then
                  m3(k)=moment(dist_n0(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
                  qfields(k,i3)=m3(k)
                end if
#if VERBOSE==1
                print*, 'WARNING: adjusting number and m3', time, params%id, n0_old, dist_n0(k,ispec), m3_old(k), m3(k)
                print*, 'new m1, m2, m3 are: ', m1(k), m2(k), m3(k)
#endif
              end if
              if (D0 < params%Dmin) then
                mu_old=dist_mu(k,ispec)
                lam_old=dist_lambda(k,ispec)
                n0_old=dist_n0(k,ispec)

                alpha=D0/params%Dmin
                dist_lambda(k,ispec)=alpha*dist_lambda(k,ispec)
                dist_n0(k,ispec)=alpha**(params%p1)*dist_n0(k,ispec)

                qfields(k,i2)=dist_n0(k,ispec)
                if (params%l_3m) then
                  m3(k)=moment(dist_n0(k,ispec),dist_lambda(k,ispec),dist_mu(k,ispec),params%p3)
                  qfields(k,i3)=m3(k)
                end if
#if VERBOSE==1
                print*, 'WARNING: adjusting number and m3', time, params%id, n0_old, dist_n0(k,ispec), m3_old(k), m3(k)
                print*, 'new m1, m2, m3 are: ', m1(k), m2(k), m3(k)
#endif

#if DEF_MODEL==MODEL_KiD
                if (nx == 1) then
                  call save_dg(k, dist_n0(k,ispec)-n0_old, 'n0_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
                  call save_dg(k, dist_mu(k,ispec)-mu_old, 'mu_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
                  call save_dg(k, dist_lambda(k,ispec)-lam_old, 'lam_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), &
                       i_dgtime)
                  call save_dg(k, m3(k)-m3_old(k), 'm3_adjust_'//trim(hydro_names(params%i_1m)//trim(chcall)), i_dgtime)
                else
                  call save_dg(k, i_here, dist_n0(k,ispec)-n0_old, 'n0_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), &
                       i_dgtime)
                  call save_dg(k, i_here, dist_mu(k,ispec)-mu_old, 'mu_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), &
                       i_dgtime)
                  call save_dg(k, i_here, dist_lambda(k,ispec)-lam_old, 'lam_adjust_'//trim(hydro_names(params%i_1m))//&
                       trim(chcall), i_dgtime)
                  call save_dg(k, i_here, m3(k)-m3_old(k), 'm3_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
                end if
#endif
              end if
            end if
          end do
        end if

#if DEF_MODEL==MODEL_KiD
        do k=1, nz
          if (icall==1) then
            if (nx > 1) then
              call save_dg(k, i_here, dist_n0(k,ispec), 'diag_n0_'//trim(hydro_names(params%i_1m)), i_dgtime)
              call save_dg(k, i_here, dist_lambda(k,ispec), 'diag_lam_'//trim(hydro_names(params%i_1m)), i_dgtime)
              call save_dg(k, i_here, dist_mu(k,ispec), 'diag_mu_'//trim(hydro_names(params%i_1m)), i_dgtime)
            else
              call save_dg(k, dist_n0(k,ispec), 'diag_n0_'//trim(hydro_names(params%i_1m)), i_dgtime)
              call save_dg(k, dist_lambda(k,ispec), 'diag_lam_'//trim(hydro_names(params%i_1m)), i_dgtime)
              call save_dg(k, dist_mu(k,ispec), 'diag_mu_'//trim(hydro_names(params%i_1m)), i_dgtime)
            end if
          else
            if (nx > 1) then
              call save_dg(k, i_here, dist_n0(k,ispec), 'diag_n0_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, i_here, dist_lambda(k,ispec), 'diag_lam_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, i_here, dist_mu(k,ispec), 'diag_mu_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
            else
              call save_dg(k, dist_n0(k,ispec), 'diag_n0_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, dist_lambda(k,ispec), 'diag_lam_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, dist_mu(k,ispec), 'diag_mu_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
            end if
          end if
        end do        
#endif
  end subroutine query_distributions
end module distributions
