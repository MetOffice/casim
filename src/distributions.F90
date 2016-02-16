module distributions
  use variable_precision, only: wp
  use mphys_parameters, only: hydro_params, nz, rain_params
  use mphys_switches, only: hydro_names, l_limit_psd, l_passive3m, max_mu_frac
  use lookup, only: get_slope_generic, max_mu, get_n0, moment, get_mu, get_lam_n0
  use special, only: GammaFunc
  use thresholds, only: thresh_tidy, thresh_sig, thresh_large

  use m3_incs, only: m3_inc_type3

#if DEF_MODEL==MODEL_KiD
  ! Kid modules
  use diagnostics, only: save_dg, i_dgtime, i_here, j_here, k_here
  use runtime, only: time
  use parameters, only: nx
#elif DEF_MODEL==MODEL_LEM
  use diaghelp_lem, only: i_here, j_here, k_here
  use com_params, only: time
#elif DEF_MODEL==MODEL_UM
  use timestep_mod, only: time => timestep_number
#elif  DEF_MODEL==MODEL_MONC
  use diaghelp_monc, only: i_here, j_here, time
#endif
  implicit none
  private

  real(wp), allocatable :: dist_lambda(:,:), dist_mu(:,:), dist_n0(:,:)

  logical :: l_verbose=.false.

  public query_distributions, dist_lambda, dist_mu, dist_n0
contains

  ! Any changes in number should be applied to the prognostic variable
  ! rather than just these parameters.  Currently this is not done.
  subroutine query_distributions(params, qfields, dist_lambda, dist_mu, dist_n0, icall)
    type(hydro_params), intent(in) :: params !< species parameters
    real(wp), intent(inout) :: qfields(:,:)
    real(wp), intent(out), target :: dist_lambda(:,:)
    real(wp), intent(out), target :: dist_mu(:,:)
    real(wp), intent(out), target :: dist_n0(:,:)
    integer, intent(in), optional :: icall

    integer :: k
    real(wp) :: mass, m1,m2, m3
    integer(wp) :: i1,i2,i3,ispec
    real(wp) :: n0_old, mu_old, lam_old
    real(wp), pointer :: lam, mu, n0
    real(wp) :: m3_old, alpha, D0, mu_pass=1.0, k1,k2,k3

    character(2) :: chcall

    ! Some diagnostic strings
    if (present(icall)) then
      write(chcall,'(A1,i1)') '_', icall
    else
      chcall=''
    end if

    ispec=params%id
    i1=params%i_1m
    i2=params%i_2m
    i3=params%i_3m

    do k=1, nz
      mass=qfields(k, i1)
      lam=>dist_lambda(k,ispec)
      mu=>dist_mu(k,ispec)
      n0=>dist_n0(k,ispec)
      lam=0.0
      mu=0.0
      n0=0.0
      if (mass <= 0.0) then
        lam=0.0
        mu=0.0
        n0=0.0
      else
        m1=mass/params%c_x
        if (params%l_2m) m2=qfields(k, i2)
        if (params%l_3m) then
          m3=qfields(k, i3)
        else
          m3=0.0
        end if
        m3_old=m3

        ! If 3rd moment has become too small then recalculate using max_mu
        ! This shouldn't happen except in small number situations - be careful
        if (params%l_3m) then
          if (m3 < spacing(m3)) then
            call m3_inc_type3(params%p1, params%p2, params%p3, m1, m2, m3, max_mu)
            if (l_verbose) then
              print*, 'WARNING: resetting negative third moment', time, params%id, m3, m3_old
              print*, 'm1 and m2 are: ', mass/params%c_x, m2
            end if
            qfields(k,i3)=m3
          else if (m3 > thresh_large(params%i_3m)) then
            call m3_inc_type3(params%p1, params%p2, params%p3, m1, m2, m3, 0.0_wp)
            if (l_verbose) then
              print*, 'WARNING: resetting large third moment', time, params%id, m3, m3_old
              print*, 'm1 and m2 are: ', mass/params%c_x, m2
            end if
            qfields(k,i3)=m3
          end if
        end if

#if DEF_MODEL==MODEL_KiD
        ! as a test calculate what mu should be
        if (l_passive3m) then
          call get_mu(m1, m2, m3, params%p1, params%p2, params%p3, mu_pass)
          if (nx>1) then
            call save_dg(k, i_here, mu_pass, 'save_mu_'//trim(hydro_names(params%i_1m)), i_dgtime)
          else
            call save_dg(k, mu_pass, 'save_mu_'//trim(hydro_names(params%i_1m)), i_dgtime)
          end if
        end if
#endif
        call get_slope_generic(k, params, n0, lam, mu, mass, m2, m3)


        ! If we diagnose a mu out of bounds, then reset m3
        if (params%l_3m) then
          if (mu<0.0 .or. mu_pass < 0.0) then
            mu_old=mu
            mu=0.0
            call get_lam_n0(m1, m2, params%p1, params%p2, mu, lam, n0)
            m3=moment(n0,lam,mu,params%p3)
            qfields(k,i3)=m3
            if (l_verbose)print*, 'WARNING: resetting negative mu', time, mu_old, m1,m2,m3,m3_old
          else if ((mu + epsilon(1.0) > max_mu .or. mu_pass +epsilon(1.0) > max_mu) .and. .not. l_limit_psd) then
            mu_old=mu
            mu=max_mu
            call get_lam_n0(m1, m2, params%p1, params%p2, mu, lam, n0)
            m3=moment(n0,lam,mu,params%p3)
            qfields(k,i3)=m3
            if (l_verbose)print*, 'WARNING: resetting large mu', time, mu_old, m1,m2,m3,m3_old
          end if
        end if

        if (l_limit_psd .and. params%l_2m) then
          if (params%l_3m .and. mu > max_mu_frac*max_mu) then
            !-----------------------
            ! Adjust mu/m3 necessary
            !-----------------------
            mu_old=mu
            mu=(mu + max_mu_frac*max_mu)*0.5
            call get_lam_n0(m1, m2, params%p1, params%p2, mu, lam, n0)
            m3_old=m3
            m3=moment(n0,lam,mu,params%p3)
            qfields(k,i3) = m3
            if (l_verbose)print*, 'WARNING: adjusting m3 with large mu', &
                 time, params%id, mu, mu_old, m1,m2,m3,m3_old
          end if

          !-----------------------
          ! Adjust D0 if necessary
          !-----------------------
          !D0 = (m1/m2)**(1./(params%p1-params%p2))
          D0=(1+mu)/lam
          if (D0 > params%Dmax) then
            mu_old=mu
            lam_old=lam
            n0_old=n0

            alpha=D0/params%Dmax
            lam=alpha*lam
            n0=alpha**(params%p1)*n0

            qfields(k,i2)=n0
            if (params%l_3m) then
              m3=moment(n0,lam,mu,params%p3)
              qfields(k,i3)=m3
            end if

            if (l_verbose) then
              print*, 'WARNING: adjusting number and m3', time, params%id, n0_old, n0, m3_old, m3
              print*, 'new m1, m2, m3 are: ', m1, m2, m3
            end if
          end if
          if (D0 < params%Dmin) then
            mu_old=mu
            lam_old=lam
            n0_old=n0

            alpha=D0/params%Dmin
            lam=alpha*lam
            n0=alpha**(params%p1)*n0

            qfields(k,i2)=n0
            if (params%l_3m) then
              m3=moment(n0,lam,mu,params%p3)
              qfields(k,i3)=m3
            end if

            if (l_verbose) then
              print*, 'WARNING: adjusting number and m3', time, params%id, n0_old, n0, m3_old, m3
              print*, 'new m1, m2, m3 are: ', m1, m2, m3
            end if

#if DEF_MODEL==MODEL_KiD
            if (nx == 1) then
              call save_dg(k, n0-n0_old, 'n0_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, mu-mu_old, 'mu_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, lam-lam_old, 'lam_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, m3-m3_old, 'm3_adjust_'//trim(hydro_names(params%i_1m)//trim(chcall)), i_dgtime)
            else
              call save_dg(k, i_here, n0-n0_old, 'n0_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, i_here, mu-mu_old, 'mu_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, i_here, lam-lam_old, 'lam_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
              call save_dg(k, i_here, m3-m3_old, 'm3_adjust_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
            end if
#endif
          end if
        end if

#if DEF_MODEL==MODEL_KiD
        if (icall==1) then
          if (nx > 1) then
            call save_dg(k, i_here, n0, 'diag_n0_'//trim(hydro_names(params%i_1m)), i_dgtime)
            call save_dg(k, i_here, lam, 'diag_lam_'//trim(hydro_names(params%i_1m)), i_dgtime)
            call save_dg(k, i_here, mu, 'diag_mu_'//trim(hydro_names(params%i_1m)), i_dgtime)
          else
            call save_dg(k, n0, 'diag_n0_'//trim(hydro_names(params%i_1m)), i_dgtime)
            call save_dg(k, lam, 'diag_lam_'//trim(hydro_names(params%i_1m)), i_dgtime)
            call save_dg(k, mu, 'diag_mu_'//trim(hydro_names(params%i_1m)), i_dgtime)
          end if
        else
          if (nx > 1) then
            call save_dg(k, i_here, n0, 'diag_n0_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
            call save_dg(k, i_here, lam, 'diag_lam_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
            call save_dg(k, i_here, mu, 'diag_mu_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
          else
            call save_dg(k, n0, 'diag_n0_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
            call save_dg(k, lam, 'diag_lam_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
            call save_dg(k, mu, 'diag_mu_'//trim(hydro_names(params%i_1m))//trim(chcall), i_dgtime)
          end if
        end if
#endif
        nullify(n0)
        nullify(mu)
        nullify(lam)
      end if
    end do
  end subroutine query_distributions
end module distributions
