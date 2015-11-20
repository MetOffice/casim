MODULE distributions

USE variable_precision, ONLY: wp
USE mphys_parameters, ONLY: hydro_params, nz, rain_params
USE mphys_switches, ONLY: hydro_names, l_limit_psd, l_passive3m, max_mu_frac
USE lookup, ONLY: get_slope_generic, max_mu, get_n0, moment, get_mu, get_lam_n0
USE special, ONLY: GammaFunc
USE thresholds, ONLY: thresh_tidy, thresh_sig, thresh_large

USE m3_incs, ONLY: m3_inc_type3

#if DEF_MODEL==MODEL_KiD
! Kid modules
USE diagnostics, ONLY: save_dg, i_dgtime, i_here, j_here, k_here
USE runtime, ONLY: time
USE parameters, ONLY: nx
#elif DEF_MODEL==MODEL_LEM
USE diaghelp_lem, ONLY: i_here, j_here, k_here
USE com_params, ONLY: time
#elif DEF_MODEL==MODEL_UM
USE timestep_mod, ONLY: time => timestep_number
#elif  DEF_MODEL==MODEL_MONC
USE diaghelp_monc, ONLY: i_here, j_here, time
#endif
IMPLICIT NONE

REAL(wp), ALLOCATABLE :: dist_lambda(:,:), dist_mu(:,:), dist_n0(:,:)

LOGICAL :: l_verbose=.FALSE.
CONTAINS

SUBROUTINE query_distributions(params, qfields, dist_lambda, dist_mu, dist_n0 &
     , icall)

  ! Any changes in number should be applied to the prognostic variable
  ! rather than just these parameters.  Currently this is not done.

TYPE(hydro_params), INTENT(IN) :: params !< species parameters
REAL(wp), INTENT(INOUT) :: qfields(:,:)
REAL(wp), INTENT(OUT), TARGET :: dist_lambda(:,:)
REAL(wp), INTENT(OUT), TARGET :: dist_mu(:,:)
REAL(wp), INTENT(OUT), TARGET :: dist_n0(:,:)

INTEGER, INTENT(IN), OPTIONAL :: icall

INTEGER :: k
REAL(wp) :: mass, m1,m2, m3
INTEGER(wp) :: i1,i2,i3,ispec
REAL(wp) :: n0_old, mu_old, lam_old
REAL(wp), POINTER :: lam, mu, n0
REAL(wp) :: m3_old, alpha, D0, mu_pass=1.0, k1,k2,k3

CHARACTER(2) :: chcall

! Some diagnostic strings
IF (PRESENT(icall)) THEN
  WRITE(chcall,'(A1,i1)') '_', icall
ELSE
  chcall=''
END IF


ispec = params%id
i1=params%i_1m
i2=params%i_2m
i3=params%i_3m

DO k=1,nz
  mass = qfields(k, i1)
  lam => dist_lambda(k,ispec)
  mu => dist_mu(k,ispec)
  n0 => dist_n0(k,ispec)
  lam=0.0
  mu=0.0
  n0=0.0
  IF (mass <= 0.0) THEN
    lam=0.0
    mu=0.0
    n0=0.0
  ELSE
    m1=mass/params%c_x
    IF (params%l_2m)m2 = qfields(k, i2)
    IF (params%l_3m) THEN
      m3 = qfields(k, i3)
    ELSE
      m3 = 0.0
    END IF


    m3_old=m3

    ! If 3rd moment has become too small then recalculate using max_mu
    ! This shouldn't happen except in small number situations - be careful
    IF (params%l_3m) THEN
      IF (m3 < SPACING(m3)) THEN
        CALL m3_inc_type3(params%p1, params%p2, params%p3,     &
             m1, m2, m3, max_mu)
        IF (l_verbose) THEN
          PRINT*, 'WARNING: resetting negative third moment',  &
             time, params%id, m3, m3_old
          PRINT*, 'm1 and m2 are: ', mass/params%c_x, m2
        END IF
        qfields(k,i3) = m3
      ELSE IF (m3 > thresh_large(params%i_3m)) THEN
        CALL m3_inc_type3(params%p1, params%p2, params%p3,     &
             m1, m2, m3, 0.0_wp)
        IF (l_verbose) THEN
          PRINT*, 'WARNING: resetting large third moment',     &
             time, params%id, m3, m3_old
          PRINT*, 'm1 and m2 are: ', mass/params%c_x, m2
        END IF
        qfields(k,i3) = m3
      END IF
    END IF

#if DEF_MODEL==MODEL_KiD
    ! as a test calculate what mu should be
    IF (l_passive3m) THEN
      CALL get_mu(m1, m2, m3, params%p1, params%p2, params%p3, mu_pass)
      IF (nx>1) THEN
        CALL save_dg(k, i_here, mu_pass,                       &
             'save_mu_'//TRIM(hydro_names(params%i_1m)), i_dgtime)
      ELSE
        CALL save_dg(k, mu_pass,                               &
             'save_mu_'//TRIM(hydro_names(params%i_1m)), i_dgtime)
      END IF
    END IF
#endif
    CALL get_slope_generic(k, params, n0, lam, mu, mass, m2, m3)


    ! If we diagnose a mu out of bounds, then reset m3
    IF (params%l_3m) THEN
      IF (mu<0.0 .OR. mu_pass < 0.0) THEN
        mu_old=mu
        mu=0.0
        CALL get_lam_n0(m1, m2, params%p1, params%p2, mu, lam, n0)
        m3=moment(n0,lam,mu,params%p3)
        qfields(k,i3) = m3
        IF (l_verbose)PRINT*, 'WARNING: resetting negative mu', &
           time, mu_old, m1,m2,m3,m3_old
      ELSE IF ((mu + EPSILON(1.0) > max_mu .OR. mu_pass +EPSILON(1.0) > max_mu)&
              .AND. .NOT. l_limit_psd) THEN
        mu_old=mu
        mu=max_mu
        CALL get_lam_n0(m1, m2, params%p1, params%p2, mu, lam, n0)
        m3=moment(n0,lam,mu,params%p3)
        qfields(k,i3) = m3
        IF (l_verbose)PRINT*, 'WARNING: resetting large mu',    &
           time, mu_old, m1,m2,m3,m3_old
      END IF
    END IF

    IF (l_limit_psd .AND. params%l_2m) THEN
      IF (params%l_3m .AND. mu > max_mu_frac*max_mu) THEN
        !-----------------------
        ! Adjust mu/m3 necessary
        !-----------------------
        mu_old=mu
        mu=(mu + max_mu_frac*max_mu)*0.5
        CALL get_lam_n0(m1, m2, params%p1, params%p2, mu, lam, n0)
        m3_old=m3
        m3=moment(n0,lam,mu,params%p3)
        qfields(k,i3) = m3
        IF (l_verbose)PRINT*, 'WARNING: adjusting m3 with large mu', &
           time, params%id, mu, mu_old, m1,m2,m3,m3_old
      END IF

      !-----------------------
      ! Adjust D0 if necessary
      !-----------------------
      !D0 = (m1/m2)**(1./(params%p1-params%p2))
      D0 = (1+mu)/lam
      IF (D0 > params%Dmax) THEN
        mu_old = mu
        lam_old=lam
        n0_old=n0

        alpha = D0/params%Dmax
        lam = alpha*lam
        n0 = alpha**(params%p1)*n0

        qfields(k,i2) = n0
        IF (params%l_3m) THEN
          m3=moment(n0,lam,mu,params%p3)
          qfields(k,i3) = m3
        END IF

        IF (l_verbose) THEN
          PRINT*, 'WARNING: adjusting number and m3', &
             time, params%id, n0_old, n0, m3_old, m3
          PRINT*, 'new m1, m2, m3 are: ', m1, m2, m3
        END IF
      END IF
      IF (D0 < params%Dmin) THEN
        mu_old = mu
        lam_old=lam
        n0_old=n0

        alpha = D0/params%Dmin
        lam = alpha*lam
        n0 = alpha**(params%p1)*n0

        qfields(k,i2) = n0
        IF (params%l_3m) THEN
          m3=moment(n0,lam,mu,params%p3)
          qfields(k,i3) = m3
        END IF

        IF (l_verbose) THEN
          PRINT*, 'WARNING: adjusting number and m3', &
             time, params%id, n0_old, n0, m3_old, m3
          PRINT*, 'new m1, m2, m3 are: ', m1, m2, m3
        END IF

#if DEF_MODEL==MODEL_KiD
        IF (nx==1) THEN
          CALL save_dg(k, n0-n0_old, 'n0_adjust_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
          CALL save_dg(k, mu-mu_old, 'mu_adjust_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
          CALL save_dg(k, lam-lam_old, 'lam_adjust_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
          CALL save_dg(k, m3-m3_old, 'm3_adjust_'//TRIM(hydro_names(params%i_1m)//TRIM(chcall)), i_dgtime)
        ELSE
          CALL save_dg(k, i_here, n0-n0_old, 'n0_adjust_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
          CALL save_dg(k, i_here, mu-mu_old, 'mu_adjust_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
          CALL save_dg(k, i_here, lam-lam_old, 'lam_adjust_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
          CALL save_dg(k, i_here, m3-m3_old, 'm3_adjust_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
        END IF
#endif
      END IF
    END IF

#if DEF_MODEL==MODEL_KiD
    IF (icall==1) THEN
      IF (nx>1) THEN
        CALL save_dg(k, i_here, n0, 'diag_n0_'//TRIM(hydro_names(params%i_1m)), i_dgtime)
        CALL save_dg(k, i_here, lam, 'diag_lam_'//TRIM(hydro_names(params%i_1m)), i_dgtime)
        CALL save_dg(k, i_here, mu, 'diag_mu_'//TRIM(hydro_names(params%i_1m)), i_dgtime)
      ELSE
        CALL save_dg(k, n0, 'diag_n0_'//TRIM(hydro_names(params%i_1m)), i_dgtime)
        CALL save_dg(k, lam, 'diag_lam_'//TRIM(hydro_names(params%i_1m)), i_dgtime)
        CALL save_dg(k, mu, 'diag_mu_'//TRIM(hydro_names(params%i_1m)), i_dgtime)
      END IF
    ELSE
      IF (nx>1) THEN
        CALL save_dg(k, i_here, n0, 'diag_n0_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
        CALL save_dg(k, i_here, lam, 'diag_lam_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
        CALL save_dg(k, i_here, mu, 'diag_mu_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
      ELSE
        CALL save_dg(k, n0, 'diag_n0_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
        CALL save_dg(k, lam, 'diag_lam_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
        CALL save_dg(k, mu, 'diag_mu_'//TRIM(hydro_names(params%i_1m))//TRIM(chcall), i_dgtime)
      END IF
    END IF
#endif

    NULLIFY(n0)
    NULLIFY(mu)
    NULLIFY(lam)

  END IF

END DO

END SUBROUTINE query_distributions

END MODULE distributions
