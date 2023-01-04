module shipway_activation_mod
  use variable_precision, only: wp

  ! A revised activation scheme
  ! See Shipway, B. J.: Revisiting Twomey's approximation for peak supersaturation,
  ! Atmos. Chem. Phys., 15, 3803-3814, doi:10.5194/acp-15-3803-2015, 2015.

  Use shipway_parameters
  Use shipway_constants
  Use shipway_lookup, only: lookup_I, xmax, ymax, xmin, ymin
  Use mphys_die, only: throw_mphys_error, incorrect_opt, mphys_message

  Implicit None
  private

  character(len=*), parameter, private :: ModuleName='SHIPWAY_ACTIVATION_MOD'

  real(wp) :: dv_flag=0 ! Determines method for calculating Dv
  real(wp) :: est_smax_time
  real(wp) :: psi1 ! placed here so it can be used for estimating smax time
  real(wp) :: C1=1.058, C2=1.904 ! See Korolev et al. paper

  real(wp) :: wdiag, Tdiag, Ndiag
  integer :: counter=0

  real(wp) :: cumulative_xmin=1.e10
  real(wp) :: cumulative_xmax=0.0
  integer :: cumulative_brent=0
  integer :: cumulative_quad=0
  integer :: cumulative_secant=0
  integer :: cumulative_mid=0

  real :: LHS_hold

  character(len=600) :: std_msg

  public solve_nccn_household,solve_nccn_brent

contains

  subroutine calc_nccn(smax, nccn, nccni_dg)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(in) :: smax
    real(wp), intent(out) :: nccn
    real(wp), intent(out), optional :: nccni_dg(max_nmodes)

    ! Local variables
    real(wp) :: nccni, erfarg
    integer :: i

    character(len=*), parameter :: RoutineName='CALC_NCCN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    nccn=0.0
    do i=1,nmodes
      if (use_mode(i)) then
        erfarg=log(smax/s0i(i))/(sqrt(2.)*log(sigmas(i)))
        nccni = 0.5*Ndi(i)*(1.0+erf(erfarg))

        nccn=nccn+nccni
        nccni_dg(i)=nccni
      end if
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine calc_nccn


  subroutine calc_LHS(w, T, pressure, alpha_c, entrain_fraction, LHSout)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(in) :: T ! Temperature (K)
    real(wp), intent(in) :: w ! vertical velocity (m/s)
    real(wp), intent(in) :: pressure ! Pressure (Pa)
    real(wp), intent(in) :: alpha_c ! Condensation coefficient
    real(wp), intent(in) :: entrain_fraction ! Entrainment_fraction
    real(wp), intent(out) :: LHSout

    ! Local variables

    real(wp) :: LHS ! intermediate value

    real(wp) :: Tc ! Temperature (C)
    real(wp) :: es ! Saturation vapour pressure

    real(wp) :: psi2, G, Dv_here

    real(wp) :: Lvt ! Temperature dependent Lv

    character(len=*), parameter :: RoutineName='CALC_LHS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    Tc = T-273.15

    es = (100.*6.1121)*exp((18.678-Tc/(234.5))*Tc/(257.14+Tc))

    LvT = Lv -(4217.4-1870.0)*Tc

    psi1 = 9.8*(LvT/(eps*cp*T)-1.0)/(T*Rd)*(1-entrain_fraction)
    psi2 = eps*pressure/es+LvT**2/(Rv*T**2*cp)

    if (int(dv_flag)==1)then
       Dv_here=Dv
    else if(dv_flag <= 0)then
       Dv_here=Dv_mean(T,alpha_c)
    else
       Dv_here=dv_flag
    end if

    G = 1/(rhow*(Rv*T/(es*Dv_here)+LvT*(LvT/(Rv*T)-1)/(ka*T)))

    LHS = sqrt(2.0*pi)*rho*(psi1*w)**1.5/(4*pi*rhow*psi2*G**1.5)

    LHSout=LHS

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine calc_LHS

  subroutine set_inputs(s,sp,sm,ds)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Routine to choose sm,s,sp and ensure they sit within (xminx, xmax)
    real(wp), intent(inout) :: s, sp, sm, ds

    ! Local variables
    logical :: adjust
    real(wp) :: s0ratioi_p, s0ratioi_m
    integer :: i

    character(len=*), parameter :: RoutineName='SET_INPUTS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ds=min(0.001, .5*abs(s))
    sp = s + ds
    sm = s - ds

    adjust=.true.
    do while (adjust)
      do i=1,nmodes
        if (use_mode(i)) then
          s0ratioi_p = sp/s0i(i)
          s0ratioi_m = sm/s0i(i)
          if (s0ratioi_p > xmax)then
            write (std_msg, *) 'Out of range X too big: ', i, xmin, s0i(i), &
                               s0ratioi_p, xmax
            call mphys_message(ModuleName//':'//RoutineName, std_msg)
            ! reset
            write (std_msg, *) 'old sm,s,sp', sm,s,sp
            call mphys_message(ModuleName//':'//RoutineName, std_msg)
            sp = xmax*s0i(i)*.99
            s  = sp - ds
            sm = s - ds
            write (std_msg, *) 'new sm,s,sp', sm,s,sp
            call mphys_message(ModuleName//':'//RoutineName, std_msg)
            exit
          end if
          if (s0ratioi_m < xmin)then
            write (std_msg, *)'Out of range X too small: ', i, xmin, s0i(i), &
                              s0ratioi_m, xmax
            call mphys_message(ModuleName//':'//RoutineName, std_msg)
            ! reset
            write (std_msg, *)'old sm,s,sp', sm,s,sp
            call mphys_message(ModuleName//':'//RoutineName, std_msg)
            sm = xmin*s0i(i)*1.01
            ds=min(0.001, sm)
            s  = sm + ds
            sp = s + ds
            write (std_msg, *)'new sm,s,sp', sm,s,sp
            call mphys_message(ModuleName//':'//RoutineName, std_msg)
            exit
          end if
        end if
        if (i==nmodes)adjust=.false.
      end do
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine set_inputs


  subroutine set_inputs_safe(sp,sm)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Routine to choose sm,sp and ensure they sit within (xminx, xmax)
    real(wp), intent(inout) :: sp, sm

    ! Local variables
    integer :: i

    character(len=*), parameter :: RoutineName='SET_INPUTS_SAFE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    sm=0.0
    sp=HUGE(sp)
    do i=1,nmodes
      if (use_mode(i)) then
        sm = max(xmin*s0i(i)*1.01, sm)
        sp = min(xmax*s0i(i)*0.99, sp)

      end if
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine set_inputs_safe

  subroutine get_extent(s, sp, sm)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Routine to choose sm,sp and ensure they sit within (xminx, xmax)
    real(wp), intent(in) :: s
    real(wp), intent(out) :: sp, sm

    ! Local variables

    integer :: i

    character(len=*), parameter :: RoutineName='GET_EXTENT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    sm=0.0
    sp=HUGE(sp)
    do i=1,nmodes
      if (use_mode(i)) then
        sm = max(s/s0i(i), sm)
        sp = min(s/s0i(i), sp)
      end if
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine get_extent

  subroutine calc_RHS(s,RHS)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(inout) :: s
    real(wp), intent(out) :: RHS

    ! Local variables

    real(wp) :: s0ratioi
    integer :: i

    real(wp) :: J1

    character(len=*), parameter :: RoutineName='CALC_RHS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    RHS=0.0
    do i=1,nmodes
      if (use_mode(i)) then
        s0ratioi = s/s0i(i)
        if (logsigmas(i) > ymax .or. logsigmas(i)<ymin)then
          write (std_msg,*) 'Out of range Y: ', ymin, logsigmas(i), ymax, &
                             s0ratioi
          call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                                 std_msg)
        end if
      end if
    end do
    do i=1,nmodes
      if (use_mode(i)) then
        call lookup_I(s0ratioi, logsigmas(i), J1)
        RHS = RHS + ai(i)*J1
      end if
    end do

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine calc_RHS

  subroutine calc_KC(T)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(in) :: T ! Temperature(K)
    real(wp) :: Ak
    integer :: i

    character(len=*), parameter :: RoutineName='CALC_KC'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    Ak = 2.*Mw*zetasa/(Ru*T*rhow)

    do i=1,nmodes
      if (use_mode(i)) then
       ! These are defined in Khvorostyanov and Curry (2006)
       s0i(i) = rdi(i)**(-(1.+betai(i)))*sqrt(4*Ak**3/(27*bi(i)))
       sigmas(i) = sigmad(i)**(1.+betai(i))

       ai(i) = Ndi(i)*s0i(i)*s0i(i)/log(sigmas(i))
       logsigmas(i) = log(sigmas(i))
     end if
    end do

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine calc_KC

  subroutine solve_nccn_household(order,niter,sa_in,  &
     w,T,pressure,alpha_c,ent_fraction,               &
     smax,nccn,nccni)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none
    ! Use householders methods (e.g. Newton, Halley) to find roots.

    integer, intent(in) :: order, niter
    real(wp), intent(in)    :: sa_in
    real(wp), intent(in)    :: w,T,pressure, alpha_c, ent_fraction
    real(wp), intent(out)   :: smax, nccn
    real(wp), intent(out)   :: nccni(3) ! diagnostic for each mode

    real(wp) :: sa
    integer :: it
    integer :: maxiter=10
    integer :: bigiter
    real(wp) :: RHSout,RHSout_p1,RHSout_m1, ds
    real(wp) :: F, dF, d2F, diff
    real(wp) :: LHS, sa_p1, sa_m1, sa_tm1

    real(wp) :: tolx=1e-6, toly ! This tolerence only necessary for very high numbers
    character(len=*), parameter :: RoutineName='SOLVE_NCCN_HOUSEHOLD'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (.not. any(use_mode)) then
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      return  ! If there's no aerosol present in any of the modes
    end if

    call calc_KC(T)
    call calc_LHS(w, T, pressure, alpha_c, ent_fraction, LHS)
    LHS_hold = LHS
    wdiag = w
    Tdiag = T
    Ndiag = sum(Ndi)

    !toly=LHS*.1 ! 10 percent error
    toly=LHS*1e5 ! Essentially not used

    sa = sa_in
    sa_tm1=sa

    diff=1000.
    RHSout=LHS+LHS
    counter=counter+1
    bigiter=0
    do while((abs(diff) > tolx .or. abs(RHSout-LHS)>toly) .and. bigiter<maxiter)
      do it=1,niter
        call set_inputs(sa, sa_p1, sa_m1, ds)
        sa_tm1=sa
        call calc_RHS(sa,RHSout)
        call calc_RHS(sa_p1,RHSout_p1)

        F=RHSout-LHS
        dF=(RHSout_p1-RHSout)/ds

        select case(order)
        case(1)
          diff=-F/dF
        case(2)
          call calc_RHS(sa_m1,RHSout_m1)
          d2F=(RHSout_p1-2*RHSout + RHSout_m1)/(ds*ds)
          diff=-2.*F*dF/(2.*dF*dF-F*d2F)
        end select
        sa = sa + diff

      end do

       bigiter=bigiter+1
    end do

    smax=sa
    est_smax_time=1./((psi1*w*C1/C2)/smax)

    call calc_nccn(smax, nccn, nccni)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

    end subroutine solve_nccn_household

  subroutine solve_nccn_brent( &
     w,T,pressure,alpha_c,ent_fraction,               &
     smax,nccn,nccni)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(in)    :: w,T,pressure, alpha_c, ent_fraction
    real(wp), intent(out)   :: smax, nccn
    real(wp), intent(out)   :: nccni(3) ! diagnostic for each mode

    real(wp) :: sa, sb
    real(wp) :: LHS, sa_p1, sa_m1
    integer :: iflag  ! flag to indicate if solution found
    integer :: ibrent ! indicates number of iteration in bre

    real(wp) :: tolx=1e-4           ! tolerence for brent
    real(wp) :: tolf=1e-1            ! tolerence for brent
    logical :: verbose=.false.  ! set brent to verbose output

    real(wp) :: ds, RHS
!    real(wp) :: J1
    integer :: nscan=101, i, iquad, isecant, imidpoint
    character(len=500) :: message

    character(len=*), parameter :: RoutineName='SOLVE_NCCN_BRENT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    if (.not. any(use_mode)) then
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      return  ! If there's no aerosol present in any of the modes
    end if

    call calc_KC(T)
    call calc_LHS(w, T, pressure, alpha_c, ent_fraction, LHS)

    LHS_hold = LHS
    wdiag = w
    Tdiag = T
    Ndiag = sum(Ndi)
    ibrent = 0
    RHS = 0.0

    call set_inputs_safe(sa_p1, sa_m1)
    sa = sa_m1
    sb = sa_p1

    call brent(smax_eq, sa, sb, smax, iflag, ibrent &
       ,tolx, tolf, verbose,iquad, isecant, imidpoint)

    if (iflag == -1) then
      if (smax_eq(sb) - LHS < 0)then
        ! In this case we've probably exceeded the bounds of the lookup table
        ! but we should be able to extrapolate using the asymptotic limit of the
        ! integral, i.e. lim x->inf I(x,ls) = x**2*I(xmax,ls)/xmax**2
        do i=1,nmodes
          if (use_mode(i)) then
!            call lookup_I(xmax, logsigmas(i), J1)
!            RHS = RHS + (ai(i)/s0i(i)/s0i(i))*(J1/xmax/xmax)
            ! NB I've left this in a long form for comparison with
            ! the integral form, but this simplifies considerable
            ! to give smax ~ w^3/4 N^-1/2
            RHS = RHS + (ai(i)/s0i(i)/s0i(i))*2*sqrt(pi)*logsigmas(i)
          end if
        end do
        smax=sqrt(LHS/RHS)
      else
        if (verbose)then
          do i=1,nmodes
            write (std_msg, *) i, smax/s0i(i)
            call mphys_message(ModuleName//':'//RoutineName, std_msg)
          end do
          write (std_msg, *) w,T, Ndi(:),logsigmas(:), LHS, smax
          call mphys_message(ModuleName//':'//RoutineName, std_msg)
          ds=(sb-sa)/nscan
          do i=1,nscan-1
            sa=sa+ds
            write(*,'(e18.4,A,e18.4,A)') sa, ',', smax_eq(sa),','
          end do
        end if
        write(message,'(A,e18.4,A,e18.4,A,e18.4,A,e18.4,A,e18.4,A,e18.4)') &
           'w=',w,'; p=',pressure,'; T=',T,'; N=',Ndiag,'; sm=',sa_m1,'; sp=',sa_p1
        call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
             'Beyond scope of numerics of the parametrization:'//trim(message))
      end if
    end if

    if (verbose)then
      write (std_msg, *) 'Number of brent calls:', ibrent
      call mphys_message(ModuleName//':'//RoutineName, std_msg)
      cumulative_brent=cumulative_brent + ibrent
      cumulative_quad=cumulative_quad + iquad
      cumulative_secant=cumulative_secant + isecant
      cumulative_mid=cumulative_mid + imidpoint
      write (std_msg, *) 'Cumulative brent calls:', cumulative_brent, &
         cumulative_quad, cumulative_secant, cumulative_mid
      call mphys_message(ModuleName//':'//RoutineName, std_msg)
      call get_extent(smax, sa, sb)
      write (std_msg, *) 'Extrema of ratios:', sa, sb
      call mphys_message(ModuleName//':'//RoutineName, std_msg)
      cumulative_xmin=min(sa, cumulative_xmin)
      cumulative_xmax=max(sb, cumulative_xmax)
      write (std_msg, *) 'Cumulative extrema', cumulative_xmin, cumulative_xmax
      call mphys_message(ModuleName//':'//RoutineName, std_msg)
    end if

    est_smax_time=1.0/((psi1*w*C1/C2)/smax)
    call calc_nccn(smax, nccn, nccni)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine solve_nccn_brent

  function smax_eq(sa)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp), intent(inout) :: sa
    real(wp) :: smax_eq

    real(wp) :: RHSout

    character(len=*), parameter :: RoutineName='SMAX_EQ'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    call calc_RHS(sa,RHSout)

    smax_eq = RHSout - LHS_hold

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end function smax_eq


  subroutine brent(f, a, b, x, iflag, ibrent, tolx, tolf, verbose, &
     iquad, isecant, imidpoint)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    real(wp) :: f
    real(wp), intent(inout) :: a, b
    real(wp), intent(out) :: x
    integer, intent(out) :: &
         iflag  & ! =-1 if f(a)f(b) > 0
                  ! =0 if successful
         , ibrent ! number of calls to function f (this does not
                  ! include calls made for verbose output.)

    real(wp), intent(inout) :: tolx,tolf
    logical, intent(inout) :: verbose
    integer, optional :: iquad, isecant, imidpoint

    real(wp) :: c, d, tmp, s
    logical :: mflag
    integer :: icount
    integer, parameter :: icountmax=50

    real(wp) :: pfs, fa, fs, fb, fc

    character(len=*), parameter :: RoutineName='BRENT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    iflag=-999
    ibrent=0
    iquad=0
    isecant=0
    imidpoint=0

    fa=f(a)
    ibrent=ibrent+1
    fb=f(b)
    ibrent=ibrent+1

    if (verbose)then
       write (std_msg, *) 'a = ', a, ', b = ', b
       call mphys_message(ModuleName//':'//RoutineName, std_msg)
       write (std_msg, *) 'F(a) = ', fa, ', F(b) = ', fb
       call mphys_message(ModuleName//':'//RoutineName, std_msg)
    end if

    if (fa == 0)then
       x=a
       iflag=0
       IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
       return
    end if

    if (fb == 0)then
       x=b
       iflag=0
       IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
       return
    end if

    if ( fa*fb > 0 )then
       x=-999
       iflag=-1
       IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
       return
    end if

    if (abs(fa) < abs(fb))then
       c=b
       fc=fb
       b=a
       fb=fa
       a=c
       fa=fc
    else
       c=a
       fc=fa
    end if

    mflag=.true.
    d=0
    fs=tolf+1.

    icount=0
    do while (fb /=0 .and. abs(b-a) > tolx .and. abs(fs) > tolf &
         .and. icount < icountmax )
       icount=icount+1
       if ((fa /= fc) .and. (fb /= fc)) then
          !inverse quadratic interpolation
          s=a*fb*fc/((fa-fb)*(fa-fc)) &
               + b*fa*fc/((fb-fa)*(fb-fc)) &
               + c*fa*fb/((fc-fa)*(fc-fb))
          if (verbose)then
             pfs=f(s)
             write (std_msg, *) 'inverse quadratic ','s=',s,'f(s)=', pfs
             call mphys_message(ModuleName//':'//RoutineName, std_msg)
          end if
          iquad=iquad+1
       else
          ! secant method
          s=b-fb*(b-a)/(fb-fa)
          if (verbose)then
             pfs=f(s)
             write (std_msg, *) 'secant ','s=',s,'f(s)=', pfs
             call mphys_message(ModuleName//':'//RoutineName, std_msg)
          end if
          isecant=isecant+1
       end if
       if ((s < b .and. s <(3*a+b)/4) &
            .or. (s > b .and. s >(3*a+b)/4) &
            .or. (mflag .and. abs(s-b)>=abs(b-c)/2.) &
            .or. ((.not.mflag) .and. (abs(s-b)>=abs(c-d)/2.)))then
          ! midpoint
          s=(a+b)/2
          if (verbose)then
             pfs=f(s)
             write (std_msg, *) 'midpoint ','s=',s,'f(s)=', pfs
             call mphys_message(ModuleName//':'//RoutineName, std_msg)
          end if
          imidpoint=imidpoint+1
          mflag=.true.
       else
          mflag=.false.
       end if
       fs=f(s)
       ibrent=ibrent+1
       d=c
       c=b
       if (fa*fs<0)then
          b=s
          fb=fs
       else
          a=s
          fa=fs
       end if


       if (abs(fa) < abs(fb))then
          tmp=b
          b=a
          a=tmp
          tmp=fb
          fb=fa
          fa=tmp
       end if

       if (verbose)then
         write (std_msg, *) '|a-b|=',abs(a-b)
         call mphys_message(ModuleName//':'//RoutineName, std_msg)
       end if
    end do


    if (icount==icountmax)then
       x=-999
       iflag=-2
    else
       x=b
       iflag=0
    end if


    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    return
  end subroutine brent

end module shipway_activation_mod
