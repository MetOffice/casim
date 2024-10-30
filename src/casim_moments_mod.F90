MODULE casim_moments_mod
!#if DEF_MODEL==UM
! Equivalent of lsp_moments on the UM side
use variable_precision, only: wp

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CASIM_MOMENTS_MOD'

CONTAINS

SUBROUTINE casim_moments(params, points, rho, t, qcf, cficei, n_out, moment_out )

  ! microphysics modules
!USE mphys_inputs_mod, ONLY: l_psd_global
USE mphys_parameters, ONLY: hydro_params

! JW: note that most of these modules below are not available to MONC
! and some thought has to be given as to how to include the code on the
! CASIM trunk. Especially note that vectlib code will have to be
! written back into something which works appropriately where those
! modules aren't available.

! General atmosphere modules
!USE conversions_mod,  ONLY: zerodegc
USE, INTRINSIC :: ISO_C_BINDING, ONLY:                   &
     C_DOUBLE
! This could be a parameter but should think about having a standardised
! value across CASIM.

! Dr Hook modules
USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim

USE vectlib_mod,      ONLY: powr_v, oneover_v, exp_v

IMPLICIT NONE

#if DEF_MODEL==MODEL_MONC
INTEGER, PARAMETER :: real_64 = C_DOUBLE
! General atmosphere modules
REAL(KIND=real_64), PARAMETER :: zerodegc            = 273.15_real_64
#else
USE conversions_mod,  ONLY: zerodegc
#endif

! Purpose:
!   Obtains a specified moment of the in-cloud  ice particle size
!   distribution given grid box mean input ice water content.

! Method:
!   Follows one of two generic ice particle size distribution
!   calculation methods dependent on the type selected.

!   a) if l_psd_global is .false., we follow Field et al, 2005.
!   Parametrization of ice-particle size distributions for mid-latitude
!   stratiform cloud. Quart. J. Royal Meterol. Soc., 131, pp 1997-2017.

!   b) if l_psd_global is .true., we follow the later method of
!   Field et al 2007. Snow Size Distribution Parameterization for
!   Midlatitude and Tropical Ice Clouds. J. Atmos. Sci., 64,
!   pp 4346-4365

!   In the midlatitudes, method b) should relax to method a), however
!   for comparisons and completeness, both methods have been included here.
!   The extra flexibility will allow selection of the appopriate routine if
!   this method is used operationally.

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   Fortran 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

!   This subroutine calculates the quantity p_moment_out = M(n_out)
!   given the quantity (rho qcf / cfice)=ai*M(bi) where M(bi) is the
!   bi'th moment of the particle size distribution (PSD), corresponding
!   to the mass diameter relationship m(D) = ai D**bi.

!   It first calculates M(2) given the input quantities rho, qcf and ai
!   and then uses M(2) to calculate the output quantity p_out*M(out),
!   which can be used directly in the transfer calculations.
!   The conversion is from Field et al and gives
!   M(n)=a(n,Tc)*M(2)**b(n,Tc) where a and b are specified parameters.

!   This subroutine can be used to convert directly from M(x) to M(y)
!   by inputting ai=1, bi=x, n_out=y, rho=[array of 1's], qcf=M(x),
!   cficei=[array of 1's] and outputting moment_out=M(y).


! Subroutine Arguments

TYPE(hydro_params), INTENT(IN) :: params

INTEGER, INTENT(IN) ::                                                &
  points
                      ! Number of points to calculate

REAL(KIND=wp), INTENT(IN) ::                                          &
  n_out,                                                              &
                      ! Moment of the PSD to be output
  rho(points),                                                        &
                      ! Air density / kg m-3
  t(points),                                                          &
                      ! Temperature / K
  qcf(points),                                                        &
                      ! Ice water content / kg kg-1
  cficei(points)
                      ! 1 / ice cloud fraction

REAL(KIND=wp), INTENT(OUT) ::                                         &
  moment_out(points)
                      ! n_out moment of the in-cloud
                      ! particle size distn

! Local Variables

INTEGER ::                                                            &
  i, j
                      ! Loop counter for points

REAL(KIND=wp) ::                                                      &
  one_over_ai,                                                        &
                      ! 1 / ai
  bi,                                                                 &
                      ! bi (m=ai*D^bi)
  tc(points),                                                         &
                      ! Temperature in degrees Celsius
  log10_abi(points),                                                  &
                      ! log10 of conversion factor a(bi,Tc)
  abi(points),                                                        &
                      ! Conversion factor a(bi,Tc)
  bbi(points),                                                        &
                      ! Conversion factor b(bi,Tc)
  m_2(points),                                                        &
                      ! Second moment of the PSD
  m_bi(points),                                                       &
                      ! bi'th moment of the PSD
  log10_an_out(points),                                               &
                      ! log10 of conversion factor a(n_out,Tc)
  an_out(points),                                                     &
                      ! Conversion factor a(n_out,Tc)
  bn_out(points),                                                     &
                      ! Conversion factor b(n_out,Tc)
  a_tot_1, b_tot_1, c_tot_1, a_tot_2, b_tot_2, c_tot_2,               &
                      ! Totals of conversion factors for global
                      ! version
  temp_in(points), tc_1_in(points) , tc_2_in(points),                 &
  temp_out(points), tc_1_out(points) , tc_2_out(points)
                      ! temporaries for vector computations



    ! The following values are from Table 1 of Field et al (2005)
    ! and represent the conversion parameters between moments.
REAL, PARAMETER:: a(10) = (/5.065339,-0.062659,-3.032362,             &
  0.029469,-0.000285,0.312550,0.000204,0.003199,0.000000,             &
  -0.015952 /)
REAL, PARAMETER:: b(10) = (/0.476221,-0.015896,0.165977,              &
  0.007468,-0.000141,0.060366,0.000079,0.000594,0.000000,             &
  -0.003577 /)


    ! The following values are from Field et al (2007)
    ! and represent the conversion parameters between moments.

REAL, PARAMETER:: gl_a(3) = (/13.6078, -7.76062, 0.478694 /)
REAL, PARAMETER:: gl_b(3) = (/-0.0360722, 0.0150830, 0.00149453 /)
REAL, PARAMETER:: gl_c(3) = (/0.806856, 0.00581022, 0.0456723 /)

REAL, PARAMETER :: smallnum = 2.2e-14
                        ! Small value used in if tests
LOGICAL :: l_psd_global

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CASIM_MOMENTS'


! Start the subroutine

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
one_over_ai = 1.0 / params%c_x
bi = params%d_x

    !--------------------------------------------------------
    ! Decide whether to use the later (global) version of the
    ! distribution or the earlier (midlat) version of the
    ! distribution as controlled by the logical l_psd_global
    !--------------------------------------------------------
l_psd_global=.TRUE.

IF (l_psd_global) THEN

  ! calculate loop invariant quantities

  a_tot_1=gl_a(1)+gl_a(2)*bi+gl_a(3)*bi**2
  b_tot_1=gl_b(1)+gl_b(2)*bi+gl_b(3)*bi**2
  c_tot_1=gl_c(1)+gl_c(2)*bi+gl_c(3)*bi**2

  ! Now do the exponential of a_tot (save variables)
  a_tot_1=EXP(a_tot_1)
  c_tot_1 = 1.0/c_tot_1

  a_tot_2=gl_a(1)+gl_a(2)*n_out+gl_a(3)*n_out**2
  b_tot_2=gl_b(1)+gl_b(2)*n_out+gl_b(3)*n_out**2
  c_tot_2=gl_c(1)+gl_c(2)*n_out+gl_c(3)*n_out**2

  a_tot_2=EXP(a_tot_2)

  !--------------------------------------------------
  ! Use the global version of the particle size distn
  ! Field et al (2007) J. Atmos. Sci.
  !--------------------------------------------------

  !-----------------------------------------------
  ! Form the bi'th moment of the ice particle size
  ! distribution from the ice water content.
  !-----------------------------------------------

  j = 0
  DO i = 1, points

    ! pack into vector
    IF (qcf(i) > smallnum) THEN
      j = j + 1
      tc_1_in(j) = (t(i) - zerodegc)*b_tot_1
      tc_2_in(j) = (t(i) - zerodegc)*b_tot_2
      temp_in(j) = rho(i) * qcf(i) * cficei(i) * one_over_ai
      temp_in(j) = MAX(temp_in(j),0.0)

    END IF

  END DO

  CALL exp_v(j, tc_1_in, tc_1_out)
  CALL exp_v(j, tc_2_in, tc_2_out)

  tc_1_in(1:j) = a_tot_1*tc_1_out(1:j)
  CALL oneover_v(j, tc_1_in, tc_1_out)
  temp_in(1:j) = temp_in(1:j)*tc_1_out(1:j)

  !-----------------------------------------------
  ! Calculate the second moment of the PSD
  !-----------------------------------------------
  CALL powr_v(j,temp_in, c_tot_1, temp_out)
  temp_in(1:j) = temp_out(1:j)
  !-----------------------------------------------
  ! Calculate the n_out moment of the PSD
  !-----------------------------------------------
  CALL powr_v(j,temp_in, c_tot_2, temp_out)

  temp_out(1:j) = a_tot_2*tc_2_out(1:j)*temp_out(1:j)
  j = 0


  DO i = 1, points
    ! unpack
    IF (qcf(i) > smallnum) THEN

      j = j + 1
      moment_out(i)=temp_out(j)

    ELSE  ! qcf > 0

          !-----------------------------------------------
          ! Set the output moment to zero
          !-----------------------------------------------
      moment_out(i) = 0.0

    END IF  ! qcf > 0

  END DO  ! i

ELSE

        !--------------------------------------------------
        ! Use the midlat version of the particle size distn
        ! Field et al (2005) Q. J. R Met Soc
        !--------------------------------------------------

  DO i = 1, points

    IF (qcf(i) > smallnum) THEN
         !-----------------------------------------------
         ! Form the bi'th moment of the ice particle size
         ! distribution from the ice water content.
         !-----------------------------------------------
      m_bi(i) = rho(i) * qcf(i) * cficei(i) * one_over_ai
      m_bi(i) = MAX(m_bi(i) , 0.0)

         !-----------------------------------------------
         ! Calculate the second moment of the PSD
         !-----------------------------------------------
      tc(i) = t(i) - zerodegc

      log10_abi(i) = a(1) + a(2)*tc(i) + a(3)*bi                      &
                + a(4)*tc(i)*bi + a(5)*tc(i)*tc(i)                    &
                + a(6)*bi*bi + a(7)*tc(i)*tc(i)*bi                    &
                + a(8)*tc(i)*bi*bi + a(9)*tc(i)*tc(i)*tc(i)           &
                + a(10)*bi*bi*bi
      abi(i) = 10.0**(log10_abi(i))

      bbi(i) =    b(1) + b(2)*tc(i) + b(3)*bi                         &
                + b(4)*tc(i)*bi + b(5)*tc(i)*tc(i)                    &
                + b(6)*bi*bi + b(7)*tc(i)*tc(i)*bi                    &
                + b(8)*tc(i)*bi*bi + b(9)*tc(i)*tc(i)*tc(i)           &
                + b(10)*bi*bi*bi

      m_2(i) = (m_bi(i) / abi(i))**(1.0/bbi(i))

         !-----------------------------------------------
         ! Calculate the n_out moment of the PSD
         !-----------------------------------------------
      log10_an_out(i) = a(1) + a(2)*tc(i) + a(3)*n_out                &
                + a(4)*tc(i)*n_out + a(5)*tc(i)*tc(i)                 &
                + a(6)*n_out*n_out + a(7)*tc(i)*tc(i)*n_out           &
                + a(8)*tc(i)*n_out*n_out + a(9)*tc(i)*tc(i)*tc(i)     &
                + a(10)*n_out*n_out*n_out
      an_out(i) = 10.0**(log10_an_out(i))

      bn_out(i) = b(1) + b(2)*tc(i) + b(3)*n_out                      &
                + b(4)*tc(i)*n_out + b(5)*tc(i)*tc(i)                 &
                + b(6)*n_out*n_out + b(7)*tc(i)*tc(i)*n_out           &
                + b(8)*tc(i)*n_out*n_out + b(9)*tc(i)*tc(i)*tc(i)     &
                + b(10)*n_out*n_out*n_out

      moment_out(i) = an_out(i) * m_2(i) ** bn_out(i)

    ELSE  ! qcf > 0
         !-----------------------------------------------
         ! Set the output moment to zero
         !-----------------------------------------------
      moment_out(i) = 0.0

    END IF  ! qcf > 0

  END DO  ! i

END IF  ! l_psd_global


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE casim_moments
!#endif
END MODULE casim_moments_mod
