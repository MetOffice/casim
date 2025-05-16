module graupel_embryo
  use variable_precision, only: wp
  use process_routines, only: process_rate, i_sacw
  use passive_fields, only: rho
  use mphys_parameters, only: snow_params, graupel_params, cloud_params
  use mphys_constants, only: rho0
  use thresholds, only: thresh_sig, cfliq_small
  use special, only: pi, Gammafunc
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use mphys_switches, only: l_prf_cfrac, i_cfs, i_cfl, mpof
  use casim_stph, only: l_rp2_casim, mpof_casim_rp, snow_a_x_rp


  implicit none
  private

  character(len=*), parameter, private :: ModuleName='GRAUPEL_EMBRYO'

  public graupel_embryos
contains

  subroutine graupel_embryos(ixy_inner, dt, nz, l_Tcold, qfields, cffields, procs)

    !< Subroutine to convert some of the small rimed snow to graupel
    !< (Ikawa & Saito 1991)
    !<    !
    !< OPTIMISATION POSSIBILITIES: See gamma functions and distribution calculations
    !<
    !< AEROSOL: NOT DONE YET - internal category transfer

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: ixy_inner
    real(wp), intent(in) :: dt
    integer, intent(in) :: nz
    logical, intent(in) :: l_Tcold(:)
    real(wp), intent(in), target :: qfields(:,:)
    real(wp), intent(in) :: cffields(:,:)
    type(process_rate), intent(inout), target :: procs(:,:)

    ! Local variables
    real(wp) :: cloud_mass, snow_mass, snow_number
    real(wp) :: dmass, dnumber
    real(wp) :: snow_n0, snow_lam, snow_mu ! distribution parameters

    real(wp) :: dnembryo ! rate of embryo creation
    real(wp) :: embryo_mass = 1.6e-10 ! mass(kg) of a new graupel embryo (should be in parameters)

    real(wp) :: pgsacw  ! Rate of mass transfer to graupel
    real(wp) :: Eff

    real(wp) :: alpha = 4.0 ! Tunable parameter see Reisner 1998, Ikawa+Saito 1991
    real(wp) :: rhogms      ! Difference between graupel density and snow density

    real(wp) :: Garg ! Argument for Gamma function
    integer  :: pid  ! Process id

    real(wp) :: cf_snow, cf_liquid, overlap_cf

    integer :: k

    character(len=*), parameter :: RoutineName='GRAUPEL_EMBRYOS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Apply RP scheme
    if ( l_rp2_casim ) then
        mpof = mpof_casim_rp
        snow_params%a_x = snow_a_x_rp
    endif

    do k = 1, nz
       if (l_Tcold(k)) then 
          if (l_prf_cfrac) then
             if (cffields(k,i_cfs) .gt. cfliq_small) then
                cf_snow=cffields(k,i_cfs)
             else
                cf_snow=cfliq_small !nonzero value - maybe move cf test higher up
             endif
             if (cffields(k,i_cfl) .gt. cfliq_small) then
                cf_liquid=cffields(k,i_cfl)
             else
                cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
             endif
          else
             cf_snow=1.0
             cf_liquid=1.0
          endif


          snow_mass=qfields(k, snow_params%i_1m)   / cf_snow
          cloud_mass=qfields(k, cloud_params%i_1m)  / cf_liquid
          
          if (snow_mass * cf_snow > thresh_sig(snow_params%i_1m) .and. &
               cloud_mass * cf_liquid > thresh_sig(cloud_params%i_1m)) then
             
             if (snow_params%l_2m) snow_number=qfields(k, snow_params%i_2m) / cf_snow

             snow_n0=dist_n0(k,snow_params%id)
             snow_mu=dist_mu(k,snow_params%id)
             snow_lam=dist_lambda(k,snow_params%id)
             !convert to reisner N0, based on Reisner et al, QJMRS, 1998 
             !( https://doi.org/10.1002/qj.49712454804 )
             snow_n0=snow_n0*(snow_lam**(snow_mu+1.0))/(GammaFunc(1.0+snow_mu))    
             
             !< This efficiency is the same as used in sacw calculation, 
             !< so that we make the assumption that collecting area is approx half 
             !< of circle - similar to operational.
             !< This efficiency should be defined consistently in the parameters
             Eff=0.5_wp
             rhogms=graupel_params%density-snow_params%density  !this is inconsistent for mass~D**2
             
             Garg=2.0+2*snow_params%b_x + snow_mu
            pgsacw = (0.75*alpha*dt*pi/rhogms)*Eff*Eff*rho(k,ixy_inner)*rho(k,ixy_inner)*cloud_mass*cloud_mass   &
           *snow_params%a_x*snow_params%a_x*snow_n0*GammaFunc(Garg)*(2*snow_params%f_x + 2*snow_lam)**(-Garg) &
           !the 2* lambda etc doesnt look the same as in reisner)
           *(rho0/rho(k,ixy_inner))**(2*snow_params%g_x) !in-graupel rate
             
             dnembryo=max(snow_params%density*pgsacw/rhogms/embryo_mass/rho(k,ixy_inner), 0.0_wp)
             dnumber=min(dnembryo, 0.95*snow_number/dt)
             
             !use mixed-phase overlap function
             overlap_cf=min(1.0,max(0.0,mpof*min(cf_snow, cf_liquid)                      &
                  + max(0.0,(1.0-mpof)*(cf_snow+cf_liquid-1.0))))
     
             pid=i_sacw%id
             pgsacw=pgsacw*overlap_cf !convert back to grid box mean
             dmass=procs(snow_params%i_1m, pid)%column_data(k)- pgsacw !grid box mean
             dnumber=dnumber*overlap_cf !convert back to grid box mean
             
             procs(snow_params%i_1m,pid)%column_data(k)=dmass           
             procs(graupel_params%i_1m,pid)%column_data(k)=pgsacw   
             if (graupel_params%l_2m) then
                procs(graupel_params%i_2m,pid)%column_data(k)=dnumber
             end if
             if (snow_params%l_2m) then
                procs(snow_params%i_2m,pid)%column_data(k)=-dnumber
             end if

          end if
       end if
    enddo
      
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine graupel_embryos
end module graupel_embryo
