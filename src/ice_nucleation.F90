module ice_nucleation
  use mphys_die, only: throw_mphys_error, bad_values, std_msg
  use variable_precision, only: wp
  use passive_fields, only: rho, pressure, w, exner
  use mphys_switches, only: i_qv, i_ql, i_qi, i_ni, i_th , hydro_complexity, i_am4, i_am6, i_an2, l_2mi, l_2ms, l_2mg, &
       i_am8, i_am9, aerosol_option, i_nl, i_ns, i_ng, iopt_inuc, i_am7, i_an6, i_an12, l_process, l_passivenumbers, &
       l_passivenumbers_ice, active_number, active_ice, isol, iinsol, l_itotsg, contact_efficiency, immersion_efficiency, &
       aero_index, l_prf_cfrac, iopt_act, i_cfl, i_cfi, l_nudge_to_cooper
  use process_routines, only: process_rate, i_inuc, i_dnuc
  use mphys_parameters, only: nucleated_ice_mass, cloud_params, ice_params
  use mphys_constants, only: Ls, cp, pi, m3_to_cm3
  use qsat_funs, only: qsaturation, qisaturation
  use thresholds, only: ql_small, w_small, ni_tidy, nl_tidy, cfliq_small
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active

  implicit none

  character(len=*), parameter, private :: ModuleName='ICE_NUCLEATION'

contains

  !> Currently this routine considers heterogeneous nucleation
  !> notionally as a combination of deposition and/or condensation freezing.
  !> Immersion (i.e. freezing through preeixisting resident IN within a cloud drop) and
  !> Contact freezing (i.e. collision between cloud drop and IN) are not
  !> yet properly concidered.  Such freezing mechanisms should consider the
  !> processing of the aerosol in different ways.
  subroutine inuc(ixy_inner, dt, nz, l_Tcold, qfields, cffields, procs, dustphys, aeroact, dustliq, &
       aerosol_procs)

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

    ! aerosol fields
    type(aerosol_phys), intent(in) :: dustphys(:)
    type(aerosol_active), intent(in) :: aeroact(:)
    type(aerosol_active), intent(in) :: dustliq(:)

    ! optional aerosol fields to be processed
    type(process_rate), intent(inout), optional, target :: aerosol_procs(:,:)


    ! Local variables
    real(wp) :: dmass, dnumber, dmad, dmac, dmadl

    ! Liquid water and ice saturation for Meyers equation
    real(wp) :: lws_meyers, is_meyers

    !coefficients for Demott parametrization
    real(wp) :: a_demott, b_demott, c_demott, d_demott, cf
    real(wp) :: Tp01 ! 273.16-Tk (or 0.01 - Tc)
    real(wp) :: Tc ! local temperature in C
    real(wp) :: Tk ! local temperature in K
    real(wp) :: th
    real(wp) :: qv
    real(wp) :: ice_number
    real(wp) :: cloud_number, cloud_mass
    real(wp) :: qs, qis, dN_imm, dN_contact, ql
    real(wp) :: Si(nz), Sw(nz)
    
    real(wp) :: cf_liquid, cf_ice

    ! parameters for Meyers et al (1992)
    ! Meyers MP, DeMott PJ, Cotton WR (1992) New primary ice-nucleation
    ! parameterizations in an explicit cloud model. J Appl Meteorol 31:708–721
    real(wp), parameter :: meyers_a = -0.639 ! Meyers eq 2.4 coeff a
    real(wp), parameter :: meyers_b = 0.1296 ! Meyers eq 2.4 coeff b

    ! parameters for Tobo et al. (2013)
    real(wp) :: a_tobo, b_tobo, c_tobo, d_tobo

    ! variables for surface site based parameterisations
    real(wp) :: n_sites, surf_area
    
    integer :: k

    logical :: l_condition

    character(len=*), parameter :: RoutineName='INUC'


    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    do k = 1, nz
      if (l_Tcold(k)) then
        th=qfields(k, i_th)
        Tk=th*exner(k,ixy_inner)
        qs=qsaturation(Tk, pressure(k,ixy_inner)/100.0)
        qis=qisaturation(Tk, pressure(k,ixy_inner)/100.0)
        qv=qfields(k, i_qv)
         
        !!    if (qs==0.0 .or. qis==0.0) then
        !!      write(std_msg, '(A)') 'Error in saturation calculation - qs or qis is zero'
        !!      call throw_mphys_error(bad_values,  ModuleName//':'//RoutineName, std_msg)
        !!    end if
        if (qs > 0.0) then
          Sw(k)=qv/qs - 1.0
        else
          Sw(k)=-999.0
        end if
        if (qis > 0.0) then
          Si(k)=qv/qis - 1.0
        else
          Si(k)=-999.0
        end if
      end if
    end do

    do k = 1, nz
       if (l_Tcold(k)) then
          if (l_prf_cfrac) then
             if (cffields(k,i_cfl) .gt. cfliq_small) then
                cf_liquid=cffields(k,i_cfl)
             else
                cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
             endif
             if (cffields(k,i_cfi) .gt. cfliq_small) then
                cf_ice=cffields(k,i_cfi)
             else
                cf_ice=cfliq_small !nonzero value - maybe move cf test higher up
             endif
          else
             cf_liquid=1.0
             cf_ice=1.0
          endif

          qv=qfields(k, i_qv)
          th=qfields(k, i_th)
          
          ql=qfields(k, i_ql)
          Tk=th*exner(k,ixy_inner)
          cloud_mass=qfields(k, i_ql) / cf_liquid
          if (cloud_params%l_2m) then
             cloud_number=qfields(k, i_nl) / cf_liquid
          else
             cloud_number=cloud_params%fix_N0
          end if
          
          Tc=Tk - 273.15

          ! What's the condition for ice nucleation...?
          ! This condition needs to be consistent with the mechanisms we're
          ! parametrized
          !l_condition=(( Sw(k) >= -1.e-8 .and. TdegC(k) < -8)) .or. Si(k) > 0.25

          select case(iopt_inuc)
          case default
             l_condition=(( Sw(k) >= -0.001 .and. Tc < -8 .and. Tc > -38) .or. &
                            Si(k) >= 0.08)
          case (2)
             ! Meyers, same condition as DeMott
             l_condition=( cloud_number >= nl_tidy .and. Tc < 0)
          case (4)
             ! DeMott Depletion of dust (contact and immersion)
             l_condition=( Sw(k) >= -0.001  .and. Tc < 0 .and. Tc > -38)
          case (6)
             l_condition=( cloud_number >= nl_tidy .and. Tc < 0)
          case (7)
             l_condition=( cloud_number >= nl_tidy .and. Tc < 0)
          case (8)
             l_condition=( cloud_number >= nl_tidy .and. Tc < 0)
          case (9)
             l_condition=( cloud_number >= nl_tidy .and. Tc < 0)
          case (10)
             l_condition=( cloud_number >= nl_tidy .and. Tc < 0)
             
          end select
          
          if (l_condition) then
             
             if (ice_params%l_2m)then
                ice_number=qfields(k, i_ni) / cf_ice
             else
                ice_number=1.e3 ! PRAGMATIC SM HACK
             end if
             
             dN_contact=0.0
             dN_imm=0.0
             select case(iopt_inuc)
             case default
                ! Cooper
                ! Cooper WA (1986) Ice Initiation in Natural Clouds. Precipitation Enhancement -
                ! A Scientific Challenge. Meteor Monogr, (Am Meteor Soc, Boston, MA), 21, pp 29-32.
                
                if (ql * cf_liquid > ql_small) then  !only make ice when liquid present
                   dN_imm=5.0*exp(-0.304*Tc)/rho(k,ixy_inner)
                   if (iopt_act .eq. 0) then !for fixed number concs adjust the rate
                                             !to nudge back to climatology- can be negative -
                                             !this just represents a nudging incr for cooper
                      if (l_nudge_to_cooper) then 
                        dN_imm =  (dN_imm-ice_number)*0.8 
                      else
                         dN_imm = MAX( dN_imm-ice_number, 0.0 )
                      end if 
                   endif
                endif
             case (2)
                ! Meyers
                ! Meyers MP, DeMott PJ, Cotton WR (1992) New primary ice-nucleation
                ! parameterizations in an explicit cloud model. J Appl Meteorol 31:708-721
                lws_meyers = 6.112 * exp(17.62*Tc/(243.12 + Tc))
                is_meyers  = 6.112 * exp(22.46*Tc/(272.62 + Tc))
                dN_imm     = 1.0e3 * exp(meyers_a + meyers_b *(100.0*(lws_meyers/is_meyers-1.0)))/rho(k,ixy_inner)
                dN_imm     = MAX( dN_imm-ice_number, 0.0 )
                ! Applied just for water saturation, deposition freezing ignored
                
             case (3)
                ! Fletcher NH (1962) The Physics of Rain Clouds (Cambridge Univ Press, Cambridge, UK)
                dN_imm=0.01*exp(-0.6*Tc)/rho(k,ixy_inner)
             case (4)
                ! DeMott Depletion of dust
                ! 'Predicting global atmospheric ice nuclei distributions and their impacts on climate',
                ! Proc. Natnl. Acad. Sci., 107 (25), 11217-11222, 2010, doi:10.1073/pnas.0910818107
                a_demott=5.94e-5
                b_demott=3.33
                c_demott=0.0264
                d_demott=0.0033
                Tp01=0.01-Tc
                
                if (dustphys(k)%N(1) > ni_tidy) then
                   dN_contact=1.0e3/rho(k,ixy_inner)*a_demott*(Tp01)**b_demott*                                &
                        (rho(k,ixy_inner) * m3_to_cm3 * contact_efficiency*dustphys(k)%N(1))**(c_demott*Tp01+d_demott)
                   dN_contact=min(.9*dustphys(k)%N(1), dN_contact)
                end if
                
                if (dustliq(k)%nact1 > ni_tidy) then
                   dN_imm=1.0e3/rho(k,ixy_inner)*a_demott*(Tp01)**b_demott*                                    &
                        (rho(k,ixy_inner) * m3_to_cm3 * dustliq(k)%nact1)**(c_demott*Tp01+d_demott)
                   dN_imm=immersion_efficiency*dN_imm
                   dN_imm=min(dustliq(k)%nact1, dN_imm)
                end if
                
             case (5)
                dN_imm=0.0
                dN_contact=max(0.0_wp, (dustphys(k)%N(1)-ice_number))
                
             case (6)
                ! DeMott Depletion of dust (2015)
                ! 'Integrating laborator and field data to quantify the immersion freezing ice nucleation
                !  activity of mineral dust particles', Atmos. Chem. Phys., 15, 393-409, doi:10.5194/acp-15-393-2015
                a_demott = 0.0
                b_demott = 1.25
                c_demott = 0.46
                d_demott = -11.6
                cf = 1.0     ! cf is the default callibration factor from Demott (eq 2, figures 5 and 6)
                Tp01 = 0.01 - Tc
                
                if (dustphys(k)%N(1) > ni_tidy) then
                   dN_contact=1.0e3/rho(k,ixy_inner)*cf*                                                        &
                        (rho(k,ixy_inner)*m3_to_cm3*contact_efficiency*dustphys(k)%N(1))**(a_demott*(273.16-Tk)+b_demott)*  &
                        exp(c_demott*(273.16-Tk)+d_demott)
                   dN_contact=min(.9*dustphys(k)%N(1), dN_contact)
                end if
                
                if (dustliq(k)%nact1 > ni_tidy) then
                   dN_imm=1.0e3/rho(k,ixy_inner)*cf*                                                            &
                        (rho(k,ixy_inner)*dustliq(k)%nact1*m3_to_cm3)**(a_demott*(273.16-Tk)+b_demott)*       &
                        exp(c_demott*(273.16-Tk)+d_demott)
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                end if
                
             case (7)
                ! Niemand et al. (2012) - using only insoluble in liquid!
                ! 'A particle-surface-area-based parameterization of immersion freezing on desert dust particles',
                ! J. Atmos. Sci., 69, 3077-3092, doi:10.1175/JAS-D-11-0249.1
                if (dustliq(k)%nact1 > ni_tidy) then
                   surf_area = 4*pi*dustliq(k)%nact1*rho(k,ixy_inner)*dustphys(k)%rd(aero_index%i_coarse_dust)**2* &
                        EXP(2*dustphys(k)%sigma(aero_index%i_coarse_dust)**2)      ! m2/m3
                   n_sites = EXP(-0.517*Tc+8.934)    ! 1/m2
                   dN_imm = n_sites*surf_area/rho(k,ixy_inner) ! 1/kg
                   ! AKM: this approximation is only valid for small particles with Sae,j*ns <<1 !
                   ! according to the paper for monodisperse aerosol this is ok for d<3mum and T < -30degC
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                end if
                
             case (8)
                ! Atkinson et al. (2013) - using only insoluble in liquid!
                ! 'The importance of feldspar for ice nucleation by mineral dust in mixed-phase clouds',
                ! Nature, 498, 355-358, doi:10.1038/nature12278
                if (dustliq(k)%nact1 > ni_tidy) then
                   surf_area=0.35*4*pi*dustliq(k)%nact1*rho(k,ixy_inner)*(dustphys(k)%rd(aero_index%i_coarse_dust))**2* &
                        EXP(2*dustphys(k)%sigma(aero_index%i_coarse_dust)**2) ! cm2/m3
                   ! AKM: assuming fraction of K-feldspar in insoluble dust is 0.35
                   n_sites = EXP(-1.038*Tk+275.26) !1/cm2
                   dN_imm = n_sites*surf_area/rho(k,ixy_inner)      ! 1/kg
                   ! AKM: this approximation is only valid for small particles ! (s. comment for case (7))
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                end if
                
             case (9)
                ! Tobo et al. (2013) - using only insoluble in liquid!
                ! 'Biological aerosol particles as a key determinant of ice nuclei populations in a forest
                !  ecosystem', J. Geophys. Res., 118, 10100-10110, doi:10.1002/jgrd.50801
                a_tobo = -0.074
                b_tobo = 3.8
                c_tobo = 0.414
                d_tobo = -9.671
                if (dustliq(k)%nact1 > ni_tidy) then
                   dN_imm=1.0e3/rho(k,ixy_inner)                                                    &
                        *(rho(k,ixy_inner)*m3_to_cm3*dustliq(k)%nact1)**(a_tobo*(273.16-Tk)+b_tobo) &
                        *EXP(c_tobo*(273.16-Tk)+d_tobo)
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                end if
                
             case (10)
                ! DeMott Depletion of dust - not distinguishing between insoluble in
                ! liquid and interstitial aerosol
                ! 'Predicting global atmospheric ice nuclei distributions and their impacts on climate',
                ! Proc. Natnl. Acad. Sci., 107 (25), 11217-11222, 2010, doi:10.1073/pnas.0910818107
                a_demott = 5.94e-5
                b_demott = 3.33
                c_demott = 0.0264
                d_demott = 0.0033
                Tp01 = 0.01 - Tc
                
                if ((dustliq(k)%nact1 > ni_tidy) .or. (dustphys(k)%N(1) > ni_tidy)) then
                   dN_imm=1.0e3/rho(k,ixy_inner)*a_demott*(Tp01)**b_demott*                               &
                        (rho(k,ixy_inner)*m3_to_cm3*(dustliq(k)%nact1+dustphys(k)%N(1)))**(c_demott*Tp01+d_demott)
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   ! distribute INP between interstital and activated dust (for budgeting
                   ! simulations with l_process > 0)
                   if ((dustliq(k)%nact1 > ni_tidy) .and. (dustphys(k)%N(1) > ni_tidy)) then
                      dN_contact = dN_imm - dustliq(k)%nact1
                      dN_imm = dN_imm - dN_contact
                      dN_contact=MIN(0.9*dustphys(k)%N(1), dN_contact)
                   else if (dustliq(k)%nact1 > ni_tidy) then
                      dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                      dN_contact=0.0
                   else if (dustphys(k)%N(1) > ni_tidy) then
                      dN_contact=MIN(0.9*dustphys(k)%N(1), dN_imm)
                      dN_imm=0.0
                   end if
                end if
                
             end select
             
             if (cloud_params%l_2m) dN_imm=min(dN_imm, cloud_number)
             
             dN_imm=dN_imm/dt
             dN_contact=dN_contact/dt
             dnumber=dN_imm + dN_contact
             dnumber=min(dnumber, cloud_number/dt)  !this is limit for condensation/imm fzg
             
             !convert back to gridbox mean
             dnumber=dnumber*cf_liquid
             cloud_number=cloud_number*cf_liquid
             cloud_mass=cloud_mass*cf_liquid
             
             if (dnumber > ni_tidy) then
                dmass=cloud_mass*dnumber/cloud_number
                procs(i_qi, i_inuc%id)%column_data(k)=dmass
                
                if (l_2mi) then
                   procs(i_ni, i_inuc%id)%column_data(k)=dnumber
                end if
                if (cloud_params%l_2m) then
                   procs(i_nl, i_inuc%id)%column_data(k)=-dnumber
                end if
                procs(i_ql, i_inuc%id)%column_data(k)=-dmass
                
                if (l_process) then
                   
                   ! New ice nuclei
                   dmad=dN_contact*dustphys(k)%M(1)/dustphys(k)%N(1)
                   
                   ! Frozen soluble aerosol
                   dmac=dnumber*aeroact(k)%mact1_mean*aeroact(k)%nratio1
                   
                   ! Dust already in the liquid phase
                   dmadl=dN_imm*dustliq(k)%mact1_mean*dustliq(k)%nratio1
                   
                   aerosol_procs(i_am8, i_dnuc%id)%column_data(k)=dmac
                   aerosol_procs(i_am4, i_dnuc%id)%column_data(k)=-dmac
                   aerosol_procs(i_am9, i_dnuc%id)%column_data(k)=-dmadl
                   
                   aerosol_procs(i_am7, i_dnuc%id)%column_data(k)=dmad+dmadl
                   aerosol_procs(i_am6, i_dnuc%id)%column_data(k)=-dmad    ! <WARNING: using coarse mode
                   aerosol_procs(i_an6, i_dnuc%id)%column_data(k)=-dN_contact ! <WARNING: using coarse mode
                   
                   if (l_passivenumbers_ice) then
                      ! we retain information on what'd been nucleated
                      aerosol_procs(i_an12, i_dnuc%id)%column_data(k)=dN_contact
                   end if
                end if
             end if
          end if
       end if ! l_Tcold
    enddo  ! k loop

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  end subroutine inuc
end module ice_nucleation
