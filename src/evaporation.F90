module evaporation
  use variable_precision, only: wp, iwp
  use special, only: pi, GammaFunc
  use passive_fields, only: rho, qws, TdegK, exner
  use mphys_switches, only: i_qr, i_nr, i_m3r, i_qv, i_ql, hydro_complexity, l_2mr, l_3mr, i_am2, &
       i_an2, i_am3, i_an3, i_am4, i_am5, l_aevp, l_process, active_cloud, active_rain, isol, aero_index, &
       l_separate_rain, i_am6, i_an6, i_am9, l_warm, i_qi, i_qs, i_qg, i_an11, i_an12, l_passivenumbers, &
       l_passivenumbers_ice, l_inhom_revp
  use mphys_constants, only: rhow, visair, rho0, Lv,ka, Dv, Rv, cp
  use mphys_parameters, only: c_r, vent_1, vent_2, a_r, b_r, f_r,p1, p2, p3, hydro_params, rain_params
  use process_routines, only: process_rate, i_prevp, i_arevp
  use thresholds, only: ql_small, qr_small, ss_small, qr_tidy
  use m3_incs, only: m3_inc_type2
  use distributions, only: dist_lambda, dist_mu, dist_n0
  use aerosol_routines, only: aerosol_phys, aerosol_chem, aerosol_active
  use ventfac, only: ventilation
  use which_mode_to_use, only : which_mode

#if DEF_MODEL==MODEL_KiD
  use diagnostics, only: save_dg, i_dgtime
#elif DEF_MODEL==MODEL_UM
  use diaghelp_um, only: i_here, j_here
#endif

  implicit none
  private

  public revp
contains

  subroutine revp(dt, k, qfields, aerofields, aerophys, aerochem, aeroact, dustliq, procs, aerosol_procs, l_sigevap)
    real(wp), intent(in) :: dt
    integer, intent(in) :: k
    real(wp), intent(in) :: qfields(:,:), aerofields(:,:)
    type(aerosol_phys), intent(in) :: aerophys(:)
    type(aerosol_chem), intent(in) :: aerochem(:)
    type(aerosol_active), intent(in) :: aeroact(:), dustliq(:)
    type(process_rate), intent(inout) :: procs(:,:)
    type(process_rate), intent(inout) :: aerosol_procs(:,:)
    logical, intent(out) :: l_sigevap ! Determines if there is significant evaporation

    real(wp) :: dmass, dnumber, dnumber_a, dnumber_d
    real(wp) :: m1, m2, m3, dm1, dm3

    real(wp) :: n0, lam, mu
    real(wp) :: V_r, AB

    real(wp) :: rain_mass
    real(wp) :: rain_number
    real(wp) :: rain_m3
    real(wp) :: qv

    logical :: l_rain_test ! conditional test on rain

    real(wp) :: dmac, dmac1, dmac2, dnac1, dnac2, dmacd

    l_sigevap=.false.

    qv=qfields(k, i_qv)
    rain_mass=qfields(k, i_qr)
    if (l_2mr)rain_number=qfields(k, i_nr)
    if (l_3mr)rain_m3=qfields(k, i_m3r)

    if (qv/qws(k) < 1.0-ss_small .and. qfields(k, i_ql) == 0.0 .and. rain_mass > qr_tidy) then

      m1=rain_mass/c_r
      if (l_2mr) m2=rain_number
      if (l_3mr) m3=rain_m3

      l_rain_test=.false.
      if (l_2mr) l_rain_test=rain_number>0
      if (rain_mass > qr_small .and. (.not. l_2mr .or. l_rain_test)) then
        n0=dist_n0(k,rain_params%id)
        mu=dist_mu(k,rain_params%id)
        lam=dist_lambda(k,rain_params%id)

        call ventilation(k, V_r, n0, lam, mu, rain_params)

        AB=1.0/(Lv**2/(Rv*ka)*rho(k)*TdegK(k)**(-2)+1.0/(Dv*qws(k)))
        dmass=(1.0-qv/qws(k))*V_r*AB
      else
        dmass=rain_mass/dt
      end if

      dm1=dmass/c_r
      if (l_2mr) then
        dnumber=0.0
        if (l_inhom_revp) dnumber=dm1*m2/m1
      end if
      if (l_3mr) dm3=dm1*m3/m1

      if (l_2mr) then
        if (dnumber*dt > rain_number .or. dmass*dt >= rain_mass-qr_tidy) then
          dmass=rain_mass/dt
          dnumber=rain_number/dt
          if (l_3mr) dm3=m3/dt
        end if
      end if

      procs(k, i_prevp%id)%source(i_qr)=-dmass
      procs(k, i_prevp%id)%source(i_qv)=dmass

      if (dmass*dt/rain_mass > .8) l_sigevap=.true.

      if (l_2mr) then
        procs(k, i_prevp%id)%source(i_nr)=-dnumber
      end if
      if (l_3mr) then
        procs(k, i_prevp%id)%source(i_m3r)=-dm3
      end if

      !============================
      ! aerosol processing
      !============================
      if (l_process .and. abs(dnumber) >0) then

        dmac=dnumber*aeroact(k)%nratio2*aeroact(k)%mact2_mean
        if (l_separate_rain) then
          aerosol_procs(k, i_arevp%id)%source(i_am5)=-dmac
        else
          aerosol_procs(k, i_arevp%id)%source(i_am4)=-dmac
        end if
        ! Return aerosol
        if (aero_index%i_accum >0 .and. aero_index%i_coarse >0) then
          ! Coarse and accumulation mode being used. Which one to return to?
          call which_mode(dmac, dnumber*aeroact(k)%nratio2, aerophys(k)%rd(aero_index%i_accum), &
               aerophys(k)%rd(aero_index%i_coarse), aerochem(k)%density(aero_index%i_accum),     &
               dmac1, dmac2, dnac1, dnac2)
          aerosol_procs(k, i_arevp%id)%source(i_am2)=dmac1
          aerosol_procs(k, i_arevp%id)%source(i_an2)=dnac1
          aerosol_procs(k, i_arevp%id)%source(i_am3)=dmac2
          aerosol_procs(k, i_arevp%id)%source(i_an3)=dnac2
        else
          if (aero_index%i_accum >0) then
            aerosol_procs(k, i_arevp%id)%source(i_am2)=dmac
            aerosol_procs(k, i_arevp%id)%source(i_an2)=dnumber
          end if
          if (aero_index%i_coarse >0) then
            aerosol_procs(k, i_arevp%id)%source(i_am3)=dmac
            aerosol_procs(k, i_arevp%id)%source(i_an3)=dnumber
          end if
        end if

        dmacd=dnumber*dustliq(k)%nratio2*dustliq(k)%mact2_mean
        if (.not. l_warm .and. dmacd /=0.0) then
          aerosol_procs(k, i_arevp%id)%source(i_am9)=-dmacd
          aerosol_procs(k, i_arevp%id)%source(i_am6)=dmacd
          aerosol_procs(k, i_arevp%id)%source(i_an6)=dnumber*dustliq(k)%nratio2
        end if

        if (l_passivenumbers) then
          dnumber_a=-dnumber*aeroact(k)%nratio2
          aerosol_procs(k, i_arevp%id)%source(i_an11)=dnumber_a
        end if

        if (l_passivenumbers_ice) then
          dnumber_d=-dnumber*dustliq(k)%nratio2
          aerosol_procs(k, i_arevp%id)%source(i_an12)=dnumber_d
        end if
      end if
    end if
  end subroutine revp
end module evaporation
