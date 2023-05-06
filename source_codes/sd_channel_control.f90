      subroutine sd_channel_control

      use sd_channel_module
      use channel_velocity_module
      use basin_module
      use hydrograph_module
      use constituent_mass_module
      use conditional_module
      use channel_data_module
      use channel_module
      use ch_pesticide_module
      use climate_module
      use water_body_module
      use time_module
    
      implicit none     
    
      !real :: rcharea                !m^2           |cross-sectional area of flow
      real :: flo_rt                  !m^3/s         |flow rate in reach for day
      integer :: isd_db               !              |
      integer :: iob                  !              |
      integer :: idb                  !none          |channel data pointer
      integer :: ihyd                 !              |
      integer :: ipest                !              |
      integer :: ihru                 !              |
      integer :: iru                  !              |
      integer :: ise                  !              |
      integer :: ielem                !              |
      integer :: id
      real :: ebtm_m                  !m             |erosion of bottom of channel
      real :: ebank_m                 !m             |meander cut on one side
      real :: erode_bank_cut          !cm            |widening caused by downcutting (both sides)
      real :: ebtm_t                  !tons          |bottom erosion
      real :: ebank_t                 !tons          |bank erosion
      real :: sedout                  !mg		     |sediment out of waterway channel
      real :: washld                  !tons          |wash load  
      real :: bedld                   !tons          |bed load
      real :: dep                     !tons          |deposition
      real :: hc_sed                  !tons          |headcut erosion
      real :: chside                  !none          |change in horizontal distance per unit
                                      !              |change in vertical distance on channel side
                                      !              |slopes; always set to 2 (slope=1/2)
      real :: a                       !m^2           |cross-sectional area of channel
      real :: b                       !m             |bottom width of channel
      real :: c                       !none          |inverse of channel side slope
      real :: p                       !m             |wetting perimeter

      real :: rh                      !m             |hydraulic radius
      real :: qman                    !m^3/s or m/s  |flow rate or flow velocity
      real :: frac                    !0-1           |fraction of hydrograph 
      real :: valint                  !              | 
      integer :: ivalint              !              |
      real :: tbase                   !none          |flow duration (fraction of 24 hr)
      real :: tb_pr                   !              |
      real :: tb                      !              |
      real :: vol_ovb                 !              |
      real :: const                   !              |
      integer :: ics                  !none          |counter
      real :: ob_const                !              |
      integer :: ii                   !none          |counter
      real :: sum_vol                 !              |
      real :: xx                      !              | 
      integer :: ic                   !              |
      real :: vol_overmx              !              |
      real :: flood_dep               !              | 
      real :: dep_e                   !              |
      real :: rto                     !none          |cloud cover factor 
      real :: e_btm                   !cm            |erosion on bottom of channel at each time step
      real :: dep_btm                 !cm            |deposition on bottom of channel
      real :: sumtime                 !              |
      real :: vc                      !m/s           |flow velocity in reach
      real :: pr_ratio                !              |
      real  :: tw                     !              |
      real :: tave                    !              |
      real :: shear_btm_cr            !              |
      real :: shear_btm_adj           !none          |take out bedld_cap adjustment
      real :: shear_btm               !              |
      real :: shear_bank_cr           !              | 
      real :: qmm                     !              | 
      real :: qh                      !              | 
      real :: hc                      !m/yr          |head cut advance
      integer :: max                  !              |
      real :: chns                    !              |
      integer :: ihval                !none          |counter 
      real :: bedld_cap               !              |
      real :: perim_bed               !              |
      real :: vol
      real :: perim_bank              !              |
      real :: s_bank                  !              |
      real :: shear_bank              !              |
      real :: shear_bank_adj          !              | 
      real :: e_bank                  !              | 
      real :: perc                    !              |
      real :: w_btm                   !m             |channel bottom width
      integer :: iaq
      integer :: iaq_ch
      real :: det                     !hr            |time step
      real :: scoef                   !none          |Storage coefficient
      integer :: dum                  !rtb gwflow
      real :: flo_ls
      real :: channel_storage         !rtb gwflow
      real :: rchvol
      real :: str_pow                 !              |stream power
      real :: str_pow_cr              !              |critical stream power
      real :: d50_m                   !m             |d50 converted to meters
      real :: vel, cohes, vel_cr, b_coef, sedcap_ppm, qcms, veg, sedcap_t, sedin_t
      real :: rad_curv, cla, pk_rto, vel_bend, vel_rch
      real :: arc_len, hyd_radius, prot_len
      
      ich = isdch
      isd_db = sd_dat(ich)%hyd
      iwst = ob(icmd)%wst
      ebtm_m = 0.
      ebank_m = 0.
      erode_bank_cut = 0.
      ebtm_t = 0.
      ebank_t = 0.
      sedout = 0.
      washld = 0.
      bedld = 0.
      dep = 0.
      hc = 0.
      hc_sed = 0.
      
      !call ch_rtmusk
      !call ch_rthr
              
      !call sd_channel_sediment (time%step)
        
      !! set ht1 to incoming hydrograph
      ht1 = ob(icmd)%hin
      !! add water transfer
      if (ob(icmd)%trans%flo > 1.e-6) then
        ht1 = ht1 + ob(icmd)%trans
        ob(icmd)%trans = hz
      end if
      
      chsd_d(ich)%flo_in = ht1%flo / 86400.     !flow for morphology output
      ch_in_d(ich) = ht1                        !set inflow om hydrograph
      ch_in_d(ich)%flo = ht1%flo / 86400.       !flow for om output
      
      !rtb hydrograph separation
      hdsep1%flo_surq = ob(icmd)%hdsep_in%flo_surq
      hdsep1%flo_latq = ob(icmd)%hdsep_in%flo_latq
      hdsep1%flo_gwsw = ob(icmd)%hdsep_in%flo_gwsw
      hdsep1%flo_swgw = ob(icmd)%hdsep_in%flo_swgw
      hdsep1%flo_satex = ob(icmd)%hdsep_in%flo_satex
      hdsep1%flo_satexsw = ob(icmd)%hdsep_in%flo_satexsw
      hdsep1%flo_tile = ob(icmd)%hdsep_in%flo_tile
      !rtb hydrograph separation
      
      !! adjust precip and temperature for elevation using lapse rates
      w = wst(iwst)%weat
      if (bsn_cc%lapse == 1) call cli_lapse
      wst(iwst)%weat = w
      ht1%temp = 5.0 + 0.75 * wst(iwst)%weat%tave
      wtemp = 5.0 + 0.75 * wst(iwst)%weat%tave

      !! if connected to aquifer - add flow
      if (sd_ch(ich)%aqu_link > 0) then
        iaq = sd_ch(ich)%aqu_link
        iaq_ch = sd_ch(ich)%aqu_link_ch
        if (aq_ch(iaq)%ch(iaq_ch)%flo_fr > 0.) then
          chsd_d(ich)%aqu_in = (aq_ch(iaq)%ch(iaq_ch)%flo_fr * aq_ch(iaq)%hd%flo) / 86400.
          chsd_d(ich)%aqu_in_mm = (aq_ch(iaq)%ch(iaq_ch)%flo_fr * aq_ch(iaq)%hd%flo) / (10. * ob(icmd)%area_ha)
          ht1 = ht1 + aq_ch(iaq)%ch(iaq_ch)%flo_fr * aq_ch(iaq)%hd
          aq_ch(iaq)%ch(iaq_ch)%flo_fr = 0.
        end if
      end if
      hcs1 = obcs(icmd)%hin
      
      !! set inflow hyds for printing
      chsd_d(ich)%flo_in = ht1%flo / 86400.     !flow for morphology output - m3/s
      chsd_d(ich)%flo_in_mm = ht1%flo / (10. * ob(icmd)%area_ha)   !flow in mm
      ch_in_d(ich) = ht1                        !set inflow om hydrograph
      ch_in_d(ich)%flo = ht1%flo / 86400.       !flow for om output - m3/s
      
      !! set outgoing flow and sediment - ht2
      ht2 = hz
      peakrate = 0.
      hyd_rad = 0.
      timeint = 0.

      !assume triangular hydrograph
      peakrate = 2. * ht1%flo / (1.5 * 60.)     !sd_chd(isd_db)%tc)
      peakrate = peakrate / 60.   !convert min to sec
      if (peakrate > 1.e-9) then
         !! compute changes in channel dimensions
          chside = sd_ch(ich)%chss
          b = sd_ch(ich)%chw
          sd_ch_vel(ich)%wid_btm = b
          sd_ch_vel(ich)%dep_bf = sd_ch(ich)%chd

          !! compute flow and travel time at bankfull depth
          p = b + 2. * sd_ch(ich)%chd * Sqrt(chside * chside + 1.)
          a = b * sd_ch(ich)%chd + chside * sd_ch(ich)%chd * sd_ch(ich)%chd
          rh = a / p
          sd_ch_vel(ich)%area = a
          sd_ch_vel(ich)%vel_bf = Qman (a, rh, sd_ch(ich)%chn, sd_ch(ich)%chs)
  
          IF (peakrate > sd_ch_vel(ich)%vel_bf) THEN
          !! OVERBANK FLOOD
            
            !rtb floodplain
            !flood_freq(ich) = 1 !flag to indicate the water is in the floodplain
            
            sd_ch(ich)%overbank = "ob"
            rcharea = sd_ch_vel(ich)%area
            rchdep = sd_ch(ich)%chd
            !calculate hydraulic radius at hydrograph time increments for degredation
            flo_rt = 0.
            rchdep = 0.
            valint = 1. / float(maxint)
            ivalint = 1
            tbase = 1.5 * 3600.     !sd_chd(isd_db)%tc * 60.   !hydrograph base time in s
            tb_pr = tbase
            DO WHILE (flo_rt < peakrate)
              rchdep = rchdep + 0.01
              rcharea = (sd_ch_vel(ich)%wid_btm + chside * rchdep) * rchdep
              p = sd_ch_vel(ich)%wid_btm + 2. * rchdep * Sqrt (1. + chside *chside)
              rh = rcharea / p
              flo_rt = Qman(rcharea, rh, sd_ch(ich)%chn, sd_ch(ich)%chs)
              !need to save hydraulic radius and time for each flow interval for downcutting and widening
              if (flo_rt > valint * peakrate) then
                hyd_rad(ivalint) = rh
                tb = (peakrate - flo_rt) * tbase / peakrate
                timeint(ivalint) = (tb_pr - tb) / 3600.   !sec to hr
                tb_pr = tb
                ivalint = ivalint + 1
                valint = float (ivalint) / float (maxint)
              end if
            END DO
            
            !! estimate overbank flow - assume a triangular hyd
            tbase = 1.5 * 3600.     !sd_chd(isd_db)%tc * 60.  !seconds
            vol_ovb = 0.5 * (peakrate - sd_ch_vel(ich)%vel_bf) * sd_ch_vel(ich)%vel_bf / peakrate * tbase
            vol_ovb = amin1(vol_ovb, ht1%flo)
            vol_ovb = peakrate - sd_ch_vel(ich)%vel_bf
            const = vol_ovb / peakrate
            ob(icmd)%hd(3) = const * ob(icmd)%hin
            
            !find current total flood volume (ht1)
            ics = ob(icmd)%props2
            if (ics > 0) then   ! flood elements are specified - link to surface elements
            ht1 = hz
            ob_const = const
            do ii = 1, ch_sur(ics)%num
              ht1 = ht1 + ch_sur(ics)%hd(ii)
            end do
            
            !add current and new flood volumes
            ht1 = ht1 + ob(icmd)%hd(3)
           if (ht1%flo > ch_sur(ics)%flood_volmx(0)) then
            !calc flood depth above channel bottom (flood_dep)
            sum_vol = 0.
            do ii = 1, ch_sur(ics)%num
              if (ht1%flo < ch_sur(ics)%flood_volmx(ii)) then
                !solve quadrative for depth above base of current element
                a = sd_ch(ich)%chl * 1000. / sd_ch(ich)%chs
                b = ch_sur(ics)%wid(ii-1) * sd_ch(ich)%chl * 1000.
                c =  ch_sur(ics)%flood_volmx(ii-1) - ht1%flo
                xx = b ** 2 - 4. * a * c
                dep = (-b + sqrt(xx)) / (2. * a)
                dep = Max(0., dep)
                ic = ii
                exit
              end if
              if (ii == ch_sur(ics)%num) then
                !flood is over the  max storage - assume linear upward wall
                ic = ch_sur(ics)%num
                vol_overmx = ht1%flo - ch_sur(ics)%flood_volmx(ii)
                dep = ch_sur(ics)%dep(ic) + (vol_overmx / (sd_ch(ich)%chl *  &
                               (sd_ch(ich)%chw + 2. * ch_sur(ics)%wid(ic))))
              end if
            end do
            !calc flood depth above channel bottom
            flood_dep = dep + ch_sur(ics)%dep(ic-1)
            !calc new flood volume for each element
            do ii = 1, ch_sur(ics)%num
              !calc depth in the element
              
              if (flood_dep < ch_sur(ics)%dep(ii-1)) then
                !no flooding on element
                ch_sur(ics)%hd(ii)%flo = 0.
              else if (flood_dep < ch_sur(ics)%dep(ii)) then
                !flood level within the element
                dep_e = flood_dep - ch_sur(ics)%dep(ii-1)
                ch_sur(ics)%hd(ii)%flo = dep_e ** 2 / sd_ch(ich)%chs * sd_ch(ich)%chl
              else
                !flood level over max element depth
                ch_sur(ics)%hd(ii)%flo = 2. * sd_ch(ich)%chw *             & 
                  (flood_dep - ch_sur(ics)%dep(ii)) + sd_ch(ich)%chw *     &
                  (ch_sur(ics)%dep(ii) - ch_sur(ics)%dep(ii-1))
              end if
              ch_sur(ics)%hd(ii)%flo = amin1 (ch_sur(ics)%hd(ii)%flo, ch_sur(ics)%flood_volmx(ii))
              sum_vol = sum_vol + ch_sur(ics)%hd(ii)%flo
            end do
            !determine fraction of total flood volume and set volume for each element
            rto = sum_vol / ht1%flo  !ensure water balance is maintained
            do ii = 1, ch_sur(ics)%num
              const = rto * ch_sur(ics)%hd(ii)%flo / ht1%flo
              ch_sur(ics)%hd(ii) = const * ht1
              !! add overbank to wetland storage of hru(s)
              if (ch_sur(ics)%obtyp(ii) == "hru") then
                ihru = ch_sur(ics)%obtypno(ii)
                wet(ihru) = wet(ihru) + ch_sur(ics)%hd(ii)
              end if
              if (ch_sur(ics)%obtyp(ii) == "ru") then
                iru = ch_sur(ics)%obtypno(ii)
                do ielem = 1, ru_def(iru)%num_tot
                  ise = ru_def(iru)%num(ielem)
                  !! assume the element is an hru
                  ihru = ru_elem(ise)%obtypno
                  wet(ihru) = wet(ihru) + ru_elem(ise)%frac * ch_sur(ics)%hd(ii)
                end do
              end if
            end do
           end if
           end if
          ELSE
            !! CHANNELIZED FLOW
            !! find the crossectional area and depth for volrt
            !! by iteration method at 1cm interval depth
            !! find the depth until the discharge rate is equal to volrt
            !zero overbank flow
            sd_ch(ich)%overbank = "ib"
            ob(icmd)%hd(3) = hz
            ob_const = 1.
            flo_rt = 0.
            rchdep = 0.
            sumtime = 0.
            valint = 1. / float(maxint)
            ivalint = 1
            tbase = 1.5 * 3600.     !sd_chd(isd_db)%tc * 60.   !hydrograph base time in s
            tb_pr = tbase
            DO WHILE (flo_rt < peakrate)
              rchdep = rchdep + 0.01
              rcharea = (sd_ch_vel(ich)%wid_btm + chside * rchdep) * rchdep
              p = sd_ch_vel(ich)%wid_btm + 2. * rchdep * Sqrt(1. + chside * chside)
              rh = rcharea / p
              !rh = Qman(rcharea, rh, sd_ch(ich)%chn, sd_ch(ich)%chs)
              flo_rt = Qman(rcharea, rh, sd_ch(ich)%chn, sd_ch(ich)%chs)
              !need to save hydraulic radius and time for each flow interval for downcutting and widening
              if (flo_rt > valint * peakrate) then
                hyd_rad(ivalint) = rh
                tb = (peakrate - flo_rt) * tbase / peakrate
                timeint(ivalint) = (tb_pr - tb) / 3600.   !sec to hr
                sumtime = sumtime + timeint(ivalint)
                tb_pr = tb
                ivalint = ivalint + 1
                valint = float (ivalint) / float(maxint)
              end if
            END DO
            timeint = timeint / sumtime
          END IF

        !! adjust peak rate for headcut advance -also adjusts CEAP gully from
        !! edge-of-field to trib (assuming rectangular shape and constant tc)
        pr_ratio = (sd_ch(ich)%chl - sd_ch(ich)%hc_len / 1000.) / sd_ch(ich)%chl
        pr_ratio = Max(pr_ratio, 0.)

        go to 777
        !! new q*qp (m3 * m3/s) equation for entire runoff event
        qmm = ht1%flo / (10. * ob(icmd)%area_ha)
        if (qmm > 3.) then
          sd_ch(ich)%hc_hgt = 0.    !gully erosion is turned off - 7/28/2021
          qh = (ht1%flo / 86400.) ** .5 * sd_ch(ich)%hc_hgt ** .225
          hc = sd_ch(ich)%hc_co * qh            !m per event
          hc = Max(hc, 0.)
          sd_ch(ich)%hc_len = sd_ch(ich)%hc_len + hc
          if (sd_ch(ich)%hc_len > sd_ch(ich)%chl * 1000.) then
            hc = hc - (sd_ch(ich)%hc_len - sd_ch(ich)%chl * 1000.)
            sd_ch(ich)%hc_len = sd_ch(ich)%chl * 1000.
          end if
            
          !! compute sediment yield from headcut- assume bd = 1.2 t/m3
          !! assume channel dimensions are same as data file
          hc_sed = hc * sd_ch(ich)%chw * sd_ch(ich)%chd * 1.2
        end if
        
        !! break hydrograph into maxint segments and compute deg at each flow increment
        do ihval = 1, maxint
          !! calc critical shear and shear on bottom of channel
          shear_btm_cr = sd_ch(ich)%d50
          shear_btm = 9800. * hyd_rad(ihval) * sd_ch(ich)%chs   !! Pa = N/m^2 * m * m/m
            
          !! degradation of the bank (widening)
          perim_bank = 2. * ((sd_ch(ich)%chd ** 2) * (1. + sd_ch(ich)%chss ** 2)) ** 0.5
          perim_bed = sd_ch(ich)%chw
          tw = perim_bed + 2. * sd_ch(ich)%chss * rchdep
          s_bank = 1.77 * (perim_bed / perim_bank + 1.5) ** (-1.4)
          !! assume bank shear is 75% of bottom shear
          shear_bank = shear_btm * 0.75     !sd_ch(ich)%shear_bnk * s_bank * (tw * perim_bed) / (2. * perim_bank)
          if (sd_ch(ich)%ch_clay >= 10.) then
            chns = .0156
          else
            chns = (sd_ch(ich)%d50 / 25.4) ** .16666 / 39.
          end if
          shear_bank_adj = shear_bank * (1. - sd_ch(ich)%cov)      !* (chns / sd_chd(isd_db)%chn) ** 2
          shear_bank_cr = 0.493 * 10. ** (.0182 * sd_ch(ich)%ch_clay)
          e_bank = 0.
          if (shear_bank_adj > shear_bank_cr) then
            e_bank = timeint(ihval) * sd_ch(ich)%cherod * (shear_bank_adj - shear_bank_cr)    !! cm = hr * cm/hr/Pa * Pa
            ebank_m = ebank_m + e_bank
            !! calc mass of sediment eroded -> t = cm * m/100cm * width (m) * length (km) * 1000 m/km * bd (t/m3)
            !! apply to only one side (perim_bank / 2.)
            ebank_t = ebank_t + 10. * e_bank * perim_bank / 2. * sd_ch(ich)%chl * sd_ch(ich)%ch_bd
          end if
              
          !! no downcutting below equilibrium slope
          e_btm = 0.
          erode_bank_cut = 0.
          if (sd_ch(ich)%chs > sd_ch(ich)%chseq) then
            !! if bottom shear > d50 -> downcut - widen to maintain width depth ratio
            if (shear_btm > shear_btm_cr) then
              e_btm = timeint(ihval) *  sd_ch(ich)%cherod * (shear_btm - shear_btm_cr)    !! cm = hr * cm/hr/Pa * Pa
              !! if downcutting - check width depth ratio to see if widens
              !if (sd_ch(ich)%chw / sd_ch(ich)%chd < sd_ch(ich)%wd_rto) then
              !  erode_bank_cut = e_btm * sd_ch(ich)%wd_rto
                !! appy to both bank sides
              !  ebank_t = ebank_t + 10. * erode_bank_cut * perim_bank * sd_ch(ich)%chl * sd_ch(ich)%ch_bd
              !end if
              ebtm_m = ebtm_m + e_btm
              !! calc mass of sediment eroded -> t = cm * m/100cm * width (m) * length (km) * 1000 m/km * bd (t/m3)
              ebtm_t = ebtm_t + 10. * e_btm * perim_bed * sd_ch(ich)%chl * sd_ch(ich)%ch_bd
            end if
          end if

        end do    ! ihval
          
        !go to 777
        !! Peters latest channel erosion model
        !!vel = 1.37 * (sd_ch(ich)%chs ** 0.31) * (12. * sd_ch(ich)%chw) ** 0.32      !annual ave for SWIFT
        !! mean daily to peak ratio developed from GARDAY - THE STUDY OF MOST PROBABLE MEAN DAILY BANKFULL RUNOFF VOLUMES 
        !! IN SMALL WATERSHEDS DOMINATED BY CONVECTIVE/FRONTAL CHANNEL FORMING EVENTS AND THE CO-INCIDENT INNER BERM CHANNELS � PART I. 
        !! Another eq from Peter - Qmax=Qmean*(1+2.66*Drainage Area^-.3)
  777   pk_rto = 0.2 + 0.5 / 250. * ob(icmd)%area_ha
        pk_rto = amin1 (1., pk_rto)
        peakrate = ht1%flo / 86400. / pk_rto    !m3/s
        
        !! interpolate rating curve using peak rate
        call rcurv_interp_flo (ich, peakrate)
        
        !! use peakrate as flow rate
        qcms = peakrate
        vel = qcms / rcurv%xsec_area
        
        cla = sd_ch(ich)%ch_clay
        
        !! check Peters calcs
        !vel = 1.42
        !sd_ch(ich)%sinu = 1.3
        !sd_ch(ich)%chw = 20.    !m  ***jga
        !rcurv%dep = 2.     !m  ***jga
        !sd_ch(ich)%chs = .004   !m/m ***jga
        !cla = 25.
        
        !! calc soil cohesion
        if (cla - 7.123 > 0.001) then
          !cohes = 177.5 * (cla - 7.123) ** 1.08
          cohes = -87. + (42.8 * cla) - (0.26 * cla ** 2) + (0.029 * cla ** 3)
        else
          cohes = 0.
        end if
        cohes = amax1 (0., cohes)
        veg = 0.    !Pa 200-10000.
        if (1.177 + 0.021 * cohes + veg > 0.001) then
          !! vel_cr = log10(8.8 * rcurv%dep / 0.004) * (0.0004 * ((1377.-1000.) * 9.81 * 0.004 + 0.021 * cohes + veg)) ** 0.5
          vel_cr = log10(8.8 * rcurv%dep / 0.004) * (0.0004 * (14.79 + 0.021 * cohes + veg)) ** 0.5
        else
          vel_cr = 0.
        end if
        
        !! calculate radius of curvature
        rad_curv = ((12. * sd_ch(ich)%chw) * sd_ch(ich)%sinu ** 1.5) / (13. * (sd_ch(ich)%sinu -1.) ** 0.5)
        vel_bend = log10(rad_curv / sd_ch(ich)%chw)
        vel_bend = vel * (1.74 - 0.52 * log10(rad_curv / sd_ch(ich)%chw))
        vel_rch = 0.33 * vel_bend + 0.66 * vel
        b_coef = 2.2    ! 8.53 / (1. + exp((-29.11 + 1.8 * cla) ** 0.0248))
        if (vel_cr > 1.e-6 .and. rcurv%dep / sd_ch(ich)%chd > 0.25) then
          ebank_m = 0.00233 * (vel_rch / vel_cr) ** b_coef    !bank erosion m/yr
        else
          ebank_m = 0.
        end if
        !! calc mass of sediment eroded -> t = bankcut (m) * depth (m) * lengthcut (m) * bd (t/m3)
        !! arc length = 0.33 * meander wavelength * sinuosity  -> protected length 
        arc_len = 0.33 *  (12. * sd_ch(ich)%chw) * sd_ch(ich)%sinu
        hyd_radius = rcurv%xsec_area / rcurv%wet_perim
        prot_len = 0.71 * (hyd_radius ** 1.1666) / sd_ch(ich)%chn
        ebank_t = ebank_m * sd_ch(ich)%chd * (arc_len + prot_len) * sd_ch(ich)%ch_bd
        ebank_t = amax1 (0., ebank_t)
        
        !! bed erosion and limiting transport capacity !!vel_cr = 0.35 * (sd_ch(ich)%d50 / 1000.) ** 0.45  !d50 in m
        !sd_ch(ich)%d50 = 10.    !mm ***jga
        !sd_ch(ich)%chw = 15.    !m  ***jga
        !sd_ch(ich)%chs = .004   !m/m ***jga
        !qcms = 45.              !cms
        d50_m = sd_ch(ich)%d50 / 1000.
        str_pow = 9810. * qcms * sd_ch(ich)%chs / sd_ch(ich)%chw !stream power
        str_pow_cr = 100. * (16.2 * d50_m) ** 1.5
        if (str_pow - str_pow_cr > 1.e-6) then
          sedcap_ppm = 0.0214 * (str_pow - str_pow_cr) ** 1.5 / d50_m * (qcms / sd_ch(ich)%chw) ** -0.833
        else
          sedcap_ppm = 0.
        end if
        sedcap_t = sedcap_ppm * ht1%flo / 1000000.  !t = ppm (g/m3) * 1 t/m3 / 1000000. g/t
        !! concentration of incoming and bank erosion
        sedin_t = ht1%sed + ebank_t
        ebtm_t = sedcap_t - sedin_t
        if (ebtm_t > 0.) then
          dep = 0.
        else
          dep = -ebtm_t
          ebtm_t = 0.
        end if
        !! check if bottom width (w_btm) is < 0
        w_btm = sd_ch(ich)%chw - 2. * sd_ch(ich)%chd * sd_ch(ich)%chss
        if (w_btm < 0.01) then
          w_btm = .5 * sd_ch(ich)%chw
          w_btm = Max(0., b)
          sd_ch(ich)%chss = (sd_ch(ich)%chw - w_btm) / (2. * sd_ch(ich)%chd)
        end if
        !! m = t / (t/m3 * m * km * 1000 m/km)
        ebtm_m = ebtm_t / (1.3 * w_btm * sd_ch(ich)%chl * 1000.)
        !erode_bank_cut = amax1 (0., erode_bank_cut)
          
          !! adjust for incoming bedload and compute deposition
          !! assume bedload is deposited
          !dep = sd_ch(ich)%bedldcoef * ht1%sed
          !dep_btm = dep / (10. * perim_bed * sd_ch(ich)%chl * sd_ch(ich)%ch_bd)
 !777     !ebtm_m = ebtm_m ! - dep_btm      !don't add in all bedload (most will be transported out)
          sd_ch(ich)%chd = sd_ch(ich)%chd + ebtm_m !/ 100.
          if (sd_ch(ich)%chd < 0.) then
            !! stream is completely filled in
            sd_ch(ich)%chd = 0.01
          end if
          
          !sd_ch(ich)%chw = sd_ch(ich)%chw + ebank_m
          !sd_ch(ich)%chs = sd_ch(ich)%chs - (ebtm_m) / (sd_ch(ich)%chl * 1000.)
          !sd_ch(ich)%chs = amax1 (sd_ch(ich)%chseq, sd_ch(ich)%chs)

      !if (ich==422) then
      !    write (888,*) time%day, time%yrs, ebank_m, sd_ch(ich)%chw, ebtm_m, sd_ch(ich)%chd
      !end if
      
        !end if

      !! compute sediment leaving the channel
      
	  washld = (1. - sd_ch(ich)%bedldcoef) * ht1%sed
      washld = amax1 (0., washld)
	  sedout = washld + hc_sed + ebtm_t + ebank_t
      dep = ht1%sed - sedout
      dep = amax1 (0., dep)
      
      !! set values for outflow hydrograph
      !! calculate flow velocity and travel time  ht2 = ht1   !***jga 
      ht2%flo = ht1%flo
      idb = ob(icmd)%props
      jrch = ich
      vc = 0.001
      if (rcharea > 1.e-4 .and. ht1%flo > 1.e-4) then
        vc = peakrate / rcharea
        !if (vc > sd_ch_vel(ich)%celerity_bf) vc = sd_ch_vel(ich)%celerity_bf
        rttime = sd_ch(ich)%chl * 1000. / (3600. * vc)
        if (time%step == 0) rt_delt = 1.
        !if (bsn_cc%wq == 1) then
          !! use modified qual-2e routines
          ht3 = ht1
          !! total incoming to output to SWIFT
          ob(icmd)%hin_tot = ob(icmd)%hin_tot + ht1
          !! convert mass to concentration
          call hyd_convert_mass_to_conc (ht3)
          jnut = sd_dat(ich)%nut
          ben_area = sd_ch(ich)%chw * sd_ch(ich)%chl
          !! convert storage hyd - mass to concentration
          if (ch_stor(ich)%flo > 0.001) then
            call hyd_convert_mass_to_conc (ch_stor(ich))    !***jga 
          else
            ch_stor(ich) = hz
          end if
          
          call ch_watqual4  !***jga 
          !! convert outflow and storage hyds - concentration to mass
          call hyd_convert_conc_to_mass (ht2)      !***jga 
          !call hyd_convert_conc_to_mass (ht3)      !***jga 
          call hyd_convert_conc_to_mass (ch_stor(ich))     !***jga
          !ht2 = ht1
         
          !! compute nutrient losses using 2-stage ditch model
          !call sd_channel_nutrients
        !end if

        !! reset sed to tons
        ht2%sed = sedout
        
        !! route constituents
        call ch_rtpest
        !! call mike winchell's new routine for pesticide routing
        ! call ch_rtpest2
        call ch_rtpath
      end if
      
      end if    ! peakrate > 0

      flo_ls = ht2%flo
      !! compute water balance - precip, evap and seep
      !! km * m * 1000 m/km * ha/10000 m2 = ha
      ch_wat_d(ich)%area_ha = sd_ch(ich)%chl * sd_ch(ich)%chw / 10.
      !! m3 = 10. * mm * ha
      ch_wat_d(ich)%precip = 10. * wst(iwst)%weat%precip * ch_wat_d(ich)%area_ha
      ch_wat_d(ich)%evap = 10. * bsn_prm%evrch * wst(iwst)%weat%pet * ch_wat_d(ich)%area_ha
      ch_wat_d(ich)%seep = 10. * sd_ch(ich)%chk * ch_wat_d(ich)%area_ha      !k units to mm/d
      !ch_wat_d(ich)%seep = 0.
      
      !! add precip
      ht2%flo = ht2%flo + ch_wat_d(ich)%precip
      
      !! subtract seepage
      if (ht2%flo < ch_wat_d(ich)%seep) then
        ch_wat_d(ich)%seep = ht2%flo
        ht2%flo = 0.
      else
        ht2%flo = ht2%flo - ch_wat_d(ich)%seep
      end if
      
      !! subtract evaporation
      if (ht2%flo < ch_wat_d(ich)%evap) then
        ch_wat_d(ich)%evap = ht2%flo
        ht2%flo = 0.
      else
        ht2%flo = ht2%flo - ch_wat_d(ich)%evap
      end if

      !! adjust entire hydrograph for losses
      if (flo_ls > 1.e-6) then
        rto = ht2%flo / flo_ls
        if (rto < 1.0) ht2 = rto * ht2
        ht2%flo = flo_ls
      end if
      
      !! total outgoing to output to SWIFT
      ob(icmd)%hout_tot = ob(icmd)%hout_tot + ht2
        
      !! calculate hydrograph leaving reach and storage in channel
      !if (time%step == 0) rt_delt = 1.
      rt_delt = 1.
      det = 24.* rt_delt
      scoef = bsn_prm%scoef * det / (rttime + det)
      scoef = amax1 (0., scoef)
      scoef = amin1 (1., scoef)
      frac = 1. - scoef
      !scoef = 1.
      !frac = 0.
      !if (rttime > det) then      ! ht1 = incoming + storage
      !  !! travel time > timestep -- then all incoming is stored and frac of stored is routed
      !  ht2 = scoef * ch_stor(ich)
      !  ch_stor(ich) = frac * ch_stor(ich) + ht1
      !  hcs2 = scoef * ch_water(ich)
      !  ch_water(ich) = frac * ch_water(ich) + hcs1
      !else
        !! travel time < timestep -- route all stored and frac of incoming
        ht3 = ht2 + ch_stor(ich)
        ht2 = scoef * ht3
        ch_stor(ich) = ht3 - ht2    !incoming + initial storage - outflow
        
        hcs2 = scoef * hcs1
        hcs2 = hcs2 + ch_water(ich)
        ch_water(ich) = frac * hcs1
      !end if
      
      !rtb hydrograph separation
      if (rttime > det) then      ! ht1 = incoming + storage
        !! travel time > timestep -- then all incoming is stored and frac of stored is routed
        hdsep2%flo_surq = scoef * ch_stor_hdsep(ich)%flo_surq
        hdsep2%flo_latq = scoef * ch_stor_hdsep(ich)%flo_latq
        hdsep2%flo_gwsw = scoef * ch_stor_hdsep(ich)%flo_gwsw
        hdsep2%flo_swgw = scoef * ch_stor_hdsep(ich)%flo_swgw
        hdsep2%flo_satex = scoef * ch_stor_hdsep(ich)%flo_satex
        hdsep2%flo_satexsw = scoef * ch_stor_hdsep(ich)%flo_satexsw
        hdsep2%flo_tile = scoef * ch_stor_hdsep(ich)%flo_tile
        ch_stor_hdsep(ich)%flo_surq = (frac*ch_stor_hdsep(ich)%flo_surq) + hdsep1%flo_surq
        ch_stor_hdsep(ich)%flo_latq = (frac*ch_stor_hdsep(ich)%flo_latq) + hdsep1%flo_latq
        ch_stor_hdsep(ich)%flo_gwsw = (frac*ch_stor_hdsep(ich)%flo_gwsw) + hdsep1%flo_gwsw
        ch_stor_hdsep(ich)%flo_swgw = (frac*ch_stor_hdsep(ich)%flo_swgw) + hdsep1%flo_swgw
        ch_stor_hdsep(ich)%flo_satex = (frac*ch_stor_hdsep(ich)%flo_satex) + hdsep1%flo_satex
        ch_stor_hdsep(ich)%flo_satexsw = (frac*ch_stor_hdsep(ich)%flo_satexsw) + hdsep1%flo_satexsw
        ch_stor_hdsep(ich)%flo_tile = (frac*ch_stor_hdsep(ich)%flo_tile) + hdsep1%flo_tile
      else
        !! travel time < timestep -- route all stored and frac of incoming
        hdsep2%flo_surq = scoef * hdsep1%flo_surq
        hdsep2%flo_latq = scoef * hdsep1%flo_latq
        hdsep2%flo_gwsw = scoef * hdsep1%flo_gwsw
        hdsep2%flo_swgw = scoef * hdsep1%flo_swgw
        hdsep2%flo_satex = scoef * hdsep1%flo_satex
        hdsep2%flo_satexsw = scoef * hdsep1%flo_satexsw
        hdsep2%flo_tile = scoef * hdsep1%flo_tile
        hdsep2%flo_surq = hdsep2%flo_surq + ch_stor_hdsep(ich)%flo_surq
        hdsep2%flo_latq = hdsep2%flo_latq + ch_stor_hdsep(ich)%flo_latq
        hdsep2%flo_gwsw = hdsep2%flo_gwsw + ch_stor_hdsep(ich)%flo_gwsw
        hdsep2%flo_swgw = hdsep2%flo_swgw + ch_stor_hdsep(ich)%flo_swgw
        hdsep2%flo_satex = hdsep2%flo_satex + ch_stor_hdsep(ich)%flo_satex
        hdsep2%flo_satexsw = hdsep2%flo_satexsw + ch_stor_hdsep(ich)%flo_satexsw
        hdsep2%flo_tile = hdsep2%flo_tile + ch_stor_hdsep(ich)%flo_tile
        ch_stor_hdsep(ich)%flo_surq = frac * hdsep1%flo_surq
        ch_stor_hdsep(ich)%flo_latq = frac * hdsep1%flo_latq
        ch_stor_hdsep(ich)%flo_gwsw = frac * hdsep1%flo_gwsw
        ch_stor_hdsep(ich)%flo_swgw = frac * hdsep1%flo_swgw
        ch_stor_hdsep(ich)%flo_satex = frac * hdsep1%flo_satex
        ch_stor_hdsep(ich)%flo_satexsw = frac * hdsep1%flo_satexsw
        ch_stor_hdsep(ich)%flo_tile = frac * hdsep1%flo_tile
      end if
      ob(icmd)%hdsep%flo_surq = hdsep2%flo_surq
      ob(icmd)%hdsep%flo_latq = hdsep2%flo_latq
      ob(icmd)%hdsep%flo_gwsw = hdsep2%flo_gwsw
      ob(icmd)%hdsep%flo_swgw = hdsep2%flo_swgw
      ob(icmd)%hdsep%flo_satex = hdsep2%flo_satex
      ob(icmd)%hdsep%flo_satexsw = hdsep2%flo_satexsw
      ob(icmd)%hdsep%flo_tile = hdsep2%flo_tile
      !store outflow components for writing (and convert from m3 --> m3/sec)
      hyd_sep_array(ich,1) = hdsep2%flo_surq / 86400.
      hyd_sep_array(ich,2) = hdsep2%flo_latq / 86400.
      hyd_sep_array(ich,3) = hdsep2%flo_gwsw / 86400.
      hyd_sep_array(ich,4) = hdsep2%flo_swgw / 86400.
      hyd_sep_array(ich,5) = hdsep2%flo_satex / 86400.
      hyd_sep_array(ich,6) = hdsep2%flo_satexsw / 86400.
      hyd_sep_array(ich,7) = 0. !hdsep2%flo_tile / 86400.
      !rtb hydrograph separation
      !end if

      ich = isdch
            
      !! check decision table for flow control - water allocation
      if (ob(icmd)%ruleset /= "null" .and. ob(icmd)%ruleset /= "0") then
        id = ob(icmd)%flo_dtbl
        d_tbl => dtbl_flo(id)
        call conditions (ich, id)
        call actions (ich, icmd, id)
      end if
 
      !! check decision table for flow control - water allocation
      if (sd_ch(isdch)%wallo > 0) then
        call wallo_control (sd_ch(isdch)%wallo)
      end if
      
      !! set outflow hyd to ht2 after diverting water
      ob(icmd)%hd(1) = ht2

      !ht2 = ob(icmd)%hd(1)  !! reset ht2 for printing
      ob(icmd)%hd(1)%temp = 5. + .75 * wst(iwst)%weat%tave
      ht2%temp = 5. + .75 * wst(iwst)%weat%tave
      ch_stor(isdch)%temp = 5. + .75 * wst(iwst)%weat%tave
      
      if (cs_db%num_pests > 0) then
        obcs(icmd)%hd(1)%pest = hcs2%pest
      end if
      
      !! output channel organic-mineral
      ch_out_d(isdch) = ob(icmd)%hd(1)                       !set outflow om hydrograph
      ch_out_d(isdch)%flo = ob(icmd)%hd(1)%flo / 86400.      !m3 -> m3/s
      
      !! output channel morphology
      chsd_d(isdch)%flo = ob(icmd)%hd(1)%flo / 86400.        !adjust if overbank flooding is moved to landscape
      chsd_d(isdch)%flo_mm = ob(icmd)%hd(1)%flo / (10. * ob(icmd)%area_ha)   !flow out in mm
      chsd_d(isdch)%peakr = peakrate 
      chsd_d(isdch)%sed_in = ob(icmd)%hin%sed
      chsd_d(isdch)%sed_out = sedout
      chsd_d(isdch)%sed_stor = ch_stor(isdch)%sed
      if (sedout > 2000.) then      !***jga
        dep = bedld
      end if
      chsd_d(isdch)%washld = washld
      chsd_d(isdch)%bedld = bedld
      chsd_d(isdch)%dep = dep
      chsd_d(isdch)%deg_btm = ebtm_t
      chsd_d(isdch)%deg_bank = ebank_t
      chsd_d(isdch)%hc_sed = hc_sed
      chsd_d(isdch)%width = sd_ch(isdch)%chw
      chsd_d(isdch)%depth = sd_ch(isdch)%chd
      chsd_d(isdch)%slope = sd_ch(isdch)%chs
      chsd_d(isdch)%deg_btm_m = ebtm_m
      chsd_d(isdch)%deg_bank_m = ebank_m
      chsd_d(isdch)%hc_m = hc
      
      !! set pesticide output variables
      do ipest = 1, cs_db%num_pests
        chpst_d(isdch)%pest(ipest)%tot_in = obcs(icmd)%hin%pest(ipest)
        chpst_d(isdch)%pest(ipest)%sol_out = frsol * obcs(icmd)%hd(1)%pest(ipest)
        chpst_d(isdch)%pest(ipest)%sor_out = frsrb * obcs(icmd)%hd(1)%pest(ipest)
        chpst_d(isdch)%pest(ipest)%react = chpst%pest(ipest)%react
        chpst_d(isdch)%pest(ipest)%volat = chpst%pest(ipest)%volat
        chpst_d(isdch)%pest(ipest)%settle = chpst%pest(ipest)%settle
        chpst_d(isdch)%pest(ipest)%resus = chpst%pest(ipest)%resus
        chpst_d(isdch)%pest(ipest)%difus = chpst%pest(ipest)%difus
        chpst_d(isdch)%pest(ipest)%react_bot = chpst%pest(ipest)%react_bot
        chpst_d(isdch)%pest(ipest)%bury = chpst%pest(ipest)%bury 
        chpst_d(isdch)%pest(ipest)%water = ch_water(ich)%pest(ipest)
        chpst_d(isdch)%pest(ipest)%benthic = ch_benthic(ich)%pest(ipest)
      end do
        
      !! set values for recharge hydrograph - should be trans losses
      !ob(icmd)%hd(2)%flo = perc  

      return
      
      end subroutine sd_channel_control