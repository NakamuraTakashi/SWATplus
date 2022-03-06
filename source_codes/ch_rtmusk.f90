      subroutine ch_rtmusk
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine routes a daily flow through a reach using the
!!    Muskingum method

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch_d(:)     |m             |average depth of main channel
!!    ch_n(2,:)   |none          |Manning"s "n" value for the main channel
!!    ch_s(2,:)   |m/m           |average slope of main channel
!!    chside(:)   |none          |change in horizontal distance per unit
!!                               |change in vertical distance on channel side
!!                               |slopes; always set to 2 (slope=1/2)
!!    flwin(:)    |m^3 H2O       |flow into reach on previous day
!!    flwout(:)   |m^3 H2O       |flow out of reach on previous day
!!    i           |none          |current day of simulation
!!    pet_day     |mm H2O        |potential evapotranspiration
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    flwin(:)    |m^3 H2O       |flow into reach on current day
!!    flwout(:)   |m^3 H2O       |flow out of reach on current day
!!    rcharea     |m^2           |cross-sectional area of flow
!!    rchdep      |m             |depth of flow on day
!!    rtevp       |m^3 H2O       |evaporation from reach on day
!!    rttime      |hr            |reach travel time
!!    rttlc       |m^3 H2O       |transmission losses from reach on day
!!    rtwtr       |m^3 H2O       |water leaving reach on day
!!    sdti        |m^3/s         |average flow on day in reach
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    c           |none          |inverse of channel side slope
!!    c1          |
!!    c2          |
!!    c3          |
!!    c4          |m^3 H2O       |
!!    det         |hr            |time step (24 hours)
!!    jrch        |none          |reach number
!!    nn          |              |number of subdaily computation points for stable 
!!                               |routing in the muskingum routing method
!!    p           |m             |wetted perimeter
!!    rh          |m             |hydraulic radius
!!    tbase       |none          |flow duration (fraction of 24 hr)
!!    topw        |m             |top width of main channel
!!    vol         |m^3 H2O       |volume of water in reach at beginning of
!!                               |day
!!    wtrin       |m^3 H2O       |water entering reach on day
!!    xkm         |hr            |storage time constant for the reach on
!!                               |current day
!!    yy          |none          |variable to hold intermediate calculation
!!                               |value
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sqrt
!!    SWAT: Qman

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!!    code provided by Dr. Valentina Krysanova, Pottsdam Institute for
!!    Climate Impact Research, Germany
!!    Modified by Balaji Narasimhan
!!    Spatial Sciences Laboratory, Texas A&M University

      use basin_module
      use channel_data_module
      use channel_module
      use hydrograph_module !, only : ob, icmd, jrch, isdch, fp_stor, ch_stor
      use time_module
      use channel_velocity_module
      use sd_channel_module
      use climate_module
      use reservoir_module
      
      implicit none
      
      integer :: nn        !none              |number of subdaily computation points for stable 
                           !                  |routing in the muskingum routing method
      integer :: ii        !none              |counter
      integer :: i         !none              |current day of simulation
      integer :: ihru
      integer :: icha
      integer :: irtstep
      integer :: isubstep
      
      real :: c1           !units             |description 
      real :: c2           !units             |description
      real :: c3           !units             |description
      real :: c4           !m^3 H2O           |
      real :: p            !m                 |wetted perimeter
      real :: vol          !m^3 H2O           |volume of water in reach at beginning of
                           !                  |day
      real :: c            !none              |inverse of channel side slope
      real :: rh           !m                 |hydraulic radius
      real :: topw         !m                 |top width of main channel
      
      real :: qinday       !units             |description 
      real :: qoutday      !units             |description  
	  real :: volrt        !units             |description 
      real :: maxrt        !units             |description 
      real :: adddep       !units             |description 
      real :: addp         !units             |description 
      real :: addarea      !units             |description 
	  real :: rttlc1       !units             |description 
      real :: rttlc2       !units             |description 
      real :: rtevp1       !units             |description 
      real :: rtevp2       !units             |description 
      real :: qman         !m^3/s or m/s      |flow rate or flow velocity
      real :: vc           !m/s               |flow velocity in reach
      real :: aaa          !units             |description 
      real :: inflo        !m^3           |inflow water volume
      real :: inflo_rate   !m^3/s         |inflow rate
      real :: xs_area      !m^2           |cross section area of channel
      real :: dep_flo      !m             |depth of flow
      real :: wet_perim    !m             |wetted perimeter
      real :: ttime        !hr            |travel time through the reach
      real :: t_inc        !hr            |time in routing step - 1/time%step
      real :: outflo       !m^3           |outflow water volume
      real :: tl           !m^3           |transmission losses during time step
      real :: trans_loss   !m^3           |transmission losses during day
      real :: rate_flo     !m^3/s         |flow rate
      real :: ev           !m^3           |evaporation during time step
      real :: evap         !m^3           |evaporation losses during day
      real :: above_prin_fr
      real :: out_decr_fr
      real :: ch_out_fr
      real :: fp_out_fr
      real :: ch_fp_fr
      real :: fp_ch_fr
      real :: rto
      real :: rto_w
      real :: outflo_rate
      real :: ch_st             !m^3        |water storage in and above channel
      real :: fp_st             !m^3        |water storage in flood plain
      real :: wet_st            !m^3        |wetland storage above principal
      real :: dts               !seconds    |time step interval for substep
      real :: dthr
      real :: scoef
      real :: vol_ch_av
      real :: vol_fp_av
      real :: vol_tot_av

      jrch = isdch
      jhyd = sd_dat(jrch)%hyd
      
      qinday = 0
      qoutday = 0
      !det = time%dtm / 60.
      
      !! volume at start of day
      vol = ch_stor(jrch)%flo + fp_stor(jrch)%flo
      
      irtstep = 1
      isubstep = 0
      dts = time%dtm / sd_ch(jrch)%msk%substeps * 60.
      dthr = dts / 3600.

      !! subdaily time step
      do ii = 1, sd_ch(jrch)%msk%nsteps

        !! water entering reach during time step - substeps for stability
        isubstep = isubstep + 1
        if (isubstep > sd_ch(jrch)%msk%substeps) then
          irtstep = irtstep + 1
          isubstep = 1
        end if
        inflo = ob(icmd)%tsin(irtstep) / sd_ch(jrch)%msk%substeps
        inflo_rate = inflo / dts 

        !! interpolate rating curve using inflow rates
        icha = jrch
        call rcurv_interp_flo (icha, inflo_rate)
        ch_rcurv(jrch)%in2 = rcurv
        
        !! save variables at each routing time step for sediment routing
        if (isubstep == 0) then
          hyd_rad(irtstep) = rcurv%xsec_area / rcurv%wet_perim
          trav_time(irtstep) = rcurv%ttime
          flo_dep(irtstep) = rcurv%dep
        end if

        !! add inflow volume of water - in channel and flood plain
        vol = vol + inflo
        vol = Max(vol, 1.e-14)
        
        !! if no water in channel - skip routing and set rating curves to zero
        if (vol < 1.e-6) then
          ch_rcurv(jrch)%in1 = rcz
          ch_rcurv(jrch)%out1 = rcz
          sd_ch(jrch)%in1_vol = 0.
          sd_ch(jrch)%out1_vol = 0.
        else
        
        !! sum flood plain wetland storage to get total flood plain storage
        fp_stor(jrch) = hz
        do ihru = 1, sd_ch(jrch)%fp%hru_tot
          if (wet(ihru)%flo > wet_ob(ihru)%pvol .and. wet(ihru)%flo > 1.e-6) then
            above_prin_fr = (wet(ihru)%flo - wet_ob(ihru)%pvol) / wet(ihru)%flo
            fp_stor(jrch) = fp_stor(jrch) + above_prin_fr * wet(ihru)
          end if
        end do
        
        !! use flow ratio (time step/total day) to get hyd for time step (assume uniform conc during day)
        if (ob(icmd)%hin%flo > 1.e-6) then
          rto = inflo / ob(icmd)%hin%flo
          rto = Min (rto, 1.)
        else
          rto = 0.
        end if
        ht1 = rto * ob(icmd)%hin
        
        !! partition inflow between channel and flood plain storage
        if (inflo > 1.e-6) then
          rto = rcurv%vol_fp / rcurv%vol
          fp_stor(jrch) = fp_stor(jrch) + rto * ht1
          do ihru = 1, sd_ch(jrch)%fp%hru_tot
            rto_w = rto * sd_ch(jrch)%fp%hru_fr(ihru)
            wet(ihru) = wet(ihru) + rto_w * ht1
          end do
          rto = rcurv%vol_ch / rcurv%vol
          ch_stor(jrch) = ch_stor(jrch) + rto * ht1
        end if

        !! Muskingum method
        outflo = sd_ch(jrch)%msk%c1 * inflo + sd_ch(jrch)%msk%c2 * sd_ch(jrch)%in1_vol +     &
                                                sd_ch(jrch)%msk%c3 * sd_ch(jrch)%out1_vol
               
        !! save inflow/outflow volumes for next day for Muskingum
        sd_ch(jrch)%in1_vol = inflo
        sd_ch(jrch)%out1_vol = outflo

        !! VSC method - sc=2*dt/(2*ttime+dt) - ttime=(in2+out1)/2
        scoef = 2. * dthr / (ch_rcurv(jrch)%in2%ttime + ch_rcurv(jrch)%out1%ttime + dthr)
        scoef = Min (scoef, 1.)
        !! outflo = scoef * (inflo + storage) --> inflow already added to storage
        !outflo = scoef * (ch_stor(jrch)%flo + fp_stor(jrch)%flo)
        
	    outflo = Min (outflo, vol)
        outflo = Max (outflo, 0.)
 
        !! compute outflow rating curve for next time step
        outflo_rate = outflo / dts      !convert to cms
        call rcurv_interp_flo (jrch, outflo_rate)
        ch_rcurv(jrch)%out2 = rcurv
        
        !! subtract outflow volume of water from channel and then adjust flood plain
        out_decr_fr = 0.
        ch_st = ch_stor(jrch)%flo - outflo
        if (ch_st < 0.) then
          fp_st = fp_stor(jrch)%flo - abs(ch_st)
          ch_out_fr = 1.
          if (fp_st < 0.) then
            fp_out_fr = 1.
            out_decr_fr = abs(fp_st) / outflo
          else
            fp_out_fr = 1. - (fp_st / fp_stor(jrch)%flo)
          end if
        else
          ch_out_fr = outflo / ch_stor(jrch)%flo
          fp_out_fr = 0.
        end if
          ht1 = ch_out_fr * ch_stor(jrch)
          ht2 = fp_out_fr * fp_stor(jrch)
          ch_stor(jrch) = ch_stor(jrch) - ht1
          fp_stor(jrch) = fp_stor(jrch) - ht2
          do ihru = 1, sd_ch(jrch)%fp%hru_tot
            wet(ihru) = wet(ihru) - sd_ch(jrch)%fp%hru_fr(ihru) * ht2
          end do
        
        !! move water from channel or flood plain to get ratio in the rating curve
        !! use average volumes from in/out rating curves and channel storage at start/end
        vol_fp_av = (rcurv%vol_fp + ch_rcurv(jrch)%out1%vol_fp) / 2.
        vol_ch_av = (rcurv%vol_ch + ch_rcurv(jrch)%out1%vol_ch) / 2.
        vol_tot_av = (rcurv%vol + ch_rcurv(jrch)%out1%vol) / 2.
        ch_st = (sd_ch(jrch)%stor + ch_stor(jrch)%flo) / 2. !! could try rating curve from ave total storage
        
        if (vol_tot_av > 1.e-6) then
          rto = vol_ch_av / vol_tot_av
          rto = Min (rto, 1.)
        else
          rto = 1.
        end if
        ch_st = rto * (ch_stor(jrch)%flo + fp_stor(jrch)%flo)
        ch_fp_fr = 0.
        fp_ch_fr = 0.
        if (ch_stor(jrch)%flo - ch_st > 1.e-6) then
          !! move fraction of water from channel to flood plain
          ch_fp_fr = (ch_stor(jrch)%flo - ch_st) / ch_stor(jrch)%flo
          fp_ch_fr = 0.
        end if
        if (ch_st - ch_stor(jrch)%flo > 1.e-6) then
          !! move fraction of water from flood plain to channel
          fp_ch_fr = (ch_st - ch_stor(jrch)%flo) / fp_stor(jrch)%flo
          ch_fp_fr = 0.
        end if
        
          if (fp_ch_fr > 1.e-6) then
            !! move water from flood plain to channel
            do ihru = 1, sd_ch(jrch)%fp%hru_tot
              ht1 = fp_ch_fr * fp_stor(jrch)
              wet_st = sd_ch(jrch)%fp%hru_fr(ihru) * ht1%flo
              if (wet(ihru)%flo - wet_ob(ihru)%pvol > wet_st) then
                !! can take all from above principal storage
                ht2 = sd_ch(jrch)%fp%hru_fr(ihru) * ht1
                wet(ihru) = wet(ihru) - ht2
                fp_stor(jrch) = fp_stor(jrch) - ht2
                ch_stor(jrch) = ch_stor(jrch) + ht2
              else
                !! just take what is above principal storage - none if below
                if (wet(ihru)%flo > wet_ob(ihru)%pvol .and. wet(ihru)%flo > 1.e-6) then
                  above_prin_fr = (wet(ihru)%flo - wet_ob(ihru)%pvol) / wet(ihru)%flo
                  ht1 = above_prin_fr * wet(ihru)
                  fp_stor(jrch) = fp_stor(jrch) - ht1
                  ch_stor(jrch) = ch_stor(jrch) + ht1
                  wet(ihru) = (1. - above_prin_fr) * wet(ihru)
                end if
              end if
            end do
          end if
          if (ch_fp_fr > 1.e-6) then
            !! move water from channel to flood plain
            ht2 = ch_stor(jrch)
            do ihru = 1, sd_ch(jrch)%fp%hru_tot
              rto = ch_fp_fr * sd_ch(jrch)%fp%hru_fr(ihru)
              ht1 =  rto * ht2
              wet(ihru) = wet(ihru) + ht1
              fp_stor(jrch) = fp_stor(jrch) + ht1
              ch_stor(jrch) = ch_stor(jrch) - ht1
            end do
          end if
        
        !! set rating curve for next time step
        ch_rcurv(jrch)%in1 = ch_rcurv(jrch)%in2
        ch_rcurv(jrch)%out1 = ch_rcurv(jrch)%out2
            
        !! add channel (ht1) and flood plain (ht2) flow to daily flow
        ob(icmd)%hd(1) = ob(icmd)%hd(1) + ht1 + ht2
          
          !! calculate transmission losses
          ttime = Min(24., rcurv%ttime)
          tl = sd_ch(jrch)%chk * sd_ch(jrch)%chl * wet_perim * ttime   !mm/hr * km * mm * hr = m3       
          tl = Min(tl, outflo)
          outflo = outflo - tl
          trans_loss = trans_loss + tl

          !! calculate evaporation
          if (outflo > 0.) then
            !! calculate width of channel at water level - flood plain evap calculated in wetlands
            if (dep_flo <= sd_ch(jrch)%chd) then
              topw = ch_rcurv(jrch)%out2%surf_area
            else
              topw = 1000. * sd_ch(jrch)%chl * sd_ch(jrch)%chw
            end if
            
            iwst = ob(icmd)%wst
            !! mm * m2 / (1000. * sd_ch(jrch)%msk%nsteps)
            ev = bsn_prm%evrch * wst(iwst)%weat%pet * topw / (1000. * sd_ch(jrch)%msk%nsteps)
            if (ev < 0.) ev = 0.
            ev = Min(ev, outflo)
            outflo = outflo - ev
            evap = evap + ev
          end if

          !! set volume of water in channel at end of hour
          vol = vol - outflo
          !! volume at start of day
          sd_ch(jrch)%stor = ch_stor(jrch)%flo + fp_stor(jrch)%flo
          ob(icmd)%hyd_flo(1,irtstep) = ob(icmd)%hyd_flo(1,irtstep) + outflo
          
        end if          !! vol < 1.e-6 loop

      end do            !! end of sub-daily loop

      return
      end subroutine ch_rtmusk