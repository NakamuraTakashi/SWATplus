      subroutine wetland_control
    
      use reservoir_data_module
      use reservoir_module
      use hru_module, only : hru, sedyld, sanyld, silyld, clayld, sagyld, lagyld, grayld, sedminps, sedminpa,   &
        surqno3, sedorgn, sedorgp, qdr, ihru, pet_day, qday, surfq, tconc, usle_cfac, cklsp
      use conditional_module
      use climate_module
      use hydrograph_module
      use time_module
      use basin_module
      use channel_module
      use water_body_module
      use soil_module
      use organic_mineral_mass_module
      use mgt_operations_module
      use constituent_mass_module
      use aquifer_module
      
      implicit none
     
      real :: bypass = 1.             !              | 
      real :: fracwet                 !              | 
      integer :: j                    !none          |counter
      integer :: iprop                !              |  
      integer :: iac                  !none          |counter
      character(len=1) :: action           !         |
      integer :: ial                  !none          |counter
      real :: b_lo                    !              |
      real :: res_h                   !              |
      real :: x1                      !              |
      real :: wet_h                   !              |
      real :: wet_h1                  !              |
      real :: flwi                    !m^3 H2O       |water entering pothole on day  
      real :: flwo                    !              |
      real :: sedi                    !metric tons   |sediment entering pothole on day
      real :: sedo                    !metric tons   |sed leaving res 
      integer :: k                    !              | 
      integer :: ii                   !none          |counter 
      integer :: jres                 !none          |reservoir number
      integer :: idat                 !              |
      integer :: ihyd                 !none          |counter
      integer :: ised                 !none          |counter
      integer :: irel                 !              |
      integer :: inut                 !none          |counter
      integer :: ipst                 !none          |counter
      integer :: ires = 0
      integer :: iob,irrop,iout,j1
      real :: wet_fr = 0.
      real :: pvol_m3
      real :: evol_m3
      real :: dep
      real :: weir_hgt,wsa1,sedppm,no3ppm,seep_rto,flo_init,qp_cms,qmm,dep_init,satk,rto
	real :: volseep, volex, swst(20)	  

      j = ihru
      ires= hru(j)%dbs%surf_stor
      ihyd = wet_dat(ires)%hyd
      ised = wet_dat(ires)%sed
      irel = wet_dat(ires)%release
      wsa1 = hru(j)%area_ha * 10. 

      !! initialize variables for wetland daily simulation
      hru(j)%water_seep = 0.
      wet_ob(j)%depth = wet(j)%flo / wsa1 / 1000. !m
      dep_init = wet_ob(j)%depth !m
      
      !! add precipitation - mm*ha*10.=m3 (used same area for infiltration and soil evap)
      wet_wat_d(j)%precip = w%precip * wsa1 !m3
      wet_ob(j)%depth = wet_ob(j)%depth + w%precip / 1000.  !m  Jaehak 2022
      wet(j)%flo =  wet(j)%flo + wet_wat_d(j)%precip !m3
      
      !add irrigated water to wetland storage
      if (wet_ob(j)%depth*1000. < hru(j)%irr_hmin) then
          
        irrig(j)%demand = max(0., hru(j)%irr_hmax - wet_ob(j)%depth*1000.) * hru(j)%area_ha * 10. !m3
          
        rto = 0.
        if (hru(j)%irr_src == 'cha') then
          if (ubound(ch_stor,1) > 0) then 
            if (ch_stor(iob)%flo > 0.001) then
              rto = min(0.99, irrig(j)%demand / ch_stor(iob)%flo)   ! ratio of water removed from channel volume
            end if
            irrig(j)%water = rto * ch_stor(iob)                     ! irrigation water
            cs_irr(iob) = rto * ch_water(iob)                       ! constituents in irrigation water
            ch_stor(iob) = (1. - rto) * ch_stor(iob)                ! remainder stays in channel
            ch_water(iob) = (1. - rto) * ch_water(iob)                      
          else
            irrig(j)%water%flo = irrig(j)%demand
          endif
            
        elseif (hru(j)%irr_src == 'res') then
          if (ubound(res,1) > 0) then 
            if (res(iob)%flo > 0.001) then
              rto = min(0.99, irrig(j)%demand / res(iob)%flo)       ! ratio of water removed from res volume
            end if
            irrig(j)%water = rto * res(iob)                         ! organics in irrigation water
            cs_irr(iob) = rto * res_water(iob)                      ! constituents in irrigation water
            res(iob) = (1. - rto) * res(iob)                        ! remainder stays in reservoir
            res_water(iob) = (1. - rto) * res_water(iob)                    
          else
            irrig(j)%water%flo = irrig(j)%demand
          endif

        elseif (hru(j)%irr_src == 'aqu') then
          if (ubound(aqu_d,1) > 0) then 
            if (aqu_d(iob)%stor > 0.001) then
              rto = min(0.99, irrig(j)%demand / aqu_d(iob)%stor)    ! ratio of water removed from aquifer volume
            end if
            irrig(j)%water%flo = rto * aqu_d(iob)%flo               ! organics in irrigation water
            cs_irr(iob) = rto * cs_aqu(iob)                         ! constituents in irrigation water
            aqu_d(iob)%stor = (1. - rto) * aqu_d(iob)%stor          ! remainder stays in aquifer
            cs_aqu(iob) = (1. - rto) * cs_aqu(iob)  
          else
            irrig(j)%water%flo = irrig(j)%demand
          endif

        else !unlimited source
          irrig(j)%water%flo = irrig(j)%demand

        end if
        
        irrig(j)%applied = irrig(j)%water%flo / wsa1 * irrig(j)%eff * (1. - irrig(j)%frac_surq) !mm
        irrig(j)%runoff = irrig(j)%water%flo / wsa1 * irrig(j)%eff * irrig(j)%frac_surq !mm
      endif
      wet(j)%flo =  wet(j)%flo + irrig(j)%applied * wsa1 !m3
      
      !! subtract evaporation  - mm*ha*10.=m3
      wet_wat_d(j)%evap = pet_day * wet_hyd(ihyd)%evrsv * wsa1 !m3
      wet_wat_d(j)%evap = min(wet_wat_d(j)%evap, wet(j)%flo)
      wet(j)%flo = wet(j)%flo - wet_wat_d(j)%evap 
      hru(j)%water_evap = wet_wat_d(j)%evap / wsa1 !mm
      !wet(j)%dep = wet(j)%dep - hru(j)%water_evap !mm
      
      wet_wat_d(j)%area_ha = 0.
      if (wet(j)%flo > 0.) then  !paddy is assumed flat
        !! update wetland surface area - solve quadratic to find new depth
        x1 = wet_hyd(ihyd)%bcoef ** 2 + 4. * wet_hyd(ihyd)%ccoef * (1. - wet(j)%flo / (wet_ob(j)%pvol + 1.e-9))
        if (x1 < 1.e-6) then
          wet_h = 0.
        else
          wet_h1 = (-wet_hyd(ihyd)%bcoef - sqrt(x1)) / (2. * wet_hyd(ihyd)%ccoef + 1.e-9)
          wet_h = wet_h1 + wet_hyd(ihyd)%bcoef
        end if
        wet_fr = (1. + wet_hyd(ihyd)%acoef * wet_h)
        wet_fr = min(wet_fr,1.)
        wet_fr = max(wet_fr,0.01)
        
        wet_wat_d(j)%area_ha = hru(j)%area_ha * wet_fr
        
        !! infiltration of the standing water to the topsoil layer. 
        !! Any excess infiltration volume estimated here is reverted (back to waterbody) in swr_satexcess.
        !! 
        wet_wat_d(j)%seep = min(wet(j)%flo, hru(j)%wet_hc * 24. * wsa1) !m3 
		
		! check potential percolation rate to refine daily seepage rate Jaehak 2022
		! actual soil moisture content is updated in percmain
		volseep = wet_wat_d(j)%seep / wsa1 !mm
		if (volseep>0.1) then
			do j1 = 1, soil(j)%nly
				swst(j1) = soil(j)%phys(j1)%st + volseep
				    if (swst(j1)>soil(j)%phys(j1)%ul*0.9) then !oversaturated Jaehak 2022
					    volex = swst(j1) - soil(j)%phys(j1)%ul*0.9  !excess water. soil is assumed to remain saturated Jaehak 2022
					volseep = min(volex, soil(j)%phys(j1)%k*24.)
					swst(j1) = swst(j1) - volseep
				else
					volseep = 0
				endif
			end do
			
			! move excess water upward to calculate daily seepage rate Jaehak 2022
			volex = 0
			do j1 = soil(j)%nly, 1, -1
				swst(j1) = swst(j1) + volex
				    if (swst(j1)>soil(j)%phys(j1)%ul*0.9) then !oversaturated
					    volex = max(0., swst(j1) - soil(j)%phys(j1)%ul*0.9)  !excess water. 
					swst(j1) = swst(j1) - volex                         !update soil water
				endif
			end do
			
			!update seepage volume
			wet_wat_d(j)%seep = max(0., wet_wat_d(j)%seep - volex * wsa1) !m3
		endif				
				
        wet(j)%flo = wet(j)%flo - wet_wat_d(j)%seep
        wet_wat_d(j)%area_ha = hru(j)%area_ha 
        hru(j)%water_seep = wet_wat_d(j)%seep / wsa1   !mm=m3/(10*ha)
      
        ! calculate dissolved nutrient infiltration Jaehak 2022
        seep_rto = wet_wat_d(j)%seep / (wet_wat_d(j)%seep + wet(j)%flo)
        soil1(j)%mn(1)%no3 = soil1(j)%mn(1)%no3 + wet(j)%no3 * seep_rto / hru(j)%area_ha !kg/ha
        soil1(j)%mn(1)%nh4 = soil1(j)%mn(1)%nh4 + wet(j)%nh3 * seep_rto / hru(j)%area_ha !kg/ha
        soil1(j)%mp(1)%act = soil1(j)%mp(1)%act + wet(j)%solp * seep_rto / hru(j)%area_ha !kg/ha
        soil1(j)%water(1)%n = soil1(j)%water(1)%n + wet(j)%orgn * seep_rto / hru(j)%area_ha !kg/ha
        soil1(j)%water(1)%p = soil1(j)%water(1)%p + wet(j)%sedp * seep_rto / hru(j)%area_ha !kg/ha
        
        ! nutrient seepage amount 
        wet_seep_day(j)%no3 = wet(j)%no3 * seep_rto !kg
        wet_seep_day(j)%nh3 = wet(j)%nh3 * seep_rto
        wet_seep_day(j)%orgn = wet(j)%orgn * seep_rto
        wet_seep_day(j)%solp = wet(j)%solp * seep_rto
        wet_seep_day(j)%sedp = wet(j)%sedp * seep_rto
        
        ! substract the seepage amount from the ponding water
        wet(j)%no3 = wet(j)%no3 - wet_seep_day(j)%no3 
        wet(j)%nh3 = wet(j)%nh3 - wet_seep_day(j)%nh3 
        wet(j)%orgn = wet(j)%orgn - wet_seep_day(j)%orgn
        wet(j)%solp = wet(j)%solp - wet_seep_day(j)%solp
        wet(j)%sedp = wet(j)%sedp - wet_seep_day(j)%sedp
      end if 
        
      !! if not a floodplain wetland
      if (hru(j)%wet_fp == "n") then
        !! calc release from decision table
        d_tbl => dtbl_res(irel)
        wbody => wet(j)
        wbody_wb => wet_wat_d(j)
        pvol_m3 = wet_ob(j)%pvol
        evol_m3 = wet_ob(j)%evol
        !if (wet_ob(j)%area_ha > 1.e-6) then
        if (hru(j)%area_ha > 1.e-6) then
          !dep = wbody%flo / wet_ob(j)%area_ha / 10000.     !m = m3 / ha / 10000m2/ha
          dep = wet(j)%flo / wsa1 / 1000.    !m 
        else
          dep = 0.
        end if
        weir_hgt = wet_ob(j)%weir_hgt !m
        wet_ob(j)%depth = dep !m


        !! use cn control for Jaehak's release algorithm - points to weir.res
        !if (bsn_cc%cn > 0) then 
        !  call res_weir_release (j, irel, ihyd, pvol_m3, evol_m3, dep, weir_hgt)
        !else
          call conditions (j, irel)
          call res_hydro (j, irel, ihyd, pvol_m3, evol_m3, dep, weir_hgt)
        !end if
          
        !! subtract outflow from storage
        wet(j)%flo =  wet(j)%flo - ht2%flo
        surfq(j) = ht2%flo / wsa1 !mm
      end if
      
      wet_ob(j)%depth = wet(j)%flo / wsa1 / 1000. !m                       
       
      !! compute sediment deposition
      call res_sediment (j, ihyd, ised)
      
      !!! subtract sediment leaving from reservoir
      !wet(j)%sed = wet(j)%sed - ht2%sed
      !wet(j)%sil = wet(j)%sil - ht2%sil
      !wet(j)%cla = wet(j)%cla - ht2%cla
          
      !! perform reservoir nutrient balance
      inut = wet_dat(ires)%nut
      call res_nutrient (ires, inut, j)

      ! calculate sediment/nutrient yield when wetlands are flushed 
      if (dep_init<0.0001 .and. ht2%flo>0.) then 
        call ero_cfactor 
        qp_cms = bsn_prm%prf / 6578.6 * hru(j)%area_ha * surfq(j) / tconc(j) / 35.3
        cklsp(j) = usle_cfac(j) * hru(j)%lumv%usle_mult
        wet(j)%sed = (10. * surfq(j) * qp_cms * hru(j)%area_ha) ** .56 * cklsp(j) !tons
        sedyld(j) = wet(j)%sed * ht2%flo / (wet(j)%flo + ht2%flo)
        wet(j)%sed = max(0., wet(j)%sed - sedyld(j)) !tons
        ht2%sed = sedyld(j) !tons sediment yield
      endif  
        
      if (wet(j)%flo>0) then
        sedppm=wet(j)%sed/wet(j)%flo*1000000.
        no3ppm=wet(j)%no3/wet(j)%flo*1000.
      else
        sedppm=0
        no3ppm=0
      endif
      

  
      !write(100100,'(3(I6,","),11(f10.1,","))') time%yrc,time%mo,time%day_mo,w%precip,irrig(j)%applied,hru(j)%water_seep,&
      !  weir_hgt*1000,wet(j)%flo/wsa1,ht2%flo/wsa1,soil(j)%sw,wet(j)%sed*1000,ht2%sed*1000,no3ppm,ht2%no3    
      !write(*,'(3(I6),11(f10.1))') time%yrc,time%mo,time%day_mo,w%precip,irrig(j)%applied,hru(j)%water_seep,&
      !  weir_hgt*1000,wet(j)%flo/wsa1,ht2%flo/wsa1,soil(j)%sw,wet(j)%sed*1000,ht2%sed*1000,wet(j)%no3,ht2%no3    
      
      !! perform reservoir pesticide transformations
      !call res_pest (ires)

      !! set values for routing variables
      ob(icmd)%hd(1)%temp = 0.                  !!undefined

      !qdr(j) = ht2%flo / (10. * hru(j)%area_ha) + ht1%flo * bypass
      sedyld(j) = ht2%sed / hru(j)%area_ha + sedyld(j) * bypass
      sanyld(j) = ht2%san / hru(j)%area_ha + sanyld(j) * bypass
      silyld(j) = ht2%sil / hru(j)%area_ha + silyld(j) * bypass
	    clayld(j) = ht2%cla / hru(j)%area_ha + clayld(j) * bypass 
	    sagyld(j) = ht2%sag / hru(j)%area_ha + sagyld(j) * bypass
	    lagyld(j) = ht2%lag / hru(j)%area_ha + lagyld(j) * bypass
	    grayld(j) = ht2%grv / hru(j)%area_ha + grayld(j) * bypass


      sedorgn(j) = ht2%orgn / hru(j)%area_ha + sedorgn(j) * bypass
      sedorgp(j) = ht2%sedp / hru(j)%area_ha + sedorgp(j) * bypass
      surqno3(j) = ht2%no3/ hru(j)%area_ha  + surqno3(j) * bypass
      !nh3 = resnh3o + 0.  !add ammonium 
      !no2  = resno2o + 0.  !add no2
      sedminps(j) = ht2%solp / hru(j)%area_ha / 2. + sedminps(j) * bypass
      sedminpa(j) = ht2%solp / hru(j)%area_ha / 2. + sedminpa(j) * bypass
      
      !! set inflow and outflow variables for reservoir_output
      if (time%yrs > pco%nyskip) then
        wet_in_d(j) = ht1 
        wet_out_d(j) = ht2
        !wet_in_d(j)%flo = wet(j)%flo / 10000.   !m^3 -> ha-m
        !wet_out_d(j)%flo = wet(j)%flo / 10000.  !m^3 -> ha-m
      end if  

      return
      end subroutine wetland_control