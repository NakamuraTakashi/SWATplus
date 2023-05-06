      subroutine swr_satexcess(j1)
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine moves water to upper layers if saturated and can't perc

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sep         |mm H2O        |micropore percolation from soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max
!!    SWAT: percmacro, percmicro

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use septic_data_module
      use hru_module, only : hru, ihru, cbodu, surfq, surqno3, surqsolp, sep_tsincefail, i_sep,  &
        isep, qday, sepday, satexq  !rtb gwflow
      use soil_module
      use hydrograph_module
      use basin_module
      use organic_mineral_mass_module
      use gwflow_module, only : gw_transfer_flag !rtb gwflow
      use reservoir_module
      
      implicit none

      integer :: j                 !none          |HRU number
      integer :: j1                !none          |counter
      integer :: ii                !none          |counter
      integer :: isp               !              | 
      real:: ul_excess             !              |
      real:: qlyr                  !              |
      real:: pormm                 !mm            |porosity in mm depth 
      real:: rtof                  !none          |weighting factor used to partition the 
                                   !              |organic N & P concentration of septic effluent
                                   !              |between the fresh organic and the stable organic pools
      real :: xx                   !              |
      real :: rto                   !              |
      integer :: jj                !              |
      integer :: l                 !              | 
      integer :: nn                !none          |number of soil layers
      integer :: ly                !none          |counter
      integer :: ires                !none          |counter

      j = ihru
      ires =  hru(j)%dbs%surf_stor  !!Jaehak 2022
      isp = sep(isep)%typ 	      !! J.Jeong 3/09/09
      rtof = 0.5

 	if (sep(isep)%opt == 2 .and. j1 == i_sep(j)) then
	  
	  ii = j1 
	  qlyr = soil(j)%phys(ii)%st
	  
	  ! distribute excess STE to upper soil layers 
	  do while (qlyr > 0 .and. ii > 1)
		 ! distribute STE to soil layers above biozone layer
        if (soil(j)%phys(ii)%st > soil(j)%phys(ii)%ul) then
	      qlyr = soil(j)%phys(ii)%st - soil(j)%phys(ii)%ul 	! excess water moving to upper layer
	      soil(j)%phys(ii)%st = soil(j)%phys(ii)%ul  ! layer saturated
	      soil(j)%phys(ii-1)%st = soil(j)%phys(ii-1)%st + qlyr ! add excess water to upper layer
	    else 
	      qlyr = 0.
	    endif
	  
	    ! Add surface ponding to the 10mm top layer when the top soil layer is saturated
		 !  and surface ponding occurs.
		 if (ii == 2) then
	     qlyr = soil(j)%phys(1)%st - soil(j)%phys(1)%ul
	     ! excess water makes surface runoff
	     if (qlyr > 0) then
            soil(j)%phys(1)%st = soil(j)%phys(1)%ul
            surfq(j) = surfq(j) + qlyr
	     endif
		 endif

	    ii = ii - 1
	  end do
	endif
	       
      if (j1 < soil(j)%nly) then
        if (soil(j)%phys(j1)%st - soil(j)%phys(j1)%ul > 1.e-4) then
          sepday = (soil(j)%phys(j1)%st - soil(j)%phys(j1)%ul)
          soil(j)%phys(j1)%st = soil(j)%phys(j1)%ul
          soil(j)%phys(j1+1)%st = soil(j)%phys(j1+1)%st + sepday
        end if
      else

        if (soil(j)%phys(j1)%st - soil(j)%phys(j1)%ul > 1.e-4) then
          ul_excess = soil(j)%phys(j1)%st - soil(j)%phys(j1)%ul
          soil(j)%phys(j1)%st = soil(j)%phys(j1)%ul
          nn = soil(j)%nly
          do ly = nn - 1, 1, -1
            soil(j)%phys(ly)%st = soil(j)%phys(ly)%st + ul_excess
            if (soil(j)%phys(ly)%st > soil(j)%phys(ly)%ul) then
              ul_excess = soil(j)%phys(ly)%st - soil(j)%phys(ly)%ul
              soil(j)%phys(ly)%st = soil(j)%phys(ly)%ul
			  soil(j)%ly(ly)%prk = max(0.,soil(j)%ly(ly)%prk - ul_excess) !jaehak 2022
            else
              ul_excess = 0.
              exit
            end if
            if (ly == 1 .and. ul_excess > 0.) then
              !! if no depressional storage (wetland), add to surface runoff
              if (ires == 0) then
                surfq(j) = surfq(j) + ul_excess
              else
                !! move water and nutrient upward and add to wetland storage Jaehak 2022
                !! this is not actual upward movement of water and nutrient, but a process computationally 
                !! rebalancing water and mass balance in the soil profile
                wet(j)%flo = wet(j)%flo + ul_excess * 10. * hru(ihru)%area_ha   !m3=mm*10*ha)
                wet_ob(j)%depth = wet(j)%flo / hru(j)%area_ha / 10000. !m
                if (hru(j)%water_seep < 1.e-6) then
                  rto = 1.
                else
                  rto = ul_excess / hru(j)%water_seep           !ratio of nutrient to be reallocated to ponding water
                end if
                hru(j)%water_seep = rto * hru(j)%water_seep     !updated infiltration volume of the standing water
                !! substract the fraction of nutrient in the top soil layer
                soil1(j)%mn(1)%no3 = soil1(j)%mn(1)%no3 - wet_seep_day(j)%no3 * rto / hru(j)%area_ha !kg/ha
                soil1(j)%mn(1)%nh4 = soil1(j)%mn(1)%nh4 - wet_seep_day(j)%nh3 * rto / hru(j)%area_ha !kg/ha
                soil1(j)%mp(1)%act = soil1(j)%mp(1)%act - wet_seep_day(j)%solp * rto / hru(j)%area_ha !kg/ha
                soil1(j)%water(1)%n = soil1(j)%water(1)%n - wet_seep_day(j)%orgn * rto / hru(j)%area_ha !kg/ha
                soil1(j)%water(1)%p = soil1(j)%water(1)%p - wet_seep_day(j)%sedp * rto / hru(j)%area_ha !kg/ha
                !! add to the wetland water nutrient storage
                wet(j)%no3 = wet(ihru)%no3 + wet_seep_day(j)%no3 * rto  !kg
                wet(j)%nh3 = wet(ihru)%nh3 + wet_seep_day(j)%nh3 * rto  !kg
                wet(j)%orgn = wet(ihru)%orgn + wet_seep_day(j)%orgn * rto !kg
                wet(j)%solp = wet(ihru)%solp + wet_seep_day(j)%solp * rto  !kg
                wet(j)%sedp = wet(ihru)%sedp + wet_seep_day(j)%sedp * rto  !kg
                
              end if
              !rtb gwflow: add ul_excess to runoff storage
              if(gw_transfer_flag.eq.1) then
                satexq(j) = satexq(j) + ul_excess !saturation excess (mm) leaving HRU soil profile on current day
              endif
            end if
          end do
          !compute tile flow again after saturation redistribution
        end if
      end if

      return
      end subroutine swr_satexcess