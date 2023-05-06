      subroutine wallo_withdraw (iwallo, idmd, isrc, dmd_m3)
      
      use water_allocation_module
      use hydrograph_module
      use aquifer_module
      use reservoir_module
      use time_module
      
      implicit none 

      integer, intent (in):: iwallo         !water allocation object number
      integer, intent (in) :: idmd          !water demand object number
      integer, intent (in) :: isrc          !source object number
      real, intent (in) :: dmd_m3           !m3     |demand
      integer :: isrc_wallo                 !source object number
      integer :: j                  !none       |hru number
      real :: res_min               !m3         |min reservoir volume for withdrawal
      real :: res_vol               !m3         |reservoir volume after withdrawal
      real :: cha_min               !m3         |minimum allowable flow in channel after withdrawal
      real :: cha_div               !m3         |maximum amount of flow that can be diverted
      real :: rto                   !none       |ratio of channel withdrawal to determine hydrograph removed
      real :: avail                 !m3         |water available to withdraw from an aquifer
        
      !! zero withdrawal hyd for the demand source
      ht5 = hz

      !! check if water is available from each source - set withdrawal and unmet
      select case (wallo(iwallo)%dmd(idmd)%src_ob(isrc)%ob_typ)
      !! divert flowing water from channel source
      case ("cha")
        j = wallo(iwallo)%dmd(idmd)%src_ob(isrc)%ob_num
        isrc_wallo = wallo(iwallo)%dmd(idmd)%src(isrc)%src
        cha_min = wallo(iwallo)%src(isrc_wallo)%limit_mon(time%mo) * 86400.  !m3 = m3/s * 86400s/d
        !! amount that can be diverted without falling below low flow limit
        cha_div = ht2%flo - cha_min
        if (dmd_m3 < cha_div) then
          rto = dmd_m3 / ht2%flo
          ht5 = rto * ht2
          ht2 = (1. - rto) * ht2
          wallod_out(iwallo)%dmd(idmd)%src(isrc)%withdr = wallod_out(iwallo)%dmd(idmd)%src(isrc)%withdr + dmd_m3
        else
          wallod_out(iwallo)%dmd(idmd)%src(isrc)%unmet = wallod_out(iwallo)%dmd(idmd)%src(isrc)%unmet + dmd_m3
        end if
            
        !! reservor source
        case ("res") 
          j = wallo(iwallo)%dmd(idmd)%src_ob(isrc)%ob_num
          isrc_wallo = wallo(iwallo)%dmd(idmd)%src(isrc)%src
          res_min = wallo(iwallo)%src(isrc_wallo)%limit_mon(time%mo) * res_ob(j)%pvol
          res_vol = res(j)%flo - dmd_m3
          if (res_vol > res_min) then
            rto = dmd_m3 / res(j)%flo
            ht5 = rto * res(j)
            res(j) = (1. - rto) * res(j)
            wallod_out(iwallo)%dmd(idmd)%src(isrc)%withdr = wallod_out(iwallo)%dmd(idmd)%src(isrc)%withdr + dmd_m3
          else
            wallod_out(iwallo)%dmd(idmd)%src(isrc)%unmet = wallod_out(iwallo)%dmd(idmd)%src(isrc)%unmet + dmd_m3
          end if
         
        !! aquifer source
        case ("aqu") 
          j = wallo(iwallo)%dmd(idmd)%src_ob(isrc)%ob_num
          isrc_wallo = wallo(iwallo)%dmd(idmd)%src(isrc)%src
          avail = (wallo(iwallo)%src(isrc_wallo)%limit_mon(time%mo) - aqu_d(j)%dep_wt)  * aqu_dat(j)%spyld
          avail = avail * 10000. * aqu_prm(j)%area_ha     !m3 = 10,000*ha*m
          if (dmd_m3 < avail) then
            !! only have flow, no3, and minp(solp) for aquifer
            ht5%flo = dmd_m3
            aqu_d(j)%stor = aqu_d(j)%stor - (dmd_m3 / (10. * aqu_prm(j)%area_ha))  !mm = m3/(10.*ha)
            ht5%no3 = rto * aqu_d(j)%no3_st
            aqu_d(j)%no3_st = (1. - rto) * aqu_d(j)%no3_st
            ht5%solp = rto * aqu_d(j)%minp
            aqu_d(j)%minp = (1. - rto) * aqu_d(j)%minp
            wallod_out(iwallo)%dmd(idmd)%src(isrc)%withdr = wallod_out(iwallo)%dmd(idmd)%src(isrc)%withdr + dmd_m3
          else
            wallod_out(iwallo)%dmd(idmd)%src(isrc)%unmet = wallod_out(iwallo)%dmd(idmd)%src(isrc)%unmet + dmd_m3
          end if
        
          !! unlimited groundwater source
          case ("unl")
            !! only have flow, no3, and minp(solp) for aquifer
            ht5%flo = dmd_m3
            wallod_out(iwallo)%dmd(idmd)%src(isrc)%withdr = wallod_out(iwallo)%dmd(idmd)%src(isrc)%withdr + dmd_m3
          end select
          
          !! add source withdrawal hyd to get total withdrawal hyd for the demand object
          wallo(iwallo)%dmd(idmd)%hd = wallo(iwallo)%dmd(idmd)%hd + ht5
          
          !! subtract withdrawal from unmet
          wallo(iwallo)%dmd(idmd)%unmet_m3 = wallo(iwallo)%dmd(idmd)%unmet_m3 - wallod_out(iwallo)%dmd(idmd)%src(isrc)%withdr
          
      return
      end subroutine wallo_withdraw