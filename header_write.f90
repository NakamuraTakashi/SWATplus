     subroutine header_write
    
     use basin_module
     use aquifer_module
     use channel_module
     use reservoir_module
     use hydrograph_module
     use sd_channel_module
     use parm
     
      if (pco%fdcout == 'y') then
        open (6000,file="flow_duration_curve.out", recl=800)
        write (9000,*) 'FDC                 flow_duration_curve.out'
        write (6000,*) fdc_hdr
      end if 

!!!!!! hru-out.cal - hru soft calibration output including soft and predicted budgets and 
!!!!!! calibration parameter adjustments
      if (sp_ob%chandeg > 0) then
	    open (4999,file="hru-out.cal", recl = 800)
	    write (9000,*)   'HRU SOFT CALIB OUT      hru-out.cal'
	    write (4999,*) calb_hdr
      end if
      
!!!!!! hru-new.cal - hru soft calibration output file.  The same format as calibration.upd and
!!!!!! can be used as input (calibration.upd) in subsequent simulations
      if (db_mx%lsu_reg > 0) then
        open (5000,file="hru-new.cal", recl = 800)
        write (5000,*) ' calibration.upd_developed_from_soft_data_calibration'
	    write (9000,*)   'HRU SOFT OUT CALIB      hru-new.cal'
        write (5000,*) calb3_hdr
      end if
      
!!!!!! hru-lte-out.cal - hru lte soft calibration output including soft and predicted budgets and 
!!!!!! calibration parameter adjustments
      !open (5003,file="hru-lte-out.cal", recl = 800)
	  !write (9000,*)   'LTE SOFT OUT CALIB      hru-lte-out.cal'
	  !write (5003,*) calb_hdr
	  
!!!!!! hru-lte-new.cal - hru lte soft calibration output file.  The same format as hru-lte.hru and
!!!!!! can be used as input (hru-lte.hru) in subsequent simulations 
      !open (5002,file="hru-lte-new.cal", recl = 800)
	  !write (9000,*)   'LTE SOFT CALIB INPUT    hru-lte-new.cal'
	  !write (5002,*) calb2_hdr
      
!! BASIN AQUIFER OUTPUT
        if (pco%aqu_bsn%d == 'y') then
          open (2090,file="aquifer_day_bsn.txt", recl = 1500)
          write (2090,*) aqu_hdr  
          write (9000,*) 'BASIN AQUIFER       aquifer_day_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2094,file="aquifer_day_bsn.csv", recl = 1500)
            write (2094,'(*(G0.3,:","))') aqu_hdr
            write (9000,*) 'BASIN AQUIFER               aquifer_day_bsn.csv'
          end if
        endif
        
      if (pco%aqu_bsn%m == 'y') then
        open (2091,file="aquifer_mon_bsn.txt",recl = 1500)      
        write (2091,*) aqu_hdr 
        write (9000,*) 'BASIN AQUIFER       aquifer_mon_bsn.txt'
         if (pco%csvout == 'y') then 
           open (2095,file="aquifer_mon_bsn.csv",recl = 1500)
           write (2095,'(*(G0.3,:","))') aqu_hdr 
           write (9000,*) 'BASIN AQUIFER               aquifer_mon_bsn.csv'
         end if
      end if 
      
      if (pco%aqu_bsn%y == 'y') then
        open (2092,file="aquifer_yr_bsn.txt",recl = 1500)      
        write (2092,*) aqu_hdr
        write (9000,*) 'BASIN AQUIFER       aquifer_yr_bsn.txt'
         if (pco%csvout == 'y') then 
           open (2096,file="aquifer_mon_yr.csv",recl = 1500)
           write (2096,'(*(G0.3,:","))') aqu_hdr 
           write (9000,*) 'BASIN AQUIFER               aquifer_yr_bsn.csv'
         end if
      end if 
      
     if (pco%aqu_bsn%a == 'y') then
        open (2093,file="aquifer_aa_bsn.txt",recl = 1500)      
        write (2093,*) aqu_hdr 
        write (9000,*) 'BASIN AQUIFER       aquifer_aa_bsn.txt'
         if (pco%csvout == 'y') then 
           open (2097,file="aquifer_aa_bsn.csv",recl = 1500)
           write (2097,'(*(G0.3,:","))') aqu_hdr 
           write (9000,*) 'BASIN AQUIFER               aquifer_aa_bsn.csv'
         end if
      end if 
!! BASIN AQUIFER OUTPUT

!! BASIN RESERVOIR OUTPUT
        if (pco%res_bsn%d == 'y') then
          open (2100,file="reservoir_day_bsn.txt", recl = 1500)
          write (2100,*) res_hdr  
          write (9000,*) 'BASIN RESERVOIR     reservoir_day_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2104,file="reservoir_day_bsn.csv", recl = 1500)
            write (2104,'(*(G0.3,:","))') res_hdr
            write (9000,*) 'BASIN RESERVOIR     reservoir_day_bsn.csv'
          end if
        endif
        
      if (pco%res_bsn%m == 'y') then
        open (2101,file="reservoir_mon_bsn.txt",recl = 1500)      
        write (2101,*) res_hdr 
        write (9000,*) 'BASIN RESERVOIR AA  reservoir_mon_bsn.txt'
       if (pco%csvout == 'y') then 
          open (2105,file="reservoir_mon_bsn.csv",recl = 1500)
          write (2105,'(*(G0.3,:","))') res_hdr 
          write (9000,*) 'BASIN RESERVOIR     reservoir_mon_bsn.csv'
       end if
      end if
       
       if (pco%res_bsn%y == 'y') then
          open (2102,file="reservoir_yr_bsn.txt", recl = 1500)
          write (2102,*) res_hdr 
          write (9000,*) 'BASIN RESERVOIR     reservoir_yr_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2106,file="reservoir_yr_bsn.csv", recl = 1500)
            write (2106,'(*(G0.3,:","))') res_hdr
            write (9000,*) 'BASIN RESERVOIR     reservoir_yr_bsn.csv'
          end if
       endif
        
      if (pco%res_bsn%a == 'y') then
       open (2103,file="reservoir_aa_bsn.txt",recl = 1500)      
        write (2103,*) res_hdr 
        write (9000,*) 'BASIN RESERVOIR AA  reservoir_aa_bsn.txt'
       if (pco%csvout == 'y') then 
          open (2107,file="reservoir_aa_bsn.csv",recl = 1500)
          write (2107,'(*(G0.3,:","))') res_hdr 
          write (9000,*) 'BASIN RESERVOIR               reservoir_aa_bsn.csv'
       end if
      end if
!! BASIN RESERVOIR OUTPUT
      
!! BASIN CHANNEL OUTPUT
        if (pco%chan_bsn%d == 'y') then
          open (2110,file="channel_day_bsn.txt", recl = 1500)
          write (2110,*) ch_hdr  
          write (9000,*) 'BASIN CHANNEL       channel_day_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2114,file="channel_day_bsn.csv", recl = 1500)
            write (2114,'(*(G0.3,:","))') ch_hdr
            write (9000,*) 'BASIN CHANNEL       channel_day_bsn.csv'
          end if
        endif
        
       if (pco%chan_bsn%m == 'y') then
        open (2111,file="channel_mon_bsn.txt",recl = 1500)      
        write (2111,*) ch_hdr
        write (9000,*) 'BASIN CHANNEL       channel_mon_bsn.txt'
         if (pco%csvout == 'y') then 
           open (2115,file="channel_mon_bsn.csv",recl = 1500)
           write (2115,'(*(G0.3,:","))') ch_hdr 
           write (9000,*) 'BASIN CHANNEL       channel_mon_bsn.txt'
         end if
        end if
       
        if (pco%chan_bsn%y == 'y') then
          open (2112,file="channel_yr_bsn.txt", recl = 1500)
          write (2112,*) ch_hdr
          write (9000,*) 'BASIN CHANNEL       channel_yr_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2116,file="channel_yr_bsn.csv", recl = 1500)
            write (2116,'(*(G0.3,:","))') ch_hdr
            write (9000,*) 'BASIN CHANNEL       channel_yr_bsn.csv'
          end if
        endif
        
        if (pco%chan_bsn%a == 'y') then
          open (2113,file="channel_aa_bsn.txt",recl = 1500)      
          write (2113,*) ch_hdr 
          write (9000,*) 'BASIN CHANNEL       channel_aa_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2117,file="channel_aa_bsn.csv",recl = 1500)
            write (2117,'(*(G0.3,:","))') ch_hdr
            write (9000,*) 'BASIN CHANNEL       channel_aa_bsn.csv'
          end if
        end if
!! BASIN CHANNEL OUTPUT

!! BASIN SWAT DEG CHANNEL OUTPUT
        if (pco%sd_chan_bsn%d == 'y') then
          open (2120,file="channel_sd_day_bsn.txt", recl = 1500)
          write (2120,*) sdch_hdr 
          write (9000,*) 'BASIN SWAT DEGCHAN  channel_sd_day_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2124,file="channel_sd_day_bsn.csv", recl = 1500)
            write (2124,'(*(G0.3,:","))') sdch_hdr
            write (9000,*) 'BASIN SWAT DEG CHAN channel_sd_day_bsn.csv'
          end if
        endif
        
       if (pco%sd_chan_bsn%m == 'y') then
        open (2121,file="channel_sd_mon_bsn.txt",recl = 1500)      
        write (2121,*) sdch_hdr 
        write (9000,*) 'BASIN SWAT DEG CHAN channel_sd_mon_bsn.txt'
         if (pco%csvout == 'y') then 
           open (2125,file="channel_sd_mon_bsn.csv",recl = 1500)
           write (2125,'(*(G0.3,:","))') sdch_hdr 
           write (9000,*) 'BASIN SWAT DEG CHANNEL               channel_sd_mon_bsn.txt'
         end if
        end if
       
        if (pco%sd_chan_bsn%y == 'y') then
          open (2122,file="channel_sd_yr_bsn.txt", recl = 1500)
          write (2122,*) sdch_hdr 
          write (9000,*) 'BASIN SWAT DEG CHAN channel_sd_yr_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2126,file="channel_sd_yr_bsn.csv", recl = 1500)
            write (2126,'(*(G0.3,:","))') sdch_hdr
            write (9000,*) 'BASIN SWAT DEG CHAN channel_sd_yr_bsn.csv'
          end if
        endif
        
        if (pco%sd_chan_bsn%a == 'y') then
          open (2123,file="channel_sd_aa_bsn.txt",recl = 1500)      
          write (2123,*) sdch_hdr 
          write (9000,*) 'BASIN SWAT DEG CHAN channel_sd_aa_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2127,file="channel_sd_aa_bsn.csv",recl = 1500)
            write (2127,'(*(G0.3,:","))') sdch_hdr 
            write (9000,*) 'BASIN SWAT DEG CHAN channel_sd_aa_bsn.csv'
          end if
        end if
!! BASIN SWAT DEG CHANNEL OUTPUT


!! BASIN RECALL OUTPUT
        if (pco%recall_bsn%d == 'y') then
          open (2130,file="pts_day_bsn.txt", recl = 1500)
          write (2130,*) hyd_hdr
          write (9000,*) 'BASIN RECALL        pts_day_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2134,file="pts_day_bsn.csv", recl = 1500)
            write (2134,'(*(G0.3,:","))') hyd_hdr
            write (9000,*) 'BASIN RECALL        pts_day_bsn.csv'
          end if
        endif
        
        if (pco%recall_bsn%m == 'y') then
        open (2131,file="pts_mon_bsn.txt",recl = 1500)      
        write (2131,*) hyd_hdr
        write (9000,*) 'BASIN RECALL        pts_mon_bsn.txt'
         if (pco%csvout == 'y') then 
            open (2135,file="pts_mon_bsn.csv",recl = 1500)
            write (2135,'(*(G0.3,:","))') hyd_hdr 
            write (9000,*) 'BASIN RECALL        pts_mon_bsn.csv'
         end if
       end if
       
        if (pco%recall_bsn%y == 'y') then
          open (2132,file="pts_yr_bsn.txt", recl = 1500)
          write (2132,*) hyd_hdr 
          write (9000,*) 'BASIN RECALL        pts_yr_bsn.txt'
          if (pco%csvout == 'y') then 
            open (2136,file="pts_yr_bsn.csv", recl = 1500)
            write (2136,'(*(G0.3,:","))') hyd_hdr
            write (9000,*) 'BASIN RECALL        pts_yr_bsn.csv'
          end if
        endif
        
        if (pco%recall_bsn%a == 'y') then 
        open (2133,file="pts_aa_bsn.txt",recl = 1500)      
        write (2133,*) hyd_hdr
        write (9000,*) 'BASIN RECALL AA     pts_aa_bsn.txt'
         if (pco%csvout == 'y') then 
            open (2137,file="pts_aa_bsn.csv",recl = 1500)
            write (2137,'(*(G0.3,:","))') hyd_hdr
            write (9000,*) 'BASIN RECALL        pts_aa_bsn.csv'
         end if
        end if
!! BASIN RECALL OUTPUT
      return
      end subroutine header_write  