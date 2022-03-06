      subroutine proc_cha
    
      use hydrograph_module
         
      implicit none
      
      integer :: irch                   !              |
      integer :: idat                   !              |
      integer :: i                      !none          |counter
         
      call ch_read_init

      call sd_hydsed_read
      call ch_read_hyd
      call ch_read_sed
      call ch_read_nut
      call ch_read
      call sd_channel_read
      call sd_hydsed_init

      call channel_allo
          
      !! intialize stream-aquifer interactions for geomorphic baseflow
      !! aquifer to channel flow
      call aqu2d_init
      
      do ich = 1, sp_ob%chan
        !! initialize flow routing variables
        call ch_ttcoef (ich)
      end do
         
      do irch = 1, sp_ob%chan
        i = sp_ob1%chan + irch - 1 
        idat = ob(i)%props
        call ch_initial (idat, irch)
      end do
      
      !set parms for sd-channel-landscape linkage
      call overbank_read
      do irch = 1, sp_ob%chandeg
        i = sp_ob1%chandeg + irch - 1
        if (ob(i)%props2 > 0) then
          call sd_channel_surf_link (ob(i)%props2)
        end if
      end do
      
      call time_conc_init

	  return
      
      end subroutine proc_cha