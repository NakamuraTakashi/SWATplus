      subroutine obj_output
      
      use time_module
      use hydrograph_module
      
      implicit none
        
      integer :: ihd           !            |
      integer :: iob           !            | 
      integer :: iunit         !            |
      integer :: itot          !none        |counter
      
      do
        do itot = 1, mobj_out
          iob = ob_out(itot)%objno
          ihd = ob_out(itot)%hydno
          iunit = ob_out(itot)%unitno          
                
          if (iob <= sp_ob%objs) then
            write (iunit+itot,*) time%yrc, time%day, ob(iob)%hd(ihd)   
          end if

        end do
        exit
      end do
      
      return

      end subroutine obj_output