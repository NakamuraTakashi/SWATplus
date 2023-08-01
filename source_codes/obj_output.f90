      subroutine obj_output
      
      use time_module
      use hydrograph_module
      use soil_module
      use hru_module, only : ihru
      use organic_mineral_mass_module
      
      implicit none
        
      integer :: ihd           !            |
      integer :: iob           !            | 
      integer :: iunit         !            |
      integer :: itot          !none        |counter
      integer :: nlyr          
      integer :: nly
      integer :: j
      integer :: ii

      j = ihru
      
      do
        do itot = 1, mobj_out
          iob = ob_out(itot)%objno
          ihd = ob_out(itot)%hydno
          iunit = ob_out(itot)%unitno 
          
          if (iob <= sp_ob%objs) then
            select case (ob_out(itot)%hydno)
            case (1-5)  ! hydrograph output
              if (ob_out(itot)%hydtyp == "subday") then
                do ii = 1, time%step 
                  !write (iunit+itot,*) iob, time%yrc,time%day, ii, ob(iob)%hyd_flo(1,ii)
                  write (iunit+itot,*) time%day, time%mo, time%day_mo, time%yrc, ob(iob)%typ, ob(iob)%name, iob, ii, ob(iob)%hyd_flo(1,ii)
	            end do
              else  
                write (iunit+itot,*) time%day, time%mo, time%day_mo, time%yrc, ob(iob)%typ, ob(iob)%name, ob(iob)%hd(ihd)
              end if
              
            case (6)    ! soil water
              if (iob == 0) then
                do j = 1, sp_ob%hru
                 write (iunit+itot,*) time%day, time%mo, time%day_mo, time%yrc, ob(j)%name, ob(j)%typ, &
                   (soil(j)%phys(nly)%st, nly = 1,soil(j)%nly)
                end do
              else
                 write (iunit+itot,*) time%day, time%mo, time%day_mo, time%yrc, ob(iob)%name, ob(iob)%typ, &
                   (soil(iob)%phys(nly)%st, nly = 1,soil(iob)%nly)
              end if
            
            case (7)    ! soil nutrients
              if (iob == 0) then
                do j = 1, sp_ob%hru
                 write (iunit+itot,*) time%day, time%mo, time%day_mo, time%yrc, ob(j)%name, ob(j)%typ,      &
                   (soil1(j)%mn(nly), soil1(j)%hact(nly)%n, soil1(j)%hsta(nly)%n, soil1(j)%hact(nly)%n,     &
                   soil1(j)%hs(nly)%n, soil1(j)%hp(nly)%n, soil1(j)%rsd(nly)%n, soil1(j)%mp(nly),           &
                   soil1(j)%hact(nly)%p, soil1(j)%hsta(nly)%p, soil1(j)%hact(nly)%p, soil1(j)%hs(nly)%p,    &
                   soil1(j)%hp(nly)%p, soil1(j)%rsd(nly)%p, nly = 1, soil(iob)%nly)
                end do
              else
                 write (iunit+itot,*) time%day, time%mo, time%day_mo, time%yrc, ob(iob)%name, ob(iob)%typ,          &
                   (soil1(iob)%mn(nly), soil1(iob)%hact(nly)%n, soil1(iob)%hsta(nly)%n, soil1(iob)%hact(nly)%n,     &
                   soil1(iob)%hs(nly)%n, soil1(iob)%hp(nly)%n, soil1(iob)%rsd(nly)%n, soil1(iob)%mp(nly),           &
                   soil1(iob)%hact(nly)%p, soil1(iob)%hsta(nly)%p, soil1(iob)%hact(nly)%p, soil1(iob)%hs(nly)%p,    &
                   soil1(iob)%hp(nly)%p, soil1(iob)%rsd(nly)%p, nly = 1, soil(iob)%nly)
              end if
            end select
            
          end if
        end do
        exit
      end do
      
      return

      end subroutine obj_output