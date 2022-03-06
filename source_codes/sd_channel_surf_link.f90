        subroutine sd_channel_surf_link (ics)
                 
        use hydrograph_module
        use sd_channel_module
        use ru_module
        use hru_module, only : hru, ihru 
        use topography_data_module
      
        implicit none 
      
        integer, intent (in) :: ics     !none          |counter
        character (len=3) :: iobtyp     !none          |object type
        integer :: ii                   !none          |counter 
        integer :: i                    !              |
        integer :: iihru                !none          |hru counter 
        integer :: ihru_tot             !none          |total number of hru in the flood plain   

        ihru_tot = 0
        
        !! determine number of hru's
        do ii = 1, sd_ch(ics)%fp%obj_tot
          iobtyp = sd_ch(ics)%fp%obtyp(ii)     !object type
          select case (iobtyp)
          case ("hru")   !hru
            ihru_tot = ihru_tot + 1
          case ("ru")   !flood plain routing unit
            iru = sd_ch(ics)%fp%obtypno(ii)
            ihru_tot = ihru_tot + ru_def(iru)%num_tot
          end select
        end do
        
        allocate (sd_ch(ics)%fp%hru(ihru_tot))
        allocate (sd_ch(ics)%fp%hru_fr(ihru_tot))
   
        !! calculate total flood plain area and set hru numbers
        ihru_tot = 0
        sd_ch(ics)%fp%ha = 0.
        do ii = 1, sd_ch(ics)%fp%obj_tot
          iobtyp = sd_ch(ics)%fp%obtyp(ii)     !object type
          select case (iobtyp)
          case ("hru")   !hru
            ihru_tot = ihru_tot + 1
            ihru = sd_ch(ics)%fp%obtypno(ii)
            sd_ch(ics)%fp%hru(ihru_tot) = ihru
            sd_ch(ics)%fp%ha = sd_ch(ics)%fp%ha + hru(ihru)%area_ha
            
          case ("ru")   !flood plain routing unit
            iru = sd_ch(ics)%fp%obtypno(ii)

            !set flood plain link and landscape element (1==closest to river)
            do iihru = 1, ru_def(iru)%num_tot
              ihru_tot = ihru_tot + 1
              ihru = ru_def(iru)%num(iihru)
              sd_ch(ics)%fp%hru(ihru_tot) = ihru
              sd_ch(ics)%fp%ha = sd_ch(ics)%fp%ha + hru(ihru)%area_ha
            end do
      
          end select
        end do
   
        !set hru flood plain area fractions
        sd_ch(ics)%fp%hru_tot = ihru_tot
        do ihru = 1, sd_ch(ics)%fp%hru_tot
          sd_ch(ics)%fp%hru_fr(ihru) = hru(ihru)%area_ha / sd_ch(ics)%fp%ha
        end do
            
            
        return

      end subroutine sd_channel_surf_link