      subroutine lsu_read_elements
   
      use input_file_module
      use maximum_data_module
      use calibration_data_module
      use hydrograph_module
      use output_landscape_module
      
      implicit none

      character (len=80) :: titldum   !           |title of file
      character (len=80) :: header    !           |header of file
      integer :: eof                  !           |end of file
      integer :: imax                 !none       |determine max number for array (imax) and total number in file
      integer :: nspu                 !           |
      integer :: i_exist              !none       |check to determine if file exists
      integer :: mcal                 !           |
      integer :: mlsu                 !none       |counter
      integer :: i                    !none       |counter
      integer :: k                    !           |
      integer :: isp                  !           |
      integer :: ielem                !none       |counter
      integer :: ii                   !none       |counter
      integer :: ie1                  !none       |counter
      integer :: ie2                  !none       |counter
      integer :: ie                   !none       |counter
      integer :: iihru                !none       |counter
            
      imax = 0
      mcal = 0
            
    !! read landscape cataloging unit definitions for output (old subbasin output file)
    inquire (file=in_regs%def_lsu, exist=i_exist)
    if (i_exist /= 0 .or. in_regs%def_lsu /= "null") then
      do
        open (107,file=in_regs%def_lsu)
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) mlsu
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        
        !allocate subbasin (landscape unit) inputs and outputs
        allocate (lsu_out(0:mlsu))
                
        !allocate subbasin (landscape unit) inputs and outputs
        allocate (ruwb_d(0:mlsu)); allocate (ruwb_m(0:mlsu)); allocate (ruwb_y(0:mlsu)); allocate (ruwb_a(0:mlsu))
        allocate (runb_d(0:mlsu)); allocate (runb_m(0:mlsu)); allocate (runb_y(0:mlsu)); allocate (runb_a(0:mlsu))
        allocate (ruls_d(0:mlsu)); allocate (ruls_m(0:mlsu)); allocate (ruls_y(0:mlsu)); allocate (ruls_a(0:mlsu))
        allocate (rupw_d(0:mlsu)); allocate (rupw_m(0:mlsu)); allocate (rupw_y(0:mlsu)); allocate (rupw_a(0:mlsu))
        
      do i = 1, mlsu

        read (107,*,iostat=eof) k, lsu_out(i)%name, lsu_out(i)%area_ha, nspu
        
        if (eof < 0) exit
        if (nspu > 0) then
          allocate (elem_cnt(nspu))
          backspace (107)
          read (107,*,iostat=eof) k, lsu_out(i)%name, lsu_out(i)%area_ha, nspu, (elem_cnt(isp), isp = 1, nspu)

          !!save the object number of each defining unit
          if (nspu == 1) then
            allocate (lsu_out(i)%num(1))
            lsu_out(i)%num_tot = 1
            lsu_out(i)%num(1) = elem_cnt(1)
          end if
          if (nspu > 1) then
          ielem = 0
          do ii = 2, nspu
            ie1 = elem_cnt(ii-1)
            if (elem_cnt(ii) > 0) then
              if (ii == nspu) then
                ielem = ielem + 1
              else
                if (elem_cnt(ii+1) > 0) then
                  ielem = ielem + 1
                end if
              end if
            else
              ielem = ielem + abs(elem_cnt(ii)) - elem_cnt(ii-1) + 1
            end if
          end do
          allocate (lsu_out(i)%num(ielem))
          lsu_out(i)%num_tot = ielem

          ielem = 0
          ii = 1
          do while (ii <= nspu)
            ie1 = elem_cnt(ii)
            if (ii == nspu) then
              ielem = ielem + 1
              ii = ii + 1
              lsu_out(i)%num(ielem) = ie1
            else
              ie2 = elem_cnt(ii+1)
              if (ie2 > 0) then
                ielem = ielem + 1
                lsu_out(i)%num(ielem) = ie1
                ielem = ielem + 1
                lsu_out(i)%num(ielem) = ie2
              else
                ie2 = abs(ie2)
                do ie = ie1, ie2
                  ielem = ielem + 1
                  lsu_out(i)%num(ielem) = ie
                end do
              end if
              ii = ii + 2
            end if
          end do
          deallocate (elem_cnt)
          end if   !nspu > 1
        else
          allocate (lsu_out(i)%num(0:0))
          !!all hrus are in region 
          !allocate (lsu_out(i)%num(sp_ob%hru))
          !lsu_out(i)%num_tot = sp_ob%hru
          !do iihru = 1, sp_ob%hru
          !  lsu_out(i)%num(iihru) = iihru
          !end do      
        end if

      end do    ! i = 1, mreg

      db_mx%lsu_out = mlsu
      
      exit
      end do 
    end if	  

      !!read data for each element in all landscape cataloging units
      inquire (file=in_regs%ele_lsu, exist=i_exist)
      if (i_exist /= 0 .or. in_regs%ele_lsu /= "null") then
      do
        open (107,file=in_regs%ele_lsu)
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        if (eof < 0) exit
        imax = 0
          do while (eof == 0)
              read (107,*,iostat=eof) i
              if (eof < 0) exit
              imax = Max(i,imax)
          end do

        allocate (lsu_elem(imax))

        rewind (107)
        read (107,*) titldum
        read (107,*) header

        db_mx%lsu_elem = imax
        
        do isp = 1, imax
          read (107,*,iostat=eof) i
          backspace (107)
          read (107,*,iostat=eof) k, lsu_elem(i)%name, lsu_elem(i)%obtyp, lsu_elem(i)%obtypno,      &
                                    lsu_elem(i)%bsn_frac, lsu_elem(i)%ru_frac
          if (eof < 0) exit
        end do
        exit
      end do
      end if

      close (107)

      return
      end subroutine lsu_read_elements