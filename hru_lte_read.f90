      subroutine hru_lte_read
      
      use jrw_datalib_module, only: pldb
      use hru_lte_module
      use climate_module
      use input_file_module
      
      character (len=80) :: titldum
      character (len=80) :: header
      integer :: eof, grow_start, grow_end
      real :: rtos, rto3
      
      eof = 0
      imax = 0

      ! soil textures
      awct(1) = .056
      awct(2) = .116
      awct(3) = .115
      awct(4) = .114
      awct(5) = .186
      awct(6) = .254
      awct(7) = .141
      awct(8) = .138
      awct(9) = .113
      awct(10) = .130
      awct(11) = .132
      awct(12) = .113  !.123
      
      port(1) = .40
      port(2) = .40
      port(3) = .40
      port(4) = .40
      port(5) = .43
      port(6) = .47
      port(7) = .43
      port(8) = .40
      port(9) = .40
      port(10) = .40
      port(11) = .48
      port(12) = .47
      
      scon(1) = 105.
      scon(2) = 60. 
      scon(3) = .26 
      scon(4) = 13.2
      scon(5) = 6.8 
      scon(6) = 5.6 
      scon(7) = 4.3 
      scon(8) = 2.3 
      scon(9) = 1.5 
      scon(10) = 1.2 
      scon(11) = 0.9 
      scon(12) = 1.5 

      a1 = .2 
      a2 = .8 
      
      inquire (file=in_hru%hru_ez, exist=i_exist)
      if (i_exist == 0 .or. in_hru%hru_ez == 'null') then
        allocate (hlt_db(0:0))
      else
      do
        open (1,file=in_hru%hru_ez)
        read (1,*,iostat=eof) titldum
        if (eof < 0) exit
        read (1,*,iostat=eof) header
        if (eof < 0) exit
         do while (eof == 0)
            read (1,*,iostat=eof) i
            if (eof < 0) exit
            imax = Max(imax,i)
         end do
         
        !assumes data for each hru -> ok since there is only one file
        allocate (hlt_db(0:imax))
        allocate (hlt(sp_ob%hru_lte))
        
        rewind (1)
        read (1,*) titldum
        read (1,*) header

      do isd_h = 1, imax
        read (1,*,iostat = eof) i
        backspace (1)
        read (1,*,iostat=eof) k, hlt_db(i)
        if (eof < 0) exit
      end do

      do i = 1, sp_ob%hru_lte
         icmd = sp_ob1%hru_lte + i - 1
         idb = ob(icmd)%props
         hlt(i)%name = ob(icmd)%name
         hlt(i)%props = idb
         hlt(i)%obj_no = icmd
         hlt(i)%km2 = hlt_db(idb)%dakm2
         hlt(i)%cn2 = hlt_db(idb)%cn2
         hlt(i)%cn2 = amax1(35., hlt(i)%cn2)
         hlt(i)%cn2 = amin1(98., hlt(i)%cn2)
         hlt(i)%etco = hlt_db(idb)%etco
         hlt(i)%perco = hlt_db(idb)%perco
         hlt(i)%tdrain = hlt_db(idb)%tdrain
         hlt(i)%revapc = hlt_db(idb)%revapc
         hlt(i)%plant = hlt_db(idb)%plant
         hlt(i)%stress = hlt_db(idb)%stress
         hlt(i)%sw = hlt_db(idb)%sw * awct(hlt_db(idb)%itext) *               &
                                       hlt_db(idb)%soildep !* 1000.
         hlt(i)%awc = awct(hlt_db(idb)%itext) * hlt_db(idb)%soildep !* 1000.
         hlt(i)%por = port(hlt_db(idb)%itext) * hlt_db(idb)%soildep !* 1000.
         hlt(i)%sc = scon(hlt_db(idb)%itext)
         hlt(i)%hk = (hlt(i)%por - hlt(i)%awc) / (scon(hlt_db(idb)%itext))
         hlt(i)%hk = Max(2., hlt(i)%hk)
         hlt_db(idb)%abf = EXP(-hlt_db(idb)%abf) 
         qn1 = hlt_db(idb)%cn2 - (20. * (100. - hlt_db(idb)%cn2)) /        &
            (100.-hlt_db(idb)%cn2 + EXP(2.533-.063*(100.-hlt_db(idb)%cn2)))
         qn1 = Max(qn1, .4 * hlt_db(idb)%cn2)
         qn3 = hlt_db(idb)%cn2 * EXP(.00673*(100.-hlt_db(idb)%cn2)) 
         hlt(i)%smx = 254. * (100. / qn1 - 1.) 
         s3 = 254. * (100. / qn3 - 1.)
         rto3 = 1. - s3 / hlt(i)%smx
         rtos = 1. - 2.54 / hlt(i)%smx
         sumul = hlt(i)%por
         sumfc = hlt(i)%awc + hlt_db(idb)%cn3_swf * (sumul - hlt(i)%awc)
         !! calculate shape parameters
         call ascrv(rto3, rtos, sumfc, sumul, hlt(i)%wrt1, hlt(i)%wrt2)
         
         xi = 30. * mo - 15. 
         xx = hlt_db(idb)%xlat / 57.3 
         hlt(i)%yls = SIN(xx) 
         hlt(i)%ylc = COS(xx) 
         hlt(i)%phu = 2000. 
         hlt(i)%dm = 0. 
         hlt(i)%alai = .15 
         hlt(i)%g = 0. 
                  
         !crosswalk plant with plants.plt
         do ipl = 1, db_mx%plantparm
            if (hlt_db(idb)%plant == pldb(ipl)%plantnm) then
              hlt(i)%iplant = ipl
              exit
            endif
          end do
         
         !compute heat units from growing season and weather generator
         iwst = ob(icmd)%wst
         iwgn = wst(iwst)%wco%wgn
         if (hlt_db(idb)%igrow2 > hlt_db(idb)%igrow1) then
           grow_start = hlt_db(idb)%igrow1
           grow_end = hlt_db(idb)%igrow2
           hu_init = .15
         else
           grow_start = hlt_db(idb)%igrow2
           grow_end = hlt_db(idb)%igrow1
           hu_init = .85
         end if
         mo = 1
         imo = 2
         phutot = 0.
         do iday = 1, 365
           if (iday > ndays(imo)) then
             imo = imo + 1
             mo = mo + 1
           end if
           if (iday > grow_start .and. iday < grow_end) then
             tave = (wgn(iwgn)%tmpmx(mo) + wgn(iwgn)%tmpmn(mo)) / 2.
             iplt = hlt(i)%iplant
             phuday = tave - pldb(iplt)%t_base
             if (phuday > 0.) then
               phutot = phutot + phuday
             end if
           end if
         end do
         ! change from growing season to time to maturity
         hlt(i)%phu = .9 * phutot
         hlt(i)%phu = Max(500., hlt(i)%phu)
         if (pldb(iplt)%idc <= 2 .or. pldb(iplt)%idc == 4 .or. pldb(iplt)%idc == 5) then
           hlt(i)%phu = Min(2000., hlt(i)%phu)
         end if

         ! compute musle factors
         ! calculate USLE slope length factor
         xm = 0.
         sin_sl = 0.
         xm = .6 * (1. - EXP(-35.835 * hlt_db(idb)%slope))
         sin_sl = SIN(Atan(hlt_db(idb)%slope))
         hlt_db(idb)%uslels = (hlt_db(idb)%slopelen/22.128)**xm *          &
                        (65.41 * sin_sl * sin_sl + 4.56 * sin_sl + .065)
!     !      calculate composite usle value
         hlt(i)%uslefac = hlt_db(idb)%uslek * hlt_db(idb)%uslep *           &
           hlt_db(idb)%uslels * hlt_db(idb)%uslec * 11.8
         
        ! compute time of concentration using Kirpich equation
        IF (hlt_db(idb)%tc < 1.e-6) THEN
         ch_len = hlt_db(idb)%dakm2
         ch_sl = hlt_db(idb)%slope
         hlt_db(idb)%tc = .0078 * (ch_len * 3210.) ** .77 * sd_sl        &   
                                                        ** (-.385)
        END IF
        hlt_db(idb)%tc = hlt_db(idb)%tc * 60.     !!min to seconds
      end do

      !! dimension swatdeg output variables
      msd_h = sp_ob%hru_lte
      allocate (hltwb_d(msd_h))
      allocate (hltwb_m(msd_h))
      allocate (hltwb_y(msd_h))
      allocate (hltwb_a(msd_h))
      allocate (hltnb_d(msd_h))
      allocate (hltnb_m(msd_h))
      allocate (hltnb_y(msd_h))
      allocate (hltnb_a(msd_h))
      allocate (hltls_d(msd_h))
      allocate (hltls_m(msd_h))
      allocate (hltls_y(msd_h))
      allocate (hltls_a(msd_h))
      allocate (hltpw_d(msd_h))
      allocate (hltpw_m(msd_h))
      allocate (hltpw_y(msd_h))
      allocate (hltpw_a(msd_h))
       
      !allocate(qday(nyrs*366),qfdc(nyrs*366)) !Jaehak Jeong for fdc
      !allocate(pr(nyrs*366))
      !allocate(sd_qday(5*366),sd_qfdc(5*366)) !Jaehak Jeong for fdc
      !allocate(pr(5*366))
      
      exit
      end do
      endif
      
      close (1)
        
      return
      end subroutine hru_lte_read   