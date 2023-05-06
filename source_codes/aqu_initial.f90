      subroutine aqu_initial 
    
      use aquifer_module  
      use hydrograph_module
      use constituent_mass_module
      use aqu_pesticide_module
       
      implicit none
      
      character (len=500) :: header    !header for output file
      character (len=80) :: titldum    !title 
      integer :: iaq                   !none      |counter
      integer :: iob                   !          | 
      integer :: iaqdb                 !          | 
      integer :: ipest                 !none      |counter
      integer :: ipath                 !          | 
      integer :: isalt                 !          | 
      integer :: i                     !none      |counter
      integer :: init                  !          | 
      integer :: idat

      !allocate objects for each aquifer
      allocate (aqu_om_init(sp_ob%aqu))
      allocate (aqu_d(sp_ob%aqu))
      allocate (aqu_dat(sp_ob%aqu))
      allocate (aqu_prm(sp_ob%aqu))
      allocate (aqu_m(sp_ob%aqu))
      allocate (aqu_y(sp_ob%aqu))
      allocate (aqu_a(sp_ob%aqu))
      allocate (cs_aqu(sp_ob%aqu))
      allocate (aqupst_d(sp_ob%aqu))
      allocate (aqupst_m(sp_ob%aqu))
      allocate (aqupst_y(sp_ob%aqu))
      allocate (aqupst_a(sp_ob%aqu))

      if (cs_db%num_pests > 0) then
        allocate (baqupst_d%pest(cs_db%num_pests))
        allocate (baqupst_m%pest(cs_db%num_pests))
        allocate (baqupst_y%pest(cs_db%num_pests))
        allocate (baqupst_a%pest(cs_db%num_pests))
      end if
      
      do iaq = 1, sp_ob%aqu
        if (cs_db%num_pests > 0) then
          !! allocate constituents
          allocate (cs_aqu(iaq)%pest(cs_db%num_pests))
          allocate (aqupst_d(iaq)%pest(cs_db%num_pests))
          allocate (aqupst_m(iaq)%pest(cs_db%num_pests))
          allocate (aqupst_y(iaq)%pest(cs_db%num_pests))
          allocate (aqupst_a(iaq)%pest(cs_db%num_pests))
          allocate (cs_aqu(iaq)%path(cs_db%num_paths))
          allocate (cs_aqu(iaq)%hmet(cs_db%num_metals))
          allocate (cs_aqu(iaq)%salt(cs_db%num_salts))
        end if
              
        iob = sp_ob1%aqu + iaq - 1
        iaqdb = ob(iob)%props

        !! initialize parameters
        aqu_dat(iaq) = aqudb(iaqdb)
        
        aqu_prm(iaq)%area_ha = ob(iob)%area_ha
        aqu_prm(iaq)%alpha_e = Exp(-aqu_dat(iaq)%alpha)
        aqu_prm(iaq)%nloss = Exp(-.693 / (aqu_dat(iaq)%hlife_n + .1))
        
        aqu_d(iaq)%flo = aqu_dat(iaq)%flo
        aqu_d(iaq)%dep_wt = aqu_dat(iaq)%dep_wt
        aqu_d(iaq)%stor = 1000. * (aqu_dat(iaq)%dep_bot - aqu_d(iaqdb)%dep_wt) * aqu_dat(iaq)%spyld
        !! convert ppm -> kg    (m3=10*mm*ha)     kg=m3*ppm/1000
        aqu_d(iaq)%no3_st = (10. * aqu_d(iaq)%flo * aqu_prm(iaq)%area_ha) * aqu_dat(iaq)%no3 / 1000.
        aqu_d(iaq)%minp = 0.
        aqu_d(iaq)%cbn = aqu_dat(iaq)%cbn
        aqu_d(iaq)%rchrg = 0.
        aqu_d(iaq)%seep = 0.
        aqu_d(iaq)%revap = 0.
        aqu_d(iaq)%no3_rchg = 0.
        aqu_d(iaq)%no3_loss = 0.
        aqu_d(iaq)%no3_lat = 0.
        aqu_d(iaq)%no3_seep = 0.
        aqu_d(iaq)%flo_cha = 0.
        aqu_d(iaq)%flo_res = 0.
        aqu_d(iaq)%flo_ls = 0
      end do
            
      ! pesticides and constituents are initialized in aqu_read_init

      return
      end subroutine aqu_initial         