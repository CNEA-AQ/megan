!---------------------------------------------------------------
! programed by: Ramiro A. Espada. April 2023.
! programed by: Hui Wang, April 2025, Change the way of interpolation
! Based on Prep_code &  MEGEFP32 (UCI-BAI-MEGAN)
!---------------------------------------------------------------
module prep_megan

  use area_mapper_grw
!, only : proj_init, discrete_frac
  use bio_types
!,       only : grid_specs, grid_cnt, grid_ndx
  use constants_module, only : rad_per_deg, earth_radius_m
  use netcdf

  implicit none

  private
  public prep
  ! Parameters
  integer, parameter :: PS = 2
  integer, parameter :: mxetype = 100
  integer, parameter :: ncantype=6, nefs=19, nldfs=4
  ! ecotype reading
  ! 5869 is the maximum of ecotype ID
  integer, parameter :: max_id = 5869, ncat = 23, nveg = 4

  ! Variables
  integer :: ids, ide, jds, jde
  integer :: ierr, astat, dimid, varid, map_proj
  !integer :: x0, y0, ncolsin, nrowsin
  real :: missing_value, scale_factor
  real :: cen_lon, cen_lat, stand_lon, truelat1, truelat2, dx

  real(8), allocatable :: xedge_megan(:), yedge_megan(:)

  logical ::  var_flag

  character(len=132) :: varname
  character(len=200) :: filespec, wrffile, megan_dir, out_dir
  character(len=200) :: inpname
  character(len=100) :: fillvalue

contains

   subroutine prep(wrffile,nlai,lai_scale_factor,                    &
                   ecotypes_file, canopy_file, growtype_file, laiv_file, GtEcoEF_file,       &
                   run_BDSNP, nitro_file, fert_file, climate_file, landtype_file,idxs)
     implicit none
     integer, intent(in) :: nlai
     real,    intent(in) :: lai_scale_factor
     integer, intent(in) :: idxs(:)
     character(200), intent(in) :: wrffile,ecotypes_file,canopy_file,growtype_file,&
                                   laiv_file,climate_file,fert_file,&
                                   landtype_file,nitro_file,GtEcoEF_file
     logical       ,intent(in)  :: run_BDSNP
   
     ! get the wrf grid definition
     call wrf_file(wrffile, ide, jde, cen_lon, cen_lat, stand_lon, truelat1, truelat2, dx)
  
     ! interpolation
     grid_ndx = 0
   
     !Static data:
     call prep_static_data(idxs,ide,jde,canopy_file,growtype_file,ecotypes_file,GtEcoEF_file,climate_file,landtype_file, run_BDSNP)
          ! `CTF` (*Canopy Type Fractions*):
          ! `EFs` (*Emission Factors*)     : (~19) VOC family, and Canopy Type (6)
          ! `LDF` (*Light Dependent EF*)   :  4 VOC families, and Canopy Type (6)
          ! `arid`     (BDSNP): arid soils mask
          ! `landtype` (BDSNP): land type classification
   
     !Time/date dependent data:
     call prep_dynamic_data(idxs,ide,jde,laiv_file,nlai,lai_scale_factor,nitro_file,fert_file,run_BDSNP) 
     !subroutine prep_dynamic_data(idxs,ide,jde,laiv_file,nlai,lai_scale_factor,nitro_file, fert_file,run_BDSNP)
          ! `LAI`     monthly.         : Leaf Area Index
          ! `N_DEP:`  monthly. (BDSNP) : Nitrogen deposition   flux
          ! `N_FERT:` daily.   (BDSNP) : Nitrogen feritization flux
   
      print*, "========================================="
      print*, " prep-megan: Completed successfully"
      print*, "========================================="
   
   end subroutine

 !----------------------------------
 !  STATIC  DATA:
 !--------------------------------
    subroutine prep_static_data(idxs,ide,jde,canopy_file,ctf_file, ecotype_file, GtEcoEF_file, climate_file,landtype_file, run_BDSNP)
      !use area_mapper_grw, only: lon, lat
      implicit none
      integer,         intent(in) :: idxs(:)
      integer,         intent(in) :: ide,jde
      character(len=*),intent(in) :: canopy_file,ctf_file, ecotype_file, GtEcoEF_file  !input  files
      character(len=*),intent(in) :: climate_file, landtype_file  !input  files
      character(len=19)           :: outfile='prep_mgn_static.nc' !output file
      logical :: run_BDSNP
 

      !Coordinates & area
      integer :: ncid,var_id
      !integer :: xt,xe,yt,ye,nx,ny
      integer :: i!,j,k
      !CTF: 
      real, allocatable :: CTF(:,:,:)          !CTF buffer
      real, allocatable :: VCF(:,:)          !CTF buffer
      real, allocatable :: NeedleFrac(:,:)          !CTF buffer
      real, allocatable :: TropFrac(:,:)          !CTF buffer
      real, allocatable :: cell_area(:,:)          !CTF buffer
      !character(len=6)  :: CTF_LIST(6)         ![Ntr, Trop, Btr, shrub, herb, crop] ! tree]
      
      !Ecotype ID & Fraction, Hui
      integer, allocatable :: ecotypeid(:,:,:)
      real, allocatable    :: ecotypefrac(:,:,:)
      real, allocatable    :: ef_growthform(:,:,:,:)
      real, allocatable    :: ef_grid(:,:,:),ef_tree(:,:,:),ef_shrub(:,:,:),ef_herb(:,:,:),ef_crop(:,:,:)
      real, allocatable    :: ldf_grid(:,:,:),ldf_tree(:,:,:),ldf_shrub(:,:,:),ldf_herb(:,:,:),ldf_crop(:,:,:)
      real, allocatable    :: landgrid(:,:,:)   !LAND buffer
      !logical              :: debug_flag

      allocate(CTF(ide,jde,ncantype+1))
      allocate(TropFrac(ide,jde))
      allocate(cell_area(ide,jde))
      allocate(NeedleFrac(ide,jde))
      allocate(VCF(ide,jde))
      allocate(ecotypeid(ide,jde,mxetype))
      allocate(ecotypefrac(ide,jde,mxetype))
      allocate(ef_growthform(ide,jde,ncat,nveg))  
      
      allocate(ef_tree(ide,jde,nefs))  
      allocate(ef_shrub(ide,jde,nefs))  
      allocate(ef_herb(ide,jde,nefs))  
      allocate(ef_crop(ide,jde,nefs))  
      allocate(ef_grid(ide,jde,nefs))     
      
      allocate(ldf_tree(ide,jde,nldfs))     
      allocate(ldf_shrub(ide,jde,nldfs))     
      allocate(ldf_herb(ide,jde,nldfs))     
      allocate(ldf_crop(ide,jde,nldfs))     
      allocate(ldf_grid(ide,jde,nldfs))     
 
      print '("prep static file: ",A19,"..")',outfile
      call create_static_file(outfile,idxs,run_BDSNP)

      !==============================================================
      !=============Creating cell area======================
      !==============================================================
      cell_area = dx*dx
      call write_2d_var(outfile,"cell_area",cell_area,idxs)


      !==============================================================
      !=============Processing Canopy Tpye Data======================
      !==============================================================

      print '(A,1X,A)', "Reading:", ctf_file
      call interpolate_area(canopy_file,"nl_tree"  ,ide,jde,NeedleFrac)
      call interpolate_area(canopy_file,"trop_tree",ide,jde,TropFrac)
      call interpolate_area(ctf_file,"shrub"    ,ide,jde,CTF(:,:,4))
      call interpolate_area(ctf_file,"grass"    ,ide,jde,CTF(:,:,5))
      call interpolate_area(ctf_file,"crop"     ,ide,jde,CTF(:,:,6))
      call interpolate_area(ctf_file,"tree"     ,ide,jde,CTF(:,:,7))
 
      !needleleaf tree
      CTF(:,:,1)=CTF(:,:,7)*(NeedleFrac/100.0)*(1.0-TropFrac/100.0)
      !tropical tree
      CTF(:,:,2)=CTF(:,:,7)*(TropFrac/100.0)
      !boradleaf tree
      CTF(:,:,3)=CTF(:,:,7)*(1.0-TropFrac/100.0)*(1.0-NeedleFrac/100.0) 


      !WRITE CTF:
      CTF=CTF*0.01 ! % to fraction.
      
      call write_3d_var(outfile,"CTF",CTF,idxs,ncantype+1 )


      !==============================================================
      !=============Processing Emission Factor Data==================
      !==============================================================
      
      print '(A,1X,A)', "Reading:", ecotype_file
      call interpolate_ecotype(ecotype_file, "ecotype" ,ide,jde, ecotypeid,ecotypefrac)
      print*,"======================================"
      print '("Processing Emission Factor: ",A19,"..")',GtEcoEF_file
      print*,"======================================"
      call compute_ef_grid(ecotypeid, ecotypefrac, GtEcoEF_file, ef_growthform)
      !'Crop', 'Herb', 'Shrub', 'Tree'
      

      ef_tree = ef_growthform(:,:,1:nefs,4)
      ef_shrub= ef_growthform(:,:,1:nefs,3)
      ef_herb = ef_growthform(:,:,1:nefs,2)
      ef_crop = ef_growthform(:,:,1:nefs,1)

      ldf_tree = ef_growthform(:,:,(nefs+1):(nefs+nldfs),4)
      ldf_shrub= ef_growthform(:,:,(nefs+1):(nefs+nldfs),3)
      ldf_herb = ef_growthform(:,:,(nefs+1):(nefs+nldfs),2)
      ldf_crop = ef_growthform(:,:,(nefs+1):(nefs+nldfs),1)


      call write_3d_var(outfile,"EFS_TREE" ,ef_tree ,idxs,nefs )
      call write_3d_var(outfile,"EFS_SHRUB",ef_shrub,idxs,nefs )
      call write_3d_var(outfile,"EFS_HERB" ,ef_herb ,idxs,nefs )
      call write_3d_var(outfile,"EFS_CROP" ,ef_crop ,idxs,nefs )
      
      call write_3d_var(outfile,"LDF_TREE" ,ldf_tree,idxs,nldfs )
      call write_3d_var(outfile,"LDF_SHRUB",ldf_shrub,idxs,nldfs )
      call write_3d_var(outfile,"LDF_HERB" ,ldf_herb,idxs,nldfs )
      call write_3d_var(outfile,"LDF_CROP" ,ldf_crop,idxs,nldfs )

      VCF = sum(CTF(:,:,4:7),dim=3)
      do i=1,nefs
        ef_growthform(:,:,i,4) = ef_growthform(:,:,i,4)*CTF(:,:,7)/VCF!Tree 
        ef_growthform(:,:,i,3) = ef_growthform(:,:,i,3)*CTF(:,:,4)/VCF!Shrub 
        ef_growthform(:,:,i,2) = ef_growthform(:,:,i,2)*CTF(:,:,5)/VCF!Herb 
        ef_growthform(:,:,i,1) = ef_growthform(:,:,i,1)*CTF(:,:,6)/VCF!Crop 
      end do
      do i=1,nldfs
        ef_growthform(:,:,nefs+i,4) = ef_growthform(:,:,nefs+i,4)*CTF(:,:,7)/VCF!Tree 
        ef_growthform(:,:,nefs+i,3) = ef_growthform(:,:,nefs+i,3)*CTF(:,:,4)/VCF!Shrub 
        ef_growthform(:,:,nefs+i,2) = ef_growthform(:,:,nefs+i,2)*CTF(:,:,5)/VCF!Herb 
        ef_growthform(:,:,nefs+i,1) = ef_growthform(:,:,nefs+i,1)*CTF(:,:,6)/VCF!Crop 
      end do

      ef_grid(:,:,:) = sum(ef_growthform(:,:,1:nefs,:),dim=4) 
      ldf_grid(:,:,:) = sum(ef_growthform(:,:,(nefs+1):(nefs+nldfs),:),dim=4) 

      call write_3d_var(outfile,"EFS" ,ef_grid,idxs,nefs )
      call write_3d_var(outfile,"LDF" ,ldf_grid,idxs,nldfs )


      !--------
      if (run_BDSNP) then
        print*,"BDSNP (LAND)"    
        allocate(landgrid(ide,jde,2))
        landgrid(:,:,1) = 0  
        call interpolate_area(climate_file,"arid",ide,jde,landgrid(:,:,1))
        landgrid(:,:,2) = 1
        call interpolate_area(landtype_file,"landtype",ide,jde,landgrid(:,:,2))
        
        call write_2d_var(outfile,"arid" ,landgrid(:,:,1),idxs)
        call write_2d_var(outfile,"landtype" ,landgrid(:,:,2),idxs)

        deallocate(landgrid)
      end if



      deallocate(CTF)
      deallocate(VCF)
      deallocate(cell_area)
      deallocate(TropFrac)
      deallocate(NeedleFrac)
      deallocate(ecotypeid)
      deallocate(ecotypefrac)
      deallocate(ef_growthform)
      deallocate(ef_grid)
      deallocate(ef_tree)
      deallocate(ef_shrub)
      deallocate(ef_herb)
      deallocate(ef_crop)
      deallocate(ldf_grid)
      deallocate(ldf_tree)
      deallocate(ldf_shrub)
      deallocate(ldf_herb)
      deallocate(ldf_crop)

    end subroutine

 !----------------------------------
 !  DYNAMIC DATA:
 !---------------------------------
    subroutine prep_dynamic_data(idxs,ide,jde,laiv_file,nlai,lai_scale_factor,nitro_file, fert_file,run_BDSNP)
      
      implicit none
      integer,         intent(in) :: idxs(:)
      integer,         intent(in) :: ide,jde
      integer,         intent(in) :: nlai
      real,            intent(in) :: lai_scale_factor
      
      character(len=200),intent(in) :: laiv_file
      character(len=200),intent(in) :: nitro_file
      character(len=200),intent(in) :: fert_file
      real, allocatable :: laiv(:,:,:)
      real, allocatable :: NDEP(:,:,:)
      real, allocatable :: NFERT(:,:,:)
      logical :: run_BDSNP
      
      character(len=200) :: out_dyn_file='prep_mgn_dynamic.nc'
      
      !local variable 
      character(len=2):: kk
      character(len=3):: kkk
      integer :: k
 
      allocate( laiv(ide,jde,nlai))  
      print '("prep dynamic file: ",A19,"..")',out_dyn_File
      print*, "========================================="
      print '(A,1X,A)', "Reading:", laiv_file
      print*, "========================================="
      call create_dynamic_file(out_dyn_file,idxs,nlai,run_BDSNP)
      
      !print*,"NLAI:",nvars
      print '(A,1X,I0)', "Number of LAI records is ", nlai
      do k=1,nlai
          write(kk,'(I0.2)') k
          call interpolate_area(laiv_file,"lai"//kk  ,ide,jde, laiv(:,:,k))
      end do
      laiv = laiv*lai_scale_factor
      where ( laiv < 0.0 )
              laiv=0.0
      endwhere
      
      call write_3d_var(out_dyn_file,"LAI",laiv,idxs,nlai)
      
      if (run_BDSNP) then
         allocate( ndep(ide,jde,12  ))  
         allocate(nfert(ide,jde,365 ))  
         !NDEP:
         do k=1,12
             write(kk,'(I0.2)') k
             !NDEP(:,:,k)  = interpolate(p,g,nitro_file, varname="nitro"//kk, method="bilinear")
             call interpolate_area(nitro_file,"nitro"//kk  ,ide,jde, ndep(:,:,k))
         enddo
         where (ndep < 0.0 )
            ndep=0.0
         endwhere
         call write_3d_var(out_dyn_file,"NDEP" ,ndep,idxs,12)

         !----
         !NFERT:
         do k=1,365
             write(kkk,'(I0.3)') k
             call interpolate_area(fert_file,"fert"//kkk  ,ide,jde, nfert(:,:,k))
             !NFERT(:,:,k)  = interpolate(p,g,fert_file, varname="fert"//kkk, method="bilinear")
         enddo
         where (NFERT< 0.0 )
            NFERT=0.0
         endwhere
         call write_3d_var(out_dyn_file,"NFERT" ,nfert,idxs,365)
         deallocate( ndep)
         deallocate(nfert)
      end if
      stop
     end subroutine prep_dynamic_data
   !===================================================
   !Other interpolation
   !===================================================
   subroutine interpolate_area(inpfile, varname, ide, jde,data_out )
       use netcdf
       character(len=*), intent(in) :: inpfile,varname
       integer, intent(in) :: ide, jde
       real, intent(inout) :: data_out(:,:)
   
       !local var
       integer :: n
       integer :: ncid
       integer :: nlon_megan, nlat_megan
       logical :: new_grid,flip_flag
       integer :: missing_value
       real, allocatable :: megan_lons(:), megan_lats(:), megan_lats_orig(:)
       character(len=256) :: message
   
   
       !==============================================================================================
       message = 'Opening MEGAN file: '//trim(inpfile)
       call handle_ncerr(nf90_open(trim(inpfile), NF90_NOWRITE, ncid), message)
       call handle_ncerr(nf90_inq_dimid(ncid, 'lon', dimid), 'Getting lon dimension')
       call handle_ncerr(nf90_inquire_dimension(ncid, dimid, len=nlon_megan), 'Inquiring lon dimension')
       call handle_ncerr(nf90_inq_dimid(ncid, 'lat', dimid), 'Getting lat dimension')
       call handle_ncerr(nf90_inquire_dimension(ncid, dimid, len=nlat_megan), 'Inquiring lat dimension')
   
       allocate(megan_lons(nlon_megan), megan_lats(nlat_megan), megan_lats_orig(nlat_megan), stat=ierr)
       call handle_ncerr(nf90_inq_varid(ncid, 'lon', varid), 'Getting lon variable')
       call handle_ncerr(nf90_get_var(ncid, varid, megan_lons), 'Reading lon data')
       call handle_ncerr(nf90_inq_varid(ncid, 'lat', varid), 'Getting lat variable')
       call handle_ncerr(nf90_get_var(ncid, varid, megan_lats_orig), 'Reading lat data')
       !==============================================================================================
   
       ! Flip latitude if decreasing
       flip_flag = .false.
       if (megan_lats_orig(2) < megan_lats_orig(1)) then
         flip_flag = .true.
         megan_lats = megan_lats_orig(nlat_megan:1:-1)
       else
         megan_lats = megan_lats_orig
       endif
       deallocate(megan_lats_orig)
   
       !==============================================================================================
       !==============================================================================================
       allocate(xedge_megan(nlon_megan+1), yedge_megan(nlat_megan+1), stat=ierr)
       xedge_megan(2:nlon_megan) = 0.5_8 * (megan_lons(1:nlon_megan-1) + megan_lons(2:nlon_megan))
       xedge_megan(1) = megan_lons(1) - 0.5_8 * (megan_lons(2) - megan_lons(1))
       xedge_megan(nlon_megan+1) = megan_lons(nlon_megan) + 0.5_8 * (megan_lons(nlon_megan) - megan_lons(nlon_megan-1))
   
       yedge_megan(2:nlat_megan) = 0.5_8 * (megan_lats(1:nlat_megan-1) + megan_lats(2:nlat_megan))
       yedge_megan(1) = megan_lats(1) - 0.5_8 * (megan_lats(2) - megan_lats(1))
       yedge_megan(nlat_megan+1) = megan_lats(nlat_megan) + 0.5_8 * (megan_lats(nlat_megan) - megan_lats(nlat_megan-1))
       !==============================================================================================
       !==============================================================================================
   
       new_grid = .true.
       do n = 1,grid_cnt
         if( grid_specs(n)%nlons /= nlon_megan .or. grid_specs(n)%nlats /= nlat_megan ) then
           cycle
         endif
         if( any( grid_specs(n)%lon(:) /= megan_lons(:) ) ) then
           cycle
         endif
         if( any( grid_specs(n)%lat(:) /= megan_lats(:) ) ) then
           cycle
         endif
         grid_ndx = n
         new_grid = .false.
         exit
       end do
   
       if(new_grid)then
         !load the new grid
         print '(A)', "This is a new grid"
         grid_ndx = grid_cnt + 1
         grid_cnt = grid_ndx
         allocate(grid_specs(grid_ndx)%lon(nlon_megan), stat=ierr)
         allocate(grid_specs(grid_ndx)%lat(nlat_megan), stat=ierr)
         grid_specs(grid_ndx)%lon = megan_lons
         grid_specs(grid_ndx)%lat = megan_lats
         grid_specs(grid_ndx)%nlons = nlon_megan
         grid_specs(grid_ndx)%nlats = nlat_megan
   
         allocate(grid_specs(grid_ndx)%model_area_type(ide, jde), stat=astat)
         grid_specs(grid_ndx)%model_area_type(:,:)%has_data = .false.
         grid_specs(grid_ndx)%model_area_type(:,:)%active_dcell_cnt = 0
         grid_specs(grid_ndx)%model_area_type(:,:)%total_dcell_cnt = 0
         grid_specs(grid_ndx)%model_area_type(:,:)%interior_dcell_cnt = 0
         grid_specs(grid_ndx)%model_area_type(:,:)%partial_dcell_cnt = 0 
       else
         print '(A)', "This is an old grid"
         print '(A,1X,I0)', "Using grid #",grid_cnt
       end if
   
       !=====================================Find missing value=======================================
       print*, "========================================="
       print '(A,1X,A)', "Area Conserving interpolation for ", varname
       print*, "========================================="
       call handle_ncerr(nf90_inq_varid(ncid, varname, varid), 'Getting var variable')
       call handle_ncerr(nf90_get_att(ncid, varid, "missing_value", missing_value),"Error reading missing_value")
       !==============================================================================================
   
       call area_interp( xedge_megan, yedge_megan, nlon_megan, nlat_megan, int(missing_value,2), &
                         data_out, ncid, varname, grid_ndx, new_grid, flip_flag)
   
       call handle_ncerr(nf90_close(ncid), 'Closing Ecotype file')
   
       deallocate(xedge_megan)
       deallocate(yedge_megan)
       deallocate(megan_lons)
       deallocate(megan_lats)
   end subroutine interpolate_area
   !===================================================
   !Ecotype interpolation
   !===================================================
   subroutine interpolate_ecotype(inpfile, varname, ide, jde, ecotypeid,ecotypefrac)
       use netcdf
       character(len=*), intent(in) :: inpfile,varname
       integer, intent(in) :: ide, jde
       integer,  intent(inout) :: ecotypeid(:,:,:)
       real,  intent(inout) :: ecotypefrac(:,:,:)
   
       !local var
       integer :: ncid
       integer :: n
       integer :: nlon_megan, nlat_megan
       logical :: new_grid,flip_flag
       integer :: missing_value
       real, allocatable :: megan_lons(:), megan_lats(:), megan_lats_orig(:)
       real, allocatable :: out_data(:,:,:,:)
       character(len=256) :: message
    
       !==============================================================================================
       message = 'Opening MEGAN file: '//trim(inpfile)
       call handle_ncerr(nf90_open(trim(inpfile), NF90_NOWRITE, ncid), message)
       call handle_ncerr(nf90_inq_dimid(ncid, 'lon', dimid), 'Getting lon dimension')
       call handle_ncerr(nf90_inquire_dimension(ncid, dimid, len=nlon_megan), 'Inquiring lon dimension')
       call handle_ncerr(nf90_inq_dimid(ncid, 'lat', dimid), 'Getting lat dimension')
       call handle_ncerr(nf90_inquire_dimension(ncid, dimid, len=nlat_megan), 'Inquiring lat dimension')
   
       allocate(megan_lons(nlon_megan), megan_lats(nlat_megan), megan_lats_orig(nlat_megan), stat=ierr)
       call handle_ncerr(nf90_inq_varid(ncid, 'lon', varid), 'Getting lon variable')
       call handle_ncerr(nf90_get_var(ncid, varid, megan_lons), 'Reading lon data')
   
       call handle_ncerr(nf90_inq_varid(ncid, 'lat', varid), 'Getting lat variable')
       call handle_ncerr(nf90_get_var(ncid, varid, megan_lats_orig), 'Reading lat data')
       !==============================================================================================
   
       ! Flip latitude if decreasing
       flip_flag = .false.
       if (megan_lats_orig(2) < megan_lats_orig(1)) then
         flip_flag = .true.
         megan_lats = megan_lats_orig(nlat_megan:1:-1)
       else
         megan_lats = megan_lats_orig
       endif
       deallocate(megan_lats_orig)
   
       !==============================================================================================
       !=========================================find the grid edge===================================
       !==============================================================================================
       allocate(xedge_megan(nlon_megan+1), yedge_megan(nlat_megan+1), stat=ierr)
       xedge_megan(2:nlon_megan) = 0.5_8 * (megan_lons(1:nlon_megan-1) + megan_lons(2:nlon_megan))
       xedge_megan(1) = megan_lons(1) - 0.5_8 * (megan_lons(2) - megan_lons(1))
       xedge_megan(nlon_megan+1) = megan_lons(nlon_megan) + 0.5_8 * (megan_lons(nlon_megan) - megan_lons(nlon_megan-1))
   
       yedge_megan(2:nlat_megan) = 0.5_8 * (megan_lats(1:nlat_megan-1) + megan_lats(2:nlat_megan))
       yedge_megan(1) = megan_lats(1) - 0.5_8 * (megan_lats(2) - megan_lats(1))
       yedge_megan(nlat_megan+1) = megan_lats(nlat_megan) + 0.5_8 * (megan_lats(nlat_megan) - megan_lats(nlat_megan-1))
       !==============================================================================================
       !==============================================================================================
   
       new_grid = .true.
       do n = 1,grid_cnt
         if( grid_specs(n)%nlons /= nlon_megan .or. grid_specs(n)%nlats /= nlat_megan ) then
           cycle
         endif
         if( any( grid_specs(n)%lon(:) /= megan_lons(:) ) ) then
           cycle
         endif
         if( any( grid_specs(n)%lat(:) /= megan_lats(:) ) ) then
           cycle
         endif
         grid_ndx = n
         new_grid = .false.
         exit
       end do
   
       if(new_grid)then
         !load the new grid
         grid_ndx = grid_cnt + 1
         grid_cnt = grid_ndx
         allocate(grid_specs(grid_ndx)%lon(nlon_megan), stat=ierr)
         allocate(grid_specs(grid_ndx)%lat(nlat_megan), stat=ierr)
         grid_specs(grid_ndx)%lon = megan_lons
         grid_specs(grid_ndx)%lat = megan_lats
   
         allocate(grid_specs(grid_ndx)%model_area_type(ide, jde), stat=astat)
         grid_specs(grid_ndx)%model_area_type(:,:)%has_data = .false.
         grid_specs(grid_ndx)%model_area_type(:,:)%active_dcell_cnt = 0
         grid_specs(grid_ndx)%model_area_type(:,:)%total_dcell_cnt = 0
         grid_specs(grid_ndx)%model_area_type(:,:)%interior_dcell_cnt = 0
         grid_specs(grid_ndx)%model_area_type(:,:)%partial_dcell_cnt = 0 
       end if
       
       !==============================================================================================
       !=========================================get the missing value================================
       !==============================================================================================
   
       call handle_ncerr(nf90_inq_varid(ncid, varname, varid), 'Getting var variable')
       call handle_ncerr(nf90_get_att(ncid, varid, "missing_value", missing_value),"Error reading missing_value")
       !==============================================================================================
   
       allocate(out_data(ide,jde,mxetype,2))
       call discrete_frac(xedge_megan, yedge_megan, nlon_megan, nlat_megan, int(missing_value, 2), &
                          out_data, ncid, varname, grid_ndx, new_grid, flip_flag)
       !ecotypeid = int(ecotypefrac)
       ecotypeid(:,:,:)   = int(out_data(:,:,:,1))
       ecotypefrac(:,:,:) = out_data(:,:,:,2)
       !==============================================================================================
       call handle_ncerr(nf90_close(ncid), 'Closing Ecotype file')
   
       deallocate(xedge_megan)
       deallocate(yedge_megan)
       deallocate(megan_lons)
       deallocate(megan_lats)
       deallocate(out_data)
   end subroutine interpolate_ecotype


  !---------------------------------------------------------------------
  !   read wrf file
  !---------------------------------------------------------------------
  subroutine wrf_file(wrffile, ide, jde, cen_lon, cen_lat, stand_lon, truelat1, truelat2, dx)
  
     use netcdf
     character(len=*), intent(in) :: wrffile
     integer, intent(out) :: ide, jde
     real, intent(out) :: cen_lon, cen_lat, stand_lon, truelat1, truelat2, dx
     character(len=80) :: message
     integer :: ncid 
   
 !---------------------------------------------------------------------
 !   open wrf input file
 !---------------------------------------------------------------------
    message = 'wrf_file: Failed to open ' // trim(wrffile)
    call handle_ncerr( nf90_open( trim(wrffile), nf90_noclobber, ncid ), message )
!---------------------------------------------------------------------
!   get wrf dimesions
!---------------------------------------------------------------------
    call handle_ncerr( nf90_inq_dimid( ncid, 'west_east', dimid ), "Failed to get lon dim. id" )
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=ide ), "Failed to get lon dim." )
    call handle_ncerr( nf90_inq_dimid( ncid, 'south_north', dimid ), "Failed to get lat dim. id" )
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=jde ), "Failed to get lat dim." )
!---------------------------------------------------------------------
!   get wrf map projection variables
!---------------------------------------------------------------------
    call handle_ncerr( nf90_get_att( ncid, nf90_global, 'MAP_PROJ', map_proj ), "Failed to get MAP_PROJ" )
    if( map_proj /= PS ) then
       write(*,*) 'wrf_file: MAP_PROJ is not polar stereographic'
    else
       write(*,*) 'wrf_file: MAP_PROJ is polar stereographic'
    endif
    call handle_ncerr( nf90_get_att( ncid, nf90_global, 'CEN_LON', cen_lon ), "Failed to get CEN_LON" )
    write(*,*) 'wrf_file: CEN_LON = ',cen_lon
    call handle_ncerr( nf90_get_att( ncid, nf90_global, 'CEN_LAT', cen_lat ), "Failed to get CEN_LAT" )
    write(*,*) 'wrf_file: CEN_LAT = ',cen_lat
    call handle_ncerr( nf90_get_att( ncid, nf90_global, 'STAND_LON', stand_lon ), "Failed to get STAND_LON" )
    write(*,*) 'wrf_file: STAND_LON = ',stand_lon
    call handle_ncerr( nf90_get_att( ncid, nf90_global, 'TRUELAT1', truelat1 ), "Failed to get TRUELAT1" )
    write(*,*) 'wrf_file: TRUELAT1 = ',truelat1
    call handle_ncerr( nf90_get_att( ncid, nf90_global, 'TRUELAT2', truelat2 ), "Failed to get TRUELAT2" )
    write(*,*) 'wrf_file: TRUELAT2 = ',truelat2
    call handle_ncerr( nf90_get_att( ncid, nf90_global, 'DX', dx ), "Failed to get DEFailed to get DEXX" )
    write(*,*) 'wrf_file: DX = ',dx
 
!---------------------------------------------------------------------
!   initialize map projection
!---------------------------------------------------------------------
      call proj_init( map_proj, cen_lon, cen_lat, truelat1, truelat2, &
                      stand_lon, dx, ide, jde )
   
      ids = 1
      jds = 1
   
   
      message = 'wrf_file: Failed to close ' // trim(wrffile)
      call handle_ncerr( nf90_close( ncid ), message )       
   
      !allocate( ecotypeid(ide,jde,mxetype),stat=astat ) 
      !if( astat /= 0 ) then
      !  write(*,*) 'wrf_file: failed to allocate ecotypeid; error = ',astat
      !  stop 'allocate failed'
      !endif
      !allocate( ecotypefrac(ide,jde,mxetype),stat=astat ) 
      !if( astat /= 0 ) then
      !  write(*,*) 'wrf_file: failed to allocate ecotypefrac; error = ',astat
      !  stop 'allocate failed'
      !endif
   
   end subroutine wrf_file
!=======================================================================
   subroutine create_static_file(outfile,idxs,run_BDSNP)
     use netcdf
     use area_mapper_grw, only: lon, lat
     implicit none
     character(len=*),intent(in) :: outfile  
     integer,         intent(in) :: idxs(:)
     logical,         intent(in) :: run_BDSNP
      
     !local var
     integer :: ncid,var_id
     integer :: x_dim_id, y_dim_id, cty_dim_id, ef_dim_id, ldf_dim_id
     integer :: xt,xe,yt,ye,nx,ny
   
     xt = idxs(1)
     yt = idxs(2)
     nx = idxs(3)
     ny = idxs(4)
     xe = idxs(1) + idxs(3) - 1
     ye = idxs(2) + idxs(4) - 1
   
     !Create File and define dimensions and variables:
     call check(nf90_create(outfile, IOR(NF90_CLOBBER, NF90_NETCDF4), ncid))
        call check(nf90_def_dim(ncid, "west_east"      , nx    , x_dim_id   ))
        call check(nf90_def_dim(ncid, "south_north"    , ny    , y_dim_id   ))
        call check(nf90_def_dim(ncid, "cantype"        , NCANTYPE+1,cty_dim_id ))
        call check(nf90_def_dim(ncid, "ef_dim"         , NEFS    ,ef_dim_id  ))
        call check(nf90_def_dim(ncid, "ldf_dim"        , NLDFS   ,ldf_dim_id  ))
        !Define variables:    
        ! Coordinates:
        call check(nf90_def_var(ncid, "lon"    , NF90_FLOAT, [x_dim_id,y_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id, "units", "degrees_east"))
        call check(nf90_put_att(ncid, var_id, "long_name", "longitude"))
        call check(nf90_def_var(ncid, "lat"    , NF90_FLOAT, [x_dim_id,y_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id, "units", "degrees_north"))
        call check(nf90_put_att(ncid, var_id, "long_name", "latitude"))
        ! AREA:
        call check(nf90_def_var(ncid, "cell_area", NF90_FLOAT, [x_dim_id,y_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "cell_area"                     ))
        call check(nf90_put_att(ncid, var_id,"units"    , "m2"                            ))
        call check(nf90_put_att(ncid, var_id,"var_desc" , "horizontal area of a gridcell" ))
   
        !ECOTYPE
        !call check(nf90_def_var(ncid, "ETY" , NF90_INT, [x_dim_id,y_dim_id],var_id)) !debug
        ! CTF:
        call check(nf90_def_var(ncid, "CTF" , NF90_FLOAT, [x_dim_id,y_dim_id,cty_dim_id],var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "CANOPY_TYPE_FRACTION"               ))
        call check(nf90_put_att(ncid, var_id,"units"    , "fraction"                                  ))
        call check(nf90_put_att(ncid, var_id,"Description",&
       "Canopy Type Fraction:1.Needleleaf Trees;2. Tropical Trees; 3.Broadleaf Tree;4. Shrub;5. Herb;6.  Crop;7. Total Tree Fraction" ))
        ! EFs:
        ! EFs for all
        call check(nf90_def_var(ncid, "EFS" , NF90_FLOAT, [x_dim_id,y_dim_id,ef_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "INTEGRATED EMISSION_FACTOR"                    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "nanomol m-2 s-1"                    )) 
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Emission Factors ISOP,MBO,MT_PINE,MT_ACYC,MT_CAMP,MT_SABI,MT_AROM,NO,SQT_HR,SQT_LR,MEOH,ACTO,ETOH,ACID,LVOC,OXPROD,STRESS,OTHER,CO" ))
        ! EFs for tree growthform
        call check(nf90_def_var(ncid, "EFS_TREE" , NF90_FLOAT, [x_dim_id,y_dim_id,ef_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "EMISSION_FACTOR FOR TREE"                    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "nanomol m-2 s-1"                    )) 
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Emission Factors ISOP,MBO,MT_PINE,MT_ACYC,MT_CAMP,MT_SABI,MT_AROM,NO,SQT_HR,SQT_LR,MEOH,ACTO,ETOH,ACID,LVOC,OXPROD,STRESS,OTHER,CO" ))
        
        ! EFs for shrub growthform
        call check(nf90_def_var(ncid, "EFS_SHRUB" , NF90_FLOAT, [x_dim_id,y_dim_id,ef_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "EMISSION_FACTOR FOR SHRUB"                    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "nanomol m-2 s-1"                    )) 
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Emission Factors ISOP,MBO,MT_PINE,MT_ACYC,MT_CAMP,MT_SABI,MT_AROM,NO,SQT_HR,SQT_LR,MEOH,ACTO,ETOH,ACID,LVOC,OXPROD,STRESS,OTHER,CO" ))
        
        ! EFs for herb growthform
        call check(nf90_def_var(ncid, "EFS_HERB" , NF90_FLOAT, [x_dim_id,y_dim_id,ef_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "EMISSION_FACTOR FOR HERB"                    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "nanomol m-2 s-1"                    )) 
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Emission Factors ISOP,MBO,MT_PINE,MT_ACYC,MT_CAMP,MT_SABI,MT_AROM,NO,SQT_HR,SQT_LR,MEOH,ACTO,ETOH,ACID,LVOC,OXPROD,STRESS,OTHER,CO" ))
        
        ! EFs for crop growthform
        call check(nf90_def_var(ncid, "EFS_CROP" , NF90_FLOAT, [x_dim_id,y_dim_id,ef_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "EMISSION_FACTOR FOR CROP"                    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "nanomol m-2 s-1"                    )) 
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Emission Factors ISOP,MBO,MT_PINE,MT_ACYC,MT_CAMP,MT_SABI,MT_AROM,NO,SQT_HR,SQT_LR,MEOH,ACTO,ETOH,ACID,LVOC,OXPROD,STRESS,OTHER,CO" ))
        ! LDF:
        ! LDF for all
        call check(nf90_def_var(ncid, "LDF" , NF90_FLOAT, [x_dim_id,y_dim_id,ldf_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "INTEGRATED LIGHT DEPENDENT EMISSION_FACTOR"    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "fraction"                    ))
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Ligth Dependent Emissions Factors: LDF01,...LDF04" ))
        ! LDF for trees
        call check(nf90_def_var(ncid, "LDF_TREE" , NF90_FLOAT, [x_dim_id,y_dim_id,ldf_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "LIGHT DEPENDENT EMISSION_FACTOR FOR TREE"    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "fraction"                    ))
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Ligth Dependent Emissions Factors: LDF01,...LDF04" ))
        ! LDF for shrub
        call check(nf90_def_var(ncid, "LDF_SHRUB" , NF90_FLOAT, [x_dim_id,y_dim_id,ldf_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "LIGHT DEPENDENT EMISSION_FACTOR FOR SHRUB"    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "fraction"                    ))
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Ligth Dependent Emissions Factors: LDF01,...LDF04" ))
        ! LDF for herb
        call check(nf90_def_var(ncid, "LDF_HERB" , NF90_FLOAT, [x_dim_id,y_dim_id,ldf_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "LIGHT DEPENDENT EMISSION_FACTOR FOR HERB"    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "fraction"                    ))
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Ligth Dependent Emissions Factors: LDF01,...LDF04" ))
        ! LDF for crop
        call check(nf90_def_var(ncid, "LDF_CROP" , NF90_FLOAT, [x_dim_id,y_dim_id,ldf_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "LIGHT DEPENDENT EMISSION_FACTOR FOR CROP"    ))
        call check(nf90_put_att(ncid, var_id,"units"    , "fraction"                    ))
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Ligth Dependent Emissions Factors: LDF01,...LDF04" ))
        if (run_BDSNP) then
           print*,"Building BDSNP_ARID, BDSNP_NONARID & BDSNP_LANDTYPE ..."
           ! LANDTYPE, ARID, NONARID (BDSNP)
           call check(nf90_def_var(ncid, "arid", NF90_INT  , [x_dim_id,y_dim_id],var_id))
           call check(nf90_put_att(ncid, var_id,"long_name", "arid"                   ))
           call check(nf90_put_att(ncid, var_id,"units"    , "1 or 0"                 ))
           call check(nf90_put_att(ncid, var_id,"var_desc" , "Arid soil mask"         ))
   
           call check(nf90_def_var(ncid,"landtype",NF90_INT, [x_dim_id,y_dim_id],var_id))
           call check(nf90_put_att(ncid, var_id,"long_name", "landtype"                ))
           call check(nf90_put_att(ncid, var_id,"units"    , "nondimension"            ))
           call check(nf90_put_att(ncid, var_id,"var_desc" , "Soil type calssification"))
        endif
        !Global Attributes
        call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "MEGAN input file"   ))
        call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
        call check(nf90_enddef(ncid))
            !Get and write variables:
        call check(nf90_open(outfile, nf90_write, ncid ))
            !Coordinates:
            call check(nf90_inq_varid(ncid,"lon" ,var_id))
            call check(nf90_put_var(ncid, var_id, lon(xt:xe,yt:ye)))
            call check(nf90_inq_varid(ncid,"lat" ,var_id))
            call check(nf90_put_var(ncid, var_id, lat(xt:xe,yt:ye)))
        call check(nf90_close(ncid))
   end subroutine create_static_file
!==============================================================
!==============================================================
   subroutine create_dynamic_file(outfile,idxs,nlai,run_BDSNP)
     use netcdf
     use area_mapper_grw, only: lon, lat
     implicit none
     character(len=*),intent(in) :: outfile  
     integer,         intent(in) :: idxs(:)
     integer,         intent(in) :: nlai
     logical,         intent(in) :: run_BDSNP
    
     !local var
     integer :: ncid,var_id
     integer :: x_dim_id, y_dim_id, nlai_dim_id, month_dim_id, day_dim_id
     integer :: xt,xe,yt,ye,nx,ny

     xt = idxs(1)
     yt = idxs(2)
     nx = idxs(3)
     ny = idxs(4)
     xe = idxs(1) + idxs(3) - 1
     ye = idxs(2) + idxs(4) - 1
     !Create File and define dimensions and variables:
     call check(nf90_create(outfile, IOR(NF90_CLOBBER, NF90_NETCDF4), ncid))
        call check(nf90_def_dim(ncid, "west_east"      , nx    , x_dim_id   ))
        call check(nf90_def_dim(ncid, "south_north"    , ny    , y_dim_id   ))
        call check(nf90_def_dim(ncid, "time"           , nlai  , nlai_dim_id ))
        if (run_BDSNP) then
          call check(nf90_def_dim(ncid, "month"          , 12    , month_dim_id))
          call check(nf90_def_dim(ncid, "day"            , 365   , day_dim_id  ))
        end if
        !Define variables:    
        ! Coordinates:
        call check(nf90_def_var(ncid, "lon"    , NF90_FLOAT, [x_dim_id,y_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id, "units", "degrees_east"))
        call check(nf90_put_att(ncid, var_id, "long_name", "longitude"))
        call check(nf90_def_var(ncid, "lat"    , NF90_FLOAT, [x_dim_id,y_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id, "units", "degrees_north"))
        call check(nf90_put_att(ncid, var_id, "long_name", "latitude"))
        
        !LAI
        call check(nf90_def_var(ncid, "LAI" , NF90_FLOAT, [x_dim_id,y_dim_id,nlai_dim_id],var_id))
        call check(nf90_put_att(ncid, var_id,"long_name", "LAI"            ))
        call check(nf90_put_att(ncid, var_id,"units"    , "m2 m-2"               ))
        call check(nf90_put_att(ncid, var_id,"var_desc" , "Leaf Area Index" ))
        if (run_BDSNP) then
           !NDEP
           call check(nf90_def_var(ncid, "NDEP", NF90_FLOAT, [x_dim_id,y_dim_id,month_dim_id],var_id))
           call check(nf90_put_att(ncid, var_id,"long_name", "NDEP"            ))
           call check(nf90_put_att(ncid, var_id,"units"    , "kg/m2/s"             ))
           call check(nf90_put_att(ncid, var_id,"var_desc" , "N Deposition   " ))
           !NFERT
           call check(nf90_def_var(ncid,"NFERT", NF90_FLOAT, [x_dim_id,y_dim_id, day_dim_id],var_id))
           call check(nf90_put_att(ncid, var_id,"long_name", "NFERT"           ))
           call check(nf90_put_att(ncid, var_id,"units"    , "mg/m3"             ))
           call check(nf90_put_att(ncid, var_id,"var_desc" , "N Fertilization" ))
        endif
        !Global Attributes
        call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "MEGAN input file"   ))
        call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
        call check(nf90_enddef(ncid))
            !Get and write variables:
        call check(nf90_open(outfile, nf90_write, ncid ))
            !Coordinates:
            call check(nf90_inq_varid(ncid,"lon" ,var_id))
            call check(nf90_put_var(ncid, var_id, lon(xt:xe,yt:ye)))
            call check(nf90_inq_varid(ncid,"lat" ,var_id))
            call check(nf90_put_var(ncid, var_id, lat(xt:xe,yt:ye)))
        call check(nf90_close(ncid))
   end subroutine create_dynamic_file

!==============================================================
!==============================================================
 
   subroutine compute_ef_grid(ecotypeid, ecotypefrac, csv_filename, EF_grid)
     implicit none
   
     ! Input arguments
     integer, intent(in) :: ecotypeid(:,:,:)
     real, intent(in)    :: ecotypefrac(:,:,:)
     character(len=*), intent(in) :: csv_filename
   
     ! Output
     real, intent(inout) :: EF_grid(:,:,:,:)
     !size(ecotypeid,1), size(ecotypeid,2), ncat)
   
     ! Local variables
     integer :: lat_size, lon_size
     integer :: i, j, k, c, id, v
     real    :: frac
     logical :: id_to_veg(0:max_id, nveg)
     real    :: ef_table(nveg, 0:max_id, ncat)
   
     lat_size = size(ecotypeid, 1)
     lon_size = size(ecotypeid, 2)
     ! Read lookup table
     call read_lookup_table(csv_filename, id_to_veg, ef_table, max_id, ncat, nveg)
   
     ! Initialize output
     EF_grid = 0.0
  
     print*,maxval(ecotypeid)
     ! Compute EF values
     do i = 1, lat_size
       do j = 1, lon_size
         do k = 1, mxetype
           id = ecotypeid(i, j, k)
           frac = ecotypefrac(i, j, k)
           if (id >= 0 .and. id <= max_id .and. frac> 0.0) then
             do v = 1, nveg
               if (id_to_veg(id, v)) then
                 do c = 1, ncat
                   EF_grid(i, j, c, v) = EF_grid(i, j, c, v) + frac * ef_table(v, id, c)
                 end do
               end if
             end do !nveg
           end if 
         end do !mxetype
       end do !lon_size
     end do !lat_size
   
   end subroutine compute_ef_grid

   subroutine read_lookup_table(filename, id_to_veg, ef_table, max_id, ncat, nveg)
     implicit none
     character(len=*), intent(in) :: filename
     integer, intent(in)          :: max_id, ncat, nveg
     logical, intent(out)         :: id_to_veg(0:max_id, nveg)
     real, intent(out)            :: ef_table(nveg, 0:max_id, ncat)

     character(len=256) :: line
     character(len=20)  :: veg
     character(len=20), dimension(nveg) :: veg_names
     real :: tmp_ef(ncat)
     integer :: id, ios, linenum, i, c, v
     integer :: unit

     veg_names = (/ 'Crop', 'Herb', 'Shrub', 'Tree' /)

     id_to_veg = .false.
     ef_table = 0.0

     !print '("Reading loo up table)'
     print*,"=============================="
     print '(A)', "Reading look up table"
     print*,"=============================="
     open(newunit=unit, file=filename, status='old', action='read')
     linenum = 0
     do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       linenum = linenum + 1
       read(line, *) veg, id, (tmp_ef(c), c=1, ncat)

       ! Map vegetation name to index
       v = 0
       do i = 1, nveg
         if (veg == veg_names(i)) then
           v = i
           exit
         end if
       end do

       if (v == 0) then
         print *, 'Warning: unknown vegetation type on line', linenum, ':', veg
       else if (id >= 0 .and. id <= max_id) then
         id_to_veg(id, v) = .true.
         do c = 1, ncat
           ef_table(v, id, c) = tmp_ef(c)
         end do
       end if
     end do
     close(unit)
     !print*,ef_table(1,1,:)
   end subroutine read_lookup_table
   
   subroutine write_2d_var(outfile,varname,data_out,idxs)
     use netcdf
     implicit none
     character(len=*),intent(in) :: outfile  
     character(len=*),intent(in) :: varname 
     integer,         intent(in) :: idxs(:)
     real,            intent(in) :: data_out(:,:)  
     real,            allocatable:: data_block(:,:)  
     !local var
     integer :: ncid,varidinp
     integer :: xt,xe,yt,ye,nx,ny
    
     
     xt = idxs(1)
     yt = idxs(2)
     nx = idxs(3)
     ny = idxs(4)
     xe = idxs(1) + idxs(3) - 1
     ye = idxs(2) + idxs(4) - 1
    
     allocate(data_block(nx,ny))
     data_block = data_out(xt:xe,yt:ye)
     !integer :: xt,xe,yt,ye,nx,ny
     call check(nf90_open(outfile, nf90_write, ncid ))
     call check(nf90_inq_varid(ncid,varname,varidinp))
     call check(nf90_put_var(ncid, varidinp, data_block ))
     call check(nf90_close(ncid))
     deallocate(data_block)
   end subroutine write_2d_var
   
   subroutine write_3d_var(outfile,varname,data_out,idxs,nt)
     use netcdf
     implicit none
     character(len=*),intent(in) :: outfile  
     character(len=*),intent(in) :: varname 
     integer,         intent(in) :: idxs(:)
     integer,         intent(in) :: nt
     real,            intent(in) :: data_out(:,:,:)  
     real,            allocatable:: data_block(:,:,:)  
     !local var
     integer :: ncid,varidinp
     integer :: xt,xe,yt,ye,nx,ny
    
     
     xt = idxs(1)
     yt = idxs(2)
     nx = idxs(3)
     ny = idxs(4)
     xe = idxs(1) + idxs(3) - 1
     ye = idxs(2) + idxs(4) - 1

     allocate(data_block(nx,ny,nt))
     data_block = data_out(xt:xe,yt:ye,:)
     !integer :: xt,xe,yt,ye,nx,ny
     call check(nf90_open(outfile, nf90_write, ncid ))
     call check(nf90_inq_varid(ncid,varname,varidinp))
     call check(nf90_put_var(ncid, varidinp, data_block ))
     call check(nf90_close(ncid))
     deallocate(data_block)
   end subroutine write_3d_var

   subroutine handle_ncerr( ret, mes )
   !---------------------------------------------------------------------
   !       ... netcdf error handling routine
   !---------------------------------------------------------------------
      integer, intent(in) :: ret
      character(len=*), intent(in) :: mes
   
      if( ret /= nf90_noerr ) then
         write(*,*) nf90_strerror( ret )
         stop 'netcdf error'
      endif
   
   end subroutine handle_ncerr
   subroutine check(status)
     integer, intent(in) :: status
     if (status /= nf90_noerr) then
       write(*,*) nf90_strerror(status)
       stop 'netcdf error'
     end if
   end subroutine check
end module 
