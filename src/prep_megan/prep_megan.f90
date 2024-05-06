!---------------------------------------------------------------
! programed by: Ramiro A. Espada. April 2023.
! Based on Prep_code &  MEGEFP32 (UCI-BAI-MEGAN)
!---------------------------------------------------------------
!program prepmegan
module prep_megan

  use netcdf

  implicit none

  private
  public prep


 INTEGER, PARAMETER :: ascii = selected_char_KIND ("ascii")
 INTEGER, PARAMETER :: ucs4  = selected_char_KIND ('ISO_10646')

 !Parameters:
 real, parameter    :: R_EARTH = 6370000.
 real, parameter    :: PI = 3.141592653589793
 real, parameter    :: RAD2DEG = 180./PI, DEG2RAD = PI/180.
 
 integer, parameter :: ncantype=6, nefs=19, nldfs=4  ! # of canopy types, # of emission factors (19 EF + 4 LDF)

 !Objects/Strucs:
 type proj_type
    character(16)    :: pName                   ! Projection name
    integer          :: typ                     ! Integer code for projection TYPE (2=lcc, 6=stere, 7=merc)
    real             :: alp,bet,gam,xcent,ycent !proj parameters.
    real             :: p1,p2,p3,p4             !extra parameters to speed up calculation once p%typ is defined.
 end type proj_type

 type grid_type
     character(12)   :: gName                   !grid-name
     integer         :: nx,ny,nz                !number of cells in x-y direction (ncols, nrows, nlevs)
     real            :: dx,dy                   !x-y cell dimension (x_cell, y_cell)
     real            :: xmin,ymin,xmax,ymax,xc,yc
     real            :: lonmin,latmin,lonmax,latmax
 end type grid_type
 
 integer :: iostat,i,j,k
contains

subroutine prep(griddesc_file, gridname,                                      &
                ecotypes_file, growtype_file, laiv_file, GtEcoEF_file,        &
                run_BDSNP, nitro_file, fert_file, climate_file, landtype_file)
 implicit none
 character(200), intent(in)   :: griddesc_file,gridname,ecotypes_file,growtype_file,laiv_file,climate_file,fert_file,landtype_file,nitro_file,GtEcoEF_file
 logical       ,intent(in)    :: run_BDSNP!=.true.

 type(proj_type)    :: proj                    !struc that describes projection
 type(grid_type)    :: grid                    !struc that describes regular grid

 real, allocatable, save  :: longitude(:,:), latitude(:,:) !coordinates

 !@namelist/control/griddesc_file,gridname,ecotypes_file,growtype_file,laiv_file,climate_file,landtype_file,nitro_file,fert_file,GtEcoEF_file,run_BDSNP 
 !@!Leo namelist:
 !@read(*,nml=control, iostat=iostat)
 !@if( iostat /= 0 ) then
 !@  write(*,*) 'prepmegan4cmaq: failed to read namelist; error = ',iostat
 !@  stop
 !@end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file,gridname, proj, grid)
   
  !Compute coordinates lat lon for all grid cells                 
  allocate(latitude(grid%nx,grid%ny));allocate(longitude(grid%nx,grid%ny))
  do j=1,grid%ny
  do i=1,grid%nx                                                                                            
     call xy2ll(proj,grid%xmin+grid%dx*i,grid%ymin+grid%dy*j,longitude(i,j),latitude(i,j)) 
  enddo
  enddo
 
print*,"prep static"
  !Static data:
  call prep_static_data(grid,proj,latitude,longitude,growtype_file,ecotypes_file,GtEcoEF_file,climate_file,landtype_file, run_BDSNP)
                                                                       ! `CTF` (*Canopy Type Fractions*):
                                                                       ! `EFs` (*Emission Factors*)     : (~19) VOC family, and Canopy Type (6)
                                                                       ! `LDF` (*Light Dependent EF*)   :  4 VOC families, and Canopy Type (6)
                                                                       ! `arid`     (BDSNP)
                                                                       ! `landtype` (BDSNP)
print*,"prep dynamic"
  !Time/date dependent data:
  call prep_dynamic_data(grid,proj,latitude,longitude,laiv_file,nitro_file,fert_file,run_BDSNP)     ! time dependent variables
                                                                       ! `LAI`     monthly.
                                                                       ! `N_DEP:`  monthly. (BDSNP)
                                                                       ! `N_FERT:` daily.   (BDSNP)
print*, "========================================="
print*, " prep-megan: Completed successfully"
print*, "========================================="
end subroutine


 !----------------------------------
 !  STATIC  DATA:
 !---------------------------------
subroutine prep_static_data(g,p,lat,lon,ctf_file, ecotype_file, GtEcoEF_file, climate_file,landtype_file, run_BDSNP)
  implicit none
  type(grid_type) ,intent(in) :: g
  type(proj_type) ,intent(in) :: p
  character(len=*),intent(in) :: ctf_file, ecotype_file, GtEcoEF_file  !input  files
  character(len=*),intent(in) :: climate_file, landtype_file  !input  files
  character(len=18)           :: outfile='mgn_static_data.nc' !output file
  logical :: run_BDSNP
  !Coordinates
  real, intent(in)  :: lat(:,:),lon(:,:)
  !netcdf indices:
  integer :: ncid,xid,yid,zid,vid,tid,var_id
  integer :: x_dim_id,y_dim_id,cty_dim_id,ef_dim_id,ldf_dim_id,i,j,k
  !CTF: 
  real, allocatable :: CTF(:,:,:)          !CTF buffer
  !EF & LDF:
  integer, allocatable :: ECOTYPE(:,:)      !este puede ser int tamb
  real   , allocatable :: GTYP(:,:,:)       !CTF (crop, tree, grass, shrub)
  integer              :: EcoID             !EcoID read in table file
  character(len=6)     :: GtID              !GtID  read in table file
  real                 :: EF(NEFS+NLDFS)    !EFs   read in table file 
  real, allocatable    :: OUTGRID(:,:,:)    !EF and LDF buffer
  character(len=6)     :: CTF_LIST(4) ! crop, tree, grass, shrub
  
  !real, allocatable :: OUTGRID(:,:,:,:)   !test v3.3
  !character(len=6) :: CTF_LIST(7) ! crop, tree, grass, shrub
  !LAND (arid, non-arid, landtype)
  real, allocatable    :: LANDGRID(:,:,:)   !LAND buffer

 !Create File and define dimensions and variables:
 call check(nf90_create(outFile, NF90_CLOBBER, ncid))
    !Define dimensions:
    call check(nf90_def_dim(ncid, "x"      , g%nx    , x_dim_id   ))
    call check(nf90_def_dim(ncid, "y"      , g%ny    , y_dim_id   ))
    call check(nf90_def_dim(ncid, "cantype", NCANTYPE, cty_dim_id ))
    call check(nf90_def_dim(ncid, "ef_dim" , NEFS    , ef_dim_id  )) !ef01, ef02, ... , ef19, ldf01, ..., ldf04
    call check(nf90_def_dim(ncid, "ldf_dim", NLDFS   ,ldf_dim_id  )) !ef01, ef02, ... , ef19, ldf01, ..., ldf04
    !Define variables:    
    ! Coordinates:
    call check(nf90_def_var(ncid, "lon"    , NF90_FLOAT, [x_dim_id,y_dim_id], var_id))
    call check(nf90_def_var(ncid, "lat"    , NF90_FLOAT, [x_dim_id,y_dim_id], var_id))
    ! CTF:
    call check(nf90_def_var(ncid, "CTF" , NF90_FLOAT, [x_dim_id,y_dim_id,cty_dim_id],var_id))
    call check(nf90_put_att(ncid, var_id,"long_name", "CANOPY_TYPE_FRACTION"               ))
    call check(nf90_put_att(ncid, var_id,"units"    , "1"                                  ))
    call check(nf90_put_att(ncid, var_id,"var_desc" , "Canopy Type Fraction"               ))
    ! EFs:
    call check(nf90_def_var(ncid, "EFS" , NF90_FLOAT, [x_dim_id,y_dim_id,ef_dim_id], var_id))
    call check(nf90_put_att(ncid, var_id,"long_name", "EMISSION_FACTOR"                    ))
    call check(nf90_put_att(ncid, var_id,"units"    , "nanomol m-2 s-1"                    )) !Q: Are this the correct units? Are the same for LDF?
    call check(nf90_put_att(ncid, var_id,"var_desc" , "Emission Factors ISOP,MBO,MT_PINE,MT_ACYC,MT_CAMP,MT_SABI,MT_AROM,NO,SQT_HR,SQT_LR,MEOH,ACTO,ETOH,ACID,LVOC,OXPROD,STRESS,OTHER,CO" ))
    ! LDF:
    call check(nf90_def_var(ncid, "LDF" , NF90_FLOAT, [x_dim_id,y_dim_id,ldf_dim_id], var_id))
    call check(nf90_put_att(ncid, var_id,"long_name", "LIGHT DEPENDENT EMISSION_FACTOR"    ))
    call check(nf90_put_att(ncid, var_id,"units"    , "nanomol m-2 s-1"                    )) !Q: Are this the correct units? Are the same for LDF?
    call check(nf90_put_att(ncid, var_id,"var_desc" , "Ligth Dependent Emissions Factors: LDF01,...LDF04" ))
    if (run_BDSNP) then
    print*,"Building BDSNP_ARID, BDSNP_NONARID & BDSNP_LANDTYPE ..."
    ! LANDTYPE, ARID, NONARID (BDSNP)
    call check(nf90_def_var(ncid, "arid", NF90_INT  , [x_dim_id,y_dim_id],var_id))
    call check(nf90_put_att(ncid, var_id,"long_name", "arid"                   ))
    call check(nf90_put_att(ncid, var_id,"units"    , "1 or 0"                 ))
    call check(nf90_put_att(ncid, var_id,"var_desc" , "Arid soils"             ))

    call check(nf90_def_var(ncid,"landtype",NF90_INT, [x_dim_id,y_dim_id],var_id))
    call check(nf90_put_att(ncid, var_id,"long_name", "landtype"                ))
    call check(nf90_put_att(ncid, var_id,"units"    , "nondimension"            ))
    call check(nf90_put_att(ncid, var_id,"var_desc" , "Canopy Type Fraction"    ))
    endif
    !Global Attributes
    call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "MEGAN input file"   ))
    call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
 call check(nf90_enddef(ncid))
 !End NetCDF define mode

 !Get and write variables:
  call check(nf90_open(outFile, nf90_write, ncid       ))
     !--------
     !Coordinates
     call check(nf90_inq_varid(ncid,"lon" ,var_id)); call check(nf90_put_var(ncid, var_id, lon ) )
     call check(nf90_inq_varid(ncid,"lat" ,var_id)); call check(nf90_put_var(ncid, var_id, lat ) )

     !--------
     ! CTF:
print*,"CTF"
     allocate(CTF(g%nx,g%ny,ncantype+1))  ! allocate(CTF(g%nx,g%ny,nvars))  

     CTF(:,:,1)=interpolate(p,g,inp_file=ctf_file, varname="nl_tree"  , method="bilinear")
     CTF(:,:,2)=interpolate(p,g,inp_file=ctf_file, varname="trop_tree", method="bilinear") !Q: is it just 100% if lat between tropics, else 0% ?
     !CTF(:,:,3) = !boradleaf tree!
     CTF(:,:,4)=interpolate(p,g,inp_file=ctf_file, varname="shrub"    , method="bilinear")
     CTF(:,:,5)=interpolate(p,g,inp_file=ctf_file, varname="grass"    , method="bilinear")
     CTF(:,:,6)=interpolate(p,g,inp_file=ctf_file, varname="crop"     , method="bilinear")
     CTF(:,:,7)=interpolate(p,g,inp_file=ctf_file, varname="tree"     , method="bilinear")

     !needleleaf tree
     CTF(:,:,1)=(    CTF(:,:,1)/100.0) * CTF(:,:,7) * (1.0-CTF(:,:,2)/100.0)
     !boradleaf tree
     CTF(:,:,3)=CTF(:,:,7) * (1.0-CTF(:,:,2)/100.0) * (1.0-CTF(:,:,1)/100.0)
     !tropical tree
     CTF(:,:,2)=CTF(:,:,7) * (    CTF(:,:,2)/100.0)

     !WRITE CTF:
     call check(nf90_inq_varid(ncid,"CTF",var_id)); call check(nf90_put_var(ncid, var_id, CTF(:,:,1:6)  ))        

     !--------
     !EFs & LDF
print*,"EFS & LDF"
     allocate( ECOTYPE(g%nx, g%ny   ))   !ecotype
     allocate(    GTYP(g%nx, g%ny,4 ))   !growthtype fracs

     ECOTYPE(:,:)=FLOOR(interpolate(p,g,ecotype_file,varname="ecotype", method="mode"))
     !
     CTF_LIST=(/'crop ','tree ','grass','shrub'/)
     GTYP(:,:,1)=interpolate(p,g,ctf_file,varname="crop" , method="bilinear")
     GTYP(:,:,2)=interpolate(p,g,ctf_file,varname="tree" , method="bilinear")
     GTYP(:,:,3)=interpolate(p,g,ctf_file,varname="grass", method="bilinear")
     GTYP(:,:,4)=interpolate(p,g,ctf_file,varname="shrub", method="bilinear")
     !CTF_LIST=["NlTree","TrTree","BrTree","shrub ","grass ","crop  ","tree  "]
     !where ( CTF < 0.0 )
     !        CTF = 0.0
     !endwhere
     !CTF=CTF/100  ! % -> (0-1) . xq las frac están en porcentaje.
     allocate(OUTGRID(g%nx, g%ny, nefs+nldfs))   !outgrids EF1,EF2,...,LDF1,LDF2,..
     !allocate(OUTGRID(g%nx, g%ny, nefs, ncantype))   !outgrids EF1,EF2,...,LDF1,LDF2,.. !test v3.3
     OUTGRID=0.0
     j=0
     open(unit=1,file=trim(GtEcoEF_file),status='unknown',action='read')  !GtEcoEF_file:
       iostat=0                                                           !cantypeId ecotypeID var1, var2, ...,var19, var20, ..., var23
       do while(iostat == 0)       !loop por cada fila                    !crop      1         EF01, EF02, ..., EF19, LDF01, ..., LDF04
          read(1,*,iostat=iostat) GtID, EcoID, EF                         !crop      2         EF01, EF02, ..., EF19, LDF01, ..., LDF04
                                                                          !....      ...       EF01, EF02, ..., EF19, LDF01, ..., LDF04
                                                                          !tree      1700      EF01, EF02, ..., EF19, LDF01, ..., LDF04
          if ( j /= FINDLOC(CTF_LIST, GtID,1) ) then
                j=FINDLOC(CTF_LIST, GtID,1)
                !if (j==7) j=1
                print*,"   Processing Growth-type: "//GtID
          endif
         !=======> (!) ACÁ está el cuello de botella <=====
         do i=1,NEFS+NLDFS  !nvars: EF/LDF
         where ( ECOTYPE == EcoID )
             OUTGRID(:,:,i) = OUTGRID(:,:,i) + GTYP(:,:,j) * EF(i)
             !OUTGRID(:,:,i,j) =  CTF(:,:,j) * EF(i)  !test v3.3
             !OUTGRID(:,:,i,j) = EF(i)                !test v3.3
             !OUTGRID(:,:,i,j) = OUTGRID(:,:,i,j) +  EF(i)                !test v3.3
         endwhere
         enddo !i (var)
         !=======> (!) ACÁ está el cuello de botella <=====
       enddo !each row.
     close(1)

   !WRITE EF & LDF:
   call check(nf90_inq_varid(ncid,"EFS",var_id)); call check(nf90_put_var(ncid, var_id, OUTGRID(:,:,1:nefs ) ))        
   call check(nf90_inq_varid(ncid,"LDF",var_id)); call check(nf90_put_var(ncid, var_id, OUTGRID(:,:,nefs+1:nefs+nldfs ) ))        
   !
   deallocate(OUTGRID)
   deallocate(ECOTYPE)
   deallocate(CTF)

   !LAND FIELDS (BDSNP)
   if (run_BDSNP) then
print*,"BDSNP (LAND)"    
    allocate(LANDGRID(g%nx,g%ny,ncantype+1))  ! allocate(CTF(g%nx,g%ny,nvars))  
     LANDGRID(:,:,1) = interpolate(p,g,climate_file , varname="arid"    , method="mode")
     LANDGRID(:,:,2) = 1  
     LANDGRID(:,:,2) = interpolate(p,g,landtype_file, varname="landtype", method="mode")

    !WRITE LAND FIELDS:
    call check(nf90_inq_varid(ncid,"arid"    ,var_id)); call check(nf90_put_var(ncid, var_id, LANDGRID(:,:,1) ))        
    call check(nf90_inq_varid(ncid,"landtype",var_id)); call check(nf90_put_var(ncid, var_id, LANDGRID(:,:,2) ))        
   end if

 !Cierro NetCDF outFile
 call check(nf90_close(ncid))

end subroutine

 !----------------------------------
 !  DYNAMIC DATA:
 !---------------------------------
 subroutine prep_dynamic_data(g,p,lat,lon,laiv_file,nitro_file, fert_file,run_BDSNP)
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(len=16) :: outfile='mgn_dynamic.nc'
    logical :: run_BDSNP
    !LAIv
    character(len=200),intent(in) :: laiv_file
    real, allocatable :: LAIv(:,:,:)
    !NDEP
    character(len=200),intent(in) :: nitro_file
    real, allocatable :: NDEP(:,:,:)
    !NFERT
    character(len=200),intent(in) :: fert_file
    real, allocatable :: NFERT(:,:,:)
    !Coordinates
    real, intent(in)  :: lat(:,:),lon(:,:)
    integer ::var_id,ncid,x_dim_id,y_dim_id,date_dim_id,day_dim_id,month_dim_id
    integer :: nvars
    character(len=10):: temporal_avg  !flag on LAIv global file that specifies temporal average of data
    character(len=2):: kk
    character(len=3):: kkk
    !real :: lat,lon
    integer :: i,j,k
 
     !----
     !LAI: 
print*,"LAI"
     !
     !Idea for implementing 8-day LAIv: just read global_var (NVARS) and depending of this read & write the output.
     call check(nf90_open(laiv_file,nf90_nowrite, ncid) )
     call check(nf90_get_att(ncid, NF90_GLOBAL, "temporal_average", temporal_avg) )
     call check(nf90_close(ncid))
     if ( trim(temporal_avg) == "monthly") then
        nvars=12
     else if ( trim(temporal_avg) == "8-day") then
        nvars=46
     else
        print*,"Couldn't get temporal_average (8-day or monthly) global attribute from LAI file",temporal_avg;stop
     endif

     allocate( LAIv(g%nx,g%ny,nvars))  
     do k=1,nvars
         write(kk,'(I0.2)') k
         LAIv(:,:,k)=interpolate(p,g,inp_file=laiv_file,varname="laiv"//kk, method="bilinear") 
     enddo
     where (LAIv < 0.0 )
             LAIv=0.0
     endwhere

    if (run_BDSNP) then
    !----
    !NDEP:
    allocate( NDEP(g%nx,g%ny,12   ))  
    do k=1,12
        write(kk,'(I0.2)') k
        NDEP(:,:,k)  = interpolate(p,g,nitro_file, varname="nitro"//kk, method="bilinear")
    enddo
    where (NDEP < 0.0 )
       NDEP=0.0
    endwhere

    !----
    !NFERT:
    allocate(NFERT(g%nx,g%ny,365  ))  
    !Levanto netcdf input files
    do k=1,365
        write(kkk,'(I0.3)') k
        NFERT(:,:,k)  = interpolate(p,g,fert_file, varname="fert"//kkk, method="bilinear")
    enddo
    where (NFERT< 0.0 )
       NFERT=0.0
    endwhere

    end if
    !Creo NetCDF file
    !Create File and define dimensions and variables:
    call check(nf90_create(outFile, NF90_CLOBBER, ncid))
       !Define dimensions:
       call check(nf90_def_dim(ncid, "x"    , g%nx      , x_dim_id    ))
       call check(nf90_def_dim(ncid, "y"    , g%ny      , y_dim_id    ))
       call check(nf90_def_dim(ncid, "time" , nvars     , date_dim_id ))
       call check(nf90_def_dim(ncid, "month", 12        , month_dim_id))
       call check(nf90_def_dim(ncid, "day"  , 365       , day_dim_id  ))
       !Define variables:    
       ! Coordinates:
       call check(nf90_def_var(ncid, "lon" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id))
       call check(nf90_def_var(ncid, "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id))
       !!Times
       !call check(nf90_def_var(ncid, "Times",NF90_INT,   [x_dim_id,y_dim_id,date_dim_id], var_id))
       !call check(nf90_put_att(ncid, var_id,"units"    , "seconds since 1970-1-1"   ))
       !LAI
       call check(nf90_def_var(ncid, "LAI" , NF90_FLOAT, [x_dim_id,y_dim_id,date_dim_id],var_id))
       call check(nf90_put_att(ncid, var_id,"long_name", "LAIv"            ))
       call check(nf90_put_att(ncid, var_id,"units"    , "1"               ))
       call check(nf90_put_att(ncid, var_id,"var_desc" , "Leaf Area Index" ))
       if (run_BDSNP) then
       !NDEP
       call check(nf90_def_var(ncid, "NDEP", NF90_FLOAT, [x_dim_id,y_dim_id,month_dim_id],var_id))
       call check(nf90_put_att(ncid, var_id,"long_name", "LAIv"            ))
       call check(nf90_put_att(ncid, var_id,"units"    , "1"               ))
       call check(nf90_put_att(ncid, var_id,"var_desc" , "Leaf Area Index" ))
       !NFERT
       call check(nf90_def_var(ncid,"NFERT", NF90_FLOAT, [x_dim_id,y_dim_id, day_dim_id],var_id))
       call check(nf90_put_att(ncid, var_id,"long_name", "LAIv"            ))
       call check(nf90_put_att(ncid, var_id,"units"    , "1"               ))
       call check(nf90_put_att(ncid, var_id,"var_desc" , "Leaf Area Index" ))
       endif
       !Global Attributes
       call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "MEGAN input file"   ))
       call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
    call check(nf90_enddef(ncid))  
    !End NetCDF define mode
    
    !Abro NetCDF outFile
    call check(nf90_open(outFile, nf90_write, ncid      ))
      !!Coordinates
      call check(nf90_inq_varid(ncid,"lon"  ,var_id)); call check(nf90_put_var(ncid, var_id, lon ) )
      call check(nf90_inq_varid(ncid,"lat"  ,var_id)); call check(nf90_put_var(ncid, var_id, lat ) )    
      !!Times:
      !call check(nf90_inq_varid(ncid,"Times",var_id)); call check(nf90_put_var(ncid, var_id, Times ))
      !LAIv
      call check(nf90_inq_varid(ncid,"LAI"  ,var_id)); call check(nf90_put_var(ncid, var_id, LAIv/1000.0 ))
      !NDEP
      call check(nf90_inq_varid(ncid,"NDEP" ,var_id)); call check(nf90_put_var(ncid, var_id, NDEP        ))
      !NFERT
      call check(nf90_inq_varid(ncid,"NFERT",var_id)); call check(nf90_put_var(ncid, var_id, NFERT       ))
    !Cierro NetCDF outFile
    call check(nf90_close( ncid ))
 end subroutine


!***********************************************************************
!UTILS *****************************************************************
 subroutine check(status)
   integer, intent(in) :: status
   if (status /= nf90_noerr) then
     write(*,*) nf90_strerror(status)
     stop 'netcdf error'
   end if
 end subroutine check

 function atoi(str)     !string -> int
   implicit none
   character(len=*), intent(in) :: str
   integer :: atoi
   read(str,*) atoi
 end function
 function itoa(i)       !int -> string
    implicit none
    integer, intent(in) :: i
    character(len=20) :: itoa
    write(itoa, '(i0)') i
    itoa = adjustl(itoa)
 end function
 function rtoa(r)       !real -> string
    implicit none
    real, intent(in) :: r
    character(len=16) :: rtoa
    write(rtoa, '(F16.3)') r
    rtoa = adjustl(rtoa)
 end function

!STATISTICS:
pure function mode(arr)      !calculate mode of a 2D array
  implicit none
  integer, intent(in) :: arr(:,:) ! Input array
  integer :: mode,modeCount       ! Most frequent value & Count of the most frequent value
  integer :: i,j,N,count
  integer, allocatable :: arr1d(:)
  !Reshape arr to 1-D
  N = product(shape(arr))
  allocate(arr1d(N))
  arr1d = reshape(arr, [N])
  ! Initialize the mode and its count to zero
  mode = 0
  modeCount = 0
  ! Find the mode
  do i = 1, N
    count = 0
    do j = 1, N
      if ( arr1d(j) == arr1d(i)) then
        count = count + 1
      endif
    end do
    if (count > modeCount) then
      mode = arr1d(i)
      modeCount = count
    endif
  end do
  deallocate(arr1d)
end function
!***********************************************************************
!GRIDDESC **************************************************************
 subroutine read_GRIDDESC(griddescFile,gridName, p, g)
  implicit none
  character(200),intent(in) :: griddescFile
  character(*) ,intent(in)  :: gridName
  type(proj_type), intent(inout) :: p
  type(grid_type), intent(inout) :: g
  character(20) :: row
  integer :: iostat
  iostat=0
  open(unit=2,file=griddescFile,status='old',action='read',access='sequential')
  do while(iostat == 0)  !loop por cada fila
     read(2,*,iostat=iostat) row
     if ( trim(row) == trim(gridname)) then
       g%gName=row
       read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%nx,g%ny !projName xorig yorig xcell ycell nrows ncols
       rewind(2)
     endif
     if (trim(row) == trim(p%pName)) then
       read(2,*) p%typ,p%alp,p%bet,p%gam,p%xcent,p%ycent  
       iostat=1
     endif
  enddo
  close(2)

  call set_additional_proj_params(p)
  call set_additional_grid_params(p,g)
 end subroutine

 subroutine set_additional_proj_params(p)
    implicit none
    type(proj_type) ,intent(inout) :: p

    if ( p%typ == 2 ) then       !lambert conformal conic:        
       if ( ABS(p%alp - p%bet) > 0.1 ) then  !secant proj case
          p%p2=     LOG( COS(p%alp           *deg2rad )/ COS(p%bet           *deg2rad)   )
          p%p2=p%p2/LOG( TAN((45.0+0.5*p%bet)*deg2rad )/ TAN((45.0+0.5*p%alp)*deg2rad)   ) !n
        else                                 !tangent proj case
          p%p2=SIN(p%alp*deg2rad) !n
       endif
       p%p3=R_EARTH*(COS(p%alp*deg2rad)*TAN((45+0.5*p%alp)*deg2rad)**p%p2)*(1/p%p2)  !F
       p%p1=p%p3/(TAN((45 + 0.5*p%ycent)*deg2rad)**p%p2)                             !rho0 

    else if ( p%typ == 6 ) then  !polar secant stereographic
       print*, "Todavia no desarrollado soporte para proyeccion polar stereografica (ptype=",p%typ,")."; stop

    else if ( p%typ == 7 ) then  !equatorial mercator
       p%p1=COS(p%alp*deg2rad)   !k0

    else if ( p%typ == 1 ) then
       print*, "sistema de coordenadas latlon!"
    else
        print*, "codigo de proyección invalido:",p%typ,"."; stop
    end if
 end subroutine

 subroutine set_additional_grid_params(p,g)
    implicit none
    type(proj_type) ,intent(inout) :: p
    type(grid_type) ,intent(inout) :: g
    real :: latmin,lonmin,latmax,lonmax
    !Obtener coordenadas del centro de la grilla, min y max:
    g%xc=0.0;g%yc=0.0;g%xmax=g%xmin+g%dx*g%nx; g%ymax=g%ymin+g%dy*g%ny

    !calculo minimos y maximos de latlon 
    !   (ojo! Dado que son transf no-lineales no corresponden necesariamente a los vertices)
    call xy2ll(p,g%xmin,g%ymin,g%lonmin,g%latmin)       !lower-left
    call xy2ll(p,g%xmax,g%ymax,g%lonmax,g%latmax)       !upper-right
 
    !latmin
    call xy2ll(p,g%xmin+g%dx*g%nx*0.5, g%ymin,lonmin,latmin)
    g%latmin=min(g%latmin,latmin)
    !latmax
    call xy2ll(p,g%xmax-g%dx*g%nx*0.5,g%ymax ,lonmax,latmax)
    g%latmax=max(g%latmax,latmax)

    !lonmin   
    call xy2ll(p,g%xmin,g%ymin+g%dy*g%ny*0.5,lonmin,latmin)
    g%lonmin=min(g%lonmin,lonmin)
    !         
    call xy2ll(p,g%xmin,g%ymax              ,lonmin,latmin)
    g%lonmin=min(g%lonmin,lonmin)

    !lonmax
    call xy2ll(p,g%xmax,g%ymax-g%dy*g%ny*0.5,lonmax,latmax)
    g%lonmax=max(g%lonmax,lonmax)
    !      
    call xy2ll(p,g%xmax,g%ymin              ,lonmax,latmax)
    g%lonmax=max(g%lonmax,lonmax)

 end subroutine
!***********************************************************************
!PROJ     **************************************************************
!COORDINATE TRANSFORMATION FUNCTIONS:======================================
 subroutine xy2ll(p,x,y,lon,lat)
     implicit none                            
     type(proj_type) ,intent(in) :: p
     real, intent(in)   :: x,y
     real, intent(inout):: lon,lat
 
     if      ( p%typ == 2 ) then  !Lambert Conformal Conic:
        call xy2ll_lcc(p,x,y,lon,lat)
     else if ( p%typ == 6 ) then  !polar secant stereographic
        call xy2ll_stere(p,x,y,lat,lon)
     else if ( p%typ == 7 ) then  !equatorial mercator
        call xy2ll_merc(p,x,y,lon,lat)
     else if ( p%typ == 1 ) then  !latlon
        lon=x;lat=y               !no transformation needed       
     else
        print*, "codigo de proyección invalido:",p%typ,"."; stop
     end if
 end subroutine
 subroutine ll2xy(p,lon,lat,x,y)
       implicit none                            
       type(proj_type) ,intent(in) :: p
       real, intent(in):: lon,lat
       real, intent(inout)   :: x,y
 
       if      ( p%typ == 2 ) then  !Lambert Conformal Conic:
          call ll2xy_lcc(p,lon,lat,x,y)
       else if ( p%typ == 6 ) then  !Polar Secant Stereographic
          call ll2xy_stere(p,lon,lat,x,y)
       else if ( p%typ == 7 ) then  !Equatorial Mercator
          call ll2xy_merc(p,lon,lat,x,y)
       else if ( p%typ == 1 ) then  !latlon
          x=lon;y=lat               !no transformation needed
       else
          print*, "codigo de proyección invalido:",p%typ,"."; stop
       end if                             
 end subroutine

 !--------------------------------------------------------------------------
 !LAMBERT CONFORMAL CONIC:
 subroutine xy2ll_lcc(p,x,y,lon,lat)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)   :: x,y
   real, intent(inout):: lon,lat
   real :: n,F,rho0,rho,theta
   
   rho0=p%p1
   n=p%p2
   F=p%p3
   
   theta=ATAN(x/(rho0-y))*rad2deg
   rho=SIGN(1.0,n) * SQRT( x*x + (rho0-y)*(rho0-y))
   
   lon=p%gam+theta/n
   lat=2.0 * ATAN( (F/rho)**(1/n) )*rad2deg - 90.0 
 end subroutine
 subroutine ll2xy_lcc(p,lon,lat,x,y)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)      :: lon,lat
   real, intent(inout)   :: x,y
   real :: n,F,rho0,rho,dlon

   !interm params:
   rho0=p%p1
   n=p%p2
   F=p%p3

   rho=F/(TAN((45.0 + 0.5*lat)*deg2rad)**n)
   dlon=lon-p%gam
   !
   x=     rho*SIN(n*dlon*deg2rad )
   y=rho0-rho*COS(n*dlon*deg2rad )
 end subroutine
 !--------------------------------------------------------------------------
 !MERCATOR                
 subroutine xy2ll_merc(p,x,y,lon,lat)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)   :: x,y
   real, intent(inout):: lon,lat
   real :: k0R,phi
 
   k0R=p%p1*R_EARTH     !es una longitud (k0*R_EARTH)
   phi=y/k0R            !es un angulo
   
   lon=p%gam + x/k0R*rad2deg
   lat=90.0-2*ATAN( EXP(-phi) )*rad2deg
 end subroutine
 subroutine ll2xy_merc(p,lon,lat,x,y)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)      :: lon,lat
   real, intent(inout)   :: x,y
   real :: k0,lon0

   k0=p%p1              !adminesional
   lon0=p%gam           !es un angulo

   x=k0*R_EARTH*(lon-lon0)*deg2rad
   y=k0*R_EARTH*LOG(TAN((45.0+0.5*lat)*deg2rad))
 end subroutine
!--------------------------------------------------------------------------
 !POLAR STEREOGRAPHIC     
 subroutine xy2ll_stere(p,x,y,lon,lat)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)   :: x,y
   real, intent(inout):: lon,lat
   real :: k,rho

   stop 'Stereographic proyection not yet tested!'
   rho = sqrt(x*x+y*y)
   k = 2.0*ATAN( rho/2.0/R_EARTH )

   lat =         ASIN(   COS(k)*SIN(p%alp*deg2rad) + y*SIN(k)*COS(p%alp*deg2rad)/rho )               * rad2deg
   lon = p%gam + ATAN( x*SIN(k)  / ( rho*COS(p%alp*deg2rad)*COS(k) - y*SIN(p%alp*deg2rad)*SIN(k) ) ) * rad2deg

 end subroutine
 subroutine ll2xy_stere(p,lon,lat,x,y)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)      :: lon,lat
   real, intent(inout)   :: x,y
   real :: k!,hemi

   stop 'Stereographic proyection not yet tested!'
   !hemi=SIGN(1.0,p%alp)
   k = 2.0*R_EARTH / (1 + SIN(p%alp*deg2rad)*SIN(lat*deg2rad) + COS(p%alp*deg2rad)*COS(lat*deg2rad)*COS( (lon-p%gam)*deg2rad ))

   x = k *   COS( lat *deg2rad) * SIN( (lon - p%gam)*deg2rad )
   y = k * ( COS(p%alp*deg2rad) * SIN(  lat         *deg2rad ) - SIN(p%alp*deg2rad)*COS(lat*deg2rad)*COS((lon-p%gam)*deg2rad) )
 end subroutine
!!END COORDINATE TRANFORMATION FUNCTIONS====================================

!***********************************************************************
!INTERPOLATION**********************************************************
!!"HOME-MADE" interpolation function:
function interpolate(p,g,inp_file,varname,method)       result(img2)
 implicit none
 type(grid_type), intent(in)  :: g  !desired grid
 type(proj_type), intent(in)  :: p  !proj of desired grid
 character(*), intent(in)     :: inp_file,varname,method
 integer                      :: methodId
 real,allocatable :: img2(:,:)      !output array

 integer :: i,j,k

 type(grid_type)  :: GG,GC          !global grid (input grid) &  global grid (CROPPED)
 real, allocatable :: img1(:,:)      !cropped img to be interpolated
 !real,allocatable :: lat(:),lon(:)

 !integer :: ncid,latid,lonid,varid  !int for netcdf handling

 integer :: is,ie,js,je     !indices that defines subarray
 real    :: px,py,x,y       !dummy variables for temporal coordinates
 real    :: w11,w12,w21,w22 !weights for bilinear interpolation
 real    :: p11,p12,p21,p22 !params. for bilinear interpolation
 real    :: x1,x2,y1,y2     !params. for interpolation
 integer :: i1,i2,j1,j2     !dummy indexes for interpolation

 integer :: scale_x,scale_y

 print*,"  Interpolando: "//trim(inp_file)//":"//trim(varname)//"..."
 ! Asumo que estoy trabajando con grillas regulares (dx/dy =cte.).    
 ! Asumo que lat y lon estan ordenados de forma creciente.

 call get_cropped_img(p,g,inp_file,varname,img1,GC)
 
 !Veo si la grilla destino es mas densa o no que la original.
 call xy2ll(p,g%xmin,g%ymin,x1,y1)  !
 call xy2ll(p,g%xmax,g%ymax,x2,y2)  !

 scale_x=CEILING((x2-x1)/(g%nx)/GC%dx)
 scale_y=CEILING((y2-y1)/(g%ny)/GC%dy)


 if ( method == "bilinear") then
         if (scale_x > 2 .or. scale_y > 2) then !Choose another method if source grid is denser than destination grid
                methodId=2
         else
                methodId=3
         endif
 else if ( method == "bicubic"  ) then
        methodId=4
 else if ( method == "avg"      ) then
        methodId=2
 else if ( method == "mode"     ) then
        methodId=1
 endif

 allocate(img2(g%nx,g%ny))  !array a interpolar:

 !REGRIDDING:
 if( methodId == 1) then
 !MODE   
    do i=1,g%nx
       do j=1,g%ny
           px=g%xmin+g%dx*i  !projected coordinate-x
           py=g%ymin+g%dy*j  !projected coordinate-y
                                                                                       
           call xy2ll(p,px,py,x,y)  !I want coordinates in same proj than global file.

           i1=MAX(     1, FLOOR( (x-GC%lonmin) / GC%dx - scale_x*0.5 ) )
           j1=MAX(     1, FLOOR( (y-GC%latmin) / GC%dy - scale_y*0.5 ) )
           i2=MIN( GC%nx, i1+scale_x                                   )
           j2=MIN( GC%ny, j1+scale_y                                   )

           if ( i1 >= 1 .and. i2 <= GC%nx .and. j1 >=1 .and. j2 <= GC%ny ) then
               !! Mode                                                                
               img2(i,j)=MODE(INT(img1(i1:i2,j1:j2)))                        !mode
           else
               img2(i,j)=0
           endif
       enddo
    enddo
 endif
 if( methodId == 2) then
 !AVERAGE (DEFAULT)
      do i=1,g%nx
         do j=1,g%ny
             px=g%xmin+g%dx*i  !projected coordinate-x
             py=g%ymin+g%dy*j  !projected coordinate-y
                                                                                         
             call xy2ll(p,px,py,x,y)  !I want coordinates in same proj than global file.
                                                                                          
             i1=MAX(     1, FLOOR( (x-GC%lonmin) / GC%dx - scale_x*0.5 ) )
             j1=MAX(     1, FLOOR( (y-GC%latmin) / GC%dy - scale_y*0.5 ) )
             i2=MIN( GC%nx, i1+scale_x                                   )
             j2=MIN( GC%ny, j1+scale_y                                   )
                                                                                          
             if ( i1 > 1 .and. i2 < GC%nx .and. j1 > 1 .and. j2 < GC%ny ) then
                 !! Average:
                 img2(i,j)=SUM(img1(i1:i2,j1:j2))/((i2-i1+1)*(j2-j1+1))  !average
             else
                 img2(i,j)=0
             endif
         enddo
      enddo
  endif
  !INTERPOLATION:  
  if ( methodId == 3 ) then
  !! Bilineal Interp:
     !print*,"Bilinear interpolation. ",scale_x,scale_y
     do i=1,g%nx
        do j=1,g%ny
            !Position where to interpolate
            px=g%xmin+g%dx*i  !projected coordinate-x
            py=g%ymin+g%dy*j  !projected coordinate-y
            call xy2ll(p,px,py,x,y)  !I want coordinates in same proj than global file.
            !indices:
            i1=FLOOR( (x-GC%lonmin) / GC%dx );  i2=i1+1
            j1=FLOOR( (y-GC%latmin) / GC%dy );  j2=j1+1
            if ( i1 > 1 .and. i2 <= GC%nx .and. j1 > 1 .and. j2 <= GC%ny ) then
               !points (coordinates)    !   p12(i1,j2)    p22(i2,j2)
               x1=GC%lonmin+GC%dx*i1    !       *- - - - - -*            
               x2=GC%lonmin+GC%dx*i2    !       |           |           
               y1=GC%latmin+GC%dy*j1    !       |           |           
               y2=GC%latmin+GC%dy*j2    !       |           |           
               !points (values)         !       |           |           
               p11=img1(i1,j1)          !       *- - - - - -*                                    
               p12=img1(i1,j2)          !   p11(i1,j1)    p21(i2,j1)
               p21=img1(i2,j1)
               p22=img1(i2,j2)
               !weights:
               w11 =(x2 - x )*(y2 - y )/(GC%dx*GC%dy)
               w12 =(x2 - x )*(y  - y1)/(GC%dx*GC%dy)
               w21 =(x  - x1)*(y2 - y )/(GC%dx*GC%dy)
               w22 =(x  - x1)*(y  - y1)/(GC%dx*GC%dy)
               !Bilineal formula:
               img2(i,j)= p11*w11 + p12*w12 + p21*w21 + p22*w22 ! DOT_PRODUCT(p,w)
            else
               img2(i,j)=0.0
            endif
        enddo
     enddo
  endif
  !if ( methodId == 4 ) then
  !!    !PROGRAMAR INTERPOLACION BICUBICA!
  !endif

 end function


subroutine get_Cropped_Img(p,g,inp_file,varname,img,GC)
   implicit none
   type(grid_type), intent(in)  :: g  !desired grid
   type(proj_type), intent(in)  :: p  !proj of desired grid
   character(*), intent(in)     :: inp_file, varname
   
   real,allocatable, intent(inout) :: img(:,:)      !output array
   type(grid_type),intent(inout)  :: GC          !global grid (input grid) &  global grid (CROPPED)

   integer :: i,j,k

   type(grid_type)  :: GG!,GC          !global grid (input grid) &  global grid (CROPPED)
   real,allocatable :: lat(:),lon(:)

   integer :: ncid,latid,lonid,varid  !int for netcdf handling

   integer :: is,ie,js,je     !indices that defines subarray

  !Leo inp_file:
  call check(nf90_open(trim(inp_file), nf90_write, ncid ))
    call check( nf90_inq_dimid(ncid, "lat",latid )             )
    call check( nf90_inquire_dimension(ncid, latid, len=GG%ny ))
    call check( nf90_inq_dimid(ncid, "lon",lonid )             )
    call check( nf90_inquire_dimension(ncid, lonid, len=GG%nx ))
    allocate(lat(GG%ny))
    allocate(lon(GG%nx))
    !lat-----------------------------------------------------     
    call check( nf90_inq_varid(ncid,trim("lat"), varid    ))
    !call check( nf90_get_var(ncid, varid , lat(GG%ny:1:-1)))     !lat viene alreves
    call check( nf90_get_var(ncid, varid , lat            ))    
    !lon-----------------------------------------------------     
    call check( nf90_inq_varid(ncid,trim("lon"), varid   ))
    call check( nf90_get_var(ncid, varid , lon           ))
    !------------------------------------------------------
    !Levanto parametros de grilla a interpolar:
    GG%latmin=lat(    1); GG%lonmin=lon(  1)     !lower-left corner?
    GG%latmax=lat(GG%ny); GG%lonmax=lon(GG%nx)   !upper-right corner?
                                                                                                                                                       
    GG%dy=ABS(GG%latmin-lat(2))  !delta lat
    GG%dx=ABS(GG%lonmin-lon(2))  !delta lon
         !Checkear que sea una grilla regular
         if( ABS(GG%dx - ABS(GG%lonmax-GG%lonmin)/(GG%nx-1)) < 1E-5  ) then; continue;else; print*,"Lon NO es regular.",GG%dx;stop;endif
         if( ABS(GG%dy - ABS(GG%latmax-GG%latmin)/(GG%ny-1)) < 1E-5  ) then; continue;else; print*,"Lat NO es regular.",GG%dy;stop;endif
    !indices de sub-array:
    is=MAX(  1   ,   FLOOR( (g%lonmin-GG%lonmin)/GG%dx) ) !calc min y max indices    
    ie=MIN(GG%nx , CEILING( (g%lonmax-GG%lonmin)/GG%dx) ) !calc min y max indices 
    js=MAX(  1   ,   FLOOR( (g%latmin-GG%latmin)/GG%dy) ) !calc min y max indices
    je=MIN(GG%ny , CEILING( (g%latmax-GG%latmin)/GG%dy) ) !calc min y max indices
    !parametros de grilla:
    GC%nx=ABS(ie-is)+1;  GC%ny=ABS(je-js)+1 
    GC%dx=GG%dx       ;  GC%dy=GG%dy       
    GC%lonmin=lon(is) ;  GC%latmin=lat(js)
    GC%lonmax=lon(ie) ;  GC%latmax=lat(je)

    deallocate(lat)
    deallocate(lon)
    allocate(img(GC%nx, GC%ny)) 
    !levanto variable a interpolat (es lo que mas tarda)-----     
    call check( nf90_inq_varid(ncid,trim(varname), varid ))
    !call check( nf90_get_var(ncid, varid , img1(1:GC%nx,GC%ny:1:-1), start=[is,GG%ny-je],count=[GC%nx,GC%ny] ) ) !Acá tarda MUCHO..!lat viene alreves
    call check( nf90_get_var(ncid, varid , img, start=[is,js],count=[GC%nx,GC%ny] ) )
    !--------------------------------------------------------     
  call check(nf90_close(ncid))
  !print*,"       File, minval, maxval",trim(inp_file),minval(img),maxval(img) !debug

end subroutine

end module 
!end program
