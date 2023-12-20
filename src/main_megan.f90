program main
   ! program:     MEGAN  v3.2
   ! description: Biogenic VOCs emission model (by Alex Guenther)
   ! coded by:    Alex Guenther, Ling Huang, Xuemei Wang, Jeff Willison, among others.
   ! adapted by:  Ramiro Espada (date: 11/2023)

   use netcdf   
   use megan_v32  !megan module: (megan_voc)
   use nox_mod    !megan module: (megan_nox)
   !use bdsnp      !megan module: (bdsnp_nox)
   
   implicit none

   type grid_type
       integer  :: nx,ny,nz,nt  !number of cells in x-y direction (ncols, nrows, nlevs, ntimes)
   end type grid_type

   type(grid_type) :: grid

   !namelist variables:
   integer           :: iostat
   integer           :: ncid,var_id
   character(len=19) :: start_date, end_date
   character(len=16) :: mechanism='CBM05' !'CBM6' 'CB6A7', 'RACM2','CRACM', 'SAPRC' 'NOCON'
   !logical           :: prep_megan=.false.
   character(250)    :: griddesc_file,gridname,ctf_file,lai_file,ef_file,ldf_file,ndep_file,fert_file,land_file
   character(250)    :: met_files(24)=""
   character(4)      :: lsm
   !date vars:
   integer      :: end_date_s,current_date_s,n_hours
   character(4) :: YYYY
   character(3) :: DDD
   character(2) :: MM,DD,HH,current_day="99",current_month="99"
   
   !input variables:
   character(19), allocatable, dimension(:) :: times                         !(t)     from wrfout
   real,    allocatable, dimension(:,:)     :: lon,lat                       !(x,y)   from wrfout
   real,    allocatable, dimension(:,:)     :: tmp,ppfd,u10,v10,pre,hum,rain !(x,y,t) from wrfout
   real,    allocatable, dimension(:,:)     :: smois,stemp                   !(x,y,t) from wrfout
   integer, allocatable, dimension(:,:)     :: stype                         !(x,y)   from wrfout
   integer, allocatable, dimension(:,:)     :: arid,non_arid,landtype        !(x,y)   from prep_megan
   real,    allocatable, dimension(:,:,:)   :: ctf,ef,ldf_in                    !(x,y,*) from prep_megan
   !real,    allocatable, dimension(:,:)    :: needl,tropi,broad,shrub,grass,crop
   real,    allocatable, dimension(:,:)     :: lai,ndep,fert                 !(x,y,t) from prep_megan

   character(250)    :: met_file
   integer           :: n_files, head
   !intermediate vars:
   integer :: t,i!,j,k
   integer, parameter :: layers=5
   real,    allocatable, dimension(:,:)     :: wind                           !(x,y,t) from wrfout
   real,    allocatable, dimension(:,:)     :: tmp_min,tmp_max,wind_max,tmp_avg,ppfd_avg !(x,y) daily meteo vars
  
   !output vars:
   real,    allocatable, dimension(:,:,:,:) :: out_buffer                     !(x,y,nclass,t)

   !---
   namelist/megan_nl/start_date,end_date,met_files,ctf_file,lai_file,ef_file,ldf_file,ndep_file,fert_file,land_file,mechanism,griddesc_file,gridname,lsm
   read(*,nml=megan_nl, iostat=iostat)
   if( iostat /= 0 ) then
     write(*,*) 'megan: failed to read namelist; error = ',iostat
     stop
   end if

   !---
   print '("Archivos meteorologicos.. ")'
   do i=1,size(met_files)
      if (met_files(i) =="") then
         n_files=i-1;
         exit
      else 
           print*,met_files(i)
      endif
   end do
   head=1
   met_file=met_files(head)

   !---
   print '("Preparo mecanísmo y especies.. ")'
   call megan_map(mechanism)
   print*,"  Mecanísmo: ",mechanism
   print*,"  Especies: ",megan_names

   !---
   print '("Read meteo params and coordinates.. ")'
   call get_meteofile_info(met_file, grid, lat, lon, times)
   
   !--- 
   print '("Get static data.. ")'
   call prep_static_data(grid) !,arid,non_arid,landtype,ctf,ef,ldf)
 
   !--- 
   print '("Allocate output buffer.. ")' 
   allocate(out_buffer(grid%nx,grid%ny,n_spca_spc,24)) !main out array

   !--- 
   print '("Comienza loop.. ")'
   current_date_s = atoi( date(start_date, "%s"))
       end_date_s = atoi( date(  end_date, "%s"))
   n_hours=floor(real((end_date_s-current_date_s)/3600))
   !do t=1,n_hours
   do while (current_date_s <= end_date_s)

      call seconds_to_date(current_date_s,YYYY,MM,DD,DDD,HH)
      !t=atoi(HH) !time index
      
      t=get_t_index(Times, current_date_s) !met t-index

      print*,"t = ",t
      !if ( t == -99 ) then
      !    head=head+1;met_file=met_files(head)
      !    do while ( t == -99 .and. head <= n_files )
      !    call prep_dynamic_data(met_file,grid,t)
      !    t=get_t_index(Times, current_date_s) !met t-index
      !    end do
      !end if
      call get_hourly_data(grid,t)

      if ( current_day /= DD .or. current_date_s == end_date_s ) then
         if ( current_day /= "99" ) then
            call write_output_file(grid,YYYY,DDD,MECHANISM)
         endif
         current_day=DD
         call get_daily_data(grid,DDD)
         if ( current_month /= MM ) then
            current_month=MM
            call get_monthly_data(grid,MM)
         endif    
      endif    

     call write_diagnostic_file(grid)

     call megan_voc(atoi(yyyy),atoi(ddd),atoi(hh),     & !MEGAN VOCs model
             grid%nx,grid%ny,layers,                   &  !dimensions
             lat,lon,                                  &
             tmp,ppfd,wind,pre,hum,                    & !INPs: 
             lai, lai,                                 &
             ctf, ef, ldf_in,                          &  !LAI (past), LAI (current), EF, LDF              &
             lsm,stype,smois,                          &  !mesea Inps: land surface model, soil type, soil moisture!
             tmp_max,tmp_min,wind_max,tmp_avg,ppfd_avg,&  !max temp, min temp, max wind, daily temp and daily ppfd
             out_buffer(:,:,:,atoi(HH))                )  !Emis array

      !soil NO model:
      !if ( bdnsp_soil_model ) then
      !    call bdsnp_nox()
      !else
      !!call megan_nox(atoi(yyyy),atoi(ddd),atoi(hh),  &
      !!         grid%nx,grid%ny,                      &
      !!         lat,                                  &
      !!         tmp,rain,                             &
      !!         lsm,stype,stemp,smois,                &
      !!         ctf, lai,                             &
      !!         out_buffer(:,:,8,atoi(HH))           )
      !endif

     !laip=laic
     !----
     !print '("  Escribiendo archivo diagnóstico.. ")'
     !call write_diagnostic_file(grid)

     current_date_s=current_date_s + 3600  !siguiente hora!
   enddo

   print '("Fin. ")'
contains


function date(date_str, fmt_str) result(output) !Interfaz a "date"
  implicit none
  integer :: iostat
  character(*), intent(in) :: date_str, fmt_str
  character(256)           :: command
  character(20)            :: output
  command="date -d '"//trim(date_str)//"' '+"//trim(fmt_str)//"'  > tmp_date_jodemequeteneselmismonombre.txt"
  call system( trim(command) )
  !print*,trim(command)
  open(9, file='tmp_date_jodemequeteneselmismonombre.txt', status='old',action='read'); read(9, '(A)', iostat=iostat) output;  close(9)
  call system('rm tmp_date_jodemequeteneselmismonombre.txt')
end function

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

subroutine seconds_to_date(date_s,yyyy,mm,dd,ddd,hh)
   implicit none
   integer     ,intent(in)    :: date_s
   character(4),intent(inout) :: YYYY
   character(3),intent(inout) :: DDD
   character(2),intent(inout) :: MM,DD,HH

    YYYY=date("@"//itoa(date_s), "%Y")!año
    MM=date("@"//itoa(date_s), "%m")  !mes
    DD=date("@"//itoa(date_s), "%d")  !dia
    DDD=date("@"//itoa(date_s), "%j") !dia juliano
    HH=date("@"//itoa(date_s), "%H")  !hora
    print '("  ",A4,"-",A2,"-",A2," (día: ",A3,"), ",2A,"hs.")',YYYY,MM,DD,DDD,HH
end subroutine


subroutine get_meteofile_info(meteo_file,g,lat,lon,times)
   implicit none
   character(250) ,intent(in)    :: meteo_file
   type(grid_type),intent(inout) :: g
   real, allocatable, dimension(:,:),intent(inout) :: lat,lon
   character(19), allocatable, intent(inout), dimension(:) :: times
   integer :: ncid,dimid,varid,time_len,i,t

   call check(nf90_open(trim(meteo_file), nf90_write, ncid ))
      !grid dimensions
      call check (nf90_inq_dimid(ncid,'Time'       ,  dimId   ))
      call check (nf90_inquire_dimension(ncid, dimId,len=g%nt ))
      call check (nf90_inq_dimid(ncid,'west_east'  ,  dimId   ))
      call check (nf90_inquire_dimension(ncid, dimId,len=g%nx ))
      call check (nf90_inq_dimid(ncid,'south_north',  dimId   ))
      call check (nf90_inquire_dimension(ncid, dimId,len=g%ny ))
      call check (nf90_inq_dimid(ncid,'bottom_top' ,  dimId   ))
      call check (nf90_inquire_dimension(ncid, dimId,len=g%nz ))
  
      g%nx=g%nx-2;g%ny=g%ny-2 !esto hay que arreglarlo! (supongo que elimina los bordes)
     
      !lat lon coordinates
      if (.not. allocated(lat)  ) allocate( lat(g%nx,g%ny))
      if (.not. allocated(lon)  ) allocate( lon(g%nx,g%ny))
      if (.not. allocated(Times)) allocate( Times(g%nt))   
      call check( nf90_inq_varid(ncId,'XLAT' , varId)); call check(nf90_get_var(ncId, varId, lat   ))
      call check( nf90_inq_varid(ncId,'XLONG', varId)); call check(nf90_get_var(ncId, varId, lon   ))
      call check( nf90_inq_varid(ncId,'Times', varId)); call check(nf90_get_var(ncId, varId, times ))
      !date and time data
   call check(nf90_close(ncid))

   do t=1,g%nt
      i=scan(Times(t),'_')
      Times(t)(i:i)=' '
   enddo

print*,"nx, ny, nt, times: ",g%nx,g%ny,g%nt,times

end subroutine

integer function get_t_index(times,current_date_s)
  implicit none
  integer :: i,secs,current_date_s!,get_t_index
  character(19), dimension(:), intent(in) :: times
  get_t_index=-99
  do i=1,size(times)
     secs=atoi(date( trim(Times(i)), "%s"))
     if ( secs == current_date_s ) then
         get_t_index=i;exit;
     endif
  enddo
endfunction

subroutine prep_static_data(g) 
  implicit none
  type(grid_type) :: g
  integer         :: i,j,k

  character(len=10),dimension(19) :: ef_vars=["EF_ISOP   ", "EF_MBO    ", "EF_MT_PINE", "EF_MT_ACYC", "EF_MT_CAMP", "EF_MT_SABI", "EF_MT_AROM", "EF_NO     ", "EF_SQT_HR ", "EF_SQT_LR ", "EF_MEOH   ", "EF_ACTO   ", "EF_ETOH   ", "EF_ACID   ", "EF_LVOC   ", "EF_OXPROD ", "EF_STRESS ", "EF_OTHER  ", "EF_CO     "]
  character(len=5),dimension(4)   :: ldf_vars=["LDF03","LDF04","LDF05","LDF06"]

print*,"NRTYP: ",nrtyp
  !Allocation of variables to use
  allocate(    arid(g%nx,g%ny               ))
  allocate(non_arid(g%nx,g%ny               ))
  allocate(landtype(g%nx,g%ny               ))
  allocate(      ef(g%nx,g%ny,size( ef_vars)))
  allocate(  ldf_in(g%nx,g%ny,size(ldf_vars)))
  allocate(     ctf(g%nx,g%ny,NRTYP         ))

  !LAND                                                                               
  call check(nf90_open(trim(land_file), nf90_write, ncid ))
     call check( nf90_inq_varid(ncid,'LANDTYPE', var_id )); call check( nf90_get_var(ncid, var_id, LANDTYPE ))
     call check( nf90_inq_varid(ncid,'ARID'    , var_id )); call check( nf90_get_var(ncid, var_id, ARID     ))
     call check( nf90_inq_varid(ncid,'NONARID' , var_id )); call check( nf90_get_var(ncid, var_id, NON_ARID ))
  call check(nf90_close(ncid))

  !CTS (canopy type fractions)
  call check(nf90_open(trim(ctf_file), nf90_write, ncid ))
      call check( nf90_inq_varid(ncid,'CTS', var_id )); call check(nf90_get_var(ncid,var_id,CTF,[1,1,1,1],[g%nx,g%ny,1,NRTYP]))
      !call check( nf90_inq_varid(ncid,'NEEDL', var_id )); call check( nf90_get_var(ncid, var_id, NEEDL ))
      !call check( nf90_inq_varid(ncid,'TROPI', var_id )); call check( nf90_get_var(ncid, var_id, TROPI ))
      !call check( nf90_inq_varid(ncid,'BROAD', var_id )); call check( nf90_get_var(ncid, var_id, BROAD ))
      !call check( nf90_inq_varid(ncid,'SHRUB', var_id )); call check( nf90_get_var(ncid, var_id, SHRUB ))
      !call check( nf90_inq_varid(ncid,'GRASS', var_id )); call check( nf90_get_var(ncid, var_id, GRASS ))
      !call check( nf90_inq_varid(ncid,'CROP' , var_id )); call check( nf90_get_var(ncid, var_id, CROP  ))
  call check(nf90_close(ncid))
  ctf=CTF*0.01  !from % to [0-1].

  !EF:
  call check(nf90_open(trim(  ef_file), nf90_write, ncid ))
  do k=1,size(ef_vars)
     call check( nf90_inq_varid(ncid,trim( ef_vars(k)), var_id )); call check( nf90_get_var(ncid, var_id , EF(:,:,k) ))
  enddo
  call check(nf90_close(ncid))

  !LDF:
  call check(nf90_open(trim( ldf_file), nf90_write, ncid ))
  do k=1,size(ldf_vars)
     call check( nf90_inq_varid(ncid,trim(ldf_vars(k)), var_id )); call check( nf90_get_var(ncid, var_id ,ldf_in(:,:,k) ))
  enddo
  call check(nf90_close(ncid))

  !From meteo:
  if (.not. allocated(stype)) then; allocate( stype(g%nx,g%ny));endif
  call check(nf90_open(trim(met_file), nf90_write, ncid ))
      call check(nf90_inq_varid(ncid,'ISLTYP', var_id));
      call check(  nf90_get_var(ncid, var_id, stype, [1,1,1], [g%nx,g%ny,1]  ))
  call check(nf90_close(ncid))
end subroutine

subroutine get_hourly_data(g,t)
  implicit none
  type(grid_type) :: g
  integer,intent(in) :: t
  integer :: ncid,time_dimid,time_len
  
  print*,"  Preparo inputs dinámicos horarios.."

  !from meteo:
  if (.not. allocated(tmp)  ) allocate(  tmp(g%nx,g%ny))! then;;endif 
  if (.not. allocated(ppfd) ) allocate( ppfd(g%nx,g%ny))! then;;endif
  if (.not. allocated(u10)  ) allocate(  u10(g%nx,g%ny))! then;;endif
  if (.not. allocated(v10)  ) allocate(  v10(g%nx,g%ny))! then;;endif
  if (.not. allocated(pre)  ) allocate(  pre(g%nx,g%ny))! then;;endif
  if (.not. allocated(hum)  ) allocate(  hum(g%nx,g%ny))! then;;endif
  if (.not. allocated(rain) ) allocate( rain(g%nx,g%ny))! then;;endif
  if (.not. allocated(stemp)) allocate(stemp(g%nx,g%ny))! then;;endif
  if (.not. allocated(smois)) allocate(smois(g%nx,g%ny))! then;;endif
  if (.not. allocated(lai)  ) allocate(  lai(g%nx,g%ny))! then;;endif
  if (.not. allocated(wind) ) allocate( wind(g%nx,g%ny))! then;;endif 
  call check(nf90_open(trim(met_file), nf90_write, ncid ))
    call check( nf90_inq_varid(ncid,'U10'   , var_id)); call check(nf90_get_var(ncid, var_id,  U10, [1,1,t]  ))!,[g%nx,g%ny,1] )) !,[1,1,1]
    call check( nf90_inq_varid(ncid,'V10'   , var_id)); call check(nf90_get_var(ncid, var_id,  V10, [1,1,t]  ))!,[g%nx,g%ny,1] )) !,[1,1,1]
    call check( nf90_inq_varid(ncid,'T2'    , var_id)); call check(nf90_get_var(ncid, var_id,  TMP, [1,1,t]  ))!,[g%nx,g%ny,1] )) !,[1,1,1]
    call check( nf90_inq_varid(ncid,'SWDOWN', var_id)); call check(nf90_get_var(ncid, var_id, PPFD, [1,1,t]  ))!,[g%nx,g%ny,1] )) !,[1,1,1]
    call check( nf90_inq_varid(ncid,'PSFC'  , var_id)); call check(nf90_get_var(ncid, var_id,  PRE, [1,1,t]  ))!,[g%nx,g%ny,1] )) !,[1,1,1]
    call check( nf90_inq_varid(ncid,'Q2'    , var_id)); call check(nf90_get_var(ncid, var_id,  HUM, [1,1,t]  ))!,[g%nx,g%ny,1] )) !,[1,1,1]
    call check( nf90_inq_varid(ncid,'RAINNC', var_id)); call check(nf90_get_var(ncid, var_id, RAIN, [1,1,t]  ))!,[g%nx,g%ny,1] )) !,[1,1,1]
    call check( nf90_inq_varid(ncid,'LAI'   , var_id)); call check(nf90_get_var(ncid, var_id,  LAI, [1,1,t]  ))!,[g%nx,g%ny,1] )) !,[1,1,1]
    call check( nf90_inq_varid(ncid,'SMOIS' , var_id)); call check(nf90_get_var(ncid, var_id,SMOIS, [1,1,1,t])) !, [g%nx,g%ny,1,t] )) !,[1,1,1]
    call check( nf90_inq_varid(ncid,'TSLB'  , var_id)); call check(nf90_get_var(ncid, var_id,STEMP, [1,1,1,t])) !, [g%nx,g%ny,1,t] )) !,[1,1,1]
  call check(nf90_close(ncid))
              
  WIND=SQRT(U10*U10 + V10*V10)
  !Ground Incident Radiation [W m-2] to PPFD (Photosynthetic Photon Flux Density [W m-2])
  ppfd=ppfd*4.5*0.45 ! ppfd = par   * 4.5    !par to ppfd
                     ! par  = rgrnd * 0.45   !total rad to Photosyntetic Active Radiation (PAR)
end subroutine
 

subroutine get_daily_data(g,DDD)
  implicit none
  type(grid_type), intent(in) :: g
  character(len=3),intent(in) :: DDD                      
  integer :: ncid,time_dimid,time_len
  real, dimension(:,:,:), allocatable :: t,u,v,r,wnd
  
     print*,"Preparo inputs dinámicos diarios.."
     !from prepmegan:
     if (.not. allocated(fert)) then; allocate( fert(g%nx,g%ny));endif
     call check(nf90_open(trim(fert_file), nf90_write, ncid ))
           call check(   nf90_inq_varid(ncid,'FERT'//DDD, var_id )); call check( nf90_get_var(ncid, var_id , FERT ))
     call check(nf90_close(ncid))
   
     time_len=size(Times)

     allocate(t(g%nx,g%ny,time_len))
     allocate(u(g%nx,g%ny,time_len))
     allocate(v(g%nx,g%ny,time_len))
     allocate(r(g%nx,g%ny,time_len))
     allocate(wnd(g%nx,g%ny,time_len)) !windspeed

     !Daily variables:
     if (.not. allocated(ppfd_avg))  allocate(ppfd_avg(g%nx,g%ny) )!then;;endif 
     if (.not. allocated(tmp_avg) )  allocate( tmp_avg(g%nx,g%ny) )!then;;endif 
     if (.not. allocated(tmp_min ))  allocate( tmp_min(g%nx,g%ny) )!then;;endif 
     if (.not. allocated(tmp_max ))  allocate( tmp_max(g%nx,g%ny) )!then;;endif 
     if (.not. allocated(wind_max))  allocate(wind_max(g%nx,g%ny) )!then;;endif 
     ppfd_avg=0;tmp_avg=0;tmp_min=0;tmp_max=0;wind_max=0;
                  
     call check(nf90_open(trim(met_file), nf90_write, ncid ))
         call check( nf90_inq_varid(ncid,'U10'   , var_id)); call check(nf90_get_var(ncid, var_id, u)) !  , [1,1,t]))!,[g%nx,g%ny,1] )) !,[1,1,1]
         call check( nf90_inq_varid(ncid,'V10'   , var_id)); call check(nf90_get_var(ncid, var_id, v)) !  , [1,1,t]))!,[g%nx,g%ny,1] )) !,[1,1,1]
         call check( nf90_inq_varid(ncid,'T2'    , var_id)); call check(nf90_get_var(ncid, var_id, t)) !  , [1,1,t]))!,[g%nx,g%ny,1] )) !,[1,1,1]
         call check( nf90_inq_varid(ncid,'SWDOWN', var_id)); call check(nf90_get_var(ncid, var_id, r)) !  , [1,1,t]))!,[g%nx,g%ny,1] )) !,[1,1,1]
     call check(nf90_close(ncid)) 
     wnd=sqrt(u*u+v*v)
     tmp_min  = minval(t, dim=3,mask=t>0)
     tmp_max  = maxval(t, dim=3)
     wind_max = maxval(wnd,dim=3)
     tmp_avg  = sum(t, dim=3)/time_len
     ppfd_avg = sum(r, dim=3)/time_len

    deallocate(t,u,v,r,wnd)
end subroutine

subroutine get_monthly_data(g,MM)!,lai,ndep)
  implicit none
  type(grid_type)   :: g
  character(len=2)  :: MM

  print*,"  Preparo inputs dinámicos mensuales.."
  if (.not. allocated(ndep)) then; allocate(ndep(g%nx,g%ny));endif
  !if (.not. allocated(lai) ) then; allocate( lai(g%nx,g%ny));endif
  !call check(nf90_open( trim(lai_file), nf90_write, ncid ))
  !   call check( nf90_inq_varid(ncid,'LAI'//MM    , var_id )); call check( nf90_get_var(ncid, var_id , LAI  ))
  !call check(nf90_close(ncid))
  call check(nf90_open( trim(ndep_file), nf90_write, ncid ))
     call check( nf90_inq_varid(ncid,'NITROGEN'//MM, var_id )); call check( nf90_get_var(ncid, var_id , NDEP ))
  call check(nf90_close(ncid))
end subroutine


subroutine write_output_file(g,YYYY,DDD,MECHANISM)
  implicit none
  type(grid_type) , intent(in) :: g
  character(len=5), intent(in) :: MECHANISM
  character(len=4), intent(in) :: YYYY
  character(len=3), intent(in) :: DDD
  !local vars
  integer           :: ncid,t_dim_id,x_dim_id,y_dim_id,z_dim_id,s_dim_id,var_dim_id
  integer           :: k!,i,j
  character(len=25) :: out_file
  !Creo NetCDF file
   out_file="out_"//YYYY//"_"//DDD//"_"//trim(MECHANISM)//".nc"
   print*,"ESCRIBIENDO: ",out_file

   !Crear NetCDF
   call check(nf90_create(trim(out_file), NF90_CLOBBER, ncid))
     !Defino dimensiones
     call check(nf90_def_dim(ncid, "Time", 24     , t_dim_id       ))
     call check(nf90_def_dim(ncid, "x"   , g%nx   , x_dim_id       ))
     call check(nf90_def_dim(ncid, "y"   , g%ny   , y_dim_id       ))
     !Defino variables
     call check(nf90_def_var(ncid,"Times",  NF90_INT    , [t_dim_id], var_id))
     call check(nf90_put_att(ncid, var_id, "units"      , "seconds from file start date" ))
     call check(nf90_put_att(ncid, var_id, "var_desc"   , "date-time variable"      ))
     !print*,"shape buffer: ",shape(out_buffer)
     !Creo variables:
     do k=1,NMGNSPC !n_scon_spc !
       print*,"Especie: ",k,n_scon_spc,trim(mech_spc(k))
        call check( nf90_def_var(ncid, trim(mech_spc(k)) , NF90_FLOAT, [x_dim_id,y_dim_id,t_dim_id], var_id)   )
     end do
   call check(nf90_enddef(ncid))   !End NetCDF define mode

   !Abro NetCDF y guardo variables de salida
   call check(nf90_open(trim(out_file), nf90_write, ncid       ))
     do k=1,NMGNSPC !n_scon_spc !,NMGNSPC
       print*,"Especie:",trim(mech_spc(k))
       call check(nf90_inq_varid(ncid,trim(mech_spc(k)),var_id)); call check(nf90_put_var(ncid, var_id, out_buffer(:,:,k,:) ))
     enddo
   call check(nf90_close( ncid ))
end subroutine write_output_file
!0000000000000000000000000000000000000000000000000000000000000000000000000000000000000

subroutine megan_map(mechanism)

  implicit none
    character( 16 )    :: mechanism     ! mechanism name
    integer, parameter :: nmechs = 16   ! dimension for number of mechanisms considered
    integer ios,indx
    integer i

    TYPE MIOG_MECH_TYPE
         CHARACTER( 32 ) :: CHEMMECH
         CHARACTER( 16 ) :: MIOGMECH
    END TYPE MIOG_MECH_TYPE

    TYPE(MIOG_MECH_TYPE) :: MIOG_MECH_MAP( NMECHS ) = [       &
       miog_mech_type('CB05E51_AE6_AQ         ','CB05    '),  &
       miog_mech_type('CB05EH51_AE6_AQ        ','CB05    '),  &
       miog_mech_type('CB05MP51_AE6_AQ        ','CB05    '),  &
       miog_mech_type('CB05TUCL51_AE6_AQ      ','CB05    '),  &
       miog_mech_type('CB6R3_AE6_AQ           ','CB6     '),  &
       miog_mech_type('CB6MP_AE6_AQ           ','CB6     '),  &
       miog_mech_type('CB6R5HAP_AE7_AQ        ','CB6_ae7 '),  &
       miog_mech_type('CB6R3_AE7_AQ           ','CB6_ae7 '),  &
       miog_mech_type('CB6R5_AE7_AQ           ','CB6_ae7 '),  &
       miog_mech_type('CB6R5M_AE7_AQ          ','CB6_ae7 '),  &
       miog_mech_type('RACM2_AE6_AQ           ','RACM2   '),  &
       miog_mech_type('SAPRC07TC_AE6_AQ       ','SAPRC07T'),  &
       miog_mech_type('SAPRC07TIC_AE7I_AQ     ','SAPRC07T'),  &
       miog_mech_type('SAPRC07TIC_AE7I_AQKMT2 ','SAPRC07T'),  &
       miog_mech_type('CRACMM1_AQ             ','CRACMM  '),  &
       miog_mech_type('CRACMM1AMORE_AQ        ','CRACMM  ')   ]

       !INDX = INDEX1( MECHNAME, NMECHS, MIOG_MECH_MAP%CHEMMECH )
       !MECHANISM = MIOG_MECH_MAP( INDX )%MIOGMECH

       SELECT CASE ( TRIM(MECHANISM) )
         CASE ('SAPRC07')
           n_scon_spc = n_saprc07
           NMGNSPC = n_saprc07_spc
         CASE ('SAPRC07T')
           n_scon_spc = n_saprc07t
           NMGNSPC = n_saprc07t_spc
         CASE ('CB05')
           n_scon_spc = n_cb05
           NMGNSPC = n_cb05_spc
         CASE ('CB6')
           n_scon_spc = n_cb6  ! 145
           NMGNSPC = n_cb6_spc ! 34
         CASE ('RACM2')
           n_scon_spc = n_racm2
           NMGNSPC = n_racm2_spc
         CASE ('CB6_ae7')
           n_scon_spc = n_cb6_ae7
           NMGNSPC = n_cb6_ae7_spc
         CASE ('CRACMM')
           n_scon_spc = n_cracmm
           NMGNSPC = n_cracmm_spc
         CASE DEFAULT
           print*,"Mechanism," // TRIM( MECHANISM) // ", is not identified."
       ENDSELECT

   allocate(MEGAN_NAMES(NMGNSPC), stat=ios)
       allocate(spmh_map(n_scon_spc), stat=ios)
       allocate(mech_map(n_scon_spc), stat=ios)
       allocate(conv_fac(n_scon_spc), stat=ios)
       allocate(mech_spc(NMGNSPC )  , stat=ios)
       allocate(mech_mwt(NMGNSPC )  , stat=ios)

       SELECT CASE ( TRIM(MECHANISM) )
         CASE ('CB05')
           spmh_map(1:n_scon_spc) = spmh_map_cb05(1:n_scon_spc)
           mech_map(1:n_scon_spc) = mech_map_cb05(1:n_scon_spc)
           conv_fac(1:n_scon_spc) = conv_fac_cb05(1:n_scon_spc)
           mech_spc(1:NMGNSPC)    = mech_spc_cb05(1:NMGNSPC)
           mech_mwt(1:NMGNSPC)    = mech_mwt_cb05(1:NMGNSPC)
           MEGAN_NAMES(1:NMGNSPC) =      mech_spc(1:NMGNSPC)
         CASE ('CB6')
           spmh_map(1:n_scon_spc) = spmh_map_cb6(1:n_scon_spc)
           mech_map(1:n_scon_spc) = mech_map_cb6(1:n_scon_spc)
           conv_fac(1:n_scon_spc) = conv_fac_cb6(1:n_scon_spc)
           mech_spc(1:NMGNSPC)    = mech_spc_cb6(1:NMGNSPC)
           mech_mwt(1:NMGNSPC)    = mech_mwt_cb6(1:NMGNSPC)
           MEGAN_NAMES(1:NMGNSPC) =     mech_spc(1:NMGNSPC)
         CASE ('RACM2')
           spmh_map(1:n_scon_spc) = spmh_map_racm2(1:n_scon_spc)
           mech_map(1:n_scon_spc) = mech_map_racm2(1:n_scon_spc)
           conv_fac(1:n_scon_spc) = conv_fac_racm2(1:n_scon_spc)
           mech_spc(1:NMGNSPC)    = mech_spc_racm2(1:NMGNSPC)
           mech_mwt(1:NMGNSPC)    = mech_mwt_racm2(1:NMGNSPC)
           MEGAN_NAMES(1:NMGNSPC) =       mech_spc(1:NMGNSPC)
         CASE ('SAPRC07')
           spmh_map(1:n_scon_spc) = spmh_map_saprc07(1:n_scon_spc)
           mech_map(1:n_scon_spc) = mech_map_saprc07(1:n_scon_spc)
           conv_fac(1:n_scon_spc) = conv_fac_saprc07(1:n_scon_spc)
           mech_spc(1:NMGNSPC)    = mech_spc_saprc07(1:NMGNSPC)
           mech_mwt(1:NMGNSPC)    = mech_mwt_saprc07(1:NMGNSPC)
           MEGAN_NAMES(1:NMGNSPC) =         mech_spc(1:NMGNSPC)
         CASE ('SAPRC07T')
           spmh_map(1:n_scon_spc) = spmh_map_saprc07t(1:n_scon_spc)
           mech_map(1:n_scon_spc) = mech_map_saprc07t(1:n_scon_spc)
           conv_fac(1:n_scon_spc) = conv_fac_saprc07t(1:n_scon_spc)
           mech_spc(1:NMGNSPC)    = mech_spc_saprc07t(1:NMGNSPC)
           mech_mwt(1:NMGNSPC)    = mech_mwt_saprc07t(1:NMGNSPC)
           MEGAN_NAMES(1:NMGNSPC) =          mech_spc(1:NMGNSPC)
         CASE ('CB6_ae7')
           spmh_map(1:n_scon_spc) = spmh_map_cb6_ae7(1:n_scon_spc)
           mech_map(1:n_scon_spc) = mech_map_cb6_ae7(1:n_scon_spc)
           conv_fac(1:n_scon_spc) = conv_fac_cb6_ae7(1:n_scon_spc)
           mech_spc(1:NMGNSPC)    = mech_spc_cb6_ae7(1:NMGNSPC)
           mech_mwt(1:NMGNSPC)    = mech_mwt_cb6_ae7(1:NMGNSPC)
           MEGAN_NAMES(1:NMGNSPC) =         mech_spc(1:NMGNSPC)
         CASE ('CRACMM')
           spmh_map(1:n_scon_spc) = spmh_map_cracmm(1:n_scon_spc)
           mech_map(1:n_scon_spc) = mech_map_cracmm(1:n_scon_spc)
           conv_fac(1:n_scon_spc) = conv_fac_cracmm(1:n_scon_spc)
           mech_spc(1:NMGNSPC)    = mech_spc_cracmm(1:NMGNSPC)
           mech_mwt(1:NMGNSPC)    = mech_mwt_cracmm(1:NMGNSPC)
           MEGAN_NAMES(1:NMGNSPC) =        mech_spc(1:NMGNSPC)
         CASE DEFAULT
           print*,"Mapping for Mechanism,"//TRIM(MECHANISM)//", is unspecified."
       ENDSELECT
end subroutine megan_map

subroutine check(status)
  integer, intent(in) :: status
  if (status /= nf90_noerr) then
    write(*,*) nf90_strerror(status)
    stop 'netcdf error'
  end if
end subroutine check


!DEBUG:
subroutine write_diagnostic_file(g) !write input and intermediate data so i can check everything is allright
  implicit none
  type(grid_type) , intent(in) :: g
  character(len=20) :: diag_file
  integer           :: ncid,t_dim_id,x_dim_id,y_dim_id,z_dim_id,s_dim_id,var_dim_id
  integer           :: k!,i,j
  character(len=5),dimension(23) :: var_names, var_types, var_dimen

  diag_file="diag_"//YYYY//"_"//DDD//"_"//HH//".nc"

  print*,'   Escrbibiendo diag_file: ',diag_file

  data var_names( 1),var_types( 1),var_dimen( 1) / 'LAT  ', 'FLOAT', 'XY   '/
  data var_names( 2),var_types( 2),var_dimen( 2) / 'LON  ', 'FLOAT', 'XY   '/
  data var_names( 3),var_types( 3),var_dimen( 3) / 'LAI  ', 'FLOAT', 'XY   '/
  data var_names( 4),var_types( 4),var_dimen( 4) / 'NDEP ', 'FLOAT', 'XY   '/
  data var_names( 5),var_types( 5),var_dimen( 5) / 'FERT ', 'FLOAT', 'XY   '/
  data var_names( 6),var_types( 6),var_dimen( 6) / 'ARID ', 'INT  ', 'XY   '/
  data var_names( 7),var_types( 7),var_dimen( 7) / 'NARID', 'INT  ', 'XY   '/
  data var_names( 8),var_types( 8),var_dimen( 8) / 'LTYPE', 'INT  ', 'XY   '/
  data var_names( 9),var_types( 9),var_dimen( 9) / 'U10  ', 'FLOAT', 'XY   '/
  data var_names(10),var_types(10),var_dimen(10) / 'V10  ', 'FLOAT', 'XY   '/
  data var_names(11),var_types(11),var_dimen(11) / 'PRE  ', 'FLOAT', 'XY   '/
  data var_names(12),var_types(12),var_dimen(12) / 'TMP  ', 'FLOAT', 'XY   '/
  data var_names(13),var_types(13),var_dimen(13) / 'PPFD ', 'FLOAT', 'XY   '/
  data var_names(14),var_types(14),var_dimen(14) / 'HUM  ', 'FLOAT', 'XY   '/
  data var_names(15),var_types(15),var_dimen(15) / 'RAIN' , 'FLOAT', 'XY   '/
  data var_names(16),var_types(16),var_dimen(16) / 'STEMP', 'FLOAT', 'XY   '/
  data var_names(17),var_types(17),var_dimen(17) / 'SMOIS', 'FLOAT', 'XY   '/
  data var_names(18),var_types(18),var_dimen(18) / 'STYPE', 'INT  ', 'XY   '/
  data var_names(19),var_types(19),var_dimen(19) / 'T_MIN', 'FLOAT', 'XY   '/
  data var_names(20),var_types(20),var_dimen(20) / 'T_MAX', 'FLOAT', 'XY   '/
  data var_names(21),var_types(21),var_dimen(21) / 'W_MAX', 'FLOAT', 'XY   '/
  data var_names(22),var_types(22),var_dimen(22) / 'T_AVG', 'FLOAT', 'XY   '/
  data var_names(23),var_types(23),var_dimen(23) / 'R_AVG', 'FLOAT', 'XY   '/
  !data var_names( 3),var_types( 3),var_dimen( 3) / 'EF   ', 'FLOAT', 'XY   '/
  !data var_names( 4),var_types( 4),var_dimen( 4) / 'LDF  ', 'FLOAT', 'XY   '/
  !data var_names( 4),var_types( 4),var_dimen( 4) / 'CTS  ', 'FLOAT', 'XY   '/

  call check(nf90_create(trim(diag_file), NF90_CLOBBER, ncid))
    !Defino dimensiones
    call check(nf90_def_dim(ncid, "Time", 1      , t_dim_id       ))
    call check(nf90_def_dim(ncid, "x"   , g%nx   , x_dim_id       ))
    call check(nf90_def_dim(ncid, "y"   , g%ny   , y_dim_id       ))
    !call check(nf90_def_dim(ncid, "z"   , layers , z_dim_id       ))
    !call check(nf90_def_dim(ncid, "c"   , nclass , s_dim_id       ))
    !Defino variables
    call check(nf90_def_var(ncid,"Times",  NF90_INT    , [t_dim_id], var_id))
    call check(nf90_put_att(ncid, var_id, "units"      , "seconds from file start date" ))
    call check(nf90_put_att(ncid, var_id, "var_desc"   , "date-time variable"      ))
    !Creo variables:
    do k=1, size(var_names)
      if ( trim(var_types(k)) == "FLOAT" ) then
         if ( trim(var_dimen(k)) == 'XY' .or. trim(var_dimen(k)) == 'XYT') then
              call check( nf90_def_var(ncid, trim(var_names(k)) , NF90_FLOAT, [x_dim_id,y_dim_id], var_id)          )
         endif
      else if ( trim(var_types(k)) == "INT") then   
              if ( trim(var_dimen(k)) == 'XY' .or. trim(var_dimen(k)) == 'XYT'  ) then
              call check( nf90_def_var(ncid, trim(var_names(k)) , NF90_INT  , [x_dim_id,y_dim_id], var_id)          )
         endif
       endif
    end do
  call check(nf90_enddef(ncid))   !End NetCDF define mode

  call check(nf90_open(trim(diag_file), nf90_write, ncid ))
      call check(nf90_inq_varid(ncid,'LAT'  ,var_id )); call check(nf90_put_var(ncid, var_id, lat     ))      !print*, " LAT'    ";
      call check(nf90_inq_varid(ncid,'LON'  ,var_id )); call check(nf90_put_var(ncid, var_id, lon     ))      !print*, " LON'    ";
      !call check(nf90_inq_varid(ncid,'EF'   ,var_id )); call check(nf90_put_var(ncid, var_id, ef      ))     !print*, " EF'     ";
      !call check(nf90_inq_varid(ncid,'LDF'  ,var_id )); call check(nf90_put_var(ncid, var_id, ldf_in     ))  !print*, " LDF'    ";
      !call check(nf90_inq_varid(ncid,'CTS'  ,var_id )); call check(nf90_put_var(ncid, var_id, CTS     ))     !print*, " CTS'    ";
      call check(nf90_inq_varid(ncid,'LAI'  ,var_id )); call check(nf90_put_var(ncid, var_id, lai     ))      !print*, " LAI'    ";
      call check(nf90_inq_varid(ncid,'NDEP' ,var_id )); call check(nf90_put_var(ncid, var_id, ndep    ))      !print*, " NDEP'   ";
      call check(nf90_inq_varid(ncid,'FERT' ,var_id )); call check(nf90_put_var(ncid, var_id, fert    ))      !print*, " FERT'   ";
      call check(nf90_inq_varid(ncid,'ARID' ,var_id )); call check(nf90_put_var(ncid, var_id, arid    ))      !print*, " ARID'   ";
      call check(nf90_inq_varid(ncid,'NARID',var_id )); call check(nf90_put_var(ncid, var_id, non_arid))      !print*, " NARID'  ";
      call check(nf90_inq_varid(ncid,'LTYPE',var_id )); call check(nf90_put_var(ncid, var_id, landtype))      !print*, " LTYPE'  ";
      call check(nf90_inq_varid(ncid,'U10'  ,var_id )); call check(nf90_put_var(ncid, var_id, u10     ))      !print*, " U10'    ";
      call check(nf90_inq_varid(ncid,'V10'  ,var_id )); call check(nf90_put_var(ncid, var_id, v10     ))      !print*, " V10'    ";
      call check(nf90_inq_varid(ncid,'PRE'  ,var_id )); call check(nf90_put_var(ncid, var_id, pre     ))      !print*, " PRE'    ";
      call check(nf90_inq_varid(ncid,'TMP'  ,var_id )); call check(nf90_put_var(ncid, var_id, tmp     ))      !print*, " TMP'    ";
      call check(nf90_inq_varid(ncid,'PPFD' ,var_id )); call check(nf90_put_var(ncid, var_id, ppfd    ))      !print*, " PPFD'   ";
      call check(nf90_inq_varid(ncid,'HUM'  ,var_id )); call check(nf90_put_var(ncid, var_id, hum     ))      !print*, " HUM'    ";
      call check(nf90_inq_varid(ncid,'STEMP',var_id )); call check(nf90_put_var(ncid, var_id, stemp   ))      !print*, " STEMP'  ";
      call check(nf90_inq_varid(ncid,'SMOIS',var_id )); call check(nf90_put_var(ncid, var_id, smois   ))      !print*, " SMOIS'  ";
      call check(nf90_inq_varid(ncid,'STYPE',var_id )); call check(nf90_put_var(ncid, var_id, stype   ))      !print*, " STYPE'  ";
      call check(nf90_inq_varid(ncid,'T_MIN',var_id )); call check(nf90_put_var(ncid, var_id, tmp_min ))      !print*, " T_MIN'  ";
      call check(nf90_inq_varid(ncid,'T_MAX',var_id )); call check(nf90_put_var(ncid, var_id, tmp_max ))      !print*, " T_MAX'  ";
      call check(nf90_inq_varid(ncid,'T_AVG',var_id )); call check(nf90_put_var(ncid, var_id, tmp_avg ))      !print*, " T_AVG'  ";
      call check(nf90_inq_varid(ncid,'W_MAX',var_id )); call check(nf90_put_var(ncid, var_id,wind_max ))      !print*, " W_MAX'  ";
      call check(nf90_inq_varid(ncid,'R_AVG',var_id )); call check(nf90_put_var(ncid, var_id, ppfd_avg ))     !print*, " R_AVG'  ";
      call check(nf90_inq_varid(ncid,'RAIN' ,var_id )); call check(nf90_put_var(ncid, var_id, rain    ))      !print*, " RAIN'   ";
  call check(nf90_close( ncid ))
end subroutine


end program main