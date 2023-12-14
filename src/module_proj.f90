module PROJ_mod
 !use utils_mod
 implicit none
 INTEGER, PARAMETER :: ascii = selected_char_KIND ("ascii")
 INTEGER, PARAMETER :: ucs4  = selected_char_KIND ('ISO_10646')

 real, parameter :: R_EARTH = 6370000.
 real, parameter :: PI = 3.141592653589793
 real, parameter :: RAD2DEG = 180./PI
 real, parameter :: DEG2RAD = PI/180.

 !Este modulo asume que existen los siguientes objetos:
 type proj_type
    character(16)    :: pName     ! Projection name
    integer          :: typ       ! Integer code for projection TYPE (2=lcc, 6=stere, 7=merc)
    real             :: alp,bet,gam,xcent,ycent !proj parameters.
    real             :: p1,p2,p3,p4 !extra parameters to speed up calculation once p%typ is defined.
 end type proj_type

 type grid_type
     character(12)   :: gName        !grid-name
     integer         :: nx,ny,nz     !number of cells in x-y direction (ncols, nrows, nlevs)
     real            :: dx,dy        !x-y cell dimension (x_cell, y_cell)
     real            :: xmin,ymin,xmax,ymax,xc,yc
     real            :: lonmin,latmin,lonmax,latmax
 end type grid_type

 contains

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

   !hemi=SIGN(1.0,p%alp)
   k = 2.0*R_EARTH / (1 + SIN(p%alp*deg2rad)*SIN(lat*deg2rad) + COS(p%alp*deg2rad)*COS(lat*deg2rad)*COS( (lon-p%gam)*deg2rad ))

   x = k *   COS( lat *deg2rad) * SIN( (lon - p%gam)*deg2rad )
   y = k * ( COS(p%alp*deg2rad) * SIN(  lat         *deg2rad ) - SIN(p%alp*deg2rad)*COS(lat*deg2rad)*COS((lon-p%gam)*deg2rad) )
 end subroutine
!!END COORDINATE TRANFORMATION FUNCTIONS====================================

end module proj_mod
