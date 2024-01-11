module bdsnp_mod
  ! Adopted from CMAQ files and modified for MEGAN3.1 by Ling Huang
  ! 2019/07/15
  ! Adopted from MEGAN3.1 and modified for CMAQ 5.4 by Jeff Willison

  use netcdf
  implicit none
  ! Parameters:
  ! Value calculated by running the 2x2.5 GEOS-Chem model
  real*8,  parameter :: tau_months   = 6. ! this is the decay time for dep. n reservoir, fert is 4 months
  real*8,  parameter :: secperday    = 86400.d0
  real*8,  parameter :: dayspermonth = 30.
  real*8,  parameter :: tau_sec      = tau_months * dayspermonth * secperday
  
  ! New soil biomes based on Steinkamp et al., 2011
  integer, parameter :: nsoil    = 24
  ! Canopy wind extinction coefficients
  ! (cf. Yienger & Levy [1995], Sec 5), now a function of the MODIS/KOPPEN biometype (J.D. Maasakkers)
  real*8, parameter :: soilexc(nsoil) = [0.10, 0.50, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 1.00, 1.00, 1.00, 1.00, 2.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 2.00, 0.10, 2.00]
  ! Steinkamp and Lawrence, 2011 A values, wet biome coefficients
  ! for each of the 24 soil biomes [ng N/m2/s] 
  real*8, parameter :: a_biome(nsoil) = [0.00, 0.00, 0.00, 0.00, 0.00, 0.06, 0.09, 0.09, 0.01, 0.84, 0.84, 0.24, 0.42, 0.62, 0.03, 0.36, 0.36, 0.35, 1.66, 0.08, 0.44, 0.57, 0.57, 0.57]
  ! Saturation values from ASX. For PX see:
  ! Jacquemin B. and Noilhan J. (1990), Bound.-Layer Meteorol., 52, 93-134.

contains

subroutine bdsnp_nox( yyyy,ddd,hh,                  &
                      ncols,nrows,                  &  
                      L_DESID_DIAG, SOILM, SOILT, RSTYP,LAI, &
                      FERT,NDEP,ARID,NONARID,LANDTYPE        &
                      CFRAC,TEMP,PPFD,                       &
                      BDSNP_NO )
     implicit none
     !input vars:
      integer, intent(in) ::   yyyy,ddd,hh
      integer, intent(in) ::   ncols,nrows
      integer, intent(in) ::   CFRAC(ncols,nrows) ! cloud fraction
      integer, intent(in) ::   RSTYP(ncols,nrows) ! soil type
      real,    intent(in) ::   SOILM(ncols,nrows) ! soil moisture [m3/m3] (PX)
      real,    intent(in) ::   SOILT(ncols,nrows) ! soil temperature [K] (PX)
      real,    intent(in) ::     LAI(ncols,nrows) ! leaf area index (m2/m2)
      real,    intent(in) ::    FERT(ncols,nrows) ! "ng N m-2" already - reservoir
      real,    intent(in) ::    NDEP(ncols,nrows) ! 
      integer, intent(in) ::    ARID(ncols,nrows) ! 
      integer, intent(in) :: NONARID(ncols,nrows) ! 
      integer, intent(in) ::LANDTYPE(ncols,nrows) ! 
      logical, intent(in) :: L_DESID_DIAG
     !output var:
      REAL,    INTENT(OUT) :: BDSNP_NO (NCOLS,NROWS) ! output NO emissions in nanomol/m^2/s

     ! Land use files for BDSNP: both time independant in CMAQ sense and absolutely - e.g. fertilizer does not vary with year
     ! Gridded Canopy NOx reduction factor for BDSNP Soil NO calculations
      REAL, ALLOCATABLE, SAVE :: CRF   ( :,: )     ! 0-1
      REAL, ALLOCATABLE, SAVE :: CFRAC ( :,: )     ! 0-1
     ! --- diagnostic variables, can be removed in final version      
      REAL, ALLOCATABLE, SAVE :: CRFAVG     ( :,: )  ! 0-1
      REAL, ALLOCATABLE, SAVE :: PULSEAVG   ( :,: )  ! 1+
      REAL, ALLOCATABLE, SAVE :: BASESUM    ( :,: )  ! used in calculating the above two averages
      REAL, ALLOCATABLE, SAVE :: THETA_DIAG ( :,: )  ! diagnositc theta
      REAL, ALLOCATABLE, SAVE :: WET_DIAG   ( :,: )  ! diagnositc wet term
      REAL, ALLOCATABLE, SAVE :: TEMP_DIAG  ( :,: )  ! diagnositc temp term
      REAL, ALLOCATABLE, SAVE :: A_DIAG     ( :,: )  ! diagnositc biome base emissions term
      REAL, ALLOCATABLE, SAVE :: AFERT_DIAG ( :,: )  ! diagnositc fert emissions term
      REAL, ALLOCATABLE, SAVE :: NRES_FERT_DIAG ( :,: )  ! diagnositc nres fert
      REAL, ALLOCATABLE, SAVE :: NRES_DEP_DIAG ( :,: )  ! diagnositc nres  dep    
     ! ---------------------------------------------------------------------------      
      REAL,              SAVE :: EMPOL,EMPOLSUM, EMPOLAVG       ! use to check reasonableness of results, g/hr
      REAL,              SAVE :: TIMECHECK ! use to output CPU_TIME(TIMECHECK) to see if this section of code is running unreasonably long

      INTEGER          SOILCAT            ! soil category
      INTEGER          NSTEPS             ! run duration (HHMMSS)
      INTEGER, SAVE :: MSTEPS             ! run no. of steps
      INTEGER          I, J, K, R, C, L   ! counters
      !ACA TOCO RAMI
      integer :: ierr, ncid, col_dim_id,row_dim_id, var_id
      !ACA TOCO RAMI
      
!      REAL,    SAVE :: EFAC
       REAL            TEMP_TERM, WET_TERM, PULSE, A_FERT
       REAL            CRF_TERM
       REAL            SOILNOX, FERTDIAG
!      REAL             CFNO               ! NO correction factor
!      REAL             CFNOGRASS          ! NO correction factor for grasslands
!      REAL             TAIR               ! surface temperature
!      REAL             TSOI               ! soil temperature
       REAL             THETA              ! water filled pore space
       REAL             THETAPREV
!      REAL             CFNOWET,  THETA
!      REAL             FAC1, FAC2, FAC3, FAC4
!-----------------------------------------------------------------------
         ! we need to initialize and allocate:
         ! pulse
         ! length of dry period
         ! soil moisture of previous time step
         ! N reservoir, deposition only
         ! These values can be provided from a restart file. The restart file is 'timeless'.
         ! This means CMAQ isn't checking to see if the restart file is actually from
         ! the immediately prior timstep.

         ! Allocate memory for data and read
         if (.not. allocated(    FERT)) allocate(    FERT(ncols,nrows))
         if (.not. allocated(   CFRAC)) allocate(   CFRAC(ncols,nrows))
         if (.not. allocated(     CRF)) allocate(     CRF(ncols,nrows))
         if (.not. allocated(  CRFAVG)) allocate(  CRFAVG(ncols,nrows))
         if (.not. allocated(PULSEAVG)) allocate(PULSEAVG(ncols,nrows))
         if (.not. allocated( BASESUM)) allocate( BASESUM(ncols,nrows))
         ! ------ Diagnostics -----------------------------------
         if (.not. allocated(    theta_diag)) allocate(    theta_diag(ncols,nrows))
         if (.not. allocated(      wet_diag)) allocate(      wet_diag(ncols,nrows))
         if (.not. allocated(     temp_diag)) allocate(     temp_diag(ncols,nrows))
         if (.not. allocated(        a_diag)) allocate(        a_diag(ncols,nrows))
         if (.not. allocated(    afert_diag)) allocate(    afert_diag(ncols,nrows))
         if (.not. allocated(nres_fert_diag)) allocate(nres_fert_diag(ncols,nrows))
         if (.not. allocated( nres_dep_diag)) allocate( nres_dep_diag(ncols,nrows))
         !-----------------------------------------------------------------------------
         ! Initial run if the model hasnt been run before, otherwise use a restart file
         ! to determine DRYPERIOD, pulse state, prev. timestep soil moisture, and N reservoir.

         ! If initial run, initialize some variables, otherwise get them from file
         PFACTOR   = 1d0   ! array
         DRYPERIOD = 0.01  ! array initialized non-zero to avoid log(0)
         CFRAC     = 0. 
         SOILMPREV = 0d0   ! array
         FERT      = 0d0   ! array
         EMPOL     = 0d0
         EMPOLSUM  = 0d0
         EMPOLAVG  = 0d0
         BASESUM   = 0.0
         CRFAVG    = 0.0
         PULSEAVG  = 0.0
         NDEPRES   = 0d0   ! array

      !attempt to use steady state condition to reduce spin up time by setting dN/dt = 0
      !or NDEPRES = Dep rate * tau, the decay time
      NDEPRES(i,j) = NDEP(i,j)*TAU_SEC

   
      ! Fertilizer N reservoir already calculated and read from file, update deposition reservoir from dep rate
      DO R = 1, NROWS
        DO C = 1, NCOLS
            CALL GET_NDEPRES( TSTEP, NDEPRES(i,j), TAU_SEC, i,j, L_DESID_DIAG)
        END DO
      END DO

      ! Calculate temporal non-speciated soil NO emissions to EMPOL
      !     If False Dont do any calculations to test - replicate 0 output

      CALL GET_CANOPY_NOX(JDATE, JTIME, Met_Data%COSZEN,              & 
           MET_DATA%TEMP2, MET_DATA%RGRND, met_data%PRSFC,            & 
           LANDTYPE, LAI, Met_Data%SNOCOV, CFRAC, Met_Data%WSPD10, CRF)

      do i=1, nrows
      do j=1, ncols
         SOILNOX  = 0d0
         FERTDIAG = 0d0

         K = LANDTYPE( i,j ) !Skip LANDTYPE not present
                             ! Temperature-dependent term of soil NOx emissions
                             ! [unitless]
                             ! Uses PX soil temperature instead of inferring from air
                             ! temperature
         TEMP_TERM = SOILTEMP( SOILT(i,j) )

         ! Use THETA instead of boolean wet/dry climate
         SOILCAT = INT( RSTYP( i,j ) )
         IF ( SOILCAT .NE. 14) THEN !not water
            THETA =  SOILM( i,j ) / Grid_Data%WSAT( i,j )
            THETAPREV = SOILMPREV( i,j ) / Grid_Data%WSAT( i,j )
            ! Soil moisture scaling of soil NOx emissions
            WET_TERM = SOILWET( THETA , BDSNP_ARID( i,j ), BDSNP_NONARID( i,j ))
         ELSE
            WET_TERM = 0d0
            THETA = 0d0
            THETAPREV = 0d0
         END IF
         PULSE = PULSING( THETA, TSTEP, THETAPREV, PFACTOR( i,j ), DRYPERIOD( i,j ) )
         A_FERT = FERTADD( FERT( i,j ) , NDEPRES( i,j ) ) !adds reservoirs returns emission rates
         ! Canopy reduction factor
         CRF_TERM  = CRF( i,j )
         !  SOILNOX includes fertilizer
         SOILNOX   = ( A_BIOME(K) + A_FERT ) * ( TEMP_TERM * WET_TERM * PULSE ) * ( 1.d0 - CRF_TERM  )  
         FERTDIAG  = ( A_FERT ) * ( TEMP_TERM * WET_TERM * PULSE ) * ( 1.d0 - CRF_TERM  )
         !scale emissions
         EMPOL = SOILNOX *  3600.0 * 10.0**-9  ![ng N/m2/s] *  s/hr * g/ng
         BDSNP_NO(i,j) = SOILNOX / 14          ![nmol/m2/s]

         ! sum various quantities for daily averaging
         EMPOLSUM = EMPOLSUM + EMPOL
         BASESUM(i,j)  = BASESUM(i,j)  + ( A_BIOME(K) + A_FERT ) * ( TEMP_TERM * WET_TERM)            *  3600.0 * 10.0**-9 ![ng N/m2/s] *  s/hr * g/ng
         PULSEAVG(i,j) = PULSEAVG(i,j) + ( A_BIOME(K) + A_FERT ) * ( TEMP_TERM * WET_TERM   * PULSE ) *  3600.0 * 10.0**-9 ![ng N/m2/s] *  s/hr * g/ng
         CRFAVG(i,j)   = CRFAVG(i,j)   + ( A_BIOME(K) + A_FERT ) * ( TEMP_TERM * WET_TERM ) * ( 1.d0 - CRF_TERM  ) *  3600.0 * 10.0**-9 ![ng N/m2/s] *  s/hr * g/ng
        !--------- MORE DIAGNOSTICS  ---------------------------------
         A_DIAG( i,j ) = A_BIOME(K)
         AFERT_DIAG( i,j ) = A_FERT
         NRES_FERT_DIAG( i,j ) = FERT( i,j )
         NRES_DEP_DIAG( i,j )  = NDEPRES( i,j )
         WET_DIAG( i,j ) = WET_TERM
         THETA_DIAG( i,j ) = THETA
         TEMP_DIAG( i,j ) = TEMP_TERM
        ! -----------------------------------------------------
        END DO ! columns
      END DO ! rows

        print*,"creando debug_diag.nc.."
        ierr=nf90_create("debug_diag.nc", NF90_CLOBBER, ncid)
           ! Defino dimensiones
           ierr=nf90_def_dim(ncid, "COL" , ncols, col_dim_id)
           ierr=nf90_def_dim(ncid, "ROW" , nrows, row_dim_id)
           !Defino variables
           ierr=nf90_def_var(ncid,'BIOME',NF90_FLOAT, [col_dim_id,row_dim_id], var_id)
           ierr=nf90_def_var(ncid,'AFERT',NF90_FLOAT, [col_dim_id,row_dim_id], var_id)
           ierr=nf90_def_var(ncid,'FERT' ,NF90_FLOAT, [col_dim_id,row_dim_id], var_id)
           ierr=nf90_def_var(ncid,'NDEP' ,NF90_FLOAT, [col_dim_id,row_dim_id], var_id)
           ierr=nf90_def_var(ncid,'WET'  ,NF90_FLOAT, [col_dim_id,row_dim_id], var_id)
           ierr=nf90_def_var(ncid,'THETA',NF90_FLOAT, [col_dim_id,row_dim_id], var_id)
           ierr=nf90_def_var(ncid,'TEMP' ,NF90_FLOAT, [col_dim_id,row_dim_id], var_id)
           ierr=nf90_def_var(ncid,'STYP' ,NF90_INT  , [col_dim_id,row_dim_id], var_id)
           ierr=nf90_def_var(ncid,'WSAT' ,NF90_FLOAT, [col_dim_id,row_dim_id], var_id)
        ierr=nf90_enddef(ncid)
        print*,"Escribiendo variables.."
        ierr=nf90_open("debug_diag.nc", NF90_WRITE, ncid  )
           ierr=nf90_inq_varid(ncid,'BIOME', var_id ); ierr=nf90_put_var(ncid, var_id, A_DIAG        )     
           ierr=nf90_inq_varid(ncid,'AFERT', var_id ); ierr=nf90_put_var(ncid, var_id, AFERT_DIAG    )     
           ierr=nf90_inq_varid(ncid,'FERT' , var_id ); ierr=nf90_put_var(ncid, var_id, NRES_FERT_DIAG)     
           ierr=nf90_inq_varid(ncid,'NDEP' , var_id ); ierr=nf90_put_var(ncid, var_id, NRES_DEP_DIAG )
           ierr=nf90_inq_varid(ncid,'WET'  , var_id ); ierr=nf90_put_var(ncid, var_id, WET_DIAG      )
           ierr=nf90_inq_varid(ncid,'THETA', var_id ); ierr=nf90_put_var(ncid, var_id, THETA_DIAG    )
           ierr=nf90_inq_varid(ncid,'TEMP' , var_id ); ierr=nf90_put_var(ncid, var_id, TEMP_DIAG     )
           ierr=nf90_inq_varid(ncid,'STYP' , var_id ); ierr=nf90_put_var(ncid, var_id, RSTYP         )
           ierr=nf90_inq_varid(ncid,'WSAT' , var_id ); ierr=nf90_put_var(ncid, var_id, Grid_Data%WSAT)
        ierr=nf90_close(ncid)

      ELSE ! add things until it dies

         do r = 1, nrows
         do c = 1, ncols
            K = BDSNP_LANDTYPE( i,j ) !Skip LANDTYPE not present
            ! Temperature-dependent term of soil NOx emissions
            ! [unitless]
            ! Uses PX soil temperature instead of inferring from air
            ! temperature
            TEMP_TERM = SOILTEMP( SOILT(i,j) )
            ! Use THETA instead of boolean wet/dry climate
            SOILCAT = INT( RSTYP( i,j ) )
            IF ( SOILCAT .NE. 14) THEN !not water
               THETA = SOILM( i,j ) / Grid_Data%WSAT(i,j)
               THETAPREV =  SOILMPREV( i,j ) / Grid_Data%WSAT(i,j)
               ! Soil moisture scaling of soil NOx emissions
               WET_TERM = SOILWET( THETA , BDSNP_ARID( i,j ), BDSNP_NONARID( i,j ))
            ELSE
               WET_TERM = 0d0
               THETA = 0d0
               THETAPREV = 0d0
            END IF
            ! Cumulative multiplication factor (over baseline emissions)
            ! that accounts for soil pulsing
            PULSE = PULSING( THETA, TSTEP, THETAPREV, PFACTOR( i,j ), DRYPERIOD( i,j ) )
         end do
         end do
      END IF ! end do nothing test if

      SOILMPREV = SOILM !save soilM array to soilMprev for next time step
      EMPOLAVG = EMPOLSUM/FLOAT(NCOLS*NROWS)
      EMPOLSUM = 0d0 !array
  
      RETURN

      ! Avoid divide by zero over water where BASESUM = 0
      WHERE ( BASESUM .eq. 0) BASESUM = 1.0  
      CRFAVG   = CRFAVG/BASESUM
      PULSEAVG = PULSEAVG/BASESUM
      WHERE ( CRFAVG   .gt. 100 ) CRFAVG   = 0.0  
      WHERE ( PULSEAVG .gt. 100 ) PULSEAVG = 0.0  

!----------------------------------------------------------------------------------
contains

REAL FUNCTION PULSING( THETA, TSTEP, THETAPREV, PFACTOR, DRYPERIOD )!_____
! !DESCRIPTION: Function PULSING calculates the increase (or "pulse") of
!  soil NOx emission that happens after preciptiation falls on dry soil.
!\\
!\\
!  According to  Yan et al., [2005] , this pulsing process is thought to
!  be due to a release of inorganic nitrogen trapped on top of the dry soil
!  and a subsequent reactivation of water-stressed bacteria, which then
!  metabolize the excess nitrogen. This can happen in seasonally dry
!  grasslands and savannahs or over freshly fertilized fields.
!  Soil NOx emissions consist of baseline emissions plus discrete "pulsing"
!  episodes.  We follow the Yan et al., [2005] algorithm, where the pulse
!  (relative to the flux pre wetting) is determined by the antecedent dry
!  period, with a simple logarithmic relationship,
!
!  PFACTOR = 13.01 ln ( DRYPERIOD ) -  53.6
!
!  ,where PFACTOR is the magnitude of peak flux relative to prewetting flux,
!  and DRYPERIOD  is the length of the antecedent dry period in hours.
!
!  The pulse decays with
!
!  PFACTOR = PFACTOR * EXP( -0.068d0 * TSTEP(HOURS) )

         IMPLICIT NONE

         INTEGER, EXTERNAL       ::   TIME2SEC

! Function arguments:
         INTEGER, INTENT( IN )    :: TSTEP( 3 )        ! time step vector (HHMMSS)
         REAL,    INTENT( IN )    :: THETA, THETAPREV  ! only avilable if PX version
         REAL,    INTENT( INOUT ) :: DRYPERIOD
         REAL,    INTENT( INOUT ) :: PFACTOR
! Local Variables
         REAL MOISTDIFF
         REAL DTHOURS
         DTHOURS = TIME2SEC(TSTEP(2))/3600.0
         ! If soil moisture less than 0.3 and no pulse is taking place
         IF ( THETA < 0.3D0 .and. PFACTOR == 1.D0) THEN

            ! Get change in soil moisture since previous timestep
            MOISTDIFF = ( THETA - THETAPREV )

            ! If change in soil moisture is > 0.01 (rains)
            IF ( MOISTDIFF > 0.01 ) THEN

               !Initialize new pulse factor (dry period hours)
               PFACTOR = 13.01 * LOG( DRYPERIOD ) - 53.6

               ! If dry period < ~3 days then no pulse
               IF ( PFACTOR < 1.0 ) PFACTOR = 1.0

                  ! Reinitialize dry period
                  DRYPERIOD = 0.001

                ! If no rain (i.e.,  change in soil moisture is < 0.01)
               ELSE
                ! Add one timestep to dry period
                DRYPERIOD = DRYPERIOD + DTHOURS

            ENDIF

         ! If box is already pulsing , then decay pulse one timestep
         ELSEIF ( PFACTOR /= 1.d0) THEN

            ! Decay pulse
            PFACTOR   = PFACTOR * EXP( -0.068d0 * DTHOURS )

            ! Update dry period
            IF ( THETA < 0.3D0 ) DRYPERIOD = DRYPERIOD + DTHOURS

            ! If end of pulse
            IF ( PFACTOR < 1.d0 ) PFACTOR = 1.d0

         ENDIF
         PULSING = PFACTOR
         RETURN

END FUNCTION PULSING!_____
!---------------------------------------------------------------------------------------------
SUBROUTINE GET_NDEPRES( TSTEP, NDEPRES, TAU_SEC,i,j,L_DESID_DIAG)
    ! Get the deposition rate of the appropriate species for the appropriate timestep, add to reservoir and decay.
    ! Return reservoir amount.
    IMPLICIT NONE
    INTEGER, EXTERNAL     ::   TIME2SEC
    ! Function arguments:
    INTEGER, INTENT( IN ) :: TSTEP( 3 )        ! time step vector (HHMMSS)
    INTEGER, INTENT( IN ) :: C
    INTEGER, INTENT( IN ) :: R
    REAL*8,  INTENT( IN ) :: TAU_SEC
    LOGICAL, INTENT( IN ) :: L_DESID_DIAG
    REAL,    INTENT( INOUT ) :: NDEPRES

    REAL*8  :: C1,C2,TS_SEC ! a factor
    real NDEPTEMP
    !           check for negatives
    ! takes the NDEP and uses it to update NDEPRES before
    ! clearing it.

    !Do mass balance (see Intro to Atm Chem Chap. 3)
    !m(t) = m(0) * exp(-t/tau) + Source * tau * (1 - exp(-t/tau))
    TS_SEC = TIME2SEC(TSTEP(2))
    C1 = EXP( - TS_SEC / TAU_SEC)
    C2 = 1.d0 - C1
    NDEPTEMP = NDEPRES
    NDEPRES = NDEPRES*C1+NDEP(i,j)*TAU_SEC*C2
    !           check for negatives
     IF( NDEPRES < 0.0 ) THEN
         print*,"ERROR! NDEP<0.0";exit;
     END IF
    ! clear NDEP for use during next time step
     
    IF (.not. L_DESID_DIAG .and. MGN_ONLN_DEP) THEN 
                                 ! need this not to be zero'd 
                                 ! out on last time step
                                 ! and don't want it zero'd if using
                                 ! offline values
        NDEP(i,j) = 0.0 
    END IF
    RETURN
END SUBROUTINE GET_NDEPRES
! -----------------------------------------------------------------------------
SUBROUTINE GET_N_DEP( SPEC,DEP,i,j )

    IMPLICIT NONE
    CHARACTER( 8 ), INTENT( IN ) :: SPEC  !  dep species
    REAL,           INTENT( IN ) :: DEP   !  deposition rate in kg/ha/s 
    INTEGER,        INTENT( IN ) :: i,j
    REAL, PARAMETER              :: HAOM2   = 1.0e-4 ! ha/m^2 conversion
    REAL, PARAMETER              :: MWNH3   = 17.031 ! molecular weight of NH3
    REAL, PARAMETER              :: MWNH4   = 18.039 ! molecular weight of NH4
    REAL, PARAMETER              :: MWHNO3  = 63.013 ! molecular weight of HNO3
    REAL, PARAMETER              :: MWNO3   = 62.005 ! molecular weight of NO3
    REAL, PARAMETER              :: MWNO2   = 46.006 ! molecular weight of NO2
    REAL, PARAMETER              :: MWPAN   = 121.05 ! molecular weight of Peroxyacyl nitrate
    REAL, PARAMETER              :: MWN     = 14.007 ! molecular weight of Nitrogen
    REAL, PARAMETER              :: NGOKG   = 1.0e12 ! ng/kg conversion
    
    ! takes Kg/hectare/s and converts to ng N / m^2/s

    IF( INDEX(TRIM( SPEC ), 'NH3') .NE. 0 ) THEN
       NDEP( i,j ) = NDEP( i,j ) + DEP*HAOM2*NGOKG*MWN/MWNH3 
    ELSE IF( INDEX(TRIM( SPEC ), 'NH4') .NE. 0 ) THEN
       NDEP( i,j ) = NDEP( i,j ) + DEP*HAOM2*NGOKG*MWN/MWNH4
    ELSE IF( INDEX(TRIM( SPEC ), 'HNO3') .NE. 0 ) THEN
       NDEP( i,j ) = NDEP( i,j ) + DEP*HAOM2*NGOKG*MWN/MWHNO3
    ELSE IF( INDEX(TRIM( SPEC ), 'NO3') .NE. 0) THEN
       NDEP( i,j ) = NDEP( i,j ) + DEP*HAOM2*NGOKG*MWN/MWNO3
    ELSE IF( INDEX(TRIM( SPEC ), 'NO2') .NE. 0 ) THEN
       NDEP( i,j ) = NDEP( i,j ) + DEP*HAOM2*NGOKG*MWN/MWNO2
    ELSE IF( INDEX(TRIM( SPEC ), 'PAN') .NE. 0 ) THEN
       NDEP( i,j ) = NDEP( i,j ) + DEP*HAOM2*NGOKG*MWN/MWPAN
    Else
       MESG = 'Invalid Species Name in Get_N_Dep: "' // SPEC // '"'
       !CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT2 )
    END IF
    
    !IF( (DEP<0.0) .OR. (NDEP( i,j )<0.0) ) then
    !    !CALL M3EXIT( 'GET_N_DEP', 0, 0, MESG, XSTAT2 )
    !END if

    RETURN
END SUBROUTINE GET_N_DEP  
! -----------------------------------------------------------------------------
REAL FUNCTION SOILTEMP( SOILT )
    ! Calculate the soil temperature factor
    IMPLICIT NONE
    ! Function arguments:
    REAL, INTENT( IN )       :: SOILT !kelvin, soil temperature
    ! Local Variables
    REAl SOILTC                          !temperature in degrees celsius
    SOILTC = SOILT - 273.16

    IF ( SOILTC <= 0d0 ) THEN        ! No soil emissions if temp below freezing
       SOILTEMP = 0d0
    ELSE IF ( SOILTC >= 30.d0 ) then ! Caps temperature response at 30C
       SOILTC = 30.d0
       SOILTEMP =  EXP( 0.103 * SOILTC )
    ENDIF
    RETURN
END FUNCTION SOILTEMP
! ---------------------------------------------------------------------------------------------------------
REAL FUNCTION FERTADD( FERT , DEPN )
   ! Add fertilizer reservoir to deposition reservoir and create N driven
   ! emission factor
   IMPLICIT NONE
   ! Function arguments:
   REAL, INTENT( IN ) :: FERT !fertilizer reservoir [ngN/m2]
   REAL, INTENT( IN ) :: DEPN !deposition reservoir [ngN/m2]
   ! Local Variables
   REAL*8,  PARAMETER :: SECPERYEAR    = 86400.d0 * 365. ! Scale factor so that fertilizer emission = 1.8 Tg N/yr (Stehfest and Bouwman, 2006)
   ! before canopy reduction
   REAL*8, PARAMETER :: FERT_SCALE = 0.0068 ! [yr -1] ! Value calculated by running the 2x2.5 GEOS-Chem model (J.D. Maasakkers)
   
   FERTADD = FERT + DEPN
   FERTADD = FERTADD / SECPERYEAR * FERT_SCALE
   RETURN
END FUNCTION FERTADD
! -------------------------------------------------------------------------------------
REAL FUNCTION SOILWET( THETA , ARID, NONARID) ! Calculate the soil moisture factor
         IMPLICIT NONE
         REAL, INTENT( IN )       :: THETA !0-1 soil moisture
         INTEGER, INTENT( IN )    :: ARID !1 indicates arid cell
         INTEGER, INTENT( IN )    :: NONARID !1 indicates nonarid cell, if both 0 then
         IF ( ARID .EQ. 1 ) THEN  !ARID, Max poison at theta = .2
             SOILWET = 8.24*THETA*EXP(-12.5*THETA*THETA)
         ELSE IF (NONARID .EQ. 1 ) THEN !NONARID Max Poisson at theta =.3
             SOILWET = 5.5*THETA*EXP(-5.55*THETA*THETA)
         ELSE !neither arid nor nonarid, water or non-emitting cell
             SOILWET = 0.0
         END IF

         IF (SOILWET > 10) THEN
             SOILWET=1.0
         ENDIF
         RETURN
END FUNCTION SOILWET
! -------------------------------------------------------------------
SUBROUTINE GET_CANOPY_NOX(JDATE, JTIME, COSZEN, TASFC, SSOLAR, PRES, LANDTYPE, LAI, SNOCOV, CFRAC, WSPD, CRF)
      ! called tmpbeis, change called BDSNP, add K argument
      IMPLICIT NONE
      INTEGER, INTENT( IN )  :: JDATE             ! current simulation date (YYYYDDD)
      INTEGER, INTENT( IN )  :: JTIME             ! current simulation time (HHMMSS)

      REAL,    INTENT( IN ) :: COSZEN( NCOLS,NROWS )        ! cosine of zenith angle
      REAL,    INTENT( IN ) :: TASFC ( NCOLS,NROWS )        ! surface air temperature [K]
      REAL,    INTENT( IN ) :: SSOLAR( NCOLS,NROWS )        ! surface radiation [w/m**2]
      REAL,    INTENT( IN ) :: PRES  ( NCOLS,NROWS )        ! surface pressure [Pa]
      INTEGER, INTENT( IN ) :: LANDTYPE( NCOLS,NROWS )     ! the biome type in each cell
      REAL,    INTENT( IN ):: LAI  ( NCOLS,NROWS )        ! leaf area index (m2/m2)
      REAL,    INTENT( IN ):: SNOCOV  ( NCOLS,NROWS )        ! snow cover
      REAL,    INTENT( IN ):: CFRAC  ( NCOLS,NROWS )        ! cloud fraction
      REAL,    INTENT( IN ):: WSPD  ( NCOLS,NROWS )        ! cloud fraction
      REAL,    INTENT( OUT ):: CRF  ( NCOLS,NROWS )        ! outputs the canopy reduction factor
      !
      CHARACTER( 16 )  :: PNAME = 'CANOPY_NOX'  ! procedure name
      INTEGER          IOS                ! IO or memory allocation status
      CHARACTER( 256 ) :: MESG            ! message buffer
      ! Scalars
     INTEGER :: i,j, K, KK, MY_NCOLS, MY_NROWS
      REAL*8  :: F0,     HSTAR, XMW              
      REAL*8  :: DTMP1,  DTMP2, DTMP3,  DTMP4, GFACT, GFACI
      REAL*8  :: RT,     RAD0,  RIX,    RIXX,  RDC,   RLUXX
      REAL*8  :: RGSX,   RCLX,  TEMPK,  TEMPC, WINDSQR
      REAL*8 :: VFNEW
      LOGICAL, SAVE          :: FIRSTCANOPY = .TRUE. 

      ! Arrays
      REAL*8  :: RI  (24)       
      REAL*8  :: RLU (24)      
      REAL*8  :: RAC (24)      
      REAL*8  :: RGSS(24)     
      REAL*8  :: RGSO(24)     
      REAL*8  :: RCLS(24)     
      REAL*8  :: RCLO(24)
     ! !DEFINED PARAMETERS:
     INTEGER, PARAMETER :: SNIRI(24)    = [9999, 200, 9999, 9999, 9999, 9999, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 400, 400, 200, 200, 200, 9999, 200]
     INTEGER, PARAMETER :: SNIRLU(24)   = [9999, 9000, 9999, 9999, 9999, 9999, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 1000, 9000, 9000, 9000, 9000, 1000, 9000, 9999, 9000]
     INTEGER, PARAMETER :: SNIRAC(24)   = [0, 300, 0, 0, 0, 0, 100, 100, 100, 100, 100, 100, 100, 100, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 200, 100, 200]
     INTEGER, PARAMETER :: SNIRGSS(24)  = [0, 0, 100, 1000, 100, 1000, 350, 350, 350, 350, 350, 350, 350, 350, 500, 200, 500, 500, 500, 500, 200, 150, 400, 150]
     INTEGER, PARAMETER :: SNIRGSO(24)  = [2000, 1000, 3500, 400, 3500, 400, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 150, 300, 150]
     INTEGER, PARAMETER :: SNIRCLS(24)  = [9999, 2500, 9999, 9999, 9999, 9999, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 9999, 2000, 2000, 2000, 2000, 9999, 2000, 9999, 2000]
     INTEGER, PARAMETER :: SNIRCLO(24)  = [9999, 1000, 1000, 9999, 1000, 9999, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 9999, 1000, 1000, 1000, 1000, 9999, 1000, 9999, 1000]
     INTEGER, PARAMETER :: SNIVSMAX(24) = [10, 100, 100, 10, 100, 10, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
     REAL, PARAMETER :: DRYCOEFF(20)    = [-0.358, 3.02, 3.85, -0.0978, -3.66, 12.0, 0.252, -7.8, 0.226, 0.274, 1.14, -2.19, 0.261, -4.62, 0.685, -0.254, 4.37, -0.266, -0.159, -0.206 ]
     ! Canopy wind extinction coefficients
     ! (cf. Yienger & Levy [1995], Sec 5), now a function of the MODIS/KOPPEN biometype (J.D. Maasakkers)
     REAL*8,  PARAMETER :: SOILEXC(24)  = [ 0.10, 0.50, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 1.00, 1.00, 1.00, 1.00, 2.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 2.00, 0.10, 2.00]

      ! Molecular weight of water [kg]
      REAL*8, PARAMETER :: XMWH2O = 18d-3
      !
      ! Ventilation velocity for NOx, day & night values [m/s]
      REAL*8, PARAMETER :: VFDAY   = 1.0d-2
      REAL*8, PARAMETER :: VFNIGHT = 0.2d-2 
      REAL*8, PARAMETER :: PRESS  = 1.5d5

      ! Set physical parameters
      HSTAR = 0.01d0              ! Henry's law constant
      F0    = 0.1d0               ! Reactivity factor for biological oxidation 
      XMW   = 46d-3               ! Molecular wt of NO2 (kg)

      CRF = 0d0 ! array
      
      ! begin calculating canopy reduction factor
      DO i=1, NCOLS
      DO j=1, NROWS
          IF(LAI(i,j) > 0.0) THEN
            TEMPC = TASFC(i,j) - 273.15d0 ! convert kelvin to Celsius
      ! Compute bulk surface resistance for gases.    
         !                                  
         !  Adjust external surface resistances for temperature; 
         !  from Wesely [1989], expression given in text on p. 1296.        
            RT = 1000.0D0 * EXP( -TEMPC - 4.0d0 )
         !--------------------------------------------------------------
         ! Get surface resistances - loop over biome types K
         !
         ! The land types within each grid square are defined using the 
         ! Olson land-type database.  Each of the Olson land types is 
         ! assigned a corresponding "deposition land type" with 
         ! characteristic values of surface resistance components.  
         ! There are 74 Olson land-types but only 11 deposition 
         ! land-types (i.e., many of the Olson land types share the 
         ! same deposition characteristics).  Surface resistance 
         ! components for the "deposition land types" are from Wesely 
         ! [1989] except for tropical forests [Jacob and Wofsy, 1990] 
         ! and for tundra [Jacob et al., 1992].  All surface resistance 
         ! components are normalized to a leaf area index of unity.
         !--------------------------------------------------------------
        !Set biometype
            K = LANDTYPE( i,j )
            
            ! Set second loop variable to K to allow snow/ice correction
             KK = K

            ! If the surface is snow or ice, then set K=3
            IF ( SNOCOV(i,j) .EQ. 1 ) KK = 3

            !USE new MODIS/KOPPEN Biometypes to read data

            ! Read the internal resistance RI (minimum stomatal resistance 
            ! for water vapor, per unit area of leaf) from the IRI array; 
            ! a '9999' value means no deposition to stomata so we impose a 
            ! very large value for RI.
            RI(K) = DBLE( SNIRI(KK) )
            IF ( RI(K) >= 9999.D0 ) RI(K)= 1.D12
            
            ! Cuticular resistances IRLU read in from 'drydep.table'
            ! are per unit area of leaf; divide them by the leaf area index 
            ! to get a cuticular resistance for the bulk canopy.  If IRLU is 
            !'9999' it means there are no cuticular surfaces on which to 
            ! deposit so we impose a very large value for RLU.
            IF ( SNIRLU(KK) >= 9999 .OR. LAI(i,j) <= 0d0 ) THEN
               RLU(K)  = 1.D6
            ELSE
               RLU(K)= DBLE( SNIRLU(KK) ) / LAI(i,j) + RT
            ENDIF

            ! The following are the remaining resistances for the Wesely
            ! resistance-in-series model for a surface canopy
            ! (see Atmos. Environ. paper, Fig.1).  
            RAC(K)  = MAX( DBLE( SNIRAC(KK)  ),      1d0 )
            RGSS(K) = MAX( DBLE( SNIRGSS(KK) ) + RT, 1d0 )
            RGSO(K) = MAX( DBLE( SNIRGSO(KK) ) + RT, 1d0 ) 
            RCLS(K) =      DBLE( SNIRCLS(KK) ) + RT           
            RCLO(K) =      DBLE( SNIRCLO(KK) ) + RT 

            IF (  RAC(K) >= 9999.D0 ) RAC(K)  = 1d12
            IF ( RGSS(K) >= 9999.D0 ) RGSS(K) = 1d12
            IF ( RGSO(K) >= 9999.D0 ) RGSO(K) = 1d12
            IF ( RCLS(K) >= 9999.D0 ) RCLS(K) = 1d12         
            IF ( RCLO(K) >= 9999.D0 ) RCLO(K) = 1d12

            !-------------------------------------------------------------
            ! Adjust stomatal resistances for insolation and temperature:  
            ! 
            ! Temperature adjustment is from Wesely [1989], equation (3).
            ! 
            ! Light adjustment by the function BIOFIT is described by Wang 
            ! [1996].  It combines:
            !
            ! - Local dependence of stomal resistance on the intensity I 
            !   of light impinging the leaf; this is expressed as a 
            !   multiplicative factor I/(I+b) to the stomatal resistance 
            !   where b = 50 W m-2
            !   (equation (7) of Baldocchi et al. [1987])
            ! - Radiative transfer of direct and diffuse radiation in the 
            !   canopy using equations (12)-(16) from Guenther et al. 
            !   [1995]
            ! - Separate accounting of sunlit and shaded leaves using
            !   equation (12) of Guenther et al. [1995]
            ! - Partitioning of the radiation at the top of the canopy 
            !   into direct and diffuse components using a 
            !   parameterization to results from an atmospheric radiative 
            !   transfer model [Wang, 1996]
            !
            ! The dependent variables of the function BIOFIT are the leaf 
            ! area index (XYLAI), the cosine of zenith angle (SUNCOS) and 
            ! the fractional cloud cover (CFRAC).  The factor GFACI 
            ! integrates the light dependence over the canopy depth; so
            ! be scaled by LAI to yield a bulk canopy value because that's 
            ! already done in the GFACI formulation.
            !-------------------------------------------------------------

            ! Radiation @ sfc [W/m2]
            RAD0 = SSOLAR(i,j)
            
            ! Internal resistance
            RIX  = RI(K)

            ! Skip the following block if the resistance RIX is high
            IF ( RIX < 9999d0 ) THEN
               GFACT = 100.0D0

               IF ( TEMPC > 0.D0 .AND. TEMPC < 40.D0) THEN
                  GFACT = 400.D0 / TEMPC / ( 40.0D0 - TEMPC )
               ENDIF

               GFACI = 100.D0

               IF ( RAD0 > 0d0 .AND. LAI(i,j) > 0d0 ) THEN
                  GFACI= 1d0 / 
     &                   BIOFIT( DRYCOEFF,       LAI(i,j),
     &                           COSZEN(i,j), CFRAC(i,j)    )
               ENDIF
            
               RIX = RIX * GFACT * GFACI
            ENDIF
            
            ! Compute aerodynamic resistance to lower elements in lower 
            ! part of the canopy or structure, assuming level terrain - 
            ! equation (5) of Wesely [1989].                     
            RDC = 100.D0*(1.0D0+1000.0D0/(RAD0 + 10.D0))

            ! Loop over species; species-dependent corrections to resistances
            ! are from equations (6)-(9) of Wesely [1989].
            !
            ! NOTE: here we only consider NO2 (bmy, 6/22/09)
            RIXX   = RIX * DIFFG( TASFC(i,j), PRESS, XMWH2O ) /
     &                     DIFFG( TASFC(i,j), PRESS, XMW    )
     &             + 1.D0 / ( HSTAR/3000.D0 + 100.D0*F0  )

            RLUXX  = 1.D12

            IF ( RLU(K) < 9999.D0 ) THEN
               RLUXX = RLU(K) / ( HSTAR / 1.0D+05 + F0 )
            ENDIF
            
            ! To prevent virtually zero resistance to species with huge HSTAR, 
            ! such as HNO3, a minimum value of RLUXX needs to be set. 
            ! The rationality of the existence of such a minimum is 
            ! demonstrated by the observed relationship between Vd(NOy-NOx) 
            ! and Ustar in Munger et al.[1996]; Vd(HNO3) never exceeds 2 cm/s 
            ! in observations. The corresponding minimum resistance is 50 s/m.
            ! was introduced by J.Y. Liang on 7/9/95.
            RGSX = 1d0 / ( HSTAR/1d5/RGSS(K) + F0/RGSO(K) )
            RCLX = 1d0 / ( HSTAR/1d5/RCLS(K) + F0/RCLO(K) )

            ! Get the bulk surface resistance of the canopy
            ! from the network of resistances in parallel and in series 
            ! (Fig. 1 of Wesely [1989])
            DTMP1 = 1.D0 / RIXX
            DTMP2 = 1.D0 / RLUXX
            DTMP3 = 1.D0 / ( RAC(K) + RGSX )
            DTMP4 = 1.D0 / ( RDC      + RCLX )

            ! Save the within canopy depvel of NOx, used in calculating 
            ! the canopy reduction factor for soil emissions [1/s]
            CRF(i,j) = DTMP1 + DTMP2 + DTMP3 + DTMP4
            
            ! Pick proper ventilation velocity for day or night
            IF ( COSZEN( i,j ) > 0d0 ) THEN
               VFNEW = VFDAY              
            ELSE 
               VFNEW = VFNIGHT            
            ENDIF

      ! If the leaf area index and the bulk surface resistance
      ! of the canopy to NOx deposition are both nonzero ...
            IF (CRF(i,j) > 0d0 ) THEN

         ! Adjust the ventilation velocity.  
         ! NOTE: SOILEXC(21) is the canopy wind extinction 
         ! coefficient for the tropical rainforest biome.
              WINDSQR=WSPD(i,j)*WSPD(i,j)
              VFNEW    = (VFNEW * SQRT( WINDSQR/9d0 * 7d0/LAI(i,j)) * ( SOILEXC(21)  / SOILEXC(K) ))

              ! Soil canopy reduction factor
              CRF(i,j) = CRF(i,j) / ( CRF(i,j) + VFNEW )
         
            ELSE ! CRF < 0.0
     
              ! Otherwise set the soil canopy reduction factor to zero
              CRF(i,j) = 0d0

            END IF

            
          ELSE
            CRF(i,j) = 0.0
          END IF !lai check

        END DO !row loop
      END DO !col loop
            
END SUBROUTINE GET_CANOPY_NOX

FUNCTION DIFFG( TK, PRESS, XM ) RESULT( DIFF_G )
! !DESCRIPTION: Function DIFFG calculates the molecular diffusivity [m2/s] in 
!  air for a gas X of molecular weight XM [kg] at temperature TK [K] and 
!  pressure PRESS [Pa].
!  We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
!  radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
!  radius of air is given in a Table on p. 479 of Levine [1988].  The Table
!  also gives radii for some other molecules.  Rather than requesting the user
!  to supply a molecular radius we specify here a generic value of 2.E-10 m for
!  all molecules, which is good enough in terms of calculating the diffusivity
!  as long as molecule is not too big.
! !REVISION HISTORY:
!     22 Jun 2009 - R. Yantosca - Copied from "drydep_mod.f"
    implicit none
    REAL,   INTENT(IN) :: TK      ! Temperature [K]
    REAL*8, INTENT(IN) :: PRESS   ! Pressure [Pa]
    REAL*8, INTENT(IN) :: XM      ! Molecular weight of gas [kg]
    REAL*8             :: DIFF_G  ! Molecular diffusivity [m2/s]
      REAL*8             :: AIRDEN, Z, DIAM, FRPATH, SPEED            
      REAL*8, PARAMETER  :: XMAIR  = 28.8d-3 
      REAL*8, PARAMETER  :: RADAIR = 1.2d-10
      REAL*8, PARAMETER  :: PI     = 3.1415926535897932d0
      REAL*8, PARAMETER  :: RADX   = 1.5d-10
      REAL*8, PARAMETER  :: RGAS   = 8.32d0
      REAL*8, PARAMETER  :: AVOGAD = 6.023d23
      !=================================================================
      ! DIFFG begins here!
      !=================================================================
      ! Air density
      AIRDEN = ( PRESS * AVOGAD ) / ( RGAS * TK )
      ! DIAM is the collision diameter for gas X with air.
      DIAM   = RADX + RADAIR
      ! Calculate the mean free path for gas X in air: ! eq. 8.5 of Seinfeld [1986];
      Z      = XM  / XMAIR
      FRPATH = 1d0 /( PI * SQRT( 1d0 + Z ) * AIRDEN*( DIAM**2 ) )
      ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      SPEED  = SQRT( 8d0 * RGAS * TK / ( PI * XM ) )
      ! Calculate diffusion coefficient of gas X in air; ! eq. 8.9 of Seinfeld [1986]
      DIFF_G = ( 3d0 * PI / 32d0 ) * ( 1d0 + Z ) * FRPATH * SPEED
      
END FUNCTION DIFFG

SUBROUTINE SUNPARAM(X)
      IMPLICIT NONE
!===============================================
! the sequence is lai,suncos,cloud fraction
!===============================================
!  NN = number of variables (lai,suncos,cloud fraction)
      INTEGER NN
      PARAMETER(NN=3)
!  ND = scaling factor for each variable
      INTEGER ND(NN),I
      DATA ND /55,20,11/
!  X0 = maximum for each variable
      REAL X(NN),X0(NN),XLOW
      DATA X0 /11.,1.,1./

      DO I=1,NN
        X(I)=MIN(X(I),X0(I))
! XLOW = minimum for each variable
        IF (I.NE.3) THEN
          XLOW=X0(I)/REAL(ND(I))
        ELSE
          XLOW= 0.
        END IF
        X(I)=MAX(X(I),XLOW)
        X(I)=X(I)/X0(I)
      END DO

      RETURN
END SUBROUTINE SUNPARAM
      
REAL*8 FUNCTION BIOFIT(COEFF1,XLAI1,SUNCOS1,CFRAC1)
      IMPLICIT NONE

!===============================================
! Calculate the light correction
!===============================================
!* BIOFIT and SUNPARAM were written by Y.H. Wang.   
!-------------------------------------------------------------
! Adjust stomatal resistances for insolation and temperature:  
! 
! Temperature adjustment is from Wesely [1989], equation (3).
! 
! Light adjustment by the function BIOFIT is described by Wang 
! [1996].  It combines:
!
! - Local dependence of stomal resistance on the intensity I 
!   of light impinging the leaf; this is expressed as a 
!   multiplicative factor I/(I+b) to the stomatal resistance 
!   where b = 50 W m-2
!   (equation (7) of Baldocchi et al. [1987])
! - Radiative transfer of direct and diffuse radiation in the 
!   canopy using equations (12)-(16) from Guenther et al. 
!   [1995]
! - Separate accounting of sunlit and shaded leaves using
!   equation (12) of Guenther et al. [1995]
! - Partitioning of the radiation at the top of the canopy 
!   into direct and diffuse components using a 
!   parameterization to results from an atmospheric radiative 
!   transfer model [Wang, 1996]
!
! The dependent variables of the function BIOFIT are the leaf 
! area index (XYLAI), the cosine of zenith angle (SUNCOS) and 
! the fractional cloud cover (CFRAC).  The factor GFACI 
! integrates the light dependence over the canopy depth; so
! be scaled by LAI to yield a bulk canopy value because that's 
! already done in the GFACI formulation.
!*************************************************************
      INTEGER KK
      PARAMETER (KK=4)
      REAL COEFF1(20),TERM(KK),REALTERM(20)
      REAL XLAI1,SUNCOS1,CFRAC1
      INTEGER K,K1,K2,K3

      TERM(1)=1.
      TERM(2)=XLAI1
      TERM(3)=SUNCOS1
      TERM(4)=CFRAC1
      CALL SUNPARAM(TERM(2))
      K=0
      DO K3=1,KK
        DO K2=K3,KK
          DO K1=K2,KK
            K=K+1
            REALTERM(K)=TERM(K1)*TERM(K2)*TERM(K3)
          END DO
        END DO
      END DO
      BIOFIT=0
      DO K=1,20
        BIOFIT=BIOFIT+COEFF1(K)*REALTERM(K)
      END DO
      IF (BIOFIT.LT.0.1) BIOFIT=0.1

      RETURN
END FUNCTION BIOFIT
      
!  References:
!  ============================================================================
!  (1 ) Baldocchi, D.D., B.B. Hicks, and P. Camara, "A canopy stomatal
!        resistance model for gaseous deposition to vegetated surfaces",
!        Atmos. Environ. 21, 91-101, 1987.
!  (2 ) Brutsaert, W., "Evaporation into the Atmosphere", Reidel, 1982.
!  (3 ) Businger, J.A., et al., "Flux-profile relationships in the atmospheric 
!        surface layer", J. Atmos. Sci., 28, 181-189, 1971.
!  (4 ) Dwight, H.B., "Tables of integrals and other mathematical data",
!        MacMillan, 1957.
!  (5 ) Guenther, A., and 15 others, A global model of natural volatile
!         organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
!  (6 ) Hicks, B.B., and P.S. Liss, "Transfer of SO2 and other reactive
!        gases across the air-sea interface", Tellus, 28, 348-354, 1976.
!  (7 ) Jacob, D.J., and S.C. Wofsy, "Budgets of reactive nitrogen,
!        hydrocarbons, and ozone over the Amazon forest during the wet season",
!        J.  Geophys. Res., 95, 16737-16754, 1990.
!  (8 ) Jacob, D.J., et al, "Deposition of ozone to tundra", J. Geophys. Res., 
!        97, 16473-16479, 1992.
!  (9 ) Levine, I.N., "Physical Chemistry, 3rd ed.", McGraw-Hill, 
!        New York, 1988.
!  (10) Munger, J.W., et al, "Atmospheric deposition of reactive nitrogen 
!        oxides and ozone in a temperate deciduous forest and a sub-arctic 
!        woodland", J. Geophys. Res., in press, 1996.
!  (11) Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, "SO2, sulfate, 
!        and HNO3 deposition velocities computed using regional landuse and
!        meteorological data", Atmos. Environ., 20, 949-964, 1986.
!  (12) Wang, Y.H., paper in preparation, 1996.
!  (13) Wesely, M.L, "Improved parameterizations for surface resistance to
!        gaseous dry deposition in regional-scale numerical models", 
!        Environmental Protection Agency Report EPA/600/3-88/025,
!        Research Triangle Park (NC), 1988.
!  (14) Wesely, M. L., Parameterization of surface resistance to gaseous dry 
!        deposition in regional-scale numerical models.  Atmos. Environ., 23
!        1293-1304, 1989. 
!  (15) Price, H., L. Jaeglé, A. Rice, P. Quay, P.C. Novelli, R. Gammon, 
!        Global Budget of Molecular Hydrogen and its Deuterium Content: 
!        Constraints from Ground Station, Cruise, and Aircraft Observations,
!        submitted to J. Geophys. Res., 2007.      

end subroutine bdsnp_nox
end module bdsnp_mod
