module megan_NOX

!include bdsnp_mod

contains

subroutine megan_nox(yyyy,ddd,hh, &
           ncols,nrows,           &
           lat,                   &
           tmp,stmp,smois,styp,   &
           ctf, lai,              &
           no_emis                )

  implicit none
  logical :: bdsnp_megan=.false.

  !input variables:
  integer, intent(in) :: yyyy,ddd,hh
  integer, intent(in) :: ncols,nrows
  real   ,intent(in),dimension(ncols,nrows)   :: lat,long, smois,stmp,lai
  real   ,intent(in),dimension(ncols,nrows,6) :: ctf 
  integer,intent(in),dimension(ncols,nrows,)  :: styp
  !local variables:

        IF (BDSNP_MEGAN) THEN
           !CALL GET_DATE(JYEAR, JDAY, MM, DD)
           !CALL HRNOBDSNP( IDATE,ITIME,TSTEP,MM, L_DESID_DIAG,SOILM1,SOILT,SLTYP,LAIc, BDSNP_NO)
           !call BDSNP ( IDATE,ITIME,TSTEP,MM, L_DESID_DIAG,SOILM1,SOILT,SLTYP,LAIc, BDSNP_NO)
        ELSE
           CALL SOILNOX(IDATE,ITIME,NCOLS,NROWS,           &
                     TEMP,LSOIL,SLTYP, SOILM1, SOILT,      &
                     LAIc, LAT, PRECADJ,                   &
                     CFNO, CFNOG )

           do i = 1,ncols
           do j = 1,nrows
               CALL GROWSEASON(yyyy,ddd,lat(I,J),GDAY,GLEN)

               IF (GDAY .EQ. 0) THEN
                  GAMNO(I,J) = CFNOG(I,J) ! non growing season ! CFNOG for everywhere

               ELSE IF (GDAY .GT. 0 .AND. GDAY .LE. 366) THEN
                  ! growing season     ! CFNOG for everywhere except crops
                  TMO1 = 0.0; TMO2 = 0.0
                  DO I_CT=1,5
                    TMO1 = TMO1 + CTF(I,J,I_CT)
                    TMO2 = TMO2 + CTF(I,J,I_CT) * CFNOG(I,J)
                  ENDDO
                  ! CFNO for crops
                  TMO1 = TMO1 + CTF(I,J,6)
                  TMO2 = TMO2 + CTF(I,J,6) * CFNO(I,J)
                  IF (TMO1 .EQ. 0.0) THEN
                     GAMNO(I,J) = 0.0
                  ELSE
                     GAMNO(I,J) = TMO2 / TMO1
                  ENDIF
               ENDIF

           enddo  !ncols
           enddo  !nrows

        END IF ! YL or BDSNP
end subroutine



SUBROUTINE SOILNOX( yyyy,ddd,hh,                 &
                nrows,ncols,                     &
                TA, LSOIL, ISLTYP, SOILM, SOILT, &
                LAIc, LAT,                       &
                PRECADJ,                         &
                CFNO, CFNOG )
        implicit none
        integer, intent (in)  :: yyyy,ddd,hh,ncols,nrows   ! 
       !input vars
        REAL, INTENT (IN)  ::  TA       ( nrows,ncols)    !  air temperature (K)
        REAL, INTENT (IN)  ::  SOILM    ( nrows,ncols)    !  soil moisture (m3/m3)
        REAL, INTENT (IN)  ::  SOILT    ( nrows,ncols)    !  soil temperature (K)
        REAL, INTENT (IN)  ::  PRECADJ  ( nrows,ncols)    !  precip adjustment
        REAL, INTENT (IN)  ::  LAIc     ( nrows,ncols)    !  soil temperature (K)
        REAL, INTENT (IN)  ::  LAT      ( nrows,ncols)    !  Latitude
       !out vars
        REAL, INTENT (IN OUT)  ::  CFNO ( nrows,ncols)    !  NO correction factor
        REAL, INTENT (IN OUT)  ::  CFNOG( nrows,ncols)    !  NO correction factor for grass

        INTEGER, INTENT (IN)  ::  ISLTYP  ( NX, NY )    !  soil type

        LOGICAL, INTENT (IN) :: LSOIL              ! true: using PX version of MCIP

        !.........  Local ARRAYS
        ! Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
        !       PLEIM-XIU LAND-SURFACE AND PBL MODEL (PX-LSM)
        ! See JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
                INTEGER, PARAMETER :: MAXSTYPES = 16
        !        REAL, PARAMETER    :: SATURATION( MAXSTYPES )     =  (/   &       
        !                              0.395, 0.410, 0.435, 0.485,         &
        !                              0.451, 0.420, 0.477, 0.476,         &
        !                              0.426, 0.482, 0.482            /)       
        
        !.........  SCRATCH LOCAL VARIABLES and their descriptions:
        INTEGER       ::   R, C, L      ! counters
        INTEGER       ::   SOILCAT      ! soil category

        REAL          ::   CF           ! NO correction factor
        REAL          ::   CFG          ! NO correction factor for grasslands
        REAL          ::  TAIR         ! surface temperature
        REAL          ::   TSOI         ! soil temperature
        REAL          ::   CFNOWET, CFNODRY, RATIO, FAC1, FAC2
        REAL, PARAMETER ::  const1 = (1. / 3.0 )  * (1.0 / 30.0)
        REAL, PARAMETER ::  const2 =EXP(-0.103 * 30.0)
        CHARACTER(256)  MESG         ! message buffer

        CHARACTER(16) :: PROGNAME = 'SOILNOX'   !  program name
!***********************************************************************

!.....  Loop through cells
        do i = 1, nrows
        do j = 1, ncols
             TAIR = TA( C, R )  ! [ºK]
             !.......  Check max bounds for temperature
             IF (TAIR > 315.0 ) THEN; TAIR = 315.0; END IF; IF( TAIR > 303.00 ) TAIR = 303.00
             IF ( TAIR > 268.8690 ) THEN
                 CFG = EXP( 0.04686 * TAIR - 14.30579 ) ! grass (from BEIS2)
             ELSE
                 CFG = 0.0
             END IF
             CFNOG(i,j) = CFG
             !   pre calculate common factors
             FAC2 = const2
             !.......  CFNO
             IF( .NOT. LSOIL ) THEN
             ! no soil
                TSOI = 0.72 * TAIR + 82.28
                IF (TSOI <= 273.16) TSOI = 273.16; IF (TSOI >= 303.16) TSOI = 303.16
                FAC1 = (TSOI- 273.16)
                CFNODRY = const1 * FAC1  ! see YL 1995 Equa 9a p. 11452
                IF (TSOI <= 283.16) THEN         ! linear cold case
                    CFNOWET =  FAC1 * FAC2 * 0.28 ! see YL 1995 Equ 7b
                ELSE                             ! exponential case
                    CFNOWET = EXP(0.103 * FAC1) *  FAC2
                END IF
                CF = 0.5 * CFNOWET + 0.5 * CFNODRY
             ELSE
             ! soil
                TSOI = SOILT( i,j )
                IF (TSOI <= 273.16) TSOI = 273.16
                IF (TSOI >= 303.16) TSOI = 303.16
                FAC1 = (TSOI- 273.16)
                CFNODRY = const1 * FAC1  ! see YL 1995 Equa 9a p. 11452
                IF (TSOI <= 283.16) THEN         ! linear cold case
                   CFNOWET = FAC1 * FAC2 * 0.28 ! see YL 1995 Equ 7b
                ELSE                             ! exponential case
                   CFNOWET = EXP(0.103 * FAC1 ) * FAC2
                END IF
                SOILCAT = INT( ISLTYP( i,j ) )
                IF( SOILCAT > 0 .AND. SOILCAT <= MAXSTYPES ) THEN
                   IF(Grid_Data%WSAT(i,j) .eq. 0) THEN
                    ! first ldesid diag call. Do nothing.
                    CF = 0.
                   ELSE
                    RATIO = SOILM( i,j ) / Grid_Data%WSAT( i,j )
                    CF = RATIO*CFNOWET + (1.0 - RATIO ) * CFNODRY
                   END IF
                ELSE
                    CF = 0.0
                END IF
             END IF  ! Endif LSOIL
             CFNO(i,j) = CF * FERTLZ_ADJ( JDATE, LAT(i,j) ) * VEG_ADJ( LAIc(i,j) ) * PRECADJ(i,j)
             if(cfno(i,j) .lt. 0) then; cfno(i,j) = 0; end if

        end do  ! loop over columns
        end do  ! loop over rows

 contains
   real function fertlz_adj( date, lat , gday, glen)
   !    This internal function computes a fertilizer adjustment factor
   !    for the given date in yyyyddd format. If it is not growing 
   !    season, the adjustment factor is 0; otherwise, it ranges from
   !    0.0 to 1.0.
        IMPLICIT NONE
   
        INTEGER, INTENT(IN) :: DATE
        REAL,    INTENT(IN) :: LAT
        INTEGER  GDAY, GLEN
   
        !CALL GROWSEASON( DATE, LAT, GDAY, GLEN )
        FERTLZ_ADJ = 0. !INITIALIZE
        IF( GDAY == 0 ) THEN
            FERTLZ_ADJ = 0.
        ELSE IF( GDAY >= 1 .AND. GDAY < 30 ) THEN
            ! first month of growing season
            FERTLZ_ADJ = 1.
        ELSE IF( GDAY >= 30 .AND. GDAY <= 366) THEN
            ! later month of growing season
            FERTLZ_ADJ = 1. + 30. / FLOAT(GLEN) - FLOAT(GDAY) / FLOAT(GLEN)
        ENDIF
        RETURN
   end function fertlz_adj
   
   real function veg_adj( lai )
   ! This internal function computes a vegetation adjustment factor based on LAIv.  See Yienger and Levy 1995
   ! VEG_ADJ = (EXP(-0.24*LAIv)+EXP(-0.0525*LAIv))*0.5 
      implicit none
      real,    intent(in) :: lai
      !veg_adj = 0.0
      veg_adj = (exp(-0.24*lai)+exp(-0.0525*lai))*0.5
      return
   end function veg_adj



subroutine growseason (year, jday,lat, gday, glen )                 
!!MEJORAR ESTA FUNCIÓN!!
   !  This internal function computes the day of the growing season
   !  corresponding to the given date in yyyyddd format.
      implicit none                                  !   NOTE: The use of "julian Day" to describe the day of tHE year is
      integer, intent(in :: year,jday,lat            !     technically incorrect. 
      integer, intent(out) :: gday, glen             ! The Julian Day Number (JDN) is the integer assigned to a whole solar 
                                                     ! day in the Julian day count starting from noon Universal time, with 
                                                     ! Julian day number 0 assigned to the day starting at noon on Monday, 
      INTEGER            :: gseason_start            ! January 1, 4713 BCE, proleptic Julian calendar (November 24, 4714 BCE, 
      INTEGER            :: gseason_end              ! in the proleptic Gregorian calendar), a date at which three 
                                                     ! multi-year cycles started (which are: Indiction, Solar, and Lunar cycles)
                                                     !  For example for January 1st, 2000 CE  at 00:00:00.0 UT1 the Julian Day 
                                                     !  is 2451544.500000 according to  the U.S. Naval Observatory.
      integer :: extra_day=0
      ! contemplo años bisiestos
      IF ( mod(year,4)   .eq. 0) extra_day=1 !LEAP = .TRUE.      
      IF ( mod(year,100) .eq. 0) extra_day=0 !LEAP = .FALSE.      
      IF ( mod(year,400) .eq. 0) extra_day=1 !LEAP = .TRUE.      

      IF ( LAT .LE. 23.0 .AND. LAT .GE. -23.0 ) THEN !tropical regions
         GSEASON_START = 001; GSEASON_END   = 355 + extra_day
         GDAY = jday        ; GLEN = 355 + extra_day
      ELSE IF ( LAT .LT. -23.0 ) THEN      !southern hemisphere (non-tropical)
         IF ( LAT .LT. -60.0 ) THEN        !  antartic/austral: no growing
            GDAY = 0;  GLEN = 0
         ELSE                              !  mid-latitude: NOV, DEC, JAN-MAY
            IF (JDAY .GE. 305 .AND. JDAY .LE. 366 ) THEN
              GSEASON_START = 305 + extra_day
              GSEASON_END   = 365 + extra_day
              GDAY = JDAY - GSEASON_START + 1
            ELSE IF (JDAY .GE. 001 .AND. JDAY .LE. 151 ) THEN
              GSEASON_START = 001;  GSEASON_END   = 151 + extra_day
              GDAY = JDAY - GSEASON_START + 1 + 61
            ELSE
              GDAY = 0
            ENDIF
            GLEN = 30 + 31 + 151 !G2J(YEAR,0531) !- G2J(YEAR,0101) + 1
         ENDIF
      ELSE IF ( LAT .GT. 23.0 ) THEN      ! northern hemisphere (non-tropical)
         IF ( LAT .GT. 65.0 ) THEN        !   arctic/boreal, no growing season
            GDAY = 0; GLEN = 0
         ELSE                              ! northern hemisphere temperate: APR-OCT
            GSEASON_START = 089+extra_day
            GSEASON_END   = 226+extra_day
            !GSEASON_START = INT( (LAT-23.0) * 4.5 )
            !GSEASON_END   = GSJULIAN_END - INT( (LAT-23.0) * 3.3 )
            GDAY=jday-GSEASON_START; GLEN = GSJULIAN_END - GSJULIAN_START + 1
         ENDIF
      ENDIF
     RETURN
END SUBROUTINE GROWSEASON

















end subroutine soilnox
