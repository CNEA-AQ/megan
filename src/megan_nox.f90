module nox_mod

   !include bdsnp_mod
   include 'tables/LSM.EXT'

contains

subroutine megan_nox(yyyy,ddd,hh,  &
           ncols,nrows,            &
           lat,                    &
           tmp,rain,               &
           lsm,styp,stemp,smois,   &
           ctf, lai,               &
           no_emis                 )

  implicit none
  !input variables:
  integer, intent(in)  :: yyyy,ddd,hh
  integer, intent(in)  :: ncols,nrows
  real   , intent(in), dimension(ncols,nrows)   :: lat,tmp,stemp,smois,rain,lai
  integer, intent(in), dimension(ncols,nrows)   :: styp
  real   , intent(in), dimension(ncols,nrows,6) :: ctf 
  character(len=4), intent(in) :: lsm          !land surface model 

  !output variables:
  real, intent(inout) :: NO_EMIS(ncols,nrows)

  !local variables:
  real    :: CFNO, CFNOG                  !  NO correction factor !  NO correction factor (for grass)
  real    :: CFNOWET,CFNODRY,CF,CFG
  real    :: TAIR,SMOI,TSOI,PREC_ADJ,LATI,LAIv
  integer :: ISTYP,MAXSTYPES
  real    :: CANF(6)
  real, allocatable :: wsat(:)            ! ver como calcular !
  real :: fac1,fac2,tmo1,tmo2,ratio
  integer :: gday,glen
  logical :: have_soil_fields=.true.

  integer :: i,j,k,i_ct
  
  select case (LSM)
         case ('NOAH' )
            allocate(wsat(size(wsat_noah))) ; wsat=wsat_noah;
            !allocate(wwlt(size(wwlt_noah))); wwlt=wwlt_noah;
         case ('JN90' )
            allocate(wsat(size(wsat_px_wrfv4p))); wsat=wsat_px_wrfv4p;
            !allocate(wwlt(size(wwlt_px_wrfv4p))); wwlt=wwlt_px_wrfv4p;
         case DEFAULT
            allocate(wsat(size(wsat_px_wrfv4p))); wsat=wsat_px_wrfv4p;
  end select

  do i = 1,ncols
  do j = 1,nrows

      !>>SOILNOX:
       TAIR  =   tmp(i,j)  !  air temperature (K)
       SMOI  = smois(i,j)  !  soil moisture (m3/m3)
       TSOI  = stemp(i,j)  !  soil temperature (K)
     PREC_ADJ= precipfact(rain(i,j))!  precip adjustment factor
       LATI  = lat(i,j)    !  latitude
       ISTYP = styp(i,j)   !  soil type
       CANF  = ctf(i,j,:)
       LAIv  = lai(i,j)

       !Check max bounds for temperature
       !IF (TAIR > 315.0 ) TAIR = 315.0  !non-sense
       IF( TAIR > 303.00 ) TAIR = 303.00 !

       !Calculate CFG:
       IF ( TAIR > 268.8690 ) THEN
           CFG = EXP( 0.04686 * TAIR - 14.30579 ) ! grass (from BEIS2)
       ELSE
           CFG = 0.0
       END IF
       CFNOG = CFG
       !   pre calculate common factors
       FAC1 = (TSOI- 273.16)
       FAC2 =  0.04550195 !const2=exp(-0.103 * 30.0)

       !Calculate CFNO
       IF( .NOT. have_soil_fields ) THEN !If soil fields aren't available (soil temp, wsat)
            ! no soil fields
            TSOI = 0.72 * TAIR + 82.28
            IF (TSOI <= 273.16) TSOI = 273.16
            IF (TSOI >= 303.16) TSOI = 303.16

            CFNODRY = 0.01111111 * FAC1        !  (see YL 1995 Eq. 9a p. 11452)
            IF (TSOI <= 283.16) THEN           ! linear cold case
                CFNOWET =  FAC1 * FAC2 * 0.28  !  (see YL 1995 Eq. 7b)
            ELSE                               ! exponential case
                CFNOWET = EXP(0.103 * FAC1) *  FAC2
            END IF
            CF = 0.5 * CFNOWET + 0.5 * CFNODRY
       ELSE
            ! soil fields available
            IF (TSOI <= 273.16) TSOI = 273.16 
            IF (TSOI >= 303.16) TSOI = 303.16
            CFNODRY = 0.01111111 * FAC1       !   (see YL 1995 Eq. 9a p. 11452)
            IF (TSOI <= 283.16) THEN          ! linear cold case
               CFNOWET = FAC1 * FAC2 * 0.28   !   (see YL 1995 Eq. 7b)
            ELSE                              ! exponential case
               CFNOWET = EXP(0.103 * FAC1 ) * FAC2
            END IF
            
            IF( ISTYP > 0 .AND. ISTYP <= MAXSTYPES ) THEN
               IF( WSAT(ISTYP) .eq. 0 ) THEN ! first ldesid diag call. Do nothing.
                 CF = 0.
               ELSE
                 RATIO = SMOI / WSAT(ISTYP)
                 CF = RATIO * CFNOWET + (1.0 - RATIO ) * CFNODRY
               END IF
            ELSE
               CF = 0.0
            END IF
       END IF  ! endif have_soil_fields
                                                                                        
       CFNO = CF * FERTLZ_ADJ(gday,glen) * VEG_ADJ(laiv) * PREC_ADJ
       if( CFNO .lt. 0 ) CFNO = 0
      !>>ENDSOILNOX
      
      call growseason(yyyy,ddd,LATI,gday,glen)

       IF (GDAY .EQ. 0) THEN                          ! non growing season ! CFNOG for everywhere
          NO_EMIS(i,j) = CFNOG 

       ELSE IF (GDAY .GT. 0 .AND. GDAY .LE. 366) THEN ! growing season     ! CFNOG for everywhere except crops
          TMO1 = 0.0; TMO2 = 0.0
          DO I_CT=1,5
            TMO1 = TMO1 + CANF(i_ct)
            TMO2 = TMO2 + CANF(i_ct) * CFNOG
          ENDDO
          ! CFNO for crops
          TMO1 = TMO1 + CANF(6)
          TMO2 = TMO2 + CANF(6) * CFNO
          IF (TMO1 .EQ. 0.0) THEN
             NO_EMIS(i,j) = 0.0
          ELSE
             NO_EMIS(i,j) = TMO2 / TMO1
          ENDIF
       ENDIF

  enddo  !ncols
  enddo  !nrows


 contains

 real function fertlz_adj(gday, glen) !    Computes fertilizer adjustment factor based ond growdate      
                                      !    If it is not growing season, the adjustment factor is 0; otherwise, it ranges from 0.0 to 1.0.
    implicit none
    integer, intent(in) ::  gday, glen
 
    if( gday == 0 ) then
        fertlz_adj = 0.
    else if( gday >= 1  .and. gday <  30 ) then      ! first month of growing season
        fertlz_adj = 1.
    else if( gday >= 30 .and. gday <= 366) then    ! later month of growing season
        fertlz_adj = 1. + 30. / float(glen) - float(gday) / float(glen)
    endif
    return
 end function fertlz_adj
   
 real function veg_adj( lai )      ! This internal function computes a vegetation adjustment factor based on LAIv.
    implicit none                  !  See Yienger and Levy 1995
    real,    intent(in) :: lai
    veg_adj = (exp(-0.24*lai)+exp(-0.0525*lai))*0.5
    return
 end function veg_adj


subroutine growseason (year,jday,lat,gday,glen)!!MEJORAR ESTA FUNCIÓN!!
   !  This internal function computes the day of the growing season
   !  corresponding to the given date in yyyyddd format.
    implicit none                          ! NOTE: The use of "julian Day" to describe the day of tHE year is
    integer, intent(in)  :: year,jday      !   technically incorrect. 
    integer, intent(out) :: gday, glen     ! The Julian Day Number (JDN) is the integer assigned to a whole solar 
    real,    intent(in)  :: lat            ! day in the Julian day count starting from noon Universal time, with 
                                           ! Julian day number 0 assigned to the day starting at noon on Monday, 
    integer            :: gseason_start    ! January 1, 4713 BCE, proleptic Julian calendar (November 24, 4714 BCE, 
    integer            :: gseason_end      ! in the proleptic Gregorian calendar), a date at which three 
                                           ! multi-year cycles started (which are: Indiction, Solar, and Lunar cycles)
                                           ! For example for January 1st, 2000 CE  at 00:00:00.0 UT1 the Julian Day 
                                           ! is 2451544.500000 according to  the U.S. Naval Observatory.
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
          GLEN = 30 + 31 + 151          !G2J(YEAR,0531) !- G2J(YEAR,0101) + 1
       ENDIF
    ELSE IF ( LAT .GT. 23.0 ) THEN      ! northern hemisphere (non-tropical)
       IF ( LAT .GT. 65.0 ) THEN        !   arctic/boreal, no growing season
          GDAY = 0; GLEN = 0
       ELSE                              ! northern hemisphere temperate: APR-OCT
          GSEASON_START = 089+extra_day  !GSEASON_START = INT( (LAT-23.0) * 4.5 )
          GSEASON_END   = 226+extra_day  !GSEASON_END   = GSJULIAN_END - INT( (LAT-23.0) * 3.3 )

          GDAY=jday-GSEASON_START; GLEN = GSEASON_END - GSEASON_START + 1
       ENDIF
    ENDIF
    RETURN
end subroutine growseason

 REAL FUNCTION PRECIPFACT(RRATE)!, JDATE, JTIME, ADATE, ATIME )
 ! This internal function computes a precipitation adjustment
 ! factor from YL 1995 based on a rain rate. The pulse type is
 ! and integer ranging from 0 to 3 indicating the type of rainfall rate.
   IMPLICIT NONE
   REAL   , intent (in) :: rrate !rainfall rate
   !...  Function arguments
   INTEGER :: PULSETYPE
   integer ::  HRDIFF = 1.0 !SECSDIFF( ADATE, ATIME, JDATE, JTIME ) / 3600.
   !pulse type
   IF( RRATE < 0.1 ) THEN
       PULSETYPE = 0
   ELSE IF( RRATE < 0.5 ) THEN
       PULSETYPE = 1
   ELSE IF( RRATE < 1.5 ) THEN
       PULSETYPE = 2
   ELSE
       PULSETYPE = 3
   ENDIF

   SELECT CASE( PULSETYPE )
   CASE( 0 )
       PRECIPFACT = 1.
   CASE( 1 )
       IF( ( HRDIFF / 24. ) < 2. ) THEN
           PRECIPFACT = 11.19 * EXP(-0.805*(HRDIFF+24)/24.)
       ELSE
           PULSETYPE = 0
           PRECIPFACT = 1.
       ENDIF
   CASE( 2 )
       IF( ( HRDIFF / 24. ) < 6. ) THEN
           PRECIPFACT = 14.68 * EXP(-0.384*(HRDIFF+24)/24.)
       ELSE
           PULSETYPE = 0
           PRECIPFACT = 1.
       ENDIF
   CASE DEFAULT
       IF( ( HRDIFF / 24. ) < 13. ) THEN
           PRECIPFACT = 18.46 * EXP(-0.208*(HRDIFF+24)/24.)
       ELSE
           PULSETYPE = 0
           PRECIPFACT = 1.
       ENDIF
   END SELECT

   RETURN

 END FUNCTION PRECIPFACT

end subroutine megan_nox


end module nox_mod
