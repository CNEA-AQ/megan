MODULE MEGAN_V32

   implicit none

   !Version dependent parameters:
   integer, parameter :: NrTyp = 6       !Number of "Canopy types": trees (needle,broad,tropical),shrub,grass,crop.
   integer, parameter :: NCLASS = 19     !MEGAN Internal Emission Categories
   !---
   integer, save :: nmgnspc           !number of megan     species
   integer, save :: n_scon_spc        !number of mechanism species
   character( 16 ), allocatable :: megan_names(:)     ! megan species names
   integer, allocatable ::  spmh_map(:),mech_map(:)   ! speciated species name

   real, allocatable :: conv_fac(:)
   real,allocatable :: mech_mwt(:)
   character( 16 ), allocatable :: mech_spc(:)

   INCLUDE 'tables/SPC_NOCONVER.EXT'
   INCLUDE 'tables/SPC_CB05.EXT'
   INCLUDE 'tables/SPC_CB6.EXT'
   INCLUDE 'tables/SPC_CB6_AE7.EXT'
   INCLUDE 'tables/SPC_RACM2.EXT'        ! new in MEGAN3
   INCLUDE 'tables/SPC_CRACMM.EXT'       ! new in CMAQ 5.4
   INCLUDE 'tables/MAP_CV2CB05.EXT'
   INCLUDE 'tables/SPC_SAPRC07.EXT'      ! new in MEGAN3
   INCLUDE 'tables/SPC_SAPRC07T.EXT'     ! new in MEGAN3
   INCLUDE 'tables/MAP_CV2CB6.EXT'
   INCLUDE 'tables/MAP_CV2CB6_AE7.EXT'
   INCLUDE 'tables/MAP_CV2RACM2.EXT'
   INCLUDE 'tables/MAP_CV2CRACMM.EXT'
   INCLUDE 'tables/MAP_CV2SAPRC07.EXT'
   INCLUDE 'tables/MAP_CV2SAPRC07T.EXT'
   
   INCLUDE 'tables/MEGAN.EXT'

contains

subroutine megan (yyyy,ddd,hh,                             & !year,jday,hour
             ncols,nrows,layers,                           & !dimensions
             lat,long,                                     & !lat,lon,canopy fractions
             temp,rad,wind,pres,qv,                        & !meteo instant variables
             laip, laic,                                   &
             ctf, efmaps, ldf_in,                          & !lai,emis factors, light emis factors
             lsm,soil_type,soil_moisture,                  & !land surface model, soil type, soil_moisture
             tmp_max, tmp_min, wind_max, tmp_avg, rad_avg, & !meteo daily
             emis                                          ) !out: Emision values

    implicit none
    ! input variables
    integer, intent(in)                           :: yyyy, ddd, hh           ! year, jday, hour
    integer, intent(in)                           :: ncols, nrows, layers    !dims x,y,levels
    real,intent(in), dimension(ncols,nrows)       :: lat, long, temp, rad, wind, pres, qv, laip,laic
    real,intent(in), dimension(ncols,nrows)       :: tmp_avg,rad_avg,tmp_min,tmp_max,wind_max

    real, intent(in) :: ctf(ncols,nrows,nrtyp) !canopy type factor array
    real, intent(in) :: efmaps(ncols,nrows,19) !only 19
    real, intent(in) :: ldf_in(ncols,nrows,4 ) !only 4 use maps

    character(len=4),intent(in)   :: LSM          !land surface model 
    integer,intent(in)     ::  soil_type(ncols,nrows)
    real,   intent(in)     ::  soil_moisture(ncols,nrows)

    ! output variables 
    real   ,intent(inout) :: emis(ncols,nrows,n_spca_spc)
    
    ! local variables
    integer :: s, t, i, j, k ! loop indices
    integer :: mm, dd

    ! megcan local variables 
     real    :: TotalCT
     real    :: month,day,hour
     real    :: SinZenith, Zenith!Sinbeta, Beta
     REAL    :: Solar, Maxsolar,Eccentricity,    &
          Difffrac, PPFDfrac, QbAbsn,            &
          Trate, Qbeamv,Qdiffv, Qbeamn, Qdiffn,  &
          QbAbsV,Ea1tCanopy, Ea1pCanopy,         &
          TairK0, HumidairPa0, Ws0, SH
     REAL,DIMENSION(LAYERS) :: VPgausWt, VPgausDis2,VPgausDis, VPslwWT, &
          QdAbsV, QsAbsV, QdAbsn,QsAbsn,                              &
          SunQv, ShadeQv, SunQn, ShadeQn,                             &
          TairK, HumidairPa, Ws, SunleafSH, sun_ppfd,shade_ppfd,      &
          SunleafLH,SunleafIR, ShadeleafSH, sun_tk,shade_tk,sun_frac, &
          ShadeleafLH,ShadeleafIR, sun_ppfd_total, shade_ppfd_total,  &
          sun_tk_total, shade_tk_total, sun_frac_total
    real, dimension(layers) :: sunt,shat,sunf,sunp,shap

    !megsea local variables:
    real, allocatable :: wwlt(:)

    !megvea local variables
    real                  :: ER !(ncols,nrows)                    !emission rate
    real                  :: non_dimgarma (ncols,nrows,nclass)  !

    logical, parameter    :: gambd_yn  = .false. !.true.!
    logical, parameter    :: gamaq_yn  = .false. !.true.!
    logical, parameter    :: gamht_yn  = .false. !.true.!
    logical, parameter    :: gamlt_yn  = .false. !.true.!
    logical, parameter    :: gamhw_yn  = .false. !.true.!
    logical, parameter    :: gamco2_yn = .false. !.true.!
    logical, parameter    :: gamsm_yn  = .false. !.true.! ! for the cmaq implementation of megan  we refer to soil moisture at layer 2, which is 1 meter for px and 0.5 m for noah.
                                                          ! Keep this in mind when enabling the GAMSM stress.

    real  :: cdea(layers)  ! Emission response to canopy depth
    real  :: gamla      ! EA leaf age response
    real  :: gamaq      ! EA response to air pollution
    real  :: gambd      ! EA bidirectional exchange LAI response
    real  :: gamht      ! EA response to high temperature
    real  :: gamlt      ! EA response to low temperature
    real  :: gamhw      ! EA response to high wind speed
    real  :: gamsm      ! EA response to soil moisture
    real  :: gamco2     ! EA response to CO2
    real  :: gamtp      ! combines GAMLD, GAMLI, GAMP to get canopy average
    real  :: ldfmap     ! light depenedent fraction map

    REAL :: VPGWT(LAYERS)
    REAL :: SUM1,SUM2,Ea1L,Ea2L
    
    !mgn2mech variables:
    integer :: nmpmg,nmpsp,nmpmc
    REAL    :: tmper(ncols, nrows, n_spca_spc)       ! Temp emission buffer
 
    ! EA response to canopy temperature/light
    IF ( Layers .EQ. 5 ) THEN
        VPGWT(1) = 0.1184635
        VPGWT(2) = 0.2393144
        VPGWT(3) = 0.284444444
        VPGWT(4) = 0.2393144
        VPGWT(5) = 0.1184635
    ELSE
        do k = 1,layers
            VPGWT(K) = 1.0 / FLOAT( Layers )
        end do
    ENDIF

    select case (LSM)
           case ('NOAH' )
              allocate(wwlt(size(wwlt_noah))); wwlt=wwlt_noah;
           case ('JN90' )
              allocate(wwlt(size(wwlt_jn90))); wwlt=wwlt_jn90;
           case DEFAULT
              allocate(wwlt(size(wwlt_jn90))); wwlt=wwlt_jn90;
    end select

   do j = 1, NROWS
      do i = 1, NCOLS! preserve stride 1 for output arrays

        !from megcan -----------
        !print*,"MEGCAN.."
        sunt(:) = temp(i,j) !default values for output
        shat(:) = temp(i,j)
        sunp(:) = rad(i,j)
        shap(:) = rad(i,j)
        sunf(:) = 1.0
        TotalCT=sum(ctf(i,j,:)) !*0.01
        if (totalCT .gt. 0.0 .AND. LAIc(i,j) .gt. 0.0 ) then

           ! Convert to "solar hour": 
           Hour  = real(HH) + long(i,j) / 15.0
           if ( hour  .lt. 0.0 ) then
             hour  = hour + 24.0; day  = real(ddd)  - 1
           elseif ( hour  .gt. 24.0 ) then
             hour  = hour - 24.0; day  = real(ddd)  + 1
           endif

            TairK0   = temp(i,j)      !temp (from meteo)
            Ws0      = wind(i,j)      !wind (from meteo)
            Solar    = rad(i,j)/2.25  !solar rad. (from meteo) [W m-2] -> [umol photons m-2 s-1]

           !(1) calc solar angle
            zenith      = CalcZenith(day,lat(i,j),hour)
            SinZenith   = sin(zenith / 57.29578) !57.29578=rad2deg
            Eccentricity= CalcEccentricity(Day)
            Maxsolar = SinZenith * SolarConstant * Eccentricity

           !(2) gaussian dist.
           call GaussianDist(VPgausDis, layers)

           !(3) determine fraction of diffuse PPFD, direct PPFD, diffuse near IR, direct near IR
           call SolarFractions(Solar, Maxsolar, Qdiffv, Qbeamv, Qdiffn, Qbeamn)

            sun_ppfd_total  = 0.0
            shade_ppfd_total= 0.0
            sun_tk_total    = 0.0
            shade_tk_total  = 0.0
            sun_frac_total  = 0.0

            do k = 1,NRTYP   !canopy types
              if ( ctf(i,j,k) .ne. 0.0 ) then
                 sun_ppfd   = 0.0
                 shade_ppfd = 0.0
                 sun_tk     = 0.0
                 shade_tk   = 0.0
                 sun_frac   = 0.0

                 !(4) canopy radiation dist (ppdf)
                 !
                 call CanopyRad(VPgausDis, layers, LAIc(i,j), SinZenith,       & !in
                       Qbeamv, Qdiffv, Qbeamn, Qdiffn, k, Canopychar, sun_frac,& !in
                       QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv,  & !in
                       ShadeQv, SunQn, ShadeQn, sun_ppfd, shade_ppfd,          & !out
                       NrCha,NrTyp)

                 HumidairPa0  =  WaterVapPres(qv(i,j), pres(i,j), waterairratio)
                 Trate        =  Stability(Canopychar, k, Solar , NrCha, NrTyp)
                 !(5) canopy energy balance (temp)
                 !
                 call CanopyEB(Trate, Layers, VPgausDis, Canopychar, k,    &
                       TairK, HumidairPa, Ws, sun_ppfd,                    &
                       shade_ppfd, SunQv, ShadeQv, SunQn, ShadeQn,         &
                       sun_tk, SunleafSH, SunleafLH, SunleafIR,            &
                       shade_tk,ShadeleafSH,ShadeleafLH,ShadeleafIR,       &
                       NrCha, NrTyp, Ws0, TairK0, HumidairPa0)
                 !(6) compute output variables
                 !
                 sun_ppfd_total(:)   = sun_ppfd_total(:)   + sun_ppfd(:)  * ctf(i,j,k)!*0.01
                 shade_ppfd_total(:) = shade_ppfd_total(:) + shade_ppfd(:)* ctf(i,j,k)!*0.01
                 sun_tk_total(:)     = sun_tk_total(:)     + sun_tk(:)    * ctf(i,j,k)!*0.01
                 shade_tk_total(:)   = shade_tk_total(:)   + shade_tk(:)  * ctf(i,j,k)!*0.01
                 sun_frac_total(:)   = sun_frac_total(:)   + sun_frac(:)  * ctf(i,j,k)!*0.01
              endif
            enddo
                sunt(:) = sun_tk_total(:)/TotalCT
                shat(:) = shade_tk_total(:)/TotalCT
                sunp(:) = sun_ppfd_total(:)/TotalCT
                shap(:) = shade_ppfd_total(:)/TotalCT
                sunf(:) = sun_frac_total(:)/TotalCT

        else if (totalCT .lt. 0) then
               print*,"Send ERROR message!"             !Send ERROR message!
        else   !totalCT == 0
               !print*,"Default values!"
        endif

        !-----------------------
        !print*,"MEGSEA.."
        !from megsea -----------
        !EA response to Soil Moisture
        gamsm=gamma_sm(soil_type(i,j),soil_moisture(i,j),wwlt(soil_type(i,j)) )  

        !from megvea -----------
        !print*,"MEGVEA.."
        ! Emission response to canopy depth
        cdea(:)=gamma_cd(layers,laic(i,j))  
        ! EA bidirectional exchange LAI response
        if ( gambd_yn )  then; gambd=gamma_laibidir(laic(i,j)); else;  gambd = 1.0; endif
        ! EA response to co2
        if ( gamco2_yn ) then; gamco2=gamma_co2(co2)          ; else; gamco2 = 1.0; endif

        do s=1,NCLASS ! Loop over all the emission classes !Now process all factors dependent on S:

            IF ( S .EQ. 3 .OR. S .EQ. 4 .OR. S .EQ. 5 .OR. S .EQ. 6 ) THEN
                LDFMAP = LDF_IN(i,j,S-2) ! only LDF 3, 4, 5, and 6 in file
            ELSE
                LDFMAP = LDF(S) !For these species,  Read LDF from previous MEGVEA.EXT 
            ENDIF
            ! EA response to leaf age 
            gamla = gamma_age(s, laip(i,j), laic(i,j), tmp_avg(i,j))
            ! EA response to air quality
            IF ( GAMAQ_YN ) THEN; GAMAQ=GAMMA_AQ(S, AQI_default)   ; ELSE; GAMAQ = 1.0; ENDIF
            ! EA response to high temperature
            IF ( GAMHT_YN ) THEN; GAMHT=GAMMA_HT(S, tmp_max(i,j))  ; ELSE; GAMHT = 1.0; ENDIF
            ! EA response to low temperature
            IF ( GAMLT_YN ) THEN; GAMLT=GAMMA_LT(S, tmp_min(i,j))  ; ELSE; GAMLT = 1.0; ENDIF
            ! EA response to high wind speed
            IF ( GAMHW_YN ) THEN; GAMHW=GAMMA_HW(S, wind_max(i,j)) ; ELSE; GAMHW = 1.0; ENDIF

            SUM1 = 0.0
            SUM2 = 0.0
            do k = 1, layers
              
              Ea1L = CDEA(K) *                                                                       &
                    GAMTLD(SunT(k),tmp_avg(i,j),S) *  GAMP(SunP(k),rad_avg(i,j)) *        SunF(k)  + &! *2.025 is the conversion to PPFD. 
                    GAMTLD(ShaT(k),tmp_avg(i,j),S) *  GAMP(ShaP(k),rad_avg(i,j)) * (1.0 - SunF(k) )   ! *2.025 is the conversion to PPFD. 
              SUM1 = SUM1 + Ea1L*VPGWT(K)

              Ea2L = GAMTLI(SunT(k),S) * SunF(k)   +    GAMTLI(ShaT(k),S) * (1.0-SunF(k))
              SUM2 = SUM2 + Ea2L*VPGWT(K)
            end do   ! end do canopy layers

            GAMTP = SUM1*LDFMAP + SUM2*( 1.0-LDFMAP )

            ! ... Calculate emission activity factors
            !ER(I,J) = LAIc(I,J) * GAMTP * GAMLA * GAMHW * GAMAQ* GAMHT * GAMLT *  GAMSM
            ER = LAIc(I,J) * GAMTP * GAMLA * GAMHW * GAMAQ* GAMHT * GAMLT *  GAMSM

            IF ( S .EQ. 1 ) THEN
                ER =ER * GAMCO2  ! GAMCO2 only applied to isoprene
            ELSE IF ( S .EQ. 13 ) THEN   
                ER = ER * GAMBD  ! GAMBD only applied to ethanol and acetaldehyde
            END IF

            !IF ( ER(I,J) .GT. 0.0 ) THEN
            IF ( ER .GT. 0.0 ) THEN
                NON_DIMGARMA (i,j,s) = ER
            ELSE                   
                NON_DIMGARMA (i,j,s) = 0.0
            END IF
        end do  ! End loop of species (S)

     end do   ! NCOLS
  end do ! NROWS


  !from mgn2mech ---------
  print*,"MGN2MECH.."
  tmper = 0.
  emis = 0.

  do s = 1, n_smap_spc
    nmpmg = mg20_map(s) !megan category
    nmpsp = spca_map(s) !megan specie
  
    IF ( nmpmg .NE. i_NO ) then !...  Not NO
       tmper(:,:,nmpsp) = non_dimgarma(:,:,nmpmg) * efmaps(:,:,nmpmg)  * effs_all(s)
    ELSEIF ( nmpmg .EQ. i_NO ) then
       tmper(:,:,nmpsp) = 0.0 ! not NO produced  by plants
      !@!!-----------------NO Stuff-----------------------
      !@IF ( .NOT. BDSNP_MEGAN ) THEN
      !@!     GAMNO is emission activity factor
      !@   tmper(:,:,nmpsp) = GAMNO(:,:) * efmaps(:,:,i_NO)        * effs_all(s)
      !@ELSE
      !@! directly use BDSNP soil NO
      !@  tmper(nmpsp,:,:) = BDSNP_NO(:,:)
      !@ENDIF
      !@!-----------------end of NO----------------------
    ENDIF     !IF ( nmpmg .NE. i_NO ) then
  enddo ! end species loop
  !.....3) Conversion from speciated species to MECHANISM species
   tmper = tmper * nmol2mol

   ! lumping to MECHANISM species
   do s = 1, n_scon_spc
     nmpsp = spmh_map(s)         ! Mapping value for SPCA
     nmpmc = mech_map(s)         ! Mapping value for MECHANISM
     if ( nmpmc .ne. 999 ) then
        emis(:,:,nmpmc) = emis(:,:,nmpmc) +  (tmper(:,:,nmpsp) * conv_fac(s))
     endif
   ENDDO ! End species loop
  !-----------------------------------------------------------------------

  RETURN
    
contains

    function gamma_cd(Layers,LAI)  result(cdea)
      implicit none
      integer, intent(in)     :: layers
      real                    :: lai 
      real, dimension(layers) :: cdea
      integer :: k
      do k=1,layers
             cdea(k)=ccd1*min(lai*(k-0.5)/float(layers),3.0)+ccd2
      enddo
      return
    end function
    !----------------------------------------------------------------
    function gamma_laibidir(lai)  result(gambd)
      real, intent(in) :: lai
      real :: gambd

      IF(LAI < 2) THEN
          GAMBD =  0.5 * LAI
      ELSEIF (LAI .LE. 6 ) THEN
          GAMBD = 1 - 0.0625 * (LAI - 2)
      ELSE
          GAMBD = 0.75
      ENDIF
    end function
    !----------------------------------------------------------------
    function gamma_co2(co2)        result(gamco2)
        implicit none
        REAL,INTENT(IN)     :: CO2
        REAL                :: GAMCO2
        ! local
        REAL    :: Ci, cxxx, cyyy

        Ci      = 0.7 * CO2
        IF ( CO2 .EQ. 400.0 ) THEN
            GAMCO2 = 1.0
        ELSE
            cxxx =  Ci**CO2h
            cyyy =  Cstar**CO2h
            GAMCO2 = ISmax- ((ISmax * cxxx) / (cyyy + cxxx))
        END IF

        RETURN
    end function gamma_co2

    !----------------------------------------------------------------
    ! EA Temperature response (light dependent emission)
    !----------------------------------------------------------------
    FUNCTION GAMTLD(T1,T24,S)
        IMPLICIT NONE
        REAL,PARAMETER :: Ct2 = 230
        INTEGER        :: S
        REAL           :: T1,T24,T240,Topt, X, Eopt, GAMTLD

        T240 = T24

        IF (T1 < 260.0) THEN
            GAMTLD = 0.0
        ELSE
            ! Temperature at which maximum emission occurs
            Topt = 312.5 + 0.6 * (T240 - 297.0)
            X    = ((1.0 / Topt) - (1.0 / T1)) / 0.00831
            ! Maximum emission (relative to emission at 30 C)
            Eopt = Cleo(S) * EXP(0.05 * (T24 - 297.0)) *          &
                  Exp(0.05*(T240-297.0))

            GAMTLD= Eopt * Ct2 * Exp(Ct1(S) * X) /                &
                  (Ct2 - Ct1(S) * (1.0 - EXP(Ct2 * X)))
        ENDIF
    END FUNCTION GAMTLD
    !----------------------------------------------------------------
    ! EA Temperature response (light independent emission)
    !----------------------------------------------------------------
    function gamtli(temp,s)
        IMPLICIT NONE

        REAL           :: temp, GAMTLI
        REAL,PARAMETER :: Ts = 303.15
        INTEGER        :: S

        GAMTLI = exp( beta(S)*(temp-Ts) )

    end function gamtli
    !----------------------------------------------------------------
    ! EA Light response
    !----------------------------------------------------------------
    function gamp(ppfd1,ppfd24)
        implicit none
        real            :: ppfd1, ppfd24, alpha, c1, gamp
        IF (PPFD24 < 0.01) THEN
            GAMP= 0.0
        ELSE
            Alpha  = 0.004
            !        C1     = 0.0468 * EXP(0.0005 * (PPFD24 - PSTD))
            !     &          * (PPFD24 ** 0.6)
            C1 = 1.03
            !        GAMP= (Alpha * C1 * PPFD1) / ((1 + Alpha**2. * PPFD1**2.)**0.5)
        !   use SQRT her for clarity and efficiency
            GAMP= (Alpha * C1 * PPFD1) / SQRT(1.0 + Alpha**2 * PPFD1**2)
        ENDIF
    end function gamp
    !----------------------------------------------------------------
    ! EA response to high temperature
    !----------------------------------------------------------------
    function GAMMA_HT(S, tmp_max)  result(gamht)
        implicit none
        ! input
        integer,intent(in)     :: s
        real,intent(in)        :: tmp_max
        ! output
        real                   :: gamht
        ! local
        REAL        :: THTK, t1
            THTK = 273.15 + THT(S)
            t1 = THTK + DTHT(S)
            IF (tmp_max <= THTK) THEN
                GAMHT = 1.0
            ELSE IF ( tmp_max < t1) THEN
                GAMHT = 1.0 + (CHT(S) - 1.0)* (tmp_max - THTK)/DTHT(S)
            ELSE
                GAMHT = CHT(S)
            ENDIF
        RETURN
    end function GAMMA_HT
    !----------------------------------------------------------------
    !  EA response to low temperature
    !----------------------------------------------------------------
    function GAMMA_LT(S, tmp_min) result(gamlt)
        implicit none
        ! input
        integer,intent(in):: s
        real,intent(in)   :: tmp_min
        ! output
        real              :: gamlt
        ! local
        real         :: tltk, t1

        TLTK = 273.15 + TLT(S)
        t1 = TLTK - DTLT(S)
        IF (tmp_min >= TLTK) THEN
            GAMLT = 1.0
        ELSE IF ( tmp_min > t1) THEN
            GAMLT = 1.0 + (CLT(S) - 1.0)* (TLTK - tmp_min)/DTLT(S)
        ELSE
            GAMLT = CLT(S)
        ENDIF

        RETURN
    end function GAMMA_LT
    !----------------------------------------------------------------
    ! EA response to high wind speed
    !----------------------------------------------------------------
    function GAMMA_HW(S, wind_max) result(gamhw)
        implicit none
        ! input
        integer,intent(in) :: S
        real,intent(in)    :: wind_max
        ! output
        real               :: GAMHW
        ! local
        real        :: t1

            t1 = THW(S) + DTHW(S)
            IF (wind_max <= THW(S)) THEN
                GAMHW = 1.0
            ELSE IF ( wind_max < t1) THEN
                GAMHW = 1.0 + (CHW(S) - 1.0)* (wind_max - THW(S))/ DTHW(S)
            ELSE
                GAMHW = CHW(S)
            ENDIF
        RETURN
    end function gamma_hw
    !----------------------------------------------------------------
    ! EA response to air quality
    !----------------------------------------------------------------
    function gamma_aq(s, aqi) result(gamaq)
        implicit none
        ! input
        integer, intent(in)   :: s
        real,    intent(in)   :: aqi
        ! output
        real                :: gamaq
        ! local
        REAL       :: t1
        t1 = TAQ(S) + DTAQ(S)
        IF (AQI <= TAQ(S)) THEN
            GAMAQ = 1.0
        ELSE IF ( AQI < t1) THEN
            GAMAQ = 1.0 + (CAQ(S) - 1.0)* (AQI - TAQ(S))/DTAQ(S)
        ELSE
            GAMAQ = CAQ(S)
        ENDIF

        RETURN
    end function gamma_aq
    !----------------------------------------------------------------
    ! EA response to leaf age
    !----------------------------------------------------------------
    function gamma_age(s, laip, laic, tt) result(gamla)

        IMPLICIT NONE
        ! input
        INTEGER,INTENT(IN):: S
        REAL,INTENT(IN)   :: Tt, LAIp, LAIc
        ! output
        REAL              :: GAMLA

        REAL :: Fnew, Fgro, Fmat, Fold
        REAL :: ti,tm
        REAL       :: TSTLEN  
        !Time step of LAI data
        !if (USE_MEGAN_LAI) THEN
        !  TSTLEN = 8.0 ! 8 daily from MEGAN file
        !else
          TSTLEN = 1.0 ! 1 Daily from soilout/metcro
        !end if

        !---------------------------------------------------
        ! local parameter arrays
        !... Calculate foliage fraction
        IF (LAIp .LT. LAIc) THEN

            !        Calculate ti and tm
            IF (Tt .LE. 303.0) THEN
                ti = 5.0 + 0.7*(300-Tt)
            ELSE
                ti = 2.9
            END IF
            tm = 2.3*ti

            !       Calculate Fnew and Fmat, then Fgro and Fold
            !       Fnew
            IF (ti .GE. TSTLEN) THEN
                Fnew = 1.0 - (LAIp/LAIc)
            ELSE
                Fnew = (ti/TSTLEN) * ( 1-(LAIp/LAIc) )
            END IF

            !       Fmat
            IF (tm .GE. TSTLEN) THEN
                Fmat = LAIp/LAIc
            ELSE
                Fmat = (LAIp/LAIc) + ( (TSTLEN-tm)/TSTLEN ) * ( 1-(LAIp/LAIc) )
            END IF
            Fgro = 1.0 - Fnew - Fmat
            Fold = 0.0

        ELSE IF (LAIp .EQ. LAIc) THEN
            Fnew = 0.0; Fgro = 0.1; Fold = 0.1            ; Fmat = 0.8
        ELSE IF (LAIp .GT. LAIc) THEN
            Fnew = 0.0; Fgro = 0.0; Fold =(LAIp-LAIc)/LAIp; Fmat = 1-Fold
        END IF
        !...  Calculate GAMLA
        GAMLA = Fnew*Anew(S) + Fgro*Agro(S) + Fmat*Amat(S) + Fold*Aold(S)

        RETURN
    end function gamma_age

    !from MEGSEA ==========================================================
    function gamma_sm(sltyp, soilm, wilt)
        implicit none
        real :: t1,soilm,wilt,gamma_sm
        integer :: sltyp  

         !wilt = wwlt(sltyp)
         t1 = wilt + d1
         if ( soilm < wilt ) then
             gamma_sm = 0
         else if ( soilm >= wilt .and. soilm < t1 ) then
             gamma_sm = (soilm - wilt)/d1
         else
             gamma_sm = 1
         end if
    end function gamma_sm
    !======================================================================






!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! MEGCAN FUNCTIONS:                                             o
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Calculates the solar zenith angle
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION CalcZenith(Day, Lat, Hour)
      IMPLICIT NONE
      REAL    :: Day
      REAL    :: Rpi,Hour,Lat,SinDelta,CosDelta,A,B,SinZenith,CalcZenith
      REAL,PARAMETER :: PI = 3.14159, Rpi180 = 57.29578

      SinDelta = -SIN(0.40907) * COS(6.28 * (Day + 10) / (365))
      CosDelta = (1 - SinDelta**2.)**0.5

      A = SIN(Lat / Rpi180) * SinDelta
      B = COS(Lat / Rpi180) * Cosdelta
      SinZenith = A + B * COS(2 * PI * (Hour - 12) / 24)
      CalcZenith= ASIN(SinZenith) * Rpi180 !57.29578
END FUNCTION CalcZenith
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   FUNCTION CalcEccentricity
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION CalcEccentricity(Day)
      IMPLICIT NONE
      REAL    :: Day
      !INTEGER :: Day
      REAL :: CalcEccentricity
      CalcEccentricity = 1 + 0.033 * COS(2*3.14*(Day-10)/365)
END FUNCTION CalcEccentricity
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   SUBROUTINE GaussianDist
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
subroutine GaussianDist(Distgauss, Layers)
      IMPLICIT NONE
      INTEGER,INTENT(IN) ::  Layers
      REAL,DIMENSION(Layers),INTENT(OUT) :: Distgauss
      ! local variables
      INTEGER ::  i
!----------------------------------------------------------------
      IF (Layers .EQ. 1) THEN
        Distgauss(1)   = 0.5
      ELSEIF (Layers .EQ. 3) THEN
        Distgauss(1)   = 0.112702
        Distgauss(2)   = 0.5
        Distgauss(3)   = 0.887298
      ELSEIF (Layers .EQ. 5) THEN
        Distgauss(1)   = 0.0469101
        Distgauss(2)   = 0.2307534
        Distgauss(3)   = 0.5
        Distgauss(4)   = 0.7692465
        Distgauss(5)   = 0.9530899
      ELSE
        DO i = 1, Layers
          Distgauss(i) = (i - 0.5) / Layers
        ENDDO
      ENDIF
      RETURN
end subroutine GaussianDist
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
subroutine SolarFractions(Solar, Maxsolar, Qdiffv,Qbeamv,Qdiffn,Qbeamn)
      IMPLICIT NONE
      REAL,INTENT(IN) :: Solar, Maxsolar
      REAL,INTENT(OUT) ::  Qdiffv,Qbeamv, Qdiffn, Qbeamn
      REAL :: FracDiff, PPFDfrac,PPFDdifFrac,Qv, Qn
      ! internal variables
      !INTEGER :: I,J
      REAL ::  Transmis

      IF (Maxsolar  <= 0) THEN
        Transmis  = 0.5
      ELSEIF (Maxsolar < Solar) THEN
        Transmis  = 1.0
      ELSE
        Transmis  = Solar  / Maxsolar
      ENDIF
      !FracDiff is based on Lizaso 2005
      FracDiff    = 0.156 + 0.86/(1 + EXP(11.1*(Transmis -0.53)))
      !PPFDfrac is based on Goudrian and Laar 1994
      PPFDfrac    = 0.55 -Transmis*0.12
      !PPFDdifFrac is based on data in Jacovides 2007
      PPFDdifFrac = FracDiff *(1.06 + Transmis*0.4)
      ! Calculate  Qdiffv,Qbeamv, Qdiffn, Qbeamn in the subroutine
      IF (PPFDdifFrac > 1.0) THEN
      PPFDdifFrac = 1.0
      ENDIF
      Qv     = PPFDfrac * Solar 
      Qdiffv = Qv * PPFDdifFrac   !diffuse PPFD
      Qbeamv = Qv - Qdiffv        !direct PPFD
      Qn     = Solar - Qv
      Qdiffn =  Qn * FracDiff     !diffuse near IR
      Qbeamn =  Qn - Qdiffn       !direct near IR
      RETURN
end subroutine SolarFractions
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!   Canopy light environment model
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!   based on Spitters et al. (1986), 
!   Goudrian and van Laar (1994), Leuning (1997)
!   Initial code 8-99, 
!   modified 7-2000, 12-2001, 1-2017
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
SUBROUTINE CanopyRad(Distgauss, Layers, LAI, SinZenith,           &
                Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype,      &
                Canopychar, Sunfrac, QbAbsV, QdAbsV, QsAbsV,  &
                QbAbsn, QdAbsn, QsAbsn, SunQv,                &
                ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD,  &
                NrCha, NrTyp)

      implicit none
      ! input
      INTEGER,INTENT(IN) :: Layers, NrCha, NrTyp, Cantype
      REAL,INTENT(IN) :: Qbeamv,Qdiffv,SinZenith,LAI,Qbeamn,Qdiffn
      REAL,DIMENSION(Layers),INTENT(IN) :: Distgauss
      ! output
      REAL,INTENT(OUT) :: QbAbsV, QbAbsn
      REAL,DIMENSION(Layers),INTENT(OUT) :: ShadePPFD, SunPPFD, &
                    QdAbsv, QsAbsv, QsAbsn, ShadeQv,  SunQn,  &
                    QdAbsn, SunQv, ShadeQn, Sunfrac 
      REAL,DIMENSION(NrCha,NrTyp),INTENT(OUT) :: Canopychar
      ! internal variables
      INTEGER :: i
      REAL :: ScatV, ScatN, RefldV, RefldN, ReflbV, ReflbN,     & 
             Kb, Kd, KbpV, KbpN, KdpV, KdpN, LAIdepth, Cluster, & 
             QdAbsVL, QsAbsVL, QdAbsNL, QsAbsNL, CANTRAN, LAIadj
      
      REAL,PARAMETER :: ConvertShadePPFD = 4.6
      REAL,PARAMETER :: ConvertSunPPFD = 4.0

      ! adjust LAI for canopy transparency
      CANTRAN = Canopychar(17,Cantype)
      LAIadj = LAI / ( 1 - CANTRAN )

     IF (((Qbeamv  + Qdiffv ) > 0.001) .AND.  (SinZenith  > 0.002) .AND. (LAIadj  > 0.001)) THEN       ! Daytime

        ! Scattering coefficients (scatV,scatN), diffuse and beam reflection 
        ! coefficients (ref..) for visible or near IR
        ScatV   = Canopychar(5,Cantype)
        ScatN   = Canopychar(6,Cantype)
        RefldV  = Canopychar(7,Cantype)
        RefldN  = Canopychar(8,Cantype)
        Cluster = Canopychar(9,Cantype)
        
        ! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
        Kb = Cluster * 0.5 / SinZenith
        ! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
        Kd = 0.8 * Cluster
        ! (0.8 assumes a spherical leaf angle distribution)

        Call CalcExtCoeff(Qbeamv,ScatV,Kb,Kd,ReflbV,KbpV,KdpV,QbAbsV)
        Call CalcExtCoeff(Qbeamn,ScatN,Kb,Kd,ReflbN,KbpN,KdpN,QbAbsn)

        DO i = 1,layers
          
          LAIdepth   = LAIadj  * Distgauss(i) ! LAI depth at this layer
          Sunfrac(i) = EXP(-Kb * LAIdepth)    !fraction of leaves that are sunlit

          Call CalcRadComponents(Qdiffv , Qbeamv , kdpV, kbpV, kb, scatV, refldV, reflbV, LAIdepth, QdAbsVL, QsAbsVL)
          Call CalcRadComponents(Qdiffn , Qbeamn , kdpN, kbpN, kb, scatN, refldN, reflbN, LAIdepth, QdAbsNL, QsAbsNL)

          ShadePPFD(i) = (QdAbsVL + QsAbsVL) * ConvertShadePPFD/(1 - scatV)
          SunPPFD(i) = ShadePPFD(i) + (QbAbsV* ConvertSunPPFD/(1 - scatV))
          QdAbsV(i) = QdAbsVL
          QsAbsV(i) = QsAbsVL
          QdAbsn(i) = QdAbsNL
          QsAbsn(i) = QsAbsNL
          ShadeQv(i) = QdAbsVL + QsAbsVL
          SunQv(i)   = ShadeQv(i) + QbAbsV
          ShadeQn(i) = QdAbsNL + QsAbsNL
          SunQn(i)   = ShadeQn(i) + QbAbsn
        ENDDO

     ELSE                           ! Night time
       QbAbsV    = 0;   QbAbsn       = 0
       Sunfrac(:)= 0.2; SunQn(:)    = 0; ShadeQn(:) = 0
       SunQv(:)  = 0  ; ShadeQv(:)  = 0; SunPPFD(:) = 0; ShadePPFD(:)= 0
       QdAbsV(:) = 0  ; QsAbsV(:)   = 0; QdAbsn(:)  = 0; QsAbsn(:)   = 0
     ENDIF
     RETURN
END SUBROUTINE CanopyRad
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Calculate light extinction coefficients
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
SUBROUTINE CalcExtCoeff(Qbeam,scat,kb,kd,reflb,kbp,kdp,QbeamAbsorb)
      IMPLICIT NONE
      REAL,INTENT(IN) :: Qbeam, scat, Kb, Kd
      REAL,INTENT(OUT) :: Reflb, Kbp, Kdp, QbeamAbsorb
      ! local variables
      REAL :: P

      P     = (1 - scat)**0.5
      Reflb = 1 - Exp((-2 * ((1 - P) / (1 + P)) * kb) / (1 + kb))
      ! Extinction coefficients
      Kbp   = Kb * P
      Kdp   = Kd * P
      ! Absorbed beam radiation
      QbeamAbsorb = kb * Qbeam * (1 - scat)
      RETURN
END SUBROUTINE CalcExtCoeff

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
SUBROUTINE CalcRadComponents(Qdiff, Qbeam, kdp, kbp, kb, scat, refld, reflb, LAIdepth, QdAbs, QsAbs)
      IMPLICIT NONE
      REAL,INTENT(IN)    :: Qdiff,Qbeam,kdp,kbp,kb,scat,refld,reflb,LAIdepth
      REAL,INTENT(OUT)   :: QdAbs, QsAbs
!-------------------------------------------------------------------
      QdAbs = Qdiff *   Kdp * (1 - Refld) * Exp(-Kdp * LAIdepth)
      QsAbs = Qbeam * ((Kbp * (1 - Reflb) * Exp(-Kbp * LAIdepth)) - (Kb * (1 - Scat) * Exp(-Kb * LAIdepth)))
      RETURN
END SUBROUTINE CalcRadComponents

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Canopy energy balance model for estimating leaf temperature
!   Coded into FORTRAN by Xuemei Wang
!   Code developed by Alex Guenther in 1990s
!   based on Goudrian and Laar (1994) and Leuning (1997)
!   Initial code 8-99, modified 7-2000 and 12-2001
!   Modified in 1-2017 by Alex Guenther and Ling Huang
!   to correct IR balance and atmos. emissivity
!   Note: i denotes an array containing a vertical profile 
!         through the canopy with 0 (above canopy conditions) 
!         plus 1 to number of canopy layers
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

SUBROUTINE CanopyEB(Trate, Layers, Distgauss, Canopychar,            &
             Cantype, TairK, HumidairPa, Ws,                         &
             SunPPFD, ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn,     &
             Sunleaftk, SunleafSH, SunleafLH,                        &
             SunleafIR, Shadeleaftk, ShadeleafSH,                    &
             ShadeleafLH, ShadeleafIR, NrCha, NrTyp, Ws0,            &
             TairK0, HumidairPa0)
     IMPLICIT NONE

! inputs
     INTEGER,INTENT(IN) :: NrCha, NrTyp, Layers, Cantype
     REAL,INTENT(IN) :: Trate, TairK0, HumidairPa0, Ws0
     REAL,DIMENSION(Layers),INTENT(IN) ::  Distgauss, SunQv,ShadeQv,SunQn, ShadeQn, SunPPFD, ShadePPFD
     REAL,DIMENSION(NrCha, NrTyp),INTENT(IN)  :: Canopychar

! outputs
      REAL,DIMENSION(Layers),INTENT(OUT) :: HumidairPa, Ws, Sunleaftk, SunleafSH, SunleafLH, SunleafIR, TairK, Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR

! local variables
      INTEGER :: i
      REAL :: Cdepth, Lwidth, Llength, Cheight, Eps, TranspireType, Deltah, EmissAtm, IRin,IRout !LeafIR, 
!     &         Deltah, UnexposedLeafIRin, ExposedLeafIRin, IRin,IRout
      REAL,DIMENSION(Layers) :: Ldepth, Wsh

      Cdepth        = Canopychar(1, Cantype)
      Lwidth        = Canopychar(2, Cantype)
      Llength       = Canopychar(3, Cantype)
      Cheight       = Canopychar(4, Cantype)
      Eps           = Canopychar(10,Cantype)
      TranspireType = Canopychar(11,Cantype)

      IF (TairK0  > 288) THEN
! Pa m-1  (humidity profile for T < 288)
        Deltah =  Canopychar(14, Cantype) / Cheight
      ELSEIF (TairK0  > 278) THEN
        Deltah =(Canopychar(14,Cantype)-((288-TairK0)/10) * (Canopychar(14,Cantype)-Canopychar(15,Cantype)))/Cheight
      ELSE
! Pa m-1  (humidity profile for T <278)
        Deltah = Canopychar(15, Cantype) / Cheight
      ENDIF

      Ldepth(:)     = Cdepth * Distgauss(:)
      TairK(:)      = TairK0  + (Trate  * Ldepth(:))      ! check this
      HumidairPa(:) = HumidairPa0  + (Deltah * Ldepth(:))

      Wsh(:) = (Cheight-Ldepth(:)) - (Canopychar(16,Cantype) * Cheight)
      Ws(:)  = (Ws0*LOG(Wsh(:))/LOG(Cheight-Canopychar(16,Cantype) * Cheight))
      WHERE (Wsh(:) < 0.001) Ws(:) = 0.05

      DO i=1,Layers

         ! REVISE - Replace UnexposedLeafIR with LeafIR
         
         !        IRin     = UnexposedLeafIRin(TairK(i), Eps)
         !        ShadeleafIR(i) = 2 * IRin
         !        SunleafIR(i) = 0.5*ExposedLeafIRin(HumidairPa0,TairK0)+1.5*IRin
         
         ! Apparent atmospheric emissivity for clear skies: 
         ! function of water vapor pressure (Pa) 
         ! and ambient Temperature (K) based on Brutsaert(1975) 
         ! referenced in Leuning (1997)
         EmissAtm        = 0.642 * (HumidairPa(i) / TairK(i))**(1./7.)   
         IRin            = LeafIR (TairK(i), EmissAtm)
         ShadeleafIR(i)  = IRin
         SunleafIR(i)    = IRin

      ! Sun
        CALL LeafEB(SunPPFD(i), SunQv(i) + SunQn(i),                    &
                   SunleafIR(i), Eps, TranspireType, Lwidth, Llength,   &
                   TairK(i), HumidairPa(i), Ws(i),                      &
                   Sunleaftk(i), SunleafSH(i),SunleafLH(i),             &
                   IRout )

         SunleafIR(i) = SunleafIR(i) - IRout

      ! Shade
        CALL LeafEB(ShadePPFD(i), ShadeQv(i)+ShadeQn(i),                &
                     ShadeleafIR(i),Eps,TranspireType, Lwidth,Llength,  &
                     TairK(i), HumidairPa(i), Ws(i),                    &
                     Shadeleaftk(i), ShadeleafSH(i),ShadeleafLH(i),     &
                     IRout)

         ShadeleafIR(i) = ShadeleafIR(i) - IRout
      ENDDO

      RETURN
END SUBROUTINE CanopyEB
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Leaf energy balance
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
SUBROUTINE LeafEB(PPFD, Q, IRin, Eps, TranspireType,         &
         Lwidth, Llength, TairK, HumidairPa, Ws, Tleaf,      &
         SH, LH, IRout)

      IMPLICIT NONE

      REAL,INTENT(IN) :: Eps, TranspireType, Lwidth, Llength,PPFD, Q, IRin, TairK, HumidairPa, Ws
      REAL,INTENT(OUT) :: IRout, Tleaf, SH, LH

! local variables

      INTEGER :: i
      REAL :: HumidAirKgm3,GHforced,StomRes,IRoutairT,LatHv,Ws1,      &
        LHairT,Tdelt,Balance,& !LeafBLC,& !LeafH,LeafLE,& !LHV,LeafIR,                 &
        GH1,SH1,LH1,E1,IRout1,GH !,ConvertHumidityPa2kgm3 !ResSC,
!     &        LHairT,Tdelt,Balance,LeafBLC,LeafH,LeafLE,LeafIRout,   

      IF (Ws <= 0) THEN
        Ws1 = 0.001
      ELSE
        Ws1 = Ws
      ENDIF

      ! Air vapor density kg m-3
      HumidAirKgm3 = ConvertHumidityPa2kgm3(HumidairPa, TairK)

      ! Heat convection coefficient (W m-2 K-1) for forced convection. 
      ! Nobel page 366
      GHforced = 0.0259 / (0.004 * ((Llength / Ws)**0.5))

      ! Stomatal resistence s m-1
      StomRes  = ResSC(PPFD)

      ! REVISE - Replace LeafIRout with LeafIR
      !      IRoutairT = LeafIROut(tairK, eps)
      !XJ      IRoutairT  = LeafIR(TairK + Tdelt, Eps)
       IRoutairT = LeafIR(TairK, Eps)

      ! Latent heat of vaporization (J Kg-1)
      LatHv = LHV(TairK)

      ! Latent heat flux
      LHairT = LeafLE(TairK,HumidAirKgm3,LatHv,GHforced,StomRes, TranspireType)

      E1 = (Q + IRin - IRoutairT - LHairT)
      IF (E1 .EQ. 0.) THEN
        E1 = -1.
      ENDIF

      Tdelt = 1.
      Balance = 10.
      DO i = 1, 10
        IF (ABS(Balance) > 2) THEN
          ! Boundary layer conductance
          GH1 = LeafBLC(GHforced, Tdelt, Llength)
          ! Convective heat flux
          SH1 = LeafH(Tdelt, GH1)
          ! Latent heat of vaporization (J Kg-1)
          LatHv = LHV(TairK + Tdelt)
          LH = LeafLE(TairK + Tdelt, HumidAirKgm3, LatHv, GH1, StomRes, TranspireType)
          LH1 = LH - LHairT
          ! REVISE - Replace LeafIROut with LeafIR
          !          IRout  = LeafIROut(TairK + Tdelt, Eps)
          IRout  = LeafIR(TairK + Tdelt, Eps)
          IRout1 = IRout - IRoutairT
          Tdelt  = E1 / ((SH1 + LH1 + IRout1) / Tdelt)
          Balance = Q + IRin - IRout - SH1 - LH
        ENDIF
      ENDDO

      If (Tdelt > 10)  Tdelt = 10
      If (Tdelt < -10) Tdelt = -10

      Tleaf = TairK + Tdelt
      GH    = LeafBLC(GHforced, Tleaf - TairK, Llength)
      SH    = LeafH(Tleaf - TairK, GH)
      LH    = LeafLE(Tleaf, HumidAirKgm3, LatHv, GH, StomRes, TranspireType)

      ! REVISE - Replace LeafIROut with LeafIR
      !      IRout = LeafIROut(Tleaf, Eps)
      IRout = LeafIR(Tleaf, Eps)

      RETURN
END SUBROUTINE LeafEB
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Convert water mixing ratio (kg/kg) to water vapor pressure 
!   (Pa or Kpa depending on units of input )
!   Mixing ratio (kg/kg), temp (C), pressure (KPa)
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION WaterVapPres(Dens, Pres, WaterAirRatio)
      IMPLICIT NONE
      REAL :: Dens, Pres, WaterVapPres, WaterAirRatio
      WaterVapPres = (Dens / (Dens + WaterAirRatio)) * Pres
END FUNCTION WaterVapPres
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   FUNCTION Stability
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION Stability(Canopychar, Cantype, Solar, NrCha, NrTyp)
      IMPLICIT NONE
      INTEGER :: Cantype, NrCha, NrTyp
      REAL :: Solar, Trateboundary, Stability
      REAL,DIMENSION(NrCha, NrTyp)  :: Canopychar

      Trateboundary = 500
      IF (Solar > Trateboundary) THEN
        ! Daytime temperature lapse rate
        Stability = Canopychar(12, Cantype)
      ELSEIF (Solar > 0) THEN
        Stability = Canopychar(12, Cantype) - ((Trateboundary - Solar) / Trateboundary) * (Canopychar(12, Cantype) - Canopychar(13, Cantype))
      ELSE
         ! Nightime temperature lapse rate
         Stability = Canopychar(13, Cantype)
      ENDIF
END FUNCTION Stability
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Saturation vapor density  (kg/m3)
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION ConvertHumidityPa2kgm3(Pa, Tk)
      implicit none
      REAL              :: ConvertHumidityPa2kgm3, Pa, Tk
      ConvertHumidityPa2kgm3 = 0.002165 * Pa / Tk
END FUNCTION ConvertHumidityPa2kgm3
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Leaf stomatal cond. resistance s m-1
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION ResSC(Par)
      IMPLICIT NONE
      REAL :: Par, SCadj, ResSC
      SCadj = ((0.0027 * 1.066 * Par) / ((1 + 0.0027 * 0.0027 * Par**2.)**0.5))
      IF (SCadj < 0.1) THEN
        ResSC = 2000
      ELSE
        ResSC = 200 / SCadj
      ENDIF
END FUNCTION ResSC
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Calculate IR transfer between leaf and air
!   Added by Alex Guenther and Ling Huang to replace previous
!   MEGAN2.1 IR balance functions
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION LeafIR(Tk, Eps)
       IMPLICIT NONE
       REAL :: Eps, Tk, LeafIR
       LeafIR = Eps * Sb * (2 * (Tk**4.))
END FUNCTION LeafIR
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Latent Heat of vaporization(J Kg-1) from Stull p641
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION LHV(Tk)
      IMPLICIT NONE
      REAL :: Tk, LHV
      ! REVISE - Replace 273 with 273.15
      !      LHV = 2501000 - (2370 * (Tk - 273))
      LHV = 2501000 - (2370 * (Tk - 273.15))
END FUNCTION LHV
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Latent energy term in Energy balance
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION LeafLE(Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType)
      IMPLICIT NONE
      REAL :: Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType, LeafRes, Vapdeficit, LeafLE, LE !SvdTk,
      LeafRes    = (1 / (1.075 * (GH / 1231))) + StomRes
      Vapdeficit = (SvdTk(Tleaf) - Ambvap)
      ! Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) / 
      !                 leaf resistence (s m-1)
      LE = TranspireType * (1 / LeafRes) * LatHv * Vapdeficit
      IF (LE < 0) THEN
        LeafLE = 0
      ELSE
        LeafLE = LE
      ENDIF
END FUNCTION  LeafLE
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Boundary layer conductance
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION LeafBLC(GHforced, Tdelta, Llength)
      IMPLICIT NONE
      REAL :: GHforced, Tdelta, Llength, Ghfree, LeafBLC
      !----------------------------------------------------------------
      ! This is based on Leuning 1995 p.1198 except using molecular 
      ! conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
      ! diffusivity so that you end up with a heat convection coefficient 
      ! (W m-2 K-1) instead of a conductance for free convection
      IF (Tdelta >= 0) THEN
         GhFree = 0.5 * 0.00253 * ((160000000 * Tdelta / (Llength**3.))**0.25) / Llength
      ELSE
        GhFree = 0
      ENDIF
      LeafBLC = GHforced + GhFree
END FUNCTION LeafBLC
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Convective energy term in Energy balance (W m-2 heat flux 
!      from both sides of leaf)
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
FUNCTION LeafH(Tdelta, GH)
      IMPLICIT NONE
      REAL :: Tdelta, GH, LeafH
      LeafH = 2 * GH * Tdelta! 2 sides X conductance X Temperature gradient
END FUNCTION LeafH
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   Saturation vapor density  (kg/m3)
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
function SvdTk(Tk)
      IMPLICIT NONE
      REAL :: Tk, Svp, SvdTk
      ! Saturation vapor pressure (millibars)
      Svp = 10**((-2937.4 / Tk) - (4.9283 * LOG10(Tk)) + 23.5518)
      SvdTk = 0.2165 * Svp / Tk
end function  SvdTk


end subroutine megan

END MODULE MEGAN_V32
