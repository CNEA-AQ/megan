# Model of Emissions of Gases and Aerosols from Nature (MEGANv3.3)

> `megan` is a model used to estimate biogenic VOCs emissions, it has been developed by Alex Guenther and his group. 

## Dependencies:

+  Fortran GNU compiler
+  NetCDF library

> &#9888; At this moment it only works on UNIX/Linux O.S. 

## Get MEGAN input data:

In order to run megan you will need meteorological and land data. 

Data required:
+ WRF output NetCDF file. <!-- with the following variables: 'XLAT', 'XLONG', 'Times', 'MAPFAC_M', 'ISLTYP', 'U10', 'V10', 'T2', 'SWDOWN', 'PSFC', 'Q2', 'RAINNC', 'LAI', 'SMOIS', 'TSLB'. -->
+ Land, vegetation and soil fields needed to run `prep_megan`. Global data is freely available from [MEGAN GLOBAL DATA](https://drive.google.com/drive/folders/1ZdohMA4f4O_Yd2HttMLjhGpgFbTQMlt0?usp=sharing).

## Build
Go to the `src` directory:

`> cd src/`

Edit the Makefile to set the compiler and path to NetCDF lib and include files.

`> make`

If the compilation is successful, the executable `meganv3.3.exe` should be created at the `./exe` directory.

> &#9888; Note that the variables must be adjusted to match the appropriate values for your system. Check your `nc-config --libdir` and `nc-config --includedir` to fill NetCDF flags on Makefile.

## Run

Edit the namelist `namelist_megan` that contains the following variables:

```fortran
&megan_nl
   start_date='2019-01-01 18:00:00',  !YYYY-MM-DD HH:MM:SS
   end_date  ='2019-01-02 07:00:00',  !YYYY-MM-DD HH:MM:SS

   met_files='wrfout_d01_<date>_<time>', !Pattern of meteo files paths

   LSM='NOAH',                         ! land surface model (LSM) used in meteo

  mechanism='CB05',                    !'RACM2','CB05',CB6_ae7','CRACMM','SAPRC07',

   static_file='prep_mgn_static.nc',   ! NetCDF file w/ CTF, EFs, LDF & LAND  layers (created by prep_megan)
  dynamic_file='prep_mgn_dynamic.nc',  ! NetCDF file w/ LAI, NDEP & NFERT leyers     (created by prep_megan)

  run_bdsnp  = .false. ,               ! run bdsnp soil model for NO emissions?
  prep_megan = .false.,                ! run prep_megan?
/
&prep_megan_nl
   griddesc='GRIDDESC',                        !GRIDDESC file (describing grid and proj)
   gridname='MERC_TEST',!'LCC_TAN_TEST',       !gridname to use in GRIDDESC file

   eco_glb='input/veg_Ecotypes.nc',            ! global ecotype
   ctf_glb='input/veg_GrowthFormFracions.nc',  ! global canopy type fraction
   lai_glb='input/veg_LAIv.nc',                ! global leaf area index

   clim_glb='input/soil_climate.nc',           ! global arid/nonarid soils
   land_glb='input/soil_landtype.nc',          ! global land type data
   ndep_glb='input/soil_nitro.nc',             ! global N-deposition flux
   fert_glb='input/soil_fert.nc',              ! global N-fertilization flux

   GtEcoEF="db/GtEFbyEcotype.csv",             !Emission factor of each GT grouped by Ecotype
/
```

> &#9888; Note that the variables must be adjusted to match the appropriate values for your run.

Then execute `megan_v3.3.exe`:

`> ./megan_v3.3.exe < namelist_megan` 

Please feel free to contact the developer if you have any issues or suggestions.


---
## Planned future improvements:

+ [ ] Portability. 
  - [ ] Find an alternative to sys calls for `date` function.
+ [ ] Input/Output. 
  - [x] Incorporate prep-megan.
  - [x] Reading multiple meteoroligcal files.
  - [ ] Support to others meteorological models.
+ [ ] Science.
  - [ ] Calculate one gamma per CANTYPE. Then use the corresponding EF for this gamma.
  - [ ] Implement bdsnp.
