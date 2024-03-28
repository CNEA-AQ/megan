# Model of Emissions of Gases and Aerosols from Nature (MEGAN)

> `megan` is a model to estimate biogenic VOCs emissions, it has been developed by Alex Guenther and his group. The code on this repository is an UN-OFFICIAL re-adaptation of that model.

## Dependencies:

+  Fortran GNU compiler
+  NetCDF library

> (!) At this moment it only works on UNIX/Linux O.S. 

## Get MEGAN input data:

In order to run megan you will need meteorological and land data. 

Data required:
+ WRF output NetCDF file. <!-- with the following variables: 'XLAT', 'XLONG', 'Times', 'MAPFAC_M', 'ISLTYP', 'U10', 'V10', 'T2', 'SWDOWN', 'PSFC', 'Q2', 'RAINNC', 'LAI', 'SMOIS', 'TSLB'. -->
+ prepmegan2cmaq output NetCDF files. <!--Land and vegetation fields: LAI, EF, LDF, CT -->
<!-- + For NO calculation you will also need: LAND, Arid/Non-Arid, Fert, NO production. -->

## Build
Go to the `src` directory:

`> cd src/`

Edit the Makefile to set the compiler and path to NetCDF lib and include files. Check your `nc-config --libdir` and `nc-config --includedir`.

`> make`

If the compilation is successful, the executable `meganv3.3.exe` should be created.

## Run

Edit the namelist `namelist_megan` that contains the following variables:

```fortran
&megan_nl
   start_date='2022-10-19 19:00:00',  !YYYY-MM-DD_HH:MM:SS
   end_date  ='2022-10-20 02:00:00',  !YYYY-MM-DD_HH:MM:SS

   met_files='wrfout_test',           !'wrfout_d01_2022-10-19_18:00:00',!
   LSM='NOAH',                        ! Land Surface Model (LSM) corresponding w/meteo 'ISLTY' variable

   ctf_file= 'CT3.nc',
   lai_file= 'LAI3.nc',
   ef_file=  'EF.nc',
   ldf_file= 'LDF.nc',
   ndep_file='NDEP.nc',
   fert_file='FERT.nc',
   land_file='LAND.nc',

   mechanism='CB6',                   !'RACM2', !'CB05', ! 'B6_ae7', !'CRACMM'  ! 'SAPRC07', !
/
```

> (!) Note that the variables must be adjusted to match the appropriate values for your system.

Then execute `megan_v3.3.exe`:

`> ./megan_v3.3.exe < namelist_megan` 

Please feel free to contact the developer if you have any issues or suggestions.

## Planned future improvements:

+ [ ] Portability. 
  - [ ] Find an alternative to sys calls for 'date' function.
+ [ ] Input/Output. 
  - [ ] Support to others meteorological models.
  - [ ] Reading multiple meteoroligcal files.
