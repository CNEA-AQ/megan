# Model of Emissions of Gases and Aerosols from Nature (MEGAN)

> `megan` is a model used to estimate biogenic VOCs emissions, it has been developed by Alex Guenther and his group. The code on this repository is an UN-OFFICIAL re-adaptation of that model.

## Dependencies:

+  Fortran GNU compiler
+  NetCDF library

> &#9888; At this moment it only works on UNIX/Linux O.S. 

## Get MEGAN input data:

In order to run megan you will need meteorological and land data. 

Data required:
+ WRF output NetCDF file. <!-- with the following variables: 'XLAT', 'XLONG', 'Times', 'MAPFAC_M', 'ISLTYP', 'U10', 'V10', 'T2', 'SWDOWN', 'PSFC', 'Q2', 'RAINNC', 'LAI', 'SMOIS', 'TSLB'. -->
+ [prepmegan4cmaq](https://github.com/CNEA-AQ/prepmegan4cmaq) output NetCDF files. <!--Land and vegetation fields: LAI, EF, LDF, CT -->
<!-- + For NO calculation you will also need: LAND, Arid/Non-Arid, Fert, NO production. -->

## Build
Go to the `src` directory:

`> cd src/`

Edit the Makefile to set the compiler and path to NetCDF lib and include files.

`> make`

If the compilation is successful, the executable `meganv3.3.exe` should be created at `./exe` directory.

> &#9888; Note that the variables must be adjusted to match the appropriate values for your system. Check your `nc-config --libdir` and `nc-config --includedir` to fill NetCDF flags on Makefile.

## Run

Edit the namelist `namelist_megan` that contains the following variables:

```fortran
&megan_nl
   start_date='2022-10-19 19:00:00',  !YYYY-MM-DD HH:MM:SS
   end_date  ='2022-10-20 02:00:00',  !YYYY-MM-DD HH:MM:SS

   met_files='wrfout_d01_<date>_<time>', !'wrfout_d01_2022-10-19_18:00:00',!
   LSM='NOAH',                  ! Land Surface Model (LSM) used in meteo

   ctf_file= 'CT3.nc',          ! canopy-type fractions
   lai_file= 'LAI3.nc',         ! leaf area index
   ef_file=  'EF.nc',           ! emission factors
   ldf_file= 'LDF.nc',          ! light-dependent factors
   ndep_file='NDEP.nc',         ! N-deposition
   fert_file='FERT.nc',         ! N-fertilization flux
   land_file='LAND.nc',         ! Land data

   mechanism='CB05',            ! mechanisms: 'RACM2','CB05','CB6_ae7','CRACMM'','SAPRC07',
/
```

> &#9888; Note that the variables must be adjusted to match the appropriate values for your run.

Then execute `megan_v3.3.exe`:

`> ./megan_v3.3.exe < namelist_megan` 

Please feel free to contact the developer if you have any issues or suggestions.

---
## Planned future improvements:

+ [ ] Portability. 
  - [ ] Find an alternative to sys calls for 'date' function.
+ [ ] Input/Output. 
  - [x] Reading multiple meteoroligcal files.
  - [ ] Support to others meteorological models.
+ [ ] Science.
  - [ ] Calculate one gamma per CANTYPE. Then use the corresponding EF for this gamma.
  - [ ] Implement bdsnp.
