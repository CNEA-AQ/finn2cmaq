# FINN Emissions

> The executable `finn2cmaq`, is a fortran program that grids FINNv2 fire emissions to the required input formats and grids for CMAQ model. 

## Dependencies
 - Fortran compiler (only `gfortran` has been tested yet)
 - NetCDF library
 - `date` command (commonly included in linux toolkits) 
 - GDAL/OGR (for projection transformations. this dependency will be removed in the future)


## Get FINN data

Emissions files are produced by Christine Wiedinmyer and are available from [NCAR ACOM FINN webpage](https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/)

A shell script for automated downloading of this data is included in this repo:
`> ./get_finn_data.sh`

## Build
`> cd src`

Edit `Makefile` to set compiler and path to NetCDF lib and include files.

`> make`

## Run
`> cd run `

Edit the 'example.inp' that contains the following variables:

```
&control

  start_date="2019-01-01",       !%Y-%m-%d (YYYY-MM-DD)
    end_date="2019-01-01",       !%Y-%m-%d (YYYY-MM-DD)
  finn_data_directory = './finn_data/'
  griddesc_directory="./GRIDDESC",   !path to the GRIDDESC file
  chemistry="GEOSchem",
  diurnal_cycle=0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0300,0.0600,0.1000,0.1400,0.1700,0.1400,0.1200,0.0900,0.0600,0.0300,0.0043,0.0043,0.0043,0.0043,0.0043
/
```

Then execute finn2cmaq:
`> ../src/finn2cmaq < example.inp` 

