# FINN Emissions

> The executable `finn2cmaq`, is a fortran program that grids FINN fire emissions to the required input formats and grids for CMAQ model. 

## Dependencies
 - Fortran GNU compiler.
 - NetCDF library.

## Get FINN data

Emissions files are produced by Christine Wiedinmyer and are available from [NCAR ACOM FINN webpage](https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/)

A shell script for automated downloading this data is included in this repository:
`> ./get_finn_data.sh`

## Build
Edit the Makefile to set the compiler and path to NetCDF lib and include files. Check your nc-config --libdir and nc-config --includedir.

`> make`

If the compilation is successful, the executable `finn2cmaq.exe` should be created.

## Run

Edit the 'example.inp' that contains the following variables:

```
&control
  start_date="2019-01-01",       		!%Y-%m-%d (YYYY-MM-DD)
    end_date="2019-01-01",       		!%Y-%m-%d (YYYY-MM-DD)
  finn_data_directory = './finn_data/'		!path to finn files directory
  griddesc_path="./GRIDDESC",   		!path to the GRIDDESC file
  chemistry="GEOSchem",
  diurnal_cycle=0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0300,0.0600,0.1000,0.1400,0.1700,0.1400,0.1200,0.0900,0.0600,0.0300,0.0043,0.0043,0.0043,0.0043,0.0043
/
```
Note that the start_date, end_date, finn_data_directory, and griddesc_path variables must be adjusted to match the appropriate values for your system.

Then execute finn2cmaq:

`> finn2cmaq.exe < example.inp` 

Please feel free to contact the developer if you have any issues or suggestions.

## Planned future improvements:
 + [ ] Optional species mapping
 + [ ] Plume rise representation

