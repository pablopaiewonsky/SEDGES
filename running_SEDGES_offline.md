# How to run:

The code is in fortran 90. The driver for the model is a bourne shell script. The processed reanalysis data that is fed into the model must be placed in subdirectories of the same directory as the one which contains the fortran code and shell script. In the case of the simulation in which CO2 is varied, the text file with the yearly CO2 data (`co2data`) must also be in a subdirectory (`co2`) of the same directory as the fortran code and driver.

Two versions of the code are included: a historically-varying CO2 simulation, and a fixed CO2 simulation. The two sets of code are essentially the same, except for the drivers. The directory containing the historically-varying CO2 also includes the file with the yearly CO2 data (`co2data`). These data come from a combination of Antarctic ice core and Mauna Loa observatory sources. There are available, respectively, at <http://cdiac.ornl.gov/trends/co2/lawdome.html> and <ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_annmean_mlo.txt>.

The fortran code files are `globalvarmod.f90`, `sedges.f90`, `hydromod.f90`, and `main.f90`, and are summarized briefly as follows:

- `globalvarmod.f90` contains most of the variable declarations and initializations, i.e. the ones with global scope. Importantly, the value for atmospheric CO2 (co2) is set in the globalvarmod module.
- `sedges.f90` contains the vegetation model, SEDGES.
- `hydromod.f90` contains a simple hydrological model that simulates runoff and soil moisture. The latter of these is required by SEDGES as input.
- `main.f90` is the main program and handles the reading of the 6-hourly reanalysis data, boundary conditions, and saved end-of-year data; it also outputs monthly and end-of-year data in ascii format. main.f90 also calculates potential evapotranspiration (PET) for use by SEDGES and actual evapotranspiration (ET) for use by the hydrological model.

The gfortran compiler is used. Compilation, linking, and object creation are as follows (speed-optimized, using the -O3 option) for experiment name of `EXPNAME`:

```
gfortran -O3 globalvarmod.f90 sedges.f90 hydromod.f90 main.f90 -o main_$\{EXPNAME\}.out
```

`EXPNAME` was originally set to `MOSES86_0.001` and `MOSES86_0.001co2var`, respectively, for the fixed CO2 and historical CO2 simulations, but the user can change these names to whatever he or she prefers. The only restriction is that `EXPNAME`, as defined in the driver script, must match what `EXPNAME` is set to for the name of the above fortran executable file, `main_$\{EXPNAME\}.out`

The names of the driver scripts are `./driver_$\{EXPNAME\}.sh` and are, again, arbitrary.

In the driver script, `WORKDIR` (the directory in which everything is found) must be renamed to whatever is appropriate for the user. Also, the lengths of the strings of the full names of the input (i.e. reanalysis-derived) data directories depends on the length of the `WORKDIR` string and must be modified accordingly in `globalvarmod` (`swdir` through to `v60dir`) to match the user-specific directory name lengths. The variables, `MONMEANSBEGINYEAR` and `MONMEANSFINALYEAR` are, respectively, the beginning and end years for the simulation over which the climatological monthly means are computed. The model is run in cycles of 30 (the number of years of reanalysis data), and so `TOTALSETNUM * 30` gives the total number of simulation years. If starting from a bare soil state, the model needs 1500 years (including acceleration of the carbon cycle up to year 260) in order to come to an almost-equilibrium state.

The fortran object file, `main_$\{EXPNAME\}.out`, is executed once every year of the simulation. The driver script outputs a data2namelist file and a yearnamelist file (this latter one is updated, yearly, with year-dependent information, including the updated CO2 in the driver version that uses historical CO2).

Climate data operators (`cdo`) is required and is used by the driver script to manipulate and process output data files as the model runs and also afterward. Output data files are in ascii format (with a `.sra` suffix), but also in `srv` (service) format, for use by `cdo`. As the driver script runs, monthly mean and end of year files are generated in `.sra` and `.srv` formats. At the end of the model run, climatological monthly means in srv format are created for the last 30 years of the simulation.

Code is included to handle `z0` (surface roughness) due to orography (`z0climo`) if the user wants to compute the latent heat transfer incorporating this value. In that case, the user would need to comment out the line
```
dz0climo(:)= 0.
```
in `main.f90`. In the current code, however, the `z0climo` data is still read in, unfortunately. Although it hasn't been tested, commenting out the line
```
read(22,*) dz0climo(i),dz0climo(i+1),dz0climo(i+2),dz0climo(i+3),dz0climo(i+4),dz0climo(i+5)
```
in the readmasks subroutine in `main.f90` should get rid of the need to have an actual orographic surface roughness ascii file.

A lot of the code is based on (i.e. some variable names are the same and the formulation for aerodynamic conductance is adapted from) that of the Planet Simulator model, to which SEDGES has already been coupled (paper in preparation).
