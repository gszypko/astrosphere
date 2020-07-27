# astrosphere.cpp
#### an astrosphere problem generator for Athena++

The goal of this project is to simulate the astrosphere of an arbitrary cool dwarf star, using parameters of that star that can be observed or inferred. As such, this problem generator allows the user to either input physical parameters of the stellar wind itself, or parameters of the host star that are then converted to wind parameters using different proposed models.

These files were created as part of my undergraduate senior thesis in Physics at Dartmouth College. This project was advised by Prof. Hans-Reinhard Müller.

## Setting up Athena++
To use this problem generator first requires a working installation of *Athena++*, which can be found at https://github.com/PrincetonUniversity/athena-public-version. More detailed installation information can be found in the accompanying wiki, but we will summarize the steps here.

1. Download the full ```athena-public-version/``` directory to your machine. All paths specified below are within this directory.
2. Ensure that your machine has a C++ compiler and Python 2.7+ or Python 3.4+. These are required for compilation and configuration of Athena++, respectively.
3. The following are optional to install on your machine (on the Discovery cluster, this involved loading the corresponding modules), but were used in the span of this project:
    1. An implementation of MPI-2. We used MPICH version 3.2.1. This allows for distributed-memory parallelization of the simulations.
4. Optionally, you can run ```tst/regression/run_tests.py``` to test the software under different a battery of different configurations.
    1. Note that if certain optional libraries are missing, parts of this test will yell at you. However, if you do not plan on using functionality related to those libraries, you can ignore this and proceed.
    2. For this project, FFTW version 3.3.7 (for self-gravity and turbulence) and HDF5 version 1.8.12 (for high-efficiency data output) were installed to pass this regression script, but were not used otherwise.

## Installing this Problem Generator and Configuring/Compiling
1. Add ```astrosphere.cpp``` from this repository to the directory ```src/pgen/```
2. Run the ```configure.py``` script, followed by the pertinent options, to configure the Makefile. ```configure.py -h``` will list all possible options.
    1. For the astrosphere problem generator, run ```configure.py --prob astrosphere -mpi``` to configure using MPI parallelism. ```configure.py --prob astrosphere -mpi -b``` will turn on magnetic fields (not yet used successfully).
3. ```make clean```, then ```make``` to compile Athena++.
    1. Note that if you do not have OpenMP (shared memory parallelism) installed and included in your configuration, compiling will flood your terminal with ```unrecognized OpenMP #pragma``` warnings. This should not cause any problems, aside from being a nuisance.
4. The compiled executable is located at ```bin/athena```.

## Input Files and Parameter Modes
1. To run without parallelization, simply execute ```athena -i athinput.inputfile``` where ```athinput.inputfile``` is the name of the input file to be used for the simulation run. This file allows you to specify the parameters of the simulation at run time.
    1. Included in this directory are several examples of these Input Files, set to solar parameters. See the Athena++ documentation to learn about the base functionality of these files.
2. To use the astrosphere problem generator, several additional parameters need to be specified under the ```<problem>``` header in the Input File. In all cases, the following need to be assigned a value:
    1. ```input_mode``` is an integer that specifies which model, if any, is used to convert stellar parameters to wind parameters. Default is 0, which is for direct wind parameter input (that is, no stellar model used). 1 is for the Wood model, 2 is for the Cohen model, 3 is for the Johnstone model, and 4 is for the Mesquita model (for M-dwarfs).
    2. ```vx_ISM```, ```vy_ISM```, ```vz_ISM```, ```d_ISM```, and ```p_ISM``` specify parameters of the interstellar medium flow: velocity components, density, and pressure, respectively. The ISM flows into the simulation domain from the lower x-boundary.
    3. ```radius``` is the outer radius of the overwriting region at the origin. ```radius_inner``` is the inner radius of the overwriting region. All wind variables are assigned to all points within ```inner_radius```, with a falloff between ```radius_inner``` and ```radius```.
    4. ```e_wind``` is the total energy density of the wind.
3. Depending on the specified ```input_mode```, additional parameters must be specified:
    1. Direct Input: ```mom_rad_wind``` and ```d_wind```, the magnitude of the radial momentum density and the density of the wind, respectively.
    2. Wood Model: ```R_star_over_solar``` and ```L_x_star_over_solar```, the radius and X-ray luminosity of the star, in solar units.
    3. Cohen Model:
    4. Johnstone Model:
    5. Mesquita Model: 
    
