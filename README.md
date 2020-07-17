# astrosphere.cpp
#### an astrosphere problem generator written by Greg Szypko for Athena++

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

## Installing our Problem Generator and Configuring/Compiling/Running
1. Add ```astrosphere.cpp``` from this repository to the directory ```src/pgen/```
2. Run the ```configure.py``` script, followed by the pertinent options, to configure the Makefile. ```configure.py -h``` will list all possible options.
    1. For the astrosphere problem generator, run ```configure.py --prob astrosphere -mpi``` to configure using MPI parallelism. ```configure.py --prob astrosphere -mpi -b``` will turn on magnetic fields (not yet used successfully).
3. ```make clean```, then ```make``` to compile Athena++.
    1. Note that if you do not have OpenMP (shared memory parallelism) installed and included in your configuration, compiling will flood your terminal with ```unrecognized OpenMP #pragma``` warnings. This should not cause any problems, aside from being a nuisance.
4. The compiled executable is located at ```bin/athena```.
