# MC3 deconvolution
This package provides a open code to perform deconvolution of teleseismic waveforms (i.e., calculation of receiver-functions) via multichannel deconvolution by rj-McMC (MC3-deconvolution). 

## How to install

Compiling soce codes can be done by simiply typing, `make` in the main directory. This process assumes that `mpif90` is available and that the `mpif90` links to GNU frotran compiler (i.e., `gfortran`). Alternatively, you can use `mpif90` that links to Intel compiler (i.e., `ifort`). In this case, some modification is necessary in Makefile (see Makefile for details).  

## How to run
Run the program by `mpirun -np (# of processes) bin/MC3deconv`. Since one of the processes are used only to control temperature swapping between McMC chains, at least two processes are needed. Also, a file, "params.in", must be exist in the current direcory, in which user-given parameters are to be specified.


## Input

## Output
---
# Samples
## Sample1: Synthetic data
## Sample2: Real OBS data
