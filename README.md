# MC3 deconvolution
This package provides a open code to perform deconvolution of teleseismic waveforms (i.e., calculation of receiver-functions) via multichannel deconvolution by rj-McMC (MC3-deconvolution). 

## How to install

Compiling soce codes can be done by simiply typing, `make` in the main directory. This process assumes that `mpif90` is available and that the `mpif90` links to GNU frotran compiler (i.e., `gfortran`). Alternatively, you can use `mpif90` that links to Intel compiler (i.e., `ifort`). In this case, some modification is necessary in Makefile (see Makefile for details).  

## How to run
Run the program by:
`mpirun -np (# of processes) bin/MC3deconv`.

## Input

## Output