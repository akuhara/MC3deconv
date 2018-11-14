# MC3deconv (Multi-Channel deconvolution by Markov-Chain Monte Calro

One of the purposes of seismology is to investigate the Earth's subsurface structure using seismic waveforms. The receiver function method extracts useful signals (i.e., P-to-S converted phases) from the teleseismic waveforms by deconvolving P component records from the corresponding SV (or SH) components. Despite the many successful applications, conventional receiver function methods often fail due to numerical instability of the deconvolution and strong multiples on the P components. 

The technique developed here, MC3deconv, nicely overcomes these issues. The method optimizes the equation of multichannel deconvolution, in which two components of the reciver-sided Green's functions are related directly without deconvolution. To regularize the inverse problem, these Green's functions are expressed in the form of successive pulses. The number of pulses, their timing, and amplitudes are inverted using Bayesian techniques, the reversible-jump Markov-chain Monte Carlo and the Parallel Tempering.

More details can be found in [the paper by T. Akuhara (currently in revision)].

### Limitations so far

* Development of MC3deconv originally aims to acquire both radial (R) and vertical (Z) components of Green's functions. However, our experiments empirically suggest that the Z-component is not estimated correctly, while it provides relatively good estimation for the R-component.
* MC3deconv assumes noise in input data is distributed according to Gaussian distribution with no correlation over the time domain. The actual noise in real seismograms, however, has temporal correlation undoubtly. 

These limitations, as well as advantages, are discussed in [the paper by T. Akuhara (currently in revision)].

---

## How to install
Use `make` command in the root directory of this package. 
* `mpif90` must be linked to the GNU fortran compiler (i.e., `gfortran`).
* If one wishes to use the Intel compiler (i.e., `ifort`), some modification is necessary in Makefile. See Makefile for more details.

## How to run

Use `mpirun -np [Nproc] bin/mc3deconv`. 
 * Nproc = Number of proccesses (must be >= 2). 
 * A parameter file named "params.in" must be exist in the working directory.
The easiest way to test is:

`cd sample1`

`mpirun -np 2 ../bin/mc3deconv`

---

## Input
A parameter file, which sets tuning parameters and input data, etc., must exist in the working directory from which `mc3deconv` is called, with the name "params.in". The format of the parameter file is as below, but you can put comment lines that start with "#" if necessary.

|Line #|parameter 1|parameter 2|Example value1| Example value2|
|:--:|:--:|:--:|:--:|:--|
|1| # of iterations in burn-in period|-|200000|-|
|2| # of iterations in sampling period|-|1200000|-|
|3| # of iterations per generating one sample|-|200|-|
|4| Random seed|-|12345678|-|
|5| # of McMC chains per process (for parallel tempering|-|10|-|
|6| # of non-tempered chains|-|2|-|
|7| Maximum temperature|-|100.0|-|
|8| Input Z component file (in SAC format) | -|input/syn_obs.z|-|
|9| Input R component file (in SAC format)|-|input/syn_obs.r|-|
|10| Sampling interval of input data (sec)|-|0.05|-|
|11| Start time of the analysis window relative to file beginning (s)| End time (s) |195.0|255.0|
|12| Lower prior limit for # of pulses | Upper limit|1|101|
|13| Lower prior limit for Z amplitude | Upper limit|-0.8|0.8|
|14| Lower prior limit for R amplitude | Upper limit|-0.8|0.8|
|15| Lower piror limit for pulse timing relative to P (s) | Upper limit|1.0|60.0|
|16| Probability of birth proposal (adding a pulse) |-|0.35|| 
|17| Probability of death proposal (removing a pulse)|-|0.35||
|18| Probability of time-shit proposal |-|0.05|-|
|19| Probability of amplitude-perturb proposal |-|0.25|-|
|20| Standard deviation to perturb Z amplitude | Same but for R amplitude| 0.05 | 0.05|
|21| Standard deviation to newly generate Z amplitude | Same but for R amplitude | 0.4 | 0.4|
|22| Standard deviation to shift timing (s) | -|0.1|-|
|23| Total time length of output (s) | -|63.0|-|
|24| Acausal time length preceding P in output (s) | - | 3.0 | -|
|25| Factor of Gaussian low-pass filter | - | 8.0 | -|
|26| Minimum amplitudes for output| Maximum | -2.2 | 2.2|
|27| Amplitude bin width for output| - | 0.02 | |


---

## Output

In the current version (ver. 0.0), five output files are created after running the program.

### dim.ppd (ascii format)

The posterior probability distribution of the number of pulses. The format is:

|1st column|2nd column|
|:--:|:--:|
|# of pulses|probability|


### Gr.ppd / Gz.ppd (ascii forat)
The posterior probability distribution of R and Z component Green's functions. The format is:

|1st column|2nd column|3rd colmun|
|:--:|:--:|:--:|
|time after P (s)|amplitude|probability|`

### Gr.mean / Gz.mean (SAC format)
The mean models of R and Z component Green's functions in [SAC](http://ds.iris.edu/files/sac-manual/) format.

---

## Samples
There are some sample data, which may be useful for testing the program.

* Sample1: Synthetic data in `sample1`
* Sample2: Real OBS data in `sample2`


## Acknowledgments

Developing this package is supported by JSPS KAKENHI Grant Number JP17H06604. OBS data in the sample2 directory are collected by K. Nakahigashi and T. Yamada, under the program "Integrated Research Project on Seismic and Tsunami Hazards Around the Sea of Japan" of the Mistry of Education, Culture, Sports, Science and Technology (MEXT), Japan. This package uses a fortran program, mt19937.f90', which is an open code to generate random numbers distributed under the GNU General Public License version 2.
