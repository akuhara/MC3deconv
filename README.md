# MC3deconv: Multi-Channel deconvolution by reversible-jump Markov-Chain Monte Carlo

One of the purposes of seismology is to investigate the Earth's subsurface structure using seismic waveforms. The receiver function method extracts useful signals (i.e., P-to-S converted phases) from the teleseismic waveforms by deconvolving P component records from the corresponding SV (or SH) components. Despite the many successful applications, conventional receiver function methods often fail due to numerical instability of the deconvolution and strong multiples on the P components. 

The technique developed here, MC3deconv, nicely overcomes these issues. The method optimizes the equation of multichannel deconvolution, in which two components of the reciver-sided Green's functions are related directly without deconvolution. To regularize the inverse problem, these Green's functions are expressed in the form of successive pulses. The number of pulses, their timing, and amplitudes are inverted using Bayesian techniques, the reversible-jump Markov-chain Monte Carlo and the Parallel Tempering.

(C) 2018 Takeshi Akuhara

Email: akuhara @ eri.u-tokyo.ac.jp

## Terms of use

* Please cite [Akuhara et al. (2019)](#Reference) when you publish an article or making presentation using this method.
* Also, make it clear that where readers or audiences can download this program package: you may put the link to the Github repository (https://github.com/akuhara/MC3deconv).
* Any bug reports are welcome! Looking forward to hearing your experience. 

## Limitations so far

* Although this method is designed to retrieve both radial (R) and vertical (Z) components of Green's functions, our experience suggest that the estimated Z-component is not so reliable as the R-component.
* This method assumes Gaussian noise without temporal correlation. This simplified treatment often leads to overfitting. 

## How to install

Use `make` command in the root directory of this package. 

* `mpifort` must be linked to the GNU fortran compiler (i.e., `gfortran`).
* If one wishes to use the Intel compiler (i.e., `ifort`), some modification is necessary in Makefile. 
* An executable file, `mc3deconv`, is created under the `bin` directory.

## How to run

`mpirun -np [N_proc] (path to the root directory of this package)/bin/mc3deconv`, for example. 
 * N_proc: Number of processes for parallel computation (must be >= 2, see the note below). 
 * A parameter file named "params.in" must exist in the current directory. 

The easiest way to test is:

`cd sample1`

`mpirun -np [N_proc] ../bin/mc3deconv`

In the `sample1` directory, all necessary data and parameter files are already prepared.

### Note on parallel computation

This program requires parallel computation. One of the processes is used to control the other processes, not performing MCMC sampling at all. Therefore, it is mandatory to use more than two processes.  

## Input files

### Parameter file (params.in)

A parameter file, which sets tuning parameters and input data, etc., must exist in the working directory from which `mc3deconv` is called, with the name "params.in". The format of the parameter file is as below, but you can put comment lines that start with "#" if necessary.

#### Format

|Line #|1st column|2nd column|
|:--:|:--:|:--:|
|1| Number of iterations in burn-in period|-|
|2| Number of iterations in sampling period|-|
|3| Number of iterations per generating one sample|-|
|4| Random number seed|-|
|5| Number of McMC chains per process (for parallel tempering)|-|
|6| Number of non-tempered chains|-|
|7| Maximum temperature|-|
|8| Input Z component file (in SAC format) | -|
|9| Input R component file (in SAC format)|-|
|10| Sampling interval of input data (sec)|-|
|11| Start time of the analysis window relative to file beginning (s)| End time of the analysis window relative to file beginning (s) |
|12| Lower prior limit for # of pulses | Upper limit for # of pulses|
|13| Lower prior limit for Z amplitude | Upper limit for Z amplitude|
|14| Lower prior limit for R amplitude | Upper limit for R amplitude|
|15| Lower prior limit for pulse timing relative to direct P arrival (s) | Upper limit for pulse timing relative to direct P arrival (s)|
|16| Probability of birth proposal (adding a pulse) |
|17| Probability of death proposal (removing a pulse)|
|18| Probability of time-shit proposal |
|19| Probability of amplitude-perturb proposal |
|20| Standard deviation to perturb Z amplitude | Standard deviation to perturb R amplitude|
|21| Standard deviation to newly generate Z amplitude | Standard deviation to newly generate R amplitude |
|22| Standard deviation to shift timing (s) ||
|23| Total time length of output (s) | -|
|24| Acausal time length preceding direct P arrival for output (s) | - |
|25| Factor of Gaussian low-pass filter | - |
|26| Minimum amplitudes for output| Maximum amplitudes for output|
|27| Amplitude bin width for output| - |

#### Example 

You can find an example of `params.in` in the `sample1` directory.

### Data file (user's given name)

#### Format

* Input waveform data should be SAC format.

#### Example

Examples of data files, `syn.r` and `syn.z` are installed in the `sample1` directory. 

## Output files

Five output files are created after running the program.

### dim.ppd 

The posterior probability distribution of the number of pulses. 

#### Format
|1st column|2nd column|
|:--:|:--:|
|# of pulses|probability|


### Gr.ppd / Gz.ppd 
The posterior probability distribution of R and Z component Green's functions. The format is:

#### Format
|1st column|2nd column|3rd colmun|
|:--:|:--:|:--:|
|time after P (s)|amplitude|probability|`

### Gr.mean / Gz.mean

The mean models of R and Z component Green's functions.

#### Format

`Gr.mean` and `Gz.mean` are written in SAC format.

## Sample datasets

There are two sample datasets, which may be useful for testing the program. These are the same data as used in [Akuhara et al. (2019)](#Reference). 

* Sample1: Synthetic data in `sample1`
* Sample2: Real OBS data in `sample2`

## Acknowledgments

Developing this package is supported by JSPS KAKENHI Grant Number JP17H06604. OBS data in the sample2 directory are collected by K. Nakahigashi and T. Yamada, under the program "Integrated Research Project on Seismic and Tsunami Hazards Around the Sea of Japan" of the Mistry of Education, Culture, Sports, Science and Technology (MEXT), Japan. This package uses a fortran program, [mt19937.f90](https://gist.github.com/ykonishi/5569005).

## Reference

* T. Akuhara, M. Bostock, A. Plourde, M. Shinohara (2019) Beyond Receiver Functions: Green's Function Estimation by Trans-Dimensional Inversion and Its Application to OBS Data, _Journal of Geophysical Research: Solid Earth_, https://doi.org/10.1029/2018JB016499  

