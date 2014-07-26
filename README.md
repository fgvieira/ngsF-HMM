

# ngsF-HMM

`ngsF-HMM` is a program to estimate per-individual inbreeding tracts using a two-state Hidden Markiv Model (HMM). Furthermore, instead of using called genotypes, it uses a probabilistic framework that takes the uncertainty of genotype's assignation into account.


### Installation

To install the entire package just download the source code:

    % git clone https://github.com/fgvieira/ngsF-HMM.git

To install these tools just run:

    % cd ngsF-HMM
    % make
    % make test

Executables are built into the main directory. If you wish to clean all binaries and intermediate files:

    % make clean

### Usage

    % ./ngsF-HMM [options] -n_ind INT -s INT -glf glf/in/file -out output/file

#### Parameters

* `-geno FILE`: Input genotype file.
* `-pos` FILE: Input site coordinates file
* `-lkl`: Input are genotype likelihoods (BEAGLE format).
* `-loglkl`: Input are genotype log-likelihoods.
* `-n_ind INT`: Sample size (number of individuals).
* `-n_sites INT`: Total number of sites.
* `-call_geno`: Call genotypes before running analyses.
* `-freqs DOUBLE or CHAR`: Initial frequency values. Can be defined by user as a DOUBLE, (r)andom, or read from a FILE.
* `-freqs_fixed`: Allele frequencies are fixed (will not be optimized).
* `-indF DOUBLE-DOUBLE or CHAR`: Initial inbreeding and transition probability values. Can be defined by user as a DOUBLE-DOUBLE, (r)andom, or read from a FILE.
* `-indF_fixed`: Inbreeding and transition probability values are fixed (will not be optimized).
* `-path DOUBLE or CHAR`: Initial state paths. Can be (r)andom, (u)niform at 0.01, or read from a FILE.
* `-path_fixed`: State paths are fixed (will not be optimized).
* `-out CHAR`: Prefix for output files.
* `-log`: Dump LOG file.
* `-log_bin`: Dump LOG file in binary.
* `-min_iters INT`: Minimum number of EM iterations. [10]
* `-max_iters INT`: Maximum number of EM iterations. [1500]
* `-min_epsilon FLOAT`: Maximum RMSD between iterations to assume convergence. [1e-5]
* `-n_threads INT`: Number of threads to use. [1]
* `-version`: Prints program version and exits.
* `-verbose INT`: Selects verbosity level. [1]
* `-seed INT`: Seed for random number generator.

### Input data
As input `ngsF-HMM` reads a Genotype Likelihood (GL) file composed of 3 genotype likelihoods, per site and individual. The file can be either in binary (as doubles) or gzipped text, in which case allows for initial non-genotype columns.

### Stopping Criteria
An issue on iterative algorithms is the stopping criteria. `ngsF-HMM` implements a dual condition threshold: relative difference in log-likelihood and estimates RMSD (F and freq). As for which threshold to use, simulations show that 1e-5 seems to be a reasonable value. However, if you're dealing with low coverage data (2x-3x), it might be worth to use lower thresholds (between 1e-6 and 1e-9).

### Thread pool
The thread pool implementation was adapted from Mathias Brossard's and is freely available from:
https://github.com/mbrossard/threadpool
