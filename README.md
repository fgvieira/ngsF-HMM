

# ngsF-HMM

`ngsF-HMM` is a program to estimate per-individual inbreeding tracts using a two-state Hidden Markov Model (HMM). Furthermore, instead of using called genotypes, it uses a probabilistic framework that takes the uncertainty of genotype's assignation into account; making it specially suited for low-quality or low-coverage datasets.


### Citation

`ngsF-HMM` was been submitted to [Bioinformatics](http://bioinformatics.oxfordjournals.org/), so please cite it if you use it in your work:

    Vieira FG, Albrechtsen A, Gilbert MT and Nielsen R
    Estimating IBD tracts from low coverage NGS data
    Bioinformatics (submitted)


### Installation

`ngsF-HMM` can be easily installed but has some external dependencies:

* `zlib`: v1.2.7 tested on Debian 7.8 (wheezy)
* `gsl` : v1.15 tested on Debian 7.8 (wheezy)
* `md5sum`: only needed for `make test`

To install the entire package just download the source code:

    % git clone https://github.com/fgvieira/ngsF-HMM.git

and run:

    % cd ngsF-HMM
    % make

To run the tests:

    % make test

Executables are built into the main directory. If you wish to clean all binaries and intermediate files:

    % make clean

### Usage

    % ./ngsF-HMM [options] --n_ind INT --n_sites INT --glf glf/in/file --out output/file

#### Parameters

* `--geno FILE`: Input genotype file.
* `--pos` FILE: Input site coordinates file
* `--lkl`: Input are genotype likelihoods (BEAGLE format).
* `--loglkl`: Input are genotype log-likelihoods.
* `--n_ind INT`: Sample size (number of individuals).
* `--n_sites INT`: Total number of sites.
* `--call_geno`: Call genotypes before running analyses.
* `--freqs DOUBLE or CHAR`: Initial frequency values. Can be defined by user as a DOUBLE, (r)andom, (e)stimated or read from a FILE.
* `--freqs_fixed`: Allele frequencies are fixed (will not be optimized).
* `--indF DOUBLE-DOUBLE or CHAR`: Initial inbreeding and transition probability values. Can be defined by user as a DOUBLE-DOUBLE, (r)andom, or read from a FILE.
* `--indF_fixed`: Inbreeding and transition probability values are fixed (will not be optimized).
* `--path DOUBLE or CHAR`: Initial state paths. Can be (r)andom, (u)niform at 0.01, or read from a FILE.
* `--path_fixed`: State paths are fixed (will not be optimized).
* `--out CHAR`: Prefix for output files.
* `--log INT`: Dump LOG file each INT iterations. [0]
* `--log_bin`: Dump LOG file in binary.
* `--min_iters INT`: Minimum number of EM iterations. [10]
* `--max_iters INT`: Maximum number of EM iterations. [100]
* `--min_epsilon FLOAT`: Maximum RMSD between iterations to assume convergence. [1e-5]
* `--n_threads INT`: Number of threads to use. [1]
* `--version`: Prints program version and exits.
* `--verbose INT`: Selects verbosity level. [1]
* `--seed INT`: Seed for random number generator.

### Input data
As input `ngsF-HMM` reads a Genotype Likelihood (GL) file composed of 3 genotype likelihoods, per site and individual. The file can be either in binary (as doubles) or gzipped text (sites on rows). The latter case allows for initial non-genotype columns.

### Stopping Criteria
An issue on iterative algorithms is the stopping criteria. `ngsF-HMM` implements a dual condition threshold: relative difference in log-likelihood and estimates RMSD (F and freq). As for which threshold to use, simulations show that 1e-5 seems to be a reasonable value. However, if you're dealing with low coverage data (2x-3x), it might be worth to use lower thresholds (between 1e-6 and 1e-9).

To avoid convergence to local maxima, ngsF-HMM should be run several times from different starting points. To make this task easier, a script (`ngsF-HMM.sh`) is provided that can be called with the exact same parameters as `ngsF-HMM`.

### Thread pool
The thread pool implementation was adapted from Mathias Brossard's and is freely available from:
https://github.com/mbrossard/threadpool
