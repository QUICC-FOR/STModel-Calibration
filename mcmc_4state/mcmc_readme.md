### Raw Data Information

* dat/GenSA\_rf\_all\_LatLon\_2\_5y\_rep1.txt: GenSA maximum likelihood parameters from Feb 24 2016 (not used in the MCMC initialization, just kept for reference)
* calibration_data.rds: raw (pre-mcmc) calibration data from Feb 24 2016
* dat/scale_info.RData: scaling for the datAll file
* dat/inits.rds: matrix of 4 sets of starting parameters. Obtained by selecting a row at random from the first MCMC chain completed in June 2015, then selecting a row from each of the next four chains so that it maximized the distance/minimized the entropy among the sets

### MCMC Datasets

* dat/mcmc_calib.txt: the transition data for calibration
* dat/mcmc_inits(N).txt: initial settings for each of (N) chains

### Scripts

1. prep_mcmc.r: This prepares the data in a text format for the mcmc run. 
2. bin/: This is the directory for holding the compiled model executable. It should be compiled using the command `make fourstate` in the `STModel-MCMC` github repo, then the exexutable should be copied here (it will be placed in `./STModel-MCMC/bin/stm4_mcmc` by default). Use the `-h` option (e.g., `./bin/stm4 -h`) to get a list of options.

### Results
All results will be put in the `res/` directory. There should be one subdirectory per chain to avoid overwriting any files. This directory should be specified via the `-o dir` option (e.g., `./bin/stm4 -o res/chain1`).

The model will write several files to this directory. The two of main interest will be `posterior.csv` containing the posterior samples (one per row in CSV format) and the `resumeData.txt` file (see **resuming** below).

### Resuming the Sampler
If a sample run ends prematurely, it can be resumed using the `-r` option at the command line (it is necessary to give it the file name to resume from; e.g., `./bin/stm4 -r res/chain1/resumeData.txt`). Use the `-h` option to get information about what other flags are required when using `-r`.

<small>Updated 2 March 2016</small>