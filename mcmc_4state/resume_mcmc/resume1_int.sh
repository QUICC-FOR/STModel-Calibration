#!/bin/bash

# resume the intercept-only model
~/mcmc_4state/bin/stm4_mcmc1.5 -r res/ch1_int/resumeData.txt -d  -t dat/mcmc_calib.txt -i 200000 2>res/ch1/log.txt &
