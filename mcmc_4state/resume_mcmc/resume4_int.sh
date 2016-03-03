#!/bin/bash

# resume the intercept-only model
~/mcmc_4state/bin/stm4_mcmc1.5 -r res/ch4_int/resumeData.txt -d  -t dat/mcmc_calib.txt -i 200000 2>res/ch4/log.txt &
