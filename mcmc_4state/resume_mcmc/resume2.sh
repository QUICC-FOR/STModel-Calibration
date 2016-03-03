#!/bin/bash

# resume the full model
~/mcmc_4state/bin/stm4_mcmc1.5 -r res/ch2/resumeData.txt -d  -t dat/mcmc_calib.txt -i 200000 2>res/ch2/log.txt &
