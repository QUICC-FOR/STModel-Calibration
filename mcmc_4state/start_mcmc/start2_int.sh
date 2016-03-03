#!/bin/bash

# start the intercept-only model
~/mcmc_4state/bin/stm4_mcmc1.5 -d -p dat/mcmc_inits_int.txt -t dat/mcmc_calib.txt -o res/ch2_int -i 200000 -c 20 -l 5 -v 2 2>res/ch2_int/log.txt &
