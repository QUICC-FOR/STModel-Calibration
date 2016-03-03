#!/bin/bash

# start the full model
~/mcmc_4state/bin/stm4_mcmc1.5 -d -p dat/mcmc_inits3.txt -t dat/mcmc_calib.txt -o res/ch3 -i 200000 -c 20 -l 5 -v 2 2>res/ch3/log.txt &

