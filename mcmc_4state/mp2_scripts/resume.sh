#/bin/bash

BIN=bin/stm4_1.5.1
DIR=~/STModel-Calibration/mcmc_4state

module load gcc/4.8.2
module load gsl64/1.16


# run the models
$BIN -r $DIR/res/full/resumeData.txt -t $DIR/dat/mcmc_calib.txt -i 10000 2>$DIR/res/full/log.txt &
$BIN -r $DIR/res/int/resumeData.txt -t $DIR/dat/mcmc_calib.txt -i 10000 2>$DIR/res/int/log.txt &

wait