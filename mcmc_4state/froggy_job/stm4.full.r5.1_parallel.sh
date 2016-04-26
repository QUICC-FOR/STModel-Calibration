#!/bin/bash
#OAR -n stm4.full.r5.1.p
#OAR -l /nodes=1/core=16,walltime=2:00:00
#OAR --project teembio
#OAR --stdout stm4.full.r5.1.p.out
#OAR --stderr stm4.full.r5.1.p.err

. /applis/ciment/v2/env.bash
module load ciment/devel_gcc-4.6.2
module load gcc/4.9.3_gcc-4.6.2
module load gsl/1.16_gcc-4.6.2

BASE=$HOME
DIR=$BASE/STModel-Calibration/mcmc_4state
BIN=$DIR/bin/stm4_1.5.1
DAT=$DIR/dat



##=================
# VARIABLES to set for each run
name=full_r5_1
calib=$DAT/mcmc_calib_r5.txt
cores=16
iter=10000
##=================

# variables set from run vars
inits=$DAT/mcmcInits_"$name".txt
outdir=$DIR/res/"$name"_par

# this is where stdout etc will get written
cd $BASE

# make sure all directories exist
mkdir -p $outdir

# run the models
$BIN -d -p $inits -t $calib -o $outdir -i $iter -c $cores -l 5 -v 2 &

wait