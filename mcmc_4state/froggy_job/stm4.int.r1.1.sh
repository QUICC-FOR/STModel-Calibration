#!/bin/bash
#OAR -n stm4.int.r1.1.p
#OAR -l /nodes=1/core=16,walltime=30:00:00
#OAR --project teembio
#OAR --stdout stm4.int.r1.1.out.log
#OAR --stderr stm4.int.r1.1.err.log

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
name=int_r1_1
calib=$DAT/mcmc_calib_r1.txt
cores=16
iter=100000
##=================

# variables set from run vars
inits=$DAT/mcmcInits_"$name".txt
outdir=$DIR/res/"$name"_par
resume="$outdir"/resumeData.txt

# this is where stdout etc will get written
cd $BASE

# make sure all directories exist
mkdir -p $outdir

# run the models
$BIN -r $resume  -t $calib -i $iter &

wait