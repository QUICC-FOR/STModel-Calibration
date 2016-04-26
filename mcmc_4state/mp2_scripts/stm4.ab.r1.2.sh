#!/bin/sh
#PBS -q qwork
#PBS -l walltime=85:00:00
#PBS -l nodes=1:ppn=1
#PBS -r n

BASE=/mnt/parallel_scratch_mp2_wipe_on_august_2016/dgravel/mtalluto
DIR=$BASE/STModel-Calibration/mcmc_4state
BIN=$DIR/bin/stm4_1.5.1
DAT=$DIR/dat
LOG=$DIR/log

module load gcc/4.8.2
module load gsl64/1.16

# some example qsubs for mp2
# this uses qwork, which will yield 24 cores
# qsub -q qwork -l walltime=1:00:00 -l nodes=1:ppn=1 myscript.sh

# qfat256 has 48 cores per node; this gives all 48
# qsub -q qfat256 -l walltime=1:00:00 -l nodes=1:ppn=1 myscript.sh

# for this script, use 2 nodes (one for int, one for full) with 48 cores each, if avail
# qsub -q qfat256 -l walltime=1:00:00 -l nodes=1:ppn=1 myscript.sh

##=================
# VARIABLES to set for each run
name=ab_r1_2
calib=$DAT/mcmc_calib_r1.txt
cores=24
iter=100000
##=================

# variables set from run vars
inits=$DAT/mcmcInits_"$name".txt
outdir=$DIR/res/$name
logfile=$LOG/"$name".log

# this is where stdout etc will get written
cd $BASE

# make sure all directories exist
mkdir -p $outdir $LOG

# run the models
$BIN -d -p $inits -t $calib -o $outdir -i $iter -c $cores -l 5 -v 2 2>$logfile &

wait