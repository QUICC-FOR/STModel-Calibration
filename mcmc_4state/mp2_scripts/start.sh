#/bin/bash

BIN=bin/stm4_1.5.1
DIR=~/STModel-Calibration/mcmc_4state

module load gcc/4.8.2
module load gsl64/1.16

# some example qsubs for mp2
# this uses qwork, which will yield 24 cores
# qsub -q qwork -l walltime=1:00:00 -l nodes=1:ppn=1 myscript.sh

# qfat256 has 48 cores per node; this gives all 48
# qsub -q qfat256 -l walltime=1:00:00 -l nodes=1:ppn=1 myscript.sh

# for this script, use 2 nodes (one for int, one for full) with 48 cores each, if avail
# qsub -q qfat256 -l walltime=1:00:00 -l nodes=1:ppn=1 myscript.sh


# run the models
$BIN -d -p $DIR/dat/inits_full.txt -t $DIR/dat/mcmc_calib.txt -o $DIR/res/full -n 25 -i 10000 -b 10000 -c 48 -l 5 -v 2 2>$DIR/res/full/log.txt &
$BIN -d -p $DIR/dat/inits_int.txt -t $DIR/dat/mcmc_calib.txt -o $DIR/res/int -n 25 -i 10000 -b 10000 -c 48 -l 5 -v 2 2>$DIR/res/int/log.txt &

wait