#/bin/bash

cat "Run the submit script interactively"
exit 1


# submit the model initially on qfat
qsub -q qfat256 -l walltime=1:00:00 -l nodes=2:ppn=1 ~/STModel-Calibration/mcmc_4state/mp2_scripts/start.sh

# resume on qfat
qsub -q qfat256 -l walltime=1:00:00 -l nodes=2:ppn=1 ~/STModel-Calibration/mcmc_4state/mp2_scripts/resume.sh


# same scripts on qwork - REMEMBER TO CHANGE THE NUM CORES IN THE SCRIPTS AND RESUMEDATA
qsub -q qwork -l walltime=1:00:00 -l nodes=2:ppn=1 ~/STModel-Calibration/mcmc_4state/mp2_scripts/start.sh

qsub -q qwork -l walltime=1:00:00 -l nodes=2:ppn=1 ~/STModel-Calibration/mcmc_4state/mp2_scripts/resume.sh
