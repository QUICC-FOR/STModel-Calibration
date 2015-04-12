#!/bin/bash
run(){
#Rscript 4-initForFit.R 0.33$1
    Rscript 5-fit_model_long.R initForFit_rf_0.33$1 &
    Rscript 5-fit_model_long.R initForFit_rf_0.33$1 _less &
    wait
}
export -f run
parallel -j 3 run ::: 1 2 3 4 5 6 7 8 9


