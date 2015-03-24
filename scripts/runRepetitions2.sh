#!/bin/bash
run(){
    Rscript 5-fit_model.R initForFit_multinom_0.3$1 &
    Rscript 5-fit_model.R initForFit_rf_0.3$1 &
}
export -f run
parallel -j 8 run ::: 1 2 3 4 5 6 7 8 9


