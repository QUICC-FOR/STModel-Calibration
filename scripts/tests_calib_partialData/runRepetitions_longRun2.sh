#!/bin/bash
run(){
    Rscript 5-fit_model.R initForFit_rf_0.33$1 3 1 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 3 5 &
    wait
}
export -f run
parallel -j  run ::: 5 6 7 8 9


