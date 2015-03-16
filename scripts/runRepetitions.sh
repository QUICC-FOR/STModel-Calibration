#!/bin/bash
run(){
    Rscript 4-init_params.R 0.003$1
    wait
    Rscript 5-fit_model.R initForFit_multinom_0.003$1 &
    Rscript 5-fit_model.R initForFit_rf_0.003$1 &
}
export -f run
parallel run ::: 1 2 3


