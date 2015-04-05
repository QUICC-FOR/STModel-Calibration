#!/bin/bash
run(){
#Rscript 4-initForFit.R 0.33$1
    Rscript 5-fit_model.R initForFit_cst_0.33$1 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 _less &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 5y &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 _less5y &
    wait
}
export -f run
parallel -j 3 run ::: 1 2 3 4 5 6 7 8 9


