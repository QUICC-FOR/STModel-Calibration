#!/bin/bash
run(){
    Rscript 5-fit_model.R initForFit_rf_0.33$1 0 1 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 0 5 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 1 1 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 1 5 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 2 1 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 2 5 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 3 1 &
    Rscript 5-fit_model.R initForFit_rf_0.33$1 3 5 &
    wait
}
export -f run
parallel -j 2 run ::: 1 2 3 4 5 6 7 8 9


