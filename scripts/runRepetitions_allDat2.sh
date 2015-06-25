#!/bin/bash
run(){
    Rscript 5-fit_model_allData_rep.R 5 $1
}
export -f run
parallel -j 6 run ::: 1 2 3 4 5 6

