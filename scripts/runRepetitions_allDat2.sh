#!/bin/bash
run(){
    Rscript 5-fit_model_allData_rep.R 5 $1
}
export -f run
parallel -j 2 run ::: 1 2 3 4

