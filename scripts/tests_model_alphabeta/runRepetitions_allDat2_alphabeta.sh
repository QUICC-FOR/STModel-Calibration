#!/bin/bash
run(){
    Rscript 5-fit_model_allData_alphabeta_rep.R 5 $1
}
export -f run
parallel -j 12 run ::: 1 2 3 4 5 6 7 8 9 10 11 12


