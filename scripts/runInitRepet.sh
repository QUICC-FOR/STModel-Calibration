#!/bin/bash
run(){
    Rscript 4-initForFit.R 0.33$1
}
export -f run
parallel -j 9 run ::: 1 2 3 4 5 6 7 8 9


