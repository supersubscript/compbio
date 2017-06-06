#!/bin/sh

## This script runs the simulations in R in parallel

run () {
  START=$1
  END=$2
  for ii in $(seq $START $END); do
    ( 
      RETURN_VAL=$(nice -n 5 Rscript --vanilla /local/data/public/hpa22/assignments/syba1/gillespie.R $ii)
    ) &
  done
}

run 1 20
wait
run 21 40
wait
run 41 60
wait
run 61 80
wait
run 81 100

exit
