#!/bin/bash

for VAR in amp_damp depol rand_x
do
  for CITYNO in 3 4
  do
    for GAMMA in 0 1 2 3 4
    do
      python energy_diff_generator.py $VAR $CITYNO 40 1 100 rand $GAMMA &
    done
  done
done
sleep 1 # to make a break before the second part

python optimization_qaoa.py 3 40
python optimization_qaoa.py 4 40
