#!/bin/sh

#  simulation_Kitas_one_day.sh
#  
#
#  Created by Roberto Moran Tovar on 03.03.21.
#  
n_testing_days=(1)
testing_days=(0 1 2 3 4)
betas=(0.18 0.36)
p_in_s=(0.11 0.011 0.0011)
taus=(0 1 2 3)
for beta in "${betas[@]}"
    do
    for p_in in "${p_in_s[@]}"
        do
        for tau in "${taus[@]}"
            do
            for d in "${testing_days[@]}"
                do
                printf "_____________________ \n"
                printf "Parameters: 1 $d $beta $p_in $tau\n"
                ./Kitas_ensemble.x 1 $d $tau $beta $p_in
                done
            done
        done
    done
