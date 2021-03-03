#!/bin/sh

#  simulation_Kitas_two_days.sh
#  
#
#  Created by Roberto Moran Tovar on 03.03.21.
#  
n_testing_days=(2)
testing_days=(02 13 24)
betas=(0.16 0.18 0.20 0.25 0.33)
p_in_s=(0.1 0.05 0.01 0.005)
taus=(1 2 3)
for beta in "${betas[@]}"
    do
    for p_in in "${p_in_s[@]}"
        do
        for tau in "${taus[@]}"
            do
            for d in "${testing_days[@]}"
                do
                printf "_____________________ \n"
                printf "Parameters: 2 $d $beta $p_in $tau\n"
                ./a.out 2 $d $tau $beta $p_in
                done
            done
        done
    done
