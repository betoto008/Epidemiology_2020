#!/bin/sh

#  simulation_Kitas.sh
#  
#
#  Created by Roberto Moran Tovar on 02.03.21.
#
n_testing_days=(0 1 2 3)
testing_days=(0 2 13 024)
betas=(0.16 0.18 0.20 0.25 0.33)
p_in_s=(0.1 0.05 0.01 0.005)
taus=(1 2 3)
for beta in "${betas[@]}"
    do
    for p_in in "${p_in_s[@]}"
        do
        for tau in "${taus[@]}"
            do
            for i in "${!n_testing_days[@]}"
                do
                printf "_____________________ \n"
                printf "Parameters: ${n_testing_days[i]} ${testing_days[i]} $beta $p_in $tau\n"
                ./a.out ${n_testing_days[i]} ${testing_days[i]} $tau $beta $p_in
                done
            done
        done
    done
