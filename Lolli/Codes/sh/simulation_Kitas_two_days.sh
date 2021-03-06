#!/bin/sh

#  simulation_Kitas_two_days.sh
#  
#
#  Created by Roberto Moran Tovar on 03.03.21.
#  
n_testing_days=(2)
testing_days=(02 13 24)
betas=(0.18 0.36)
p_in_s=(0.011 0.0011)
taus=(0, 1, 2, 3, 4)
c++ Kitas_ensemble.cpp -lgsl -o Kitas_ensemble.x
chmod +x Kitas_ensemble.x
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
                ./Kitas_ensemble.x 2 $d $tau $beta $p_in
                done
            done
        done
    done
