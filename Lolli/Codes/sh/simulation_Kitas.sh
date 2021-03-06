#!/bin/sh

#  simulation_Kitas.sh
#  
#
#  Created by Roberto Moran Tovar on 02.03.21.
#
n_testing_days=(0 1 2 3)
testing_days=(0 2 13 024)
#betas=(0.18, 0.36,  0.54)
betas=(0.36, 0.54, 0.72, 0.90)
p_in_s=(0.001, 0.01, 0.1)
taus=(4, 5)
c++ Kitas_ensemble.cpp -lgsl -o Kitas_ensemble.x
chmod +x Kitas_ensemble.x
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
                ./Kitas_ensemble.x ${n_testing_days[i]} ${testing_days[i]} $tau $beta $p_in
                done
            done
        done
    done
python analysis.py
