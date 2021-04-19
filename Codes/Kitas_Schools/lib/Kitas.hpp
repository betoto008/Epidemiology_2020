//
//  Kitas.hpp
//  
//
//  Created by Roberto Moran Tovar on 27.02.21.
//

#ifndef Kitas_hpp
#define Kitas_hpp

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <time.h>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "./random.cpp"

#endif /* Kitas_hpp */

void update_state_kids(int N, int *H, int *Inc, int *Inf, int *Rec, int *Det, std::vector<int> *kids){
    
    *H = std::count(kids->begin(), kids->end(), 0);
    *Inc = std::count(kids->begin(), kids->end(), 1);
    *Inf = std::count(kids->begin(), kids->end(), 2);
    *Rec = std::count(kids->begin(), kids->end(), 3);
    *Det = std::count(kids->begin(), kids->end(), 4);
    
    
}

double sensitivity(int infected_days, int t_inc, int tau){
    double a = 0.0;
    if(tau<4){
        if(infected_days > (t_inc-tau)){
            a = 1.0;
        }
    }
    else{
        a = (1.0/(1.0+pow(log10(exp(3.0*infected_days))/(2.5),-4)));
    }
    
    return a;
}
