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
#include <time.h>
#include <algorithm>

#include "./random.cpp"

#endif /* Kitas_hpp */

void update_state_kids(int N, int *H, int *Inc, int *Inf, int *Rec, int *Det, std::vector<int> *kids){
    
    *H = std::count(kids->begin(), kids->end(), 0);
    *Inc = std::count(kids->begin(), kids->end(), 1);
    *Inf = std::count(kids->begin(), kids->end(), 2);
    *Rec = std::count(kids->begin(), kids->end(), 3);
    *Det = std::count(kids->begin(), kids->end(), 4);
    
    
}
