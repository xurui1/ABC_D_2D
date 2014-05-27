//
//  main.cpp
//  SCFT: Crank-Nicolson in 1D linear system
//
//  Created by Rui Xu on 2014-05-22.
//  Copyright (c) 2014 McMaster. All rights reserved.
//
// Basic idea:
/* 1. initialize w fields with random values
   2. solve diffusion equation for the two end integrated propagators (q, q+)
   3. determine monomer densities (Qalpha, Phialpha)
   4. Find new omega field, determine free energy
   5. Repeat steps 2-4 until free energy converges to a single value.
*/




#include <iostream>
#include "NumMeth.h"
#include "global.h"


using namespace std;




int main() {
    
    
}
