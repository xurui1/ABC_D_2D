//
//  global.h
//  SCFT
//
//  Created by Rui Xu on 2014-05-26.
//  Copyright (c) 2014 McMaster. All rights reserved.
//

#ifndef SCFT_global_h
#define SCFT_global_h
#include <iostream>
#include <cstdlib>

using namespace std;



class global {
    public:
        int N,M;
        const double pi = 3.14195265;
        double delr,delz,delt,sig,sig2;
        double Q_AB,Q_C,Q_ED,fE,fE_old;
        double xAB,xAC,xBC,xAE,xAD,xBE,xBD,xCE,xCD,xED;
        double Vol;
        double Conv_w,Conv_p,dfffE,fE_hom;
        double kappaC,kappaED,fracA,fracE,muAB,muC,muED,r_0;
        double XM[5][5];
        float phiAB_hom,phiC_hom;
    
    int bilayer,disk,ten_find;
    int Mtip,Ntip,iter;
    
        double OP; // Order parameter
    
        int n1, n2, n3;
    
        double* eta = new double;
        double* dpp = new double;
        double* eta2 = new double;
    
    
   // need to figure out how to dynamically allocate memory
    
   /* eta = new double [2];
    dpp = new double [2];
    eta2 = new double [2];*/
    
    /***************************A-Block*******************************/
    int Na;
    double*  bAx = new double; double*  bAy = new double; double*  wA = new double; double*  pA = new double;
    double*  qA = new double; double*  qA_0 = new double; double*  qAD = new double; double*  qAD_0 = new double;
    double*  DiagAx = new double; double*  DiagAUx = new double; double*  DiagALx = new double;
    double*  DiagAdx = new double; double*  DiagAUdx = new double; double*  DiagALdx = new double;
    double*  DiagAy = new double; double*  DiagAUy = new double; double*  DiagALy = new double;
    double*  DiagAdy = new double; double*  DiagAUdy = new double; double*  DiagALdy = new double;
   
    // repeat for B-Block, C-Block
    
    /***************************D-Chain*******************************/
    
};
#endif
