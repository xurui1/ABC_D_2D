

#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <stdlib.h>    //Include standard fucntion libraries
#include <math.h>      //Use the math function libraries
#include <time.h>      //Call system time libraries to define the integer seed for random numbers
#include <cmath>
#include <vector>

using namespace std;

        int Nr,Nz,Ns;
        int D_r=5.0, D_z=11.0;
        int NA=50.0,NB=50.0,NC=50.0,ND=50.0;
        const double pi = 3.14195265;
        double delr,delz,delt,sig=0.05,sig2=0.05;
        double Q_ABC,Q_D,fE,fE_old;
        double xAB=30.0,xAC=30.0,xAD=0.0,xBC=0.0,xBD=30.0,xCD=30.0;
        double Vol;
        double Conv_w,Conv_p,dfffE,fE_hom;
        double kappaD,fracA,fracB,fracC,muABC=0.0,muD=-5.0,r_0;
        double XM[4][4];
        float phiABC_hom,phiD_hom;


    
    int bilayer,disk,ten_find;
    int Mtip,Ntip,iter;

/********************************Parameters*********************************/

double OP; // Order parameter
    
int n1, n2, n3;
    
double **eta;
double **dpp;
double **eta2;

/********************************A-Block*********************************/

double *bAr;
double *bAz;
double **wA;
double **dwA;
double **pA;
double ***qA;
double **qA_0;
double ***qdagA;
double **qdagA_0;
double *DiagAr;
double *DiagAUr;
double *DiagALr;
double *DiagArdr;
double *DiagAUrdr;
double *DiagALrdr;
double *DiagAz;
double *DiagAUz;
double *DiagALz;
double *DiagAzdz;
double *DiagAUzdz;
double *DiagALzdz;


/********************************B-Block*********************************/

double *bBr;
double *bBz;
double **wB;
double **dwB;
double **pB;
double ***qB;
double **qB_0;
double ***qdagB;
double **qdagB_0;
double *DiagBr;
double *DiagBUr;
double *DiagBLr;
double *DiagBrdr;
double *DiagBUrdr;
double *DiagBLrdr;
double *DiagBz;
double *DiagBUz;
double *DiagBLz;
double *DiagBzdz;
double *DiagBUzdz;
double *DiagBLzdz;
 


/********************************C-Block*********************************/

double *bCr;
double *bCz;
double **wC;
double **dwC;
double **pC;
double ***qC;
double **qC_0;
double ***qdagC;
double **qdagC_0;
double *DiagCr;
double *DiagCUr;
double *DiagCLr;
double *DiagCrdr;
double *DiagCUrdr;
double *DiagCLrdr;
double *DiagCz;
double *DiagCUz;
double *DiagCLz;
double *DiagCzdz;
double *DiagCUzdz;
double *DiagCLzdz;


/********************************D-chain***********************************/

double *bDr;
double *bDz;
double **wD;
double **dwD;
double **pD;
double ***qD;
double **qD_0;
double *DiagDr;
double *DiagDUr;
double *DiagDLr;
double *DiagDrdr;
double *DiagDUrdr;
double *DiagDLrdr;
double *DiagDz;
double *DiagDUz;
double *DiagDLz;
double *DiagDzdz;
double *DiagDUzdz;
double *DiagDLzdz;

