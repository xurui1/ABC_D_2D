

#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <stdlib.h>    //Include standard function libraries
#include <math.h>      //Use the math function libraries
#include <cmath>
#include <ctgmath>
#include <ctime>
#include <vector>

using namespace std;

        int Nr=50,Nz=50,Ns=100;
        int D_r=5, D_z=11;                    //Box size in the r and z direction is Rg^2
        int NA=50.0,NB=50.0,NC=50.0,ND=100.0;
        const double pi = 3.14195265;
        double delr,delz,delt,sig=0.005,sig2=0.005;
        double Q_ABC,Q_D,fE,fE_old;
        double xAB=30.0,xAC=30.0,xAD=0.0,xBC=0.0,xBD=30.0,xCD=30.0;
        double Vol;
        double Conv_w,Conv_p,dfffE,fE_hom;
        double kappaD,fracA,fracB,fracC,muABC=0.0,muD=-4.596,r_0;
        double XM[4][4];
        float phiABC_hom,phiD_hom;
        double phistar;                           //concentration of star copolymer
        double ptot_in,phiA,phiB,phiC,phiD;


    int bilayer,disk,ten_find;
    int Mtip,Ntip,iter;

/********************************Parameters*********************************/

double OP; // Order parameter
    
int n1, n2, n3;
    
double **eta;
double **dpp;
double **eta2;

/********************************A-Block************************************/
double *bAr;
double *bAz;
double **wA;
double **dwA;
double **pA;
double ***qA;
double **qA_0;
double ***qdagA;
double **qdagA_0;
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
/********************************D-chain***********************************/
double *bDr;
double *bDz;
double **wD;
double **dwD;
double **pD;
double ***qD;
double **qD_0;


