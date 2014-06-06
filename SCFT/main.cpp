// Basic idea:
/* 1. initialize w fields with random values
   2. solve diffusion equation for the two end integrated propagators (q, q+)
   3. determine monomer densities (Qalpha, Phialpha)
   4. Find new omega field, determine free energy
   5. Repeat steps 2-4 until free energy converges to a single value.
*/

#include <iostream>                 //Standard functions
#include <fstream>                  //For file input/output
#include "smemory.h"                //Ash's memory allocation class
#include "global.h"                 //definition of global variables, functions
#include <cstring>                  //For manipulation of strings and arrays
#include <cstdlib>                  //General utilities library (has random num)
#include "xmatrix.h"                //The interaction matrix
#include "phi.h"                    //Concentrations of chains
#include "etas.h"                   //incompressibility and pinning conditions
#include "Q_partition.h"            //Chain partition function
#include "rand_field.h"             //Generate random field or input from file
#include "profile.h"                //output profile of system
#include "cleanme.h"                //clean old data files
#include "new_fields.h"             //generate new potential field based on concentration
#include "fE_homo.h"                //Homogeneous free energy
#include "FreeEnergy.h"             //Free energy
#include "read.h"                   //extract parameters from file
#include "TDMA.h"                   //Tridiagonal Matrix Algorithm
#include "modA.h"                   //Solve Diffusion Equation for A block
#include "modB.h"                   //Solve Diffusion Equation for B block
#include "modC.h"                   //Solve Diffusion Equation for C block
#include "modD.h"                   //Solve Diffusion Equation for D chain
#include "modphi.h"                 //concentrations
#include <string>                   //manipulation of strings (not sure if needed)
#include <sstream>                  //stringstream (not sure if needed)


using namespace std;

int secant();                                   //function prototype, defined below

int main() {
    
    int once;
    int msg1,msg2;
    int s,s2,initial,i,j,jj;                               //counting indices
    int muD_up,muD_down;                                   //For stepping up and down.
    double ptot,phiA,phiB,phiC,phiD;                       //concentrations of chains
    double fracA, fracB, fracC;                            //chain fractions
    double fEnow,fEold;                                    // Used for calculating the dfE
    double D_r,D_z;                                        //Box size in the r and z direction is Rg^2
    double Area = 0.0,Tip_R;                               //Area and perimeter?? of system
    double f_int_0, ABC_0;
    
    eta =create_2d_double_array(Nr,Nz, "eta");             //incompressibility condition
    dpp =create_2d_double_array(Nr,Nz, "dpp");             //not sure what dpp is used for (update w fields?)
    eta2 =create_2d_double_array(Nr,Nz, "eta2");           //pinning condition
    
    wA=create_2d_double_array(Nr,Nz, "wA");                //interaction potential of chain A
    dwA=create_2d_double_array(Nr,Nz, "dwA");              //differential of wA
    pA=create_2d_double_array(Nr,Nz, "pA");                //probability amplitude of chain A
    
    wB=create_2d_double_array(Nr,Nz, "wB");                //interaction potential of chain B
    dwB=create_2d_double_array(Nr,Nz, "dwB");              //differential of wB
    pB=create_2d_double_array(Nr,Nz, "pB");                //probability amplitude of chain B
    
    wC=create_2d_double_array(Nr,Nz, "wC");                //interaction potential of chain C
    dwC=create_2d_double_array(Nr,Nz, "dwC");              //differential of wC
    pC=create_2d_double_array(Nr,Nz, "pC");                //probability amplitude of chain C
    
    wD=create_2d_double_array(Nr,Nz, "wD");                //interaction potential of chain D
    dwD=create_2d_double_array(Nr,Nz, "dwD");              //differential of wD
    pD=create_2d_double_array(Nr,Nz, "pD");                //probability amplitude of chain D
    
    
    f_int_0=0.0;
    ABC_0=0.0;
    
    initial=1;
    
    read();
  
    
    delt=1.0/(NA+NB+NC);
    delr=D_r/double(Nr-1);
    delz=D_z/double(Nz-1);
    
    /***************************************Initial settings********************************/
    
    
    
    r_0=10.0;                            // Distance from the center of cylinder
    initial=1;                           // 0=Random 1=Read from file 2=make your own condition
                                         //if intial=1, then can choose config. 1=on, 0=off
    
    ten_find=1;                          // If ten_find=1 then it will find the mu for the tensionless bilayer
    
    muD_up=0;                            // increase mu of homopolymer D
    muD_down=0;                          // decrease mu of homopolymer D
    
    bilayer=1;                          //If ten_find is on, turn off bilayer
    once=1;
    disk=0;
    
    
    
    /******************************** Clean the data files*****************************/
    
    if (disk==1) {
        ofstream myfile;
        myfile.open("./results/fE_disk.dat", std::ofstream::trunc);
        myfile.close();
    }
    if (bilayer==1) {
        ofstream myfile;
        myfile.open("./results/fE_bilayer.dat", std::ofstream::trunc);
        myfile.close();
        }

    
    XMatrix();                                //define the interaction matrix
    
    
    /*********************Define the pinning condition*************************************/
    
    if (disk==1){
        Ntip=20;
        Mtip=1;
        Tip_R=(Ntip-1)*delr;
    }
    else {
        Tip_R=0;
    }
    
    /*********************Set Kappa and fA**************************************************/
    
    kappaD=ND/(NA+NB+NC);
    fracA=NA/(NA+NB+NC);
    fracB=NB/(NA+NB+NC);
    fracC=NC/(NA+NB+NC);
    
    for (;;)
        if (ten_find==1){
            secant();
            
            rand_field(initial);
            iter=iter+1;
            **eta2=0;
            **eta=0;
            
            fE_homo();
            
            Vol=2.0*pi*(0.5*(Nz-1)*delz*pow(((Nr-1)*delr+r_0),2)-pow((r_0),2));
            if (bilayer==1){Area=pi*(pow(((Nr-1)*delr+r_0),2)-pow((r_0),2));}
            if (disk==1) {Area=pi*(pow((Tip_R+r_0),2)-pow((r_0),2));}
            
            cleanme();
            s2=0;
            fE_old=0.0;
            
            for (s=1;s<=100000;s++){
                qA_forward();
                qB_forward();
                qC_forward();
                qD_forward();
                
                qdagA_forward();
                qdagB_forward();
                qdagC_forward();
                
                Q_partition();
                phi();
                pressure();
                if (disk==1){pressure2();}
                FreeEnergy();
                totalphi(ptot,phiA,phiB,phiC,phiD);
                new_fields();
                double phistar= phiA+phiB+phiC;
                
                
                /********************print during run***********************/
                cout<< fE-fE_hom<< phistar << phiD<<endl;
                //write(4,*) real(s),fE,dfffE
                //write temporary data for plotting
                
                if (Conv_p<1.0e-4 and Conv_w<1.0e-4 and dfffE<1.0e-4)
                    break;
            }
            //OP=((phiA+phiB)-(phiE+phiD))/((phiA+phiB)+(phiE+phiD));
            cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
            if (bilayer==1){
                cout<<"Shape chosen is: Bilayer"<< endl;
                //cout<<"order_parameter"<<OP<< endl;
                cout<<"muABC="<<muABC<<"   "<<"muD="<< muD<< endl;
                cout<<"fE="<<((fE-fE_hom)*Vol)/Area;
            }
            else if (disk==1){
                cout<<"Shape chosen is: Disk"<< endl;
                //cout<<"order_parameter"<<OP<< endl;
                cout<<"muABC="<<muABC<<"   "<<"muD="<< muD<< endl;
                cout<<"fE="<<((fE-fE_hom)*Vol)/Area;
                }
            cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
            
            if (disk==1){
                fstream myfile;
                myfile.open("./results/fE_disk.dat");
                myfile<<ND<<((fE-fE_hom)*Vol)/Area<<((phiA+phiB+phiD)-phiABC_hom)*Vol/Area<< (r_0+Tip_R)<<
                (phiA+phiB+phiC)<<(phiD)<<Area<<muABC<<muD<< endl;
                myfile.close();
            }
            if (bilayer==1){
                fstream myfile;
                myfile.open("./results/fE_bilayer.dat");
                myfile<<ND<<((fE-fE_hom)*Vol)/Area<<((phiA+phiB+phiD)-phiABC_hom)*Vol/Area<< (r_0+Tip_R)<<
                (phiA+phiB+phiC)<<(phiD)<<Area<<muABC<<muD<< endl;
                myfile.close();
            }
            profile(2);
            
            if ((disk==1) and (bilayer==1)) {break;}
            
            if (disk==1){
                if (muD_up==1){
                    muD=muD+0.1;
                    if (OP<(-0.99)){break;}}
                else if (muD_down==1){
                    muD=muD-0.1;
                    if (OP>0.99){break;}}
                }
            if ((disk==0) and (bilayer==0)){break;}
           
            
            
        }
    
    
    
    
    
}

int secant(){
    
    
    int s,s2,ii,msg,msg1;
    double mu1,mu2,muD;
    int initial;
     /*int i,j;
    double mu3
     double fE1,fE2,fE3;*/
    
    
    mu1=muD;
    mu2=muD+0.01;
    
    bilayer=1;
    msg1=1;
    
    rand_field(initial);
    ii=1;
    msg=0;
    
    for (;;){
        **eta2=0.0;
        **eta=0.0;
        fE_homo( );
        Vol=2.0*pi*(0.5*(Nz-1)*delz*(pow(((Nr-1)*delr+r_0),2)-pow((r_0),2)));
        
        s2=0;
        fE_old=0.0;
        cleanme();
        
        for (s=1;s<=10000;s++){
            qA_forward();
            qB_forward();
            qC_forward();
            qD_forward();
            qdagA_forward();
            qdagB_forward();
            qdagC_forward();
            Q_partition( );
            phi( );
            pressure( );
            FreeEnergy( );
            new_fields( );
            
            
            cout<< Conv_w<< fE-fE_hom<< muD << ii << endl;
            
        }
        
    }
    return secant();
}
