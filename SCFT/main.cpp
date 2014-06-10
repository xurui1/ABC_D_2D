// Basic idea:
/* 1. initialize w fields with random values
   2. solve diffusion equation for the two end integrated propagators (q, q+)
   3. determine monomer densities (Qalpha, Phialpha)
   4. Find new omega field, determine free energy
   5. Repeat steps 2-4 until free energy converges to a single value.
*/

#include <iostream>                 //Standard functions
#include <fstream>                  //For file input/output
#include "global.h"                 //definition of global variables, functions
#include "smemory.h"                //Ash's memory allocation class
#include <cstring>                  //For manipulation of strings and arrays
#include <cstdlib>                  //General utilities library (has random num)
#include "xmatrix.h"                //The interaction matrix
#include "phi.h"                    //Concentrations of chains
#include "etas.h"                   //incompressibility and pinning conditions
#include "Q_partition.h"            //Chain partition function
#include "rand_field.h"             //Generate random field or input from file
#include "profile.h"                //output profile of system
#include "new_fields.h"             //generate new potential field based on concentration
#include "fE_homo.h"                //Homogeneous free energy
#include "FreeEnergy.h"             //Free energy
#include "read_write.h"                   //extract parameters from file
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
    int s,s2,initial,i,j;                               //counting indices
    int muD_up,muD_down;                                   //For stepping up and down.
    double ptot,phiA,phiB,phiC,phiD;                       //concentrations of chains
    
                                                           //Box size in the r and z direction is Rg^2
    double Area = 0.0,Tip_R;                               //Area and perimeter?? of system
    double f_int_0, ABC_0;                                 //not sure of purpose
    double phistar;                                        //concentration of star copolymer
    
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
    
    /*********************Set Kappa and fA**************************************************/
    kappaD=(double)ND/((double)(NA+NB+NC));
    fracA=((double)NA/((double)(NA+NB+NC)));
    fracB=(double)NB/(double)(NA+NB+NC);
    fracC=(double)NC/(double)(NA+NB+NC);
    
    f_int_0=0.0;
    ABC_0=0.0;
    
    initial=1;
    
    //read(); not working right now for some reason, cannot figure out why
  
    
    delt=1.0/(NA+NB+NC);
    delr=(double)D_r/(double(Nr-1));
    delz=(double)D_z/(double(Nz-1));
    
    /***************************************Initial settings********************************/
    
    r_0=10.0;                            // Distance from the center of cylinder
    initial=0;                           // 0=Random 1=Read from file 2=make your own condition
                                         //if intial=1, then can choose config. 1=on, 0=off
    
    ten_find=0;                          // If ten_find=1 then it will find the mu for the tensionless bilayer
    muD_up=0;                            // increase mu of homopolymer D
    muD_down=0;                          // decrease mu of homopolymer D
    
    bilayer=1;                          //If ten_find is on, turn off bilayer
    once=1;
    disk=0;
    
    XMatrix();                           //call the interaction matrix
    clean_data();                        // Clean the data files
    

    /*********************Define the pinning condition*************************************/
    
    if (disk==1){
        Ntip=20;
        Mtip=1;
        Tip_R=(Ntip-1)*delr;
    }
    else {
        Tip_R=0;
    }
    
   
    
    for (iter=0;;iter++){
        if (ten_find==1)
            {secant();}
        rand_field(initial);
        
        for (i=0;i<=Nr-1;i++){
            for (j=0;j<=Nz-1;j++){
                eta2[i][j]=0.0;
                eta[i][j]=0.0;
            }
        }
            
       //fE_homo();
            
        Vol=2.0*pi*(0.5*(Nz-1)*delz*pow(((Nr-1)*delr+r_0),2)-pow((r_0),2));
        if (bilayer==1){Area=pi*(pow(((Nr-1)*delr+r_0),2)-pow((r_0),2));}
        if (disk==1) {Area=pi*(pow((Tip_R+r_0),2)-pow((r_0),2));}
            
        cleanme();
    
        s2=0;
        fE_old=0.0;
            
        for (s=1;s<=100000;s++){
            qA_forward();                       //solve diffusion equation for A block
            qB_forward();                       //solve diffusion equation for B block
            qC_forward();                       //solve diffusion equation for C block
            qD_forward();                       //solve diffusion equation for D block
                
            qdagA_forward();                    //solve diffusion equation for complementary A propagator
            qdagB_forward();                    //solve diffusion equation for complementary B propagator
            qdagC_forward();                    //solve diffusion equation for complementary C propagator
                
            Q_partition();                      //Solve chain partition function
            phi();                              //Determine new monomer concentrations
            pressure();                         //Determine/apply incompressibility condition
            if (disk==1){pressure2();}          //Determine/apply pinning condition
            FreeEnergy();                       //Solve for free energy
            totalphi(ptot,phiA,phiB,phiC,phiD); //Determine total monomer concentrations
            new_fields();   //Set new interaction fields
            
            phistar= phiA+phiB+phiC;            //Concentration of star copolymer
                
                
            
            cout<< fE-fE_hom<< phistar<< phiD<<endl;
            
            //write(4,*) real(s),fE,dfffE
            
            //write temporary data for plotting
                
            if (Conv_p<1.0e-4 and Conv_w<1.0e-4 and dfffE<1.0e-4){
                    break;
            }
            
            OP=((phiA+phiB+phiC)-(phiD))/((phiA+phiB+phiC)+(phiD));
            
            show_data(Area);
        }
        
        save_data(Area,phiA,phiB,phiC,phiD,Tip_R);
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
    
    destroy_2d_double_array(eta);
    destroy_2d_double_array(dpp);
    destroy_2d_double_array(eta2);
    
    destroy_2d_double_array(wA);
    destroy_2d_double_array(dwA);
    destroy_2d_double_array(pA);
    
    destroy_2d_double_array(wB);
    destroy_2d_double_array(dwB);
    destroy_2d_double_array(pB);

    destroy_2d_double_array(wC);
    destroy_2d_double_array(dwC);
    destroy_2d_double_array(pC);

    destroy_2d_double_array(wD);
    destroy_2d_double_array(dwD);
    destroy_2d_double_array(pD);
    
    destroy_3d_double_array(qA);
    destroy_3d_double_array(qB);
    destroy_3d_double_array(qC);
    destroy_3d_double_array(qD);
    
    destroy_3d_double_array(qdagA);
    destroy_3d_double_array(qdagB);
    destroy_3d_double_array(qdagC);
}
/******************************************************************************************/
/**************************************The Secant method***********************************/

int secant(){
    
    
    int s,s2,ii,msg,msg1,i,j;
    double mu1,mu2,muD;
    int initial;
    
    
    
    mu1=muD;
    mu2=muD+0.01;
    
    bilayer=1;
    msg1=1;
    
    rand_field(initial);
    ii=1;
    msg=0;
    
    for (;;){
        for (i=0;i<=Nr;i++){
            for (j=0;j<=Nz;j++){
                eta2[i][j]=0;
                eta[i][j]=0;
            }
        }
        
        
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
