// Basic idea:
/* 1. initialize w fields with random values
   2. solve diffusion equation for the two end integrated propagators (q, q+)
   3. determine monomer densities (Qalpha, Phialpha)
   4. Find new omega field, determine free energy
   5. Repeat steps 2-4 until free energy converges to a single value.
*/

#include <iostream>                 //Standard functions
#include <fstream>                  //For file input/output
#include <string>                   //manipulation of strings (not sure if needed)
#include <sstream>                  //stringstream (not sure if needed)
#include <cstring>                  //For manipulation of strings and arrays
#include <cstdlib>                  //General utilities library (has random num)
#include "global.h"                 //definition of global variables, functions
#include "smemory.h"                //Ash's memory allocation class
#include "xmatrix.h"                //The interaction matrix
#include "Q_partition.h"            //Chain partition function
#include "fields.h"                 //initial and update fields (omega, eta1, eta2)
#include "energy.h"                 //Free energy
#include "read_write.h"             //file input/ouput parameters from file
#include "TDMA.h"                   //Tridiagonal Matrix Algorithm
#include "modA.h"                   //Solve Diffusion Equation for A block
#include "modB.h"                   //Solve Diffusion Equation for B block
#include "modC.h"                   //Solve Diffusion Equation for C block
#include "modD.h"                   //Solve Diffusion Equation for D chain
#include "modphi.h"                 //concentrations
#include "Create_Destroy.h"         //define and destroy arrays
#include "modsecant.h"              //The secant method

using namespace std;

int main() {
    
    int once;
    int s,s2,initial,i,j;                                  //counting indices
    int muD_up,muD_down;                                   //For stepping up and down potential E.
    double Area = 0.0,Tip_R;                               //Area and perimeter?? of system
    double f_int_0, ABC_0;                                 //not sure of purpose
    double nan(const char* tagp);
    
    /*********************Set Kappa and fA**************************************************/
    kappaD=(double)ND/((double)(NA+NB+NC));
    fracA=((double)NA/((double)(NA+NB+NC)));
    fracB=(double)NB/(double)(NA+NB+NC);
    fracC=(double)NC/(double)(NA+NB+NC);
    
    create();
    
    f_int_0=0.0;
    ABC_0=0.0;
    
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
    clean_data();                        // Clean the result files
    
    /*********************Define the pinning condition*************************************/
    
    if (disk==1){
        Ntip=20;Mtip=1;
        Tip_R=(Ntip-1)*delr;}
    else {Tip_R=0;}
    
    /***********************************for LOOP*******************************************/
    
    for (iter=0;;iter++){
        if (ten_find==1)
            {secant();}
        
        rand_field(initial);                    //Found in rand_field
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0;j<=int(Nz-1);j++){
                eta2[i][j]=0.0;
                eta[i][j]=0.0;}}

        fE_homo();                               //Found in Energy
            
        Vol=2.0*pi*(0.5*(Nz-1)*delz*pow(((Nr-1)*delr+r_0),2)-pow((r_0),2));
        if (bilayer==1){Area=pi*(pow(((Nr-1)*delr+r_0),2)-pow((r_0),2));}
        if (disk==1) {Area=pi*(pow((Tip_R+r_0),2)-pow((r_0),2));}
        
        cleanme();                              //Found in read_write
    
        s2=0;fE_old=0.0;
        
        for (s=1;s<=100000;s++){
            qA_forward();qB_forward();qC_forward();qD_forward();  //solve forward propagators
            qdagA_forward();  qdagB_forward();  qdagC_forward();  //solve complementary propagators

            Q_partition();                      //Solve chain partition function
            phi();                              //Determine new monomer concentrations
            pressure();                         //Determine/apply incompressibility condition
            if (disk==1){pressure2();}          //Determine/apply pinning condition
            FreeEnergy();                       //Solve for free energy
            totalphi();                         //Determine total monomer concentrations
            new_fields();                       //Set new interaction fields
            
            
            cout<< "fE: "<<(fE-fE_hom)*Vol/Area<< " Phi-star: "<<phiA+phiB+phiC<<" phi-D: "<< phiD<<" "<<iter<<endl;
            
            if (abs((fE-fE_hom)*Vol/Area)>=(1.0e5) or isnan(fE)==true){cout<<"Free energy out of bounds";return 0;}

            if (s>s2){
                profile(1);
                write_data();
                s2=s2+10;}
                
            if (Conv_p<1.0e-4 and Conv_w<3.0e-4 and dfffE<1.0e-4){break;}
            
        }
            
        //OP=((phiA+phiB+phiC)-(phiD))/((phiA+phiB+phiC)+(phiD)); //irrelevant, only one type of star copolymer
        show_data(Area);
        save_data(Area,phiA,phiB,phiC,phiD,Tip_R);
        profile(2);
            
        if (bilayer==1 and once==1){break;}
        if ((disk==1) and (bilayer==1)) {break;}
        if ((disk==0) and (bilayer==0)){break;}
        if (disk==1){
                if (muD_up==1){
                    muD=muD+0.1;
                    if (OP<(-0.99)){break;}}
                else if (muD_down==1){
                    muD=muD-0.1;
                    if (OP>0.99){break;}}
                }
        if (bilayer==1){
            if (muD<-5.0){break;}
            muD=muD-0.1;}
 
        
    /********************************end of iteration loop*********************************/
    }

    Destroy();                      //Destroy all allocated arrays
    
}




