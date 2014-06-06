void XMatrix(){
#include "global.h"
    
 XM[0][0]=0.0;
 XM[1][0]=xAB;
 XM[2][0]=xAC;
 XM[3][0]=xAD;
 
 XM[0][1]=xAB;
 XM[1][1]=0.0;
 XM[2][1]=xBC;
 XM[3][1]=xBD;
 
 XM[0][2]=xAC;
 XM[1][2]=xBC;
 XM[2][2]=0.0;
 XM[3][2]=xCD;
 
 XM[0][3]=xAD;
 XM[1][3]=xBD;
 XM[2][3]=xCD;
 XM[3][3]=0.0;
 
 
 }