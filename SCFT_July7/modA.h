

void A_Matrix_r(int ii, double *DiagAz, double *DiagALz, double *DiagAUz){


    int i;
    
    for (i=0;i<=int(Nr-1);i++){
        DiagAz[i]=0;}
    for (i=0;i<=int(Nr-2);i++){
        DiagAUz[i]=0;
        DiagALz[i]=0;}
   
    
    for (i=0; i<=int(Nr-1); i++){
        DiagAz[i]=1.0+(delt/(pow((double)delr,(int)2)))+((delt/2.0)*wA[i][ii]);}
    for (i=1; i<=int(Nr-2);i++){
        DiagAUz[i]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(i-1))+r_0)*4.0*delr));}
    for (i=0; i<=int(Nr-3);i++){
        DiagALz[i]=-(delt/(2.0*pow((double) delr,(int)2)))+(delt/(((delr*(i+1-1))+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagAUz[0]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(-1))+r_0)*4.0*delr));
    
    DiagALz[int(Nr-2)]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(int(Nr-1)))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(int(Nr-1)))+r_0)*4.0*delr));
    
}

void A_Matrix_z(int ii, double *DiagAr, double *DiagALr, double *DiagAUr){

    int j;
    
    for (j=0;j<=int(Nz-1);j++){
       DiagAr[j]=0;}
    for (j=0;j<(int)(Nz-1);j++){
        DiagAUr[j]=0;
        DiagALr[j]=0;}
   
    
    for (j=0; j<=int(Nz-1); j++){
        DiagAr[j]=1.0+delt/(pow((double)delz,(int)2));}
    for (j=1; j<=int(Nz-2); j++){
        DiagAUr[j]=-delt/(2.0*pow((double)delz,(int)2));}
    for (j=0; j<=(Nz-3); j++){
        DiagALr[j]=-delt/(2.0*pow((double)delz,(int)2));}
    
    DiagAUr[0]=2.0*DiagAUr[0];
    DiagALr[int(Nz-2)]=2.0*DiagALr[int(Nz-2)];
    
}


/**************************Finally time to define some propagators**********************/
void qA_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
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

    bAr=create_1d_double_array(Nz, "bAr");
    bAz=create_1d_double_array(Nr, "bAz");
    DiagAr=create_1d_double_array(Nz, "DiagAr");
    DiagAUr=create_1d_double_array(((int)Nz-1), "DiagAUr");
    DiagALr=create_1d_double_array(((int)Nz-1), "DiagALr");
    DiagAz=create_1d_double_array(Nr, "DiagAz");
    DiagAUz=create_1d_double_array(((int)Nr-1), "DiagAUz");
    DiagALz=create_1d_double_array(((int)Nr-1), "DiagALz");
    DiagArdr=create_1d_double_array(Nz, "DiagArdr");
    DiagAUrdr=create_1d_double_array(((int)Nz-1), "DiagAUrdr");
    DiagALrdr=create_1d_double_array(((int)Nz-1), "DiagALrdr");
    DiagAzdz=create_1d_double_array(Nr, "DiagAzdz");
    DiagAUzdz=create_1d_double_array(((int)Nr-1), "DiagAUzdz");
    DiagALzdz=create_1d_double_array(((int)Nr-1), "DiagALzdz");
    
    for (i=0; i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            qA_0[i][j]=0;
            for (s=0;s<=int(Ns-1);s++){
                 qA[i][j][s]=0;}}}
   
    
//Initialize the qs
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=1;j<=int(Nz-1);j++){
            qA_0[i][j]=1.0;
            qA[i][j][0]=1.0;}}
        
    for (s=0;s<=int(NA-1);s++){
        
        for (j=0;j<=int(Nz-1);j++){
            bAr[j]=0;}
        
/********************************scan over z***********************************************/
        for (i=0;i<=int(Nr-1);i++){
            A_Matrix_z(i, DiagAr, DiagALr, DiagALr);
            if (i==0){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wA[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bAr[j]=gamma*qA_0[i][j]+betaU*qA_0[int(i+1)][j]+betaL*qA_0[int(i+1)][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wA[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bAr[j]=gamma*qA_0[i][j]+betaU*qA_0[int(i-1)][j]+betaL*qA_0[int(i-1)][j];}
            }
            else {
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wA[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bAr[j]=gamma*qA_0[i][j]+betaU*qA_0[int(i+1)][j]+betaL*qA_0[int(i-1)][j];}
            }
            for (j=0;j<=int(Nz-1);j++){
                DiagArdr[j]=DiagAr[j];}
            for (j=0;j<=int(Nz-2);j++){
                DiagAUrdr[j]=DiagAUr[j];
                DiagALrdr[j]=DiagALr[j];}
            
            TDMA(Nz,DiagALrdr,DiagArdr,DiagAUrdr,bAr);
            
            for (j=0; j<=int(Nz-1); j++){
                qA[i][j][s]=bAr[j];}
            
            }
        for (j=0;j<=int(Nz-1);j++){
            bAr[j]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qA_0[i][j]=qA[i][j][s];}
            }
        
        /***********************************scan over r***************************************************/

        for (i=0;i<=int(Nr-1);i++){
            bAz[i]=0;
        }
        for (j=0;j<=int(Nz-1);j++){
            A_Matrix_r(j,DiagAz, DiagALz,DiagAUz);
            if (j==0){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bAz[i]=gamma*qA_0[i][j]+beta*qA_0[i][int(j+1)]+beta*qA_0[i][int(j+1)];}
            }
            else if (j==int(Nz-1)){
                for (i=0;i<=int(Nz-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bAz[i]=gamma*qA_0[i][j]+beta*qA_0[i][int(j-1)]+beta*qA_0[i][int(j-1)];}
            }
            else {
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bAz[i]=gamma*qA_0[i][j]+beta*qA_0[i][int(j+1)]+beta*qA_0[i][int(j-1)];}
            }
            for (i=0;i<=int(Nr-1);i++){
                DiagAzdz[i]=DiagAz[i];}
            for (i=0;i<=int(Nr-2);i++){
                DiagAUzdz[i]=DiagAUz[i];
                DiagALzdz[i]=DiagALz[i];}
            
            TDMA(Nr,DiagALzdz,DiagAzdz,DiagAUzdz,bAz);
            
            for (i=0; i<=int(Nr-1); i++){
                qA[i][j][s]=bAz[i];}
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bAz[i]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qA_0[i][j]=qA[i][j][s];}
        }
    }
    
    
    destroy_1d_double_array(DiagALr);
    destroy_1d_double_array(DiagAUr);
    destroy_1d_double_array(DiagAr);
    destroy_1d_double_array(DiagALrdr);
    destroy_1d_double_array(DiagAUrdr);
    destroy_1d_double_array(DiagArdr);
    destroy_1d_double_array(DiagALz);
    destroy_1d_double_array(DiagAUz);
    destroy_1d_double_array(DiagAz);
    destroy_1d_double_array(DiagALzdz);
    destroy_1d_double_array(DiagAUzdz);
    destroy_1d_double_array(DiagAzdz);
    destroy_1d_double_array(bAr);
    destroy_1d_double_array(bAz);
}

/***********************************Define the complementary propagator***********************************/

  void qdagA_forward(){
        
        int s,i,j;
        double gamma,betaU,betaL,beta;
      
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
      
      bAr=create_1d_double_array(Nz, "bAr");
      bAz=create_1d_double_array(Nr, "bAz");
      DiagAr=create_1d_double_array(Nz, "DiagAr");
      DiagAUr=create_1d_double_array(((int)Nz-1), "DiagAUr");
      DiagALr=create_1d_double_array(((int)Nz-1), "DiagALr");
      DiagAz=create_1d_double_array(Nr, "DiagAz");
      DiagAUz=create_1d_double_array(((int)Nr-1), "DiagAUz");
      DiagALz=create_1d_double_array(((int)Nr-1), "DiagALz");
      DiagArdr=create_1d_double_array(Nz, "DiagArdr");
      DiagAUrdr=create_1d_double_array(((int)Nz-1), "DiagAUrdr");
      DiagALrdr=create_1d_double_array(((int)Nz-1), "DiagALrdr");
      DiagAzdz=create_1d_double_array(Nr, "DiagAzdz");
      DiagAUzdz=create_1d_double_array(((int)Nr-1), "DiagAUzdz");
      DiagALzdz=create_1d_double_array(((int)Nr-1), "DiagALzdz");
        
        for (i=0; i<=int(Nr-1);i++){
            for (j=0;j<=int(Nz-1);j++){
                qdagA_0[i][j]=0;
                for (s=0;s<=int(Ns-1);s++){
                    qdagA[i][j][s]=0;}}}
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0;j<=int(Nz-1);j++){
                qdagA[i][j][0]=qB[i][j][int(NB-1)]*qC[i][j][int(NC-1)];
                qdagA_0[i][j]=qB[i][j][int(NB-1)]*qC[i][j][int(NC-1)];
            }
        }
        
        for (s=0;s<=int(NA-1);s++){
            
            for (j=0;j<=int(Nz-1);j++){
                bAr[j]=0;}
            
            /********************************scan over z***********************************************/
            for (i=0;i<=int(Nr-1);i++){
                A_Matrix_z(i,DiagAr, DiagALr, DiagALr);
                if (i==0){
                    for (j=0;j<=int(Nz-1);j++){
                        gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wA[i][j]);
                        betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        bAr[j]=gamma*qdagA_0[i][j]+betaU*qdagA_0[int(i+1)][j]+betaL*qdagA_0[int(i+1)][j];}
                }
                else if (i==int(Nr-1)){
                    for (j=0;j<=int(Nz-1);j++){
                        gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wA[i][j]);
                        betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        bAr[j]=gamma*qdagA_0[i][j]+betaU*qdagA_0[int(i-1)][j]+betaL*qdagA_0[int(i-1)][j];}
                }
                
                else {
                    for (j=0;j<=int(Nz-1);j++){
                        gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wA[i][j]);
                        betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        bAr[j]=gamma*qdagA_0[i][j]+betaU*qdagA_0[int(i+1)][j]+betaL*qdagA_0[int(i-1)][j];}}
                
                for (j=0;j<=int(Nz-1);j++){
                    DiagArdr[j]=DiagAr[j];}
                for (j=0;j<=int(Nz-2);j++){
                    DiagAUrdr[j]=DiagAUr[j];
                    DiagALrdr[j]=DiagALr[j];}
                
                TDMA(Nz,DiagALrdr,DiagArdr,DiagAUrdr,bAr);
                
                for (j=0; j<=int(Nz-1); j++){
                    qdagA[i][j][s]=bAr[j];}
            }
            
            for (j=0;j<=int(Nz-1);j++){
                bAr[j]=0;}
            
            for (i=0;i<=int(Nr-1);i++){
                for (j=0; j<=int(Nz-1); j++){
                    qdagA_0[i][j]=qdagA[i][j][s];}}
            
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bAz[i]=0;
        }
        /***********************************scan over r***************************************************/
        for (j=0;j<=int(Nz-1);j++){
            A_Matrix_r(j,DiagAz, DiagALz, DiagALz);
            if (j==0){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bAz[i]=gamma*qdagA_0[i][j]+beta*qdagA_0[i][int(j+1)]+beta*qdagA_0[i][int(j+1)];}
            }
            else if (j==int(Nz-1)){
                for (i=0;i<=int(Nz-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bAz[i]=gamma*qdagA_0[i][j]+beta*qdagA_0[i][int(j-1)]+beta*qdagA_0[i][int(j-1)];}
            }
            else {
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bAz[i]=gamma*qdagA_0[i][j]+beta*qdagA_0[i][int(j+1)]+beta*qdagA_0[i][int(j-1)];}
            }
            
            for (i=0;i<=int(Nr-1);i++){
                DiagAzdz[i]=DiagAz[i];}
            for (i=0;i<=int(Nr-2);i++){
                DiagAUzdz[i]=DiagAUz[i];
                DiagALzdz[i]=DiagALz[i];}
            
            TDMA(Nr,DiagALzdz,DiagAzdz,DiagAUzdz,bAz);
            
            for (i=0; i<=int(Nr-1); i++){
                qdagA[i][j][s]=bAz[i];}
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bAz[i]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qdagA_0[i][j]=qdagA[i][j][s];}
        }
        
      destroy_1d_double_array(DiagALr);
      destroy_1d_double_array(DiagAUr);
      destroy_1d_double_array(DiagAr);
      destroy_1d_double_array(DiagALrdr);
      destroy_1d_double_array(DiagAUrdr);
      destroy_1d_double_array(DiagArdr);
      destroy_1d_double_array(DiagALz);
      destroy_1d_double_array(DiagAUz);
      destroy_1d_double_array(DiagAz);
      destroy_1d_double_array(DiagALzdz);
      destroy_1d_double_array(DiagAUzdz);
      destroy_1d_double_array(DiagAzdz);
      destroy_1d_double_array(bAr);
      destroy_1d_double_array(bAz);
    }

        
        
        
        








