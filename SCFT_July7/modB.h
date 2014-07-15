
void B_Matrix_r(int ii, double *DiagBz, double *DiagBLz,double *DiagBUz){
    
    int i;
    
    for (i=0;i<=int(Nr-1);i++){
        DiagBz[i]=0;}
    for (i=0;i<=int(Nr-2);i++){
        DiagBUz[i]=0;
        DiagBLz[i]=0;}
    
    
    for (i=0; i<=int(Nr-1); i++){
        DiagBz[i]=1.0+(delt/(pow((double)delr,(int)2)))+((delt/2.0)*wB[i][ii]);}
    for (i=1; i<=int(Nr-2);i++){
        DiagBUz[i]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(i-1))+r_0)*4.0*delr));}
    for (i=0; i<=int(Nr-3);i++){
        DiagBLz[i]=-(delt/(2.0*pow((double) delr,(int)2)))+(delt/(((delr*(i+1-1))+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagBUz[0]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(-1))+r_0)*4.0*delr));
    
    DiagBLz[int(Nr-2)]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(int(Nr-2)))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(int(Nr-2)))+r_0)*4.0*delr));
    
}

void B_Matrix_z(int ii, double *DiagBr, double *DiagBLr, double *DiagBUr){
    
    int j;
    
    for (j=0;j<=int(Nz-1);j++){
        DiagBr[j]=0;}
    for (j=0;j<=int(Nz-2);j++){
        DiagBUr[j]=0;
        DiagBLr[j]=0;}
    
    
    for (j=0; j<=int(Nz-1); j++){
        DiagBr[j]=1.0+delt/(pow((double)delz,(int)2));}
    for (j=1; j<=int(Nz-2); j++){
        DiagBUr[j]=-delt/(2.0*pow((double)delz,(int)2));}
    for (j=0; j<=int(Nz-3); j++){
        DiagBLr[j]=-delt/(2.0*pow((double)delz,(int)2));}
    
    DiagBUr[0]=2.0*DiagBUr[0];
    DiagBLr[int(Nz-2)]=2.0*DiagBLr[int(Nz-2)];
    
}


/**************************Finally time to define some propagators**********************/
void qB_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
    
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
    
    
    bBr=create_1d_double_array(Nz, "bBr");
    bBz=create_1d_double_array(Nr, "bBz");
    DiagBr=create_1d_double_array(Nz, "DiagBr");
    DiagBUr=create_1d_double_array(((int)Nz-1), "DiagBUr");
    DiagBLr=create_1d_double_array(((int)Nz-1), "DiagBLr");
    DiagBz=create_1d_double_array(Nr, "DiagBz");
    DiagBUz=create_1d_double_array(((int)Nr-1), "DiagBUz");
    DiagBLz=create_1d_double_array(((int)Nr-1), "DiagBLz");
    DiagBrdr=create_1d_double_array(Nz, "DiagBrdr");
    DiagBUrdr=create_1d_double_array(((int)Nz-1), "DiagBUrdr");
    DiagBLrdr=create_1d_double_array(((int)Nz-1), "DiagBLrdr");
    DiagBzdz=create_1d_double_array(Nr, "DiagBzdz");
    DiagBUzdz=create_1d_double_array(((int)Nr-1), "DiagBUzdz");
    DiagBLzdz=create_1d_double_array(((int)Nr-1), "DiagBLzdz");
    
    for (i=0; i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            qB_0[i][j]=0;
            for (s=0;s<=int(Ns-1);s++){
                qB[i][j][s]=0;}}}
    
    
    //Initialize the qs
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=1;j<=int(Nz-1);j++){
            qB_0[i][j]=1.0;
            qB[i][j][0]=1.0;}}
    
    for (s=0;s<=int(NB-1);s++){
        
        for (j=0;j<=int(Nz-1);j++){
            bBr[j]=0;}
        
        /********************************scan over z***********************************************/
        for (i=0;i<=int(Nr-1);i++){
            B_Matrix_z(i,DiagBr,DiagBLr,DiagBUr);
            if (i==0){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qB_0[i][j]+betaU*qB_0[int(i+1)][j]+betaL*qB_0[int(i+1)][j];}
            }
            else if (i==int(Nr-1)){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qB_0[i][j]+betaU*qB_0[int(i-1)][j]+betaL*qB_0[int(i-1)][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qB_0[i][j]+betaU*qB_0[int(i+1)][j]+betaL*qB_0[int(i-1)][j];}
            }
            for (j=0;j<=int(Nz-1);j++){
                DiagBrdr[j]=DiagBr[j];}
            for (j=0;j<=int(Nz-2);j++){
                DiagBUrdr[j]=DiagBUr[j];
                DiagBLrdr[j]=DiagBLr[j];}
            
            TDMA(Nz,DiagBLrdr,DiagBrdr,DiagBUrdr,bBr);
            
            for (j=0; j<=int(Nz-1); j++){
                qB[i][j][s]=bBr[j];}
            
        }
        for (j=0;j<=int(Nz-1);j++){
            bBr[j]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qB_0[i][j]=qB[i][j][s];}
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bBz[i]=0;
        }
        /***********************************scan over r***************************************************/
        for (j=0;j<=int(Nz-1);j++){
            B_Matrix_r(j,DiagBz,DiagBLz, DiagBUz);
            if (j==0){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bBz[i]=gamma*qB_0[i][j]+beta*qB_0[i][int(j+1)]+beta*qB_0[i][int(j+1)];}
            }
            else if (j==(Nz-1)){
                for (i=0;i<=(Nz-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bBz[i]=gamma*qB_0[i][j]+beta*qB_0[i][int(j-1)]+beta*qB_0[i][int(j-1)];}
            }
            else {
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bBz[i]=gamma*qB_0[i][j]+beta*qB_0[i][int(j+1)]+beta*qB_0[i][int(j-1)];}
            }
            for (i=0;i<=int(Nr-1);i++){
                DiagBzdz[i]=DiagBz[i];}
            for (i=0;i<=int(Nr-2);i++){
                DiagBUzdz[i]=DiagBUz[i];
                DiagBLzdz[i]=DiagBLz[i];}
            
            TDMA(Nr,DiagBLzdz,DiagBzdz,DiagBUzdz,bBz);
            
            for (i=0; i<=int(Nr-1); i++){
                qB[i][j][s]=bBz[i];}
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bBz[i]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qB_0[i][j]=qB[i][j][s];}
        }
    }
    
    destroy_1d_double_array(DiagBLr);
    destroy_1d_double_array(DiagBUr);
    destroy_1d_double_array(DiagBr);
    destroy_1d_double_array(DiagBLrdr);
    destroy_1d_double_array(DiagBUrdr);
    destroy_1d_double_array(DiagBrdr);
    destroy_1d_double_array(DiagBLz);
    destroy_1d_double_array(DiagBUz);
    destroy_1d_double_array(DiagBz);
    destroy_1d_double_array(DiagBLzdz);
    destroy_1d_double_array(DiagBUzdz);
    destroy_1d_double_array(DiagBzdz);
    destroy_1d_double_array(bBr);
    destroy_1d_double_array(bBz);
}

/***********************************Define the complementary propagator***********************************/

void qdagB_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
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
    
    bBr=create_1d_double_array(Nz, "bBr");
    bBz=create_1d_double_array(Nr, "bBz");
    DiagBr=create_1d_double_array(Nz, "DiagBr");
    DiagBUr=create_1d_double_array(((int)Nz-1), "DiagBUr");
    DiagBLr=create_1d_double_array(((int)Nz-1), "DiagBLr");
    DiagBz=create_1d_double_array(Nr, "DiagBz");
    DiagBUz=create_1d_double_array(((int)Nr-1), "DiagBUz");
    DiagBLz=create_1d_double_array(((int)Nr-1), "DiagBLz");
    DiagBrdr=create_1d_double_array(Nz, "DiagBrdr");
    DiagBUrdr=create_1d_double_array(((int)Nz-1), "DiagBUrdr");
    DiagBLrdr=create_1d_double_array(((int)Nz-1), "DiagBLrdr");
    DiagBzdz=create_1d_double_array(Nr, "DiagBzdz");
    DiagBUzdz=create_1d_double_array(((int)Nr-1), "DiagBUzdz");
    DiagBLzdz=create_1d_double_array(((int)Nr-1), "DiagBLzdz");
    
    for (i=0; i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            qdagB_0[i][j]=0;
            for (s=0;s<=int(Ns-1);s++){
                qdagB[i][j][s]=0;
            }
        }
    }
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            qdagB[i][j][0]=qA[i][j][int(NA-1)]*qC[i][j][int(NC-1)];
            qdagB_0[i][j]=(qA[i][j][int(NA-1)])*(qC[i][j][int(NC-1)]);
        }
    }
    
    for (s=0;s<=int(NB-1);s++){
        
        for (j=0;j<=int(Nz-1);j++){
            bBr[j]=0;
        }
        
        /********************************scan over z***********************************************/
        for (i=0;i<=int(Nr-1);i++){
            B_Matrix_z(i,DiagBr,DiagBLr,DiagBUr);
            if (i==0){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qdagB_0[i][j]+betaU*qdagB_0[int(i+1)][j]+betaL*qdagB_0[int(i+1)][j];}
            }
            else if (i==int(Nr-1)){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qdagB_0[i][j]+betaU*qdagB_0[int(i-1)][j]+betaL*qdagB_0[int(i-1)][j];}
            }
            else {
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qdagB_0[i][j]+betaU*qdagB_0[int(i+1)][j]+betaL*qdagB_0[int(i-1)][j];}
            }
            for (j=0;j<=int(Nz-1);j++){
                DiagBrdr[j]=DiagBr[j];}
            for (j=0;j<=int(Nz-2);j++){
                DiagBUrdr[j]=DiagBUr[j];
                DiagBLrdr[j]=DiagBLr[j];}
            
            TDMA(Nz,DiagBLrdr,DiagBrdr,DiagBUrdr,bBr);
            
            for (j=0; j<=int(Nz-1); j++){
                qdagB[i][j][s]=bBr[j];}
            
        }
        
        for (j=0;j<=int(Nz-1);j++){
            bBr[j]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qdagB_0[i][j]=qdagB[i][j][s];}
        }
        
        
    }
    
    for (i=0;i<=int(Nr-1);i++){
        bBz[i]=0;
    }
    /***********************************scan over r***************************************************/
    for (j=0;j<=int(Nz-1);j++){
        B_Matrix_r(j,DiagBz,DiagBLz, DiagBUz);

        if (j==0){
            for (i=0;i<=int(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bBz[i]=gamma*qdagB_0[i][j]+beta*qdagB_0[i][int(j+1)]+beta*qdagB_0[i][int(j+1)];}
        }
        else if (j==int(Nz-1)){
            for (i=0;i<=int(Nz-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bBz[i]=gamma*qdagB_0[i][j]+beta*qdagB_0[i][int(j-1)]+beta*qdagB_0[i][int(j-1)];}
        }
        else {
            for (i=0;i<=int(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bBz[i]=gamma*qdagB_0[i][j]+beta*qdagB_0[i][int(j+1)]+beta*qdagB_0[i][int(j-1)];}
        }
        for (i=0;i<=int(Nr-1);i++){
            DiagBzdz[i]=DiagBz[i];}
        for (i=0;i<=int(Nr-2);i++){
            DiagBUzdz[i]=DiagBUz[i];
            DiagBLzdz[i]=DiagBLz[i];}
        
        TDMA(Nr,DiagBLzdz,DiagBzdz,DiagBUzdz,bBz);
        
        for (i=0; i<=int(Nr-1); i++){
            qdagB[i][j][s]=bBz[i];}
    }
    
    for (i=0;i<=int(Nr-1);i++){
        bBz[i]=0;
    }
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=0; j<=int(Nz-1); j++){
            qdagB_0[i][j]=qdagB[i][j][s];}
    }
    
    destroy_1d_double_array(DiagBLr);
    destroy_1d_double_array(DiagBUr);
    destroy_1d_double_array(DiagBr);
    destroy_1d_double_array(DiagBLrdr);
    destroy_1d_double_array(DiagBUrdr);
    destroy_1d_double_array(DiagBrdr);
    destroy_1d_double_array(DiagBLz);
    destroy_1d_double_array(DiagBUz);
    destroy_1d_double_array(DiagBz);
    destroy_1d_double_array(DiagBLzdz);
    destroy_1d_double_array(DiagBUzdz);
    destroy_1d_double_array(DiagBzdz);
    destroy_1d_double_array(bBr);
    destroy_1d_double_array(bBz);
    
   
}









