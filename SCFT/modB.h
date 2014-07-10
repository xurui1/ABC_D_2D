
void B_Matrix_z(int ii, double *DiagBr, double *DiagBLr,double *DiagBUr){
    
    int i;
    
    for (i=0;i<=int(Nr-1);i++){
        DiagBr[i]=0;}
    for (i=0;i<=int(Nr-2);i++){
        DiagBUr[i]=0;
        DiagBLr[i]=0;}
    
    
    for (i=0; i<=int(Nr-1); i++){
        DiagBr[i]=1.0+(delt/(pow((double)delr,(int)2)))+((delt/4.0)*wB[i][ii]);}
    for (i=1; i<=int(Nr-2);i++){
        DiagBUr[i]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*i)+r_0)*4.0*delr));}
    for (i=0; i<=int(Nr-3);i++){
        DiagBLr[i]=-(delt/(2.0*pow((double) delr,(int)2)))+(delt/(((delr*i)+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagBUr[0]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(0))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(0))+r_0)*4.0*delr));
    
    DiagBLr[int(Nr-2)]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(int(Nr-2)))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(int(Nr-2)))+r_0)*4.0*delr));
    
}

void B_Matrix_r(int ii, double *DiagBz, double *DiagBLz, double *DiagBUz){
    
    int j;
    
    for (j=0;j<=int(Nz-1);j++){
        DiagBz[j]=0;}
    for (j=0;j<=int(Nz-2);j++){
        DiagBUz[j]=0;
        DiagBLz[j]=0;}
    
    
    for (j=0; j<=int(Nz-1); j++){
        DiagBz[j]=1.0+delt/(pow((double)delz,(int)2));}
    for (j=1; j<=int(Nz-2); j++){
        DiagBUz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    for (j=0; j<=int(Nz-3); j++){
        DiagBLz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    
    DiagBUz[0]=2.0*DiagBUz[0];
    DiagBLz[int(Nz-2)]=2.0*DiagBLz[int(Nz-2)];
    
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
            for (s=0;s<=Ns;s++){
                qB[i][j][s]=0;}}}
    
    
    //Initialize the qs
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=1;j<=int(Nz-1);j++){
            qB_0[i][j]=1.0;
            qB[i][j][0]=1.0;}}
    
    for (s=0;s<=NB;s++){
        
        for (j=0;j<=int(Nz-1);j++){
            bBz[j]=0;}
        
        /********************************scan over z******************************/
        for (i=0;i<=int(Nr-1);i++){
            B_Matrix_r(i,DiagBz,DiagBLz,DiagBUz);
            if (i==0){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bBz[j]=gamma*qB_0[i][j]+betaU*qB_0[int(i+1)][j]+betaL*qB_0[int(i+1)][j];}
            }
            else if (i==int(Nr-1)){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bBz[j]=gamma*qB_0[i][j]+betaU*qB_0[int(i-1)][j]+betaL*qB_0[int(i-1)][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bBz[j]=gamma*qB_0[i][j]+betaU*qB_0[int(i+1)][j]+betaL*qB_0[int(i-1)][j];}
            }
            for (j=0;j<=int(Nz-1);j++){
                DiagBzdz[j]=DiagBz[j];}
            for (j=0;j<=int(Nz-2);j++){
                DiagBUzdz[j]=DiagBUz[j];
                DiagBLzdz[j]=DiagBLz[j];}
            
            TDMA(Nz,DiagBLzdz,DiagBzdz,DiagBUzdz,bBz);
            
            for (j=0; j<=int(Nz-1); j++){
                qB[i][j][s]=bBz[j];}
            
        }
        for (j=0;j<=int(Nz-1);j++){
            bBz[j]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qB_0[i][j]=qB[i][j][s];}
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bBr[i]=0;
        }
        /***********************************scan over r***************************/
        for (j=0;j<=int(Nz-1);j++){
            B_Matrix_z(j,DiagBr,DiagBLr, DiagBUr);
            if (j==0){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bBr[i]=gamma*qB_0[i][j]+beta*qB_0[i][int(j+1)]+beta*qB_0[i][int(j+1)];}
            }
            else if (j==(Nz-1)){
                for (i=0;i<=(Nz-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bBr[i]=gamma*qB_0[i][j]+beta*qB_0[i][int(j-1)]+beta*qB_0[i][int(j-1)];}
            }
            else {
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bBr[i]=gamma*qB_0[i][j]+beta*qB_0[i][int(j+1)]+beta*qB_0[i][int(j-1)];}
            }
            for (i=0;i<=int(Nr-1);i++){
                DiagBrdr[i]=DiagBr[i];}
            for (i=0;i<=int(Nr-2);i++){
                DiagBUrdr[i]=DiagBUr[i];
                DiagBLrdr[i]=DiagBLr[i];}
            
            TDMA(Nr,DiagBLrdr,DiagBrdr,DiagBUrdr,bBr);
            
            for (i=0; i<=int(Nr-1); i++){
                qB[i][j][s]=bBr[i];}
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bBr[i]=0;
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
            for (s=0;s<=Ns;s++){
                qdagB[i][j][s]=0;
            }
        }
    }
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            qdagB[i][j][0]=qA[i][j][NA]*qC[i][j][NC];
            qdagB_0[i][j]=(qA[i][j][NA])*(qC[i][j][NC]);
        }
    }
    
    for (s=0;s<=NB;s++){
        
        for (j=0;j<=int(Nz-1);j++){
            bBz[j]=0;
        }
        
        /********************************scan over z***********************************************/
        for (i=0;i<=int(Nr-1);i++){
            B_Matrix_r(i,DiagBz,DiagBLz,DiagBUz);
            if (i==0){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bBz[j]=gamma*qdagB_0[i][j]+betaU*qdagB_0[int(i+1)][j]+betaL*qdagB_0[int(i+1)][j];}
            }
            else if (i==int(Nr-1)){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bBz[j]=gamma*qdagB_0[i][j]+betaU*qdagB_0[int(i-1)][j]+betaL*qdagB_0[int(i-1)][j];}
            }
            else {
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bBz[j]=gamma*qdagB_0[i][j]+betaU*qdagB_0[int(i+1)][j]+betaL*qdagB_0[int(i-1)][j];}
            }
            for (j=0;j<=int(Nz-1);j++){
                DiagBzdz[j]=DiagBz[j];}
            for (j=0;j<=int(Nz-2);j++){
                DiagBUzdz[j]=DiagBUz[j];
                DiagBLzdz[j]=DiagBLz[j];}
            
            TDMA(Nz,DiagBLzdz,DiagBzdz,DiagBUzdz,bBz);
            
            for (j=0; j<=int(Nz-1); j++){
                qdagB[i][j][s]=bBz[j];}
            
        }
        
        for (j=0;j<=int(Nz-1);j++){
            bBz[j]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qdagB_0[i][j]=qdagB[i][j][s];}
        }
        
        
    }
    
    for (i=0;i<=int(Nr-1);i++){
        bBr[i]=0;
    }
    /***********************************scan over r*****************************/
    for (j=0;j<=int(Nz-1);j++){
        B_Matrix_z(j,DiagBr,DiagBLr, DiagBUr);

        if (j==0){
            for (i=0;i<=int(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bBr[i]=gamma*qdagB_0[i][j]+beta*qdagB_0[i][int(j+1)]+beta*qdagB_0[i][int(j+1)];}
        }
        else if (j==int(Nz-1)){
            for (i=0;i<=int(Nz-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bBr[i]=gamma*qdagB_0[i][j]+beta*qdagB_0[i][int(j-1)]+beta*qdagB_0[i][int(j-1)];}
        }
        else {
            for (i=0;i<=int(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bBr[i]=gamma*qdagB_0[i][j]+beta*qdagB_0[i][int(j+1)]+beta*qdagB_0[i][int(j-1)];}
        }
        for (i=0;i<=int(Nr-1);i++){
            DiagBrdr[i]=DiagBr[i];}
        for (i=0;i<=int(Nr-2);i++){
            DiagBUrdr[i]=DiagBUr[i];
            DiagBLrdr[i]=DiagBLr[i];}
        
        TDMA(Nr,DiagBLrdr,DiagBrdr,DiagBUrdr,bBr);
        
        for (i=0; i<=int(Nr-1); i++){
            qdagB[i][j][s]=bBr[i];}
    }
    
    for (i=0;i<=int(Nr-1);i++){
        bBr[i]=0;
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









