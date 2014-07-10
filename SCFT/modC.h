


void C_Matrix_r(int ii,double *DiagCz, double *DiagCLz, double *DiagCUz){
    
    
    int i;
    
    for (i=0;i<=int(Nr-1);i++){
        DiagCz[i]=0;}
    for (i=0;i<=int(Nr-2);i++){
        DiagCUz[i]=0;
        DiagCLz[i]=0;}
    
    
    for (i=0; i<=int(Nr-1); i++){
        DiagCz[i]=1.0+(delt/(pow((double)delr,(int)2)))+((delt/2.0)*wC[i][ii]);}
    for (i=1; i<=int(Nr-2);i++){
        DiagCUz[i]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*i)+r_0)*4.0*delr));}
    for (i=0; i<=int(Nr-3);i++){
        DiagCLz[i]=-(delt/(2.0*pow((double) delr,(int)2)))+(delt/(((delr*(i+1))+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagCUz[0]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(0))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(0))+r_0)*4.0*delr));
    
    DiagCLz[Nr-2]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(int(Nr-2)))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(int(Nr-2)))+r_0)*4.0*delr));

}

void C_Matrix_z(int ii,double *DiagCr, double *DiagCLr, double *DiagCUr){
    
    int j;
    
    for (j=0;j<=int(Nz-1);j++){
        DiagCr[j]=0;}
    for (j=0;j<=int(Nz-2);j++){
        DiagCUr[j]=0;
        DiagCLr[j]=0;}
    
    
    for (j=0; j<=int(Nz-1); j++){
        DiagCr[j]=1.0+delt/(pow((double)delz,(int)2));}
    for (j=1; j<=int(Nz-2); j++){
        DiagCUr[j]=-delt/(2.0*pow((double)delz,(int)2));}
    for (j=0; j<=int(Nz-3); j++){
        DiagCLr[j]=-delt/(2.0*pow((double)delz,(int)2));}
    
    DiagCUr[0]=2.0*DiagCUr[0];
    DiagCLr[int(Nz-2)]=2.0*DiagCLr[int(Nz-2)];
    
}


/**************************Finally time to define some propagators**********************/
void qC_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
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
    
    bCr=create_1d_double_array(Nz, "bCr");
    bCz=create_1d_double_array(Nr, "bCz");
    DiagCr=create_1d_double_array(Nz, "DiagCr");
    DiagCUr=create_1d_double_array(((int)Nz-1), "DiagCUr");
    DiagCLr=create_1d_double_array(((int)Nz-1), "DiagCLr");
    DiagCz=create_1d_double_array(Nr, "DiagCz");
    DiagCUz=create_1d_double_array(((int)Nr-1), "DiagCUz");
    DiagCLz=create_1d_double_array(((int)Nr-1), "DiagCLz");
    DiagCrdr=create_1d_double_array(Nz, "DiagCrdr");
    DiagCUrdr=create_1d_double_array(((int)Nz-1), "DiagCUrdr");
    DiagCLrdr=create_1d_double_array(((int)Nz-1), "DiagCLrdr");
    DiagCzdz=create_1d_double_array(Nr, "DiagCzdz");
    DiagCUzdz=create_1d_double_array(((int)Nr-1), "DiagCUzdz");
    DiagCLzdz=create_1d_double_array(((int)Nr-1), "DiagCLzdz");

    
    for (i=0; i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            qC_0[i][j]=0;
            for (s=0;s<=int(Ns-1);s++){
                qC[i][j][s]=0;}}}
    
    
    //Initialize the qs
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=1;j<=int(Nz-1);j++){
            qC_0[i][j]=1.0;
            qC[i][j][0]=1.0;}}
    
    for (s=0;s<=NC;s++){
        
        for (j=0;j<=int(Nz-1);j++){
            bCr[j]=0;}
        
        /********************************scan over z***********************************************/
        for (i=0;i<=int(Nr-1);i++){
            C_Matrix_z(i,DiagCr,DiagCLr,DiagCUr);
            if (i==0){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qC_0[i][j]+betaU*qC_0[int(i+1)][j]+betaL*qC_0[int(i+1)][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qC_0[i][j]+betaU*qC_0[int(i-1)][j]+betaL*qC_0[int(i-1)][j];}
            }
            else {
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qC_0[i][j]+betaU*qC_0[int(i+1)][j]+betaL*qC_0[int(i-1)][j];}
            }
            for (j=0;j<=int(Nz-1);j++){
                DiagCrdr[j]=DiagCr[j];}
            for (j=0;j<=int(Nz-2);j++){
                DiagCUrdr[j]=DiagCUr[j];
                DiagCLrdr[j]=DiagCLr[j];}
            
            TDMA(Nz,DiagCLrdr,DiagCrdr,DiagCUrdr,bCr);
            
            for (j=0; j<=int(Nz-1); j++){
                qC[i][j][s]=bCr[j];}
            
        }
        for (j=0;j<=int(Nz-1);j++){
            bCr[j]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qC_0[i][j]=qC[i][j][s];}
        }
        
        for (i=0;j<=int(Nr-1);i++){
            bCz[i]=0;
        }
        /***********************************scan over r***************************************************/
        for (j=0;j<=int(Nz-1);j++){
            C_Matrix_r(j,DiagCz,DiagCLz,DiagCUz);
            if (j==0){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bCz[i]=gamma*qC_0[i][j]+beta*qC_0[i][int(j+1)]+beta*qC_0[i][int(j+1)];}
            }
            else if (j==int(Nz-1)){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bCz[i]=gamma*qC_0[i][j]+beta*qC_0[i][int(j-1)]+beta*qC_0[i][int(j-1)];}
            }
            else {
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bCz[i]=gamma*qC_0[i][j]+beta*qC_0[i][int(j+1)]+beta*qC_0[i][int(j-1)];}
            }
            for (i=0;i<=int(Nr-1);i++){
                DiagCzdz[i]=DiagCz[i];}
            for (i=0;i<=int(Nr-2);i++){
                DiagCUzdz[i]=DiagCUz[i];
                DiagCLzdz[i]=DiagCLz[i];}
            
            TDMA(Nr,DiagCLzdz,DiagCzdz,DiagCUzdz,bCz);
            
            for (i=0; i<=int(Nr-1); i++){
                qC[i][j][s]=bCz[i];}
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bCz[i]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qC_0[i][j]=qC[i][j][s];}
        }
    }
    
    destroy_1d_double_array(DiagCLr);
    destroy_1d_double_array(DiagCUr);
    destroy_1d_double_array(DiagCr);
    destroy_1d_double_array(DiagCLrdr);
    destroy_1d_double_array(DiagCUrdr);
    destroy_1d_double_array(DiagCrdr);
    destroy_1d_double_array(DiagCLz);
    destroy_1d_double_array(DiagCUz);
    destroy_1d_double_array(DiagCz);
    destroy_1d_double_array(DiagCLzdz);
    destroy_1d_double_array(DiagCUzdz);
    destroy_1d_double_array(DiagCzdz);
    destroy_1d_double_array(bCr);
    destroy_1d_double_array(bCz);
}

/***********************************Define the complementary propagator***********************************/

void qdagC_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
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
    
    bCr=create_1d_double_array(Nz, "bCr");
    bCz=create_1d_double_array(Nr, "bCz");
    DiagCr=create_1d_double_array(Nz, "DiagCr");
    DiagCUr=create_1d_double_array(((int)Nz-1), "DiagCUr");
    DiagCLr=create_1d_double_array(((int)Nz-1), "DiagCLr");
    DiagCz=create_1d_double_array(Nr, "DiagCz");
    DiagCUz=create_1d_double_array(((int)Nr-1), "DiagCUz");
    DiagCLz=create_1d_double_array(((int)Nr-1), "DiagCLz");
    DiagCrdr=create_1d_double_array(Nz, "DiagCrdr");
    DiagCUrdr=create_1d_double_array(((int)Nz-1), "DiagCUrdr");
    DiagCLrdr=create_1d_double_array(((int)Nz-1), "DiagCLrdr");
    DiagCzdz=create_1d_double_array(Nr, "DiagCzdz");
    DiagCUzdz=create_1d_double_array(((int)Nr-1), "DiagCUzdz");
    DiagCLzdz=create_1d_double_array(((int)Nr-1), "DiagCLzdz");

    for (i=0; i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            qdagC_0[i][j]=0;
            for (s=0;s<=int(Ns-1);s++){
                qdagC[i][j][s]=0;
            }
        }
    }
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            qdagC[i][j][0]=qA[i][j][NA]*qB[i][j][NB];
            qdagC_0[i][j]=(qA[i][j][NA])*(qB[i][j][NB]);
        }
    }
    
    for (s=0;s<=NC;s++){
        
        for (j=0;j<=int(Nz-1);j++){
            bCr[j]=0;
            
        }
        
        /********************************scan over z***********************************************/
        for (i=0;i<=int(Nr-1);i++){
            C_Matrix_z(i,DiagCr,DiagCLr,DiagCUr);
            if (i==0){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qdagC_0[i][j]+betaU*qdagC_0[int(i+1)][j]+betaL*qdagC_0[int(i+1)][j];}
            }
            else if (i==int(Nr-1)){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qdagC_0[i][j]+betaU*qdagC_0[int(i-1)][j]+betaL*qdagC_0[int(i-1)][j];}
            }
            else {
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qdagC_0[i][j]+betaU*qdagC_0[int(i+1)][j]+betaL*qdagC_0[int(i-1)][j];}
            }
            for (j=0;j<=int(Nz-1);j++){
                DiagCrdr[j]=DiagCr[j];}
            for (j=0;j<=int(Nz-2);j++){
                DiagCUrdr[j]=DiagCUr[j];
                DiagCLrdr[j]=DiagCLr[j];}
            
            TDMA(Nz,DiagCLrdr,DiagCrdr,DiagCUrdr,bCr);
            
            for (j=0; j<=int(Nz-1); j++){
                qdagC[i][j][s]=bCr[j];}
            
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bCr[j]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qdagC_0[i][j]=qdagC[i][j][s];}
        }
    }
    
    for (i=0;i<=int(Nr-1);i++){
        bCz[i]=0;
    }
    /***********************************scan over r***************************************************/
    for (j=0;j<=int(Nz-1);j++){
        C_Matrix_r(j,DiagCz,DiagCLz,DiagCUz);
        if (j==0){
            for (i=0;i<=int(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bCz[i]=gamma*qdagC_0[i][j]+beta*qdagC_0[i][int(j+1)]+beta*qdagC_0[i][int(j+1)];}
        }
        else if (j==int(Nz-1)){
            for (i=0;i<=int(Nz-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bCz[i]=gamma*qdagC_0[i][j]+beta*qdagC_0[i][int(j-1)]+beta*qdagC_0[i][int(j-1)];}
        }
        else {
            for (i=0;i<=int(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bCz[i]=gamma*qdagC_0[i][j]+beta*qdagC_0[i][int(j+1)]+beta*qdagC_0[i][int(j-1)];}
        }
        for (i=0;i<=int(Nr-1);i++){
            DiagCzdz[i]=DiagCz[i];}
        for (i=0;i<=int(Nr-2);i++){
            DiagCUzdz[i]=DiagCUz[i];
            DiagCLzdz[i]=DiagCLz[i];}
        
        TDMA(Nr,DiagCLzdz,DiagCzdz,DiagCUzdz,bCz);
        
        for (i=0; i<=int(Nr-1); i++){
            qdagC[i][j][s]=bCz[i];}
    }
    
    for (i=0;i<=int(Nr-1);i++){
        bCz[i]=0;
    }
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=0; j<=int(Nz-1); j++){
            qdagC_0[i][j]=qdagC[i][j][s];}
    }
    destroy_1d_double_array(DiagCLr);
    destroy_1d_double_array(DiagCUr);
    destroy_1d_double_array(DiagCr);
    destroy_1d_double_array(DiagCLrdr);
    destroy_1d_double_array(DiagCUrdr);
    destroy_1d_double_array(DiagCrdr);
    destroy_1d_double_array(DiagCLz);
    destroy_1d_double_array(DiagCUz);
    destroy_1d_double_array(DiagCz);
    destroy_1d_double_array(DiagCLzdz);
    destroy_1d_double_array(DiagCUzdz);
    destroy_1d_double_array(DiagCzdz);
    destroy_1d_double_array(bCr);
    destroy_1d_double_array(bCz);

}









