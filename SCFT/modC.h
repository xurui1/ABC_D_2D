



void C_Matrix_r(int ii){
    
    
    int i;
    
    for (i=0;i<=Nr-1;i++){
        DiagCr[i]=0;}
    for (i=0;i<=Nr-2;i++){
        DiagCUr[i]=0;
        DiagCLr[i]=0;}
    
    
    for (i=0; i<=(Nr-1); i++){
        DiagCr[i]=1.0+(delt/(pow((double)delr,(int)2)))+((delt/2.0)*wC[i][ii]);}
    for (i=1; i<=(Nr-2);i++){
        DiagCUr[i]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(i-1))+r_0)*4.0*delr));}
    for (i=0; i<=(Nr-2);i++){
        DiagCLr[i]=-(delt/(2.0*pow((double) delr,(int)2)))+(delt/(((delr*(i+1-1))+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagCUr[0]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(1-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(1-1))+r_0)*4.0*delr));
    
    DiagCLr[Nr-1]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(Nr-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(Nr-1))+r_0)*4.0*delr));
    
}

void C_Matrix_z(int ii){
    
    int j;
    
    for (j=0;j<Nz;j++){
        DiagCz[j]=0;}
    for (j=0;j<(int)(Nz-1);j++){
        DiagCUz[j]=0;
        DiagCLz[j]=0;}
    
    
    for (j=0; j<=(Nz-1); j++){
        DiagCz[j]=1.0+delt/(pow((double)delz,(int)2));}
    for (j=0; j<=(Nz-2); j++){
        DiagCUz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    for (j=0; j<=(Nz-2); j++){
        DiagCLz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    
    DiagCUz[0]=2.0*DiagCUz[1];
    DiagCLz[Nz-1]=2.0*DiagCLz[Nz-1];
    
}


/**************************Finally time to define some propagators**********************/
void qC_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
    
    qC=create_3d_double_array(Nr,Nz,Ns,"qC");
    qC_0=create_2d_double_array(Nr,Nz, "qC_0");
    bCr=create_1d_double_array(Nr, "bCr");
    bCz=create_1d_double_array(Nz, "bCz");
    DiagCr=create_1d_double_array(Nr, "DiagCr");
    DiagCUr=create_1d_double_array(((int)Nr-1), "DiagCUr");
    DiagCLr=create_1d_double_array(((int)Nr-1), "DiagCLr");
    DiagCz=create_1d_double_array(Nr, "DiagCz");
    DiagCUz=create_1d_double_array(((int)Nr-1), "DiagCUz");
    DiagCLz=create_1d_double_array(((int)Nr-1), "DiagCLz");
    DiagCrdr=create_1d_double_array(Nr, "DiagCrdr");
    DiagCUrdr=create_1d_double_array(((int)Nr-1), "DiagCUrdr");
    DiagCLrdr=create_1d_double_array(((int)Nr-1), "DiagCLrdr");
    DiagCzdz=create_1d_double_array(Nr, "DiagCzdz");
    DiagCUzdz=create_1d_double_array(((int)Nr-1), "DiagCUzdz");
    DiagCLzdz=create_1d_double_array(((int)Nr-1), "DiagCLzdz");
    
    for (i=0; i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            qC_0[i][j]=0;
            for (s=0;s<=Ns-1;s++){
                qC[i][j][s]=0;}}}
    
    
    //Initialize the qs
    
    for (i=0;i<=(Nr-1);i++){
        for (j=1;j<=(Nz-1);j++){
            qC_0[i][j]=1.0;
            qC[i][j][0]=1.0;}}
    
    for (s=0;s<=(NC-1);s++){
        
        for (i=0;i<=Nr-1;i++){
            bCr[i]=0;}
        
        /********************************scan over z***********************************************/
        for (i=0;i<=(Nr-1);i++){
            B_Matrix_z(i);
            if (i==0){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qC_0[i][j]+betaU*qC_0[i+1][j]+betaL*qC_0[i+1][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qC_0[i][j]+betaU*qC_0[i-1][j]+betaL*qC_0[i-1][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qC_0[i][j]+betaU*qC_0[i+1][j]+betaL*qC_0[i-1][j];}
            }
            DiagCrdr[i]=DiagCr[i];
            DiagCUrdr[i]=DiagCUr[i];
            DiagCLrdr[i]=DiagCLr[i];
            
            TDMA(Nz,DiagCLrdr,DiagCrdr,DiagCUrdr,bCr);
            
            for (j=0; j<=Nz-1; j++){
                qC[i][j][s]=bCr[j];}
            
        }
        for (i=0;i<=Nr-1;i++){
            bCr[i]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qC_0[i][j]=qC[i][j][s];}
        }
        
        for (j=0;j<=Nz-1;j++){
            bCz[j]=0;
        }
        /***********************************scan over r***************************************************/
        for (j=0;j<=(Nz-1);j++){
            B_Matrix_r(j);
            if (j==0){
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bCz[i]=gamma*qC_0[i][j]+beta*qC_0[j+1][i]+beta*qC_0[j+1][i];}
            }
            else if (j==(Nz-1)){
                for (i=0;i<=(Nz-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bCz[i]=gamma*qC_0[i][j]+beta*qC_0[j-1][i]+beta*qC_0[j-1][i];}
            }
            else {
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bCz[i]=gamma*qC_0[i][j]+beta*qC_0[j+1][i]+beta*qC_0[j-1][i];}
            }
            DiagCzdz[j]=DiagCz[j];
            DiagCUzdz[j]=DiagCUz[j];
            DiagCLzdz[j]=DiagCLz[j];
            
            TDMA(Nr,DiagCLzdz,DiagCzdz,DiagCUzdz,bCz);
            
            for (i=0; i<=Nr-1; i++){
                qC[i][j][s]=bCz[i];}
        }
        
        for (i=0;i<=Nr-1;i++){
            bCr[i]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qC_0[i][j]=qC[i][j][s];}
        }
    }
    
    destroy_2d_double_array(qC_0);
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
    
    qdagC=create_3d_double_array(Nr,Nz,Ns,"qdagC");
    qdagC_0=create_2d_double_array(Nr,Nz, "qdagC_0");
    bCr=create_1d_double_array(Nr, "bCr");
    bCz=create_1d_double_array(Nz, "bCz");
    DiagCr=create_1d_double_array(Nr, "DiagCr");
    DiagCUr=create_1d_double_array(((int)Nr-1), "DiagCUr");
    DiagCLr=create_1d_double_array(((int)Nr-1), "DiagCLr");
    DiagCz=create_1d_double_array(Nr, "DiagCz");
    DiagCUz=create_1d_double_array(((int)Nr-1), "DiagCUz");
    DiagCLz=create_1d_double_array(((int)Nr-1), "DiagCLz");
    DiagCrdr=create_1d_double_array(Nr, "DiagCrdr");
    DiagCUrdr=create_1d_double_array(((int)Nr-1), "DiagCUrdr");
    DiagCLrdr=create_1d_double_array(((int)Nr-1), "DiagCLrdr");
    DiagCzdz=create_1d_double_array(Nr, "DiagCzdz");
    DiagCUzdz=create_1d_double_array(((int)Nr-1), "DiagCUzdz");
    DiagCLzdz=create_1d_double_array(((int)Nr-1), "DiagCLzdz");
    
    
    for (i=0; i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            qdagC_0[i][j]=0;
            for (s=0;s<=Ns-1;s++){
                qdagC[i][j][s]=0;
            }
        }
    }
    
    for (i=0;i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            qdagC[i][j][0]=qA[i][j][NA]*qB[i][j][NB];
            qdagC_0[i][j]=(qA[i][j][NA])*(qB[i][j][NB]);
        }
    }
    
    for (s=0;s<=NC-1;s++){
        
        for (i=0;i<=Nr-1;i++){
            bCr[i]=0;
        }
        
        /********************************scan over z***********************************************/
        for (i=0;i<=(Nr-1);i++){
            B_Matrix_z(i);
            if (i==0){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qdagC_0[i][j]+betaU*qdagC_0[i+1][j]+betaL*qdagC_0[i+1][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qdagC_0[i][j]+betaU*qdagC_0[i-1][j]+betaL*qdagC_0[i-1][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qdagC_0[i][j]+betaU*qdagC_0[i+1][j]+betaL*qdagC_0[i-1][j];}
            }
            DiagCrdr[i]=DiagCr[i];
            DiagCUrdr[i]=DiagCUr[i];
            DiagCLrdr[i]=DiagCLr[i];
            
            TDMA(Nz,DiagCLrdr,DiagCrdr,DiagCUrdr,bCr);
            
            for (j=0; j<=Nz-1; j++){
                qdagC[i][j][s]=bCr[j];}
            
        }
        
        for (i=0;i<=Nr-1;i++){
            bCr[i]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qdagC_0[i][j]=qdagC[i][j][s];}
        }
    }
    
    for (j=0;j<=Nz-1;j++){
        bCz[j]=0;
    }
    /***********************************scan over r***************************************************/
    for (j=0;j<=(Nz-1);j++){
        B_Matrix_r(j);
        if (j==0){
            for (i=0;i<=(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bCz[i]=gamma*qdagC_0[i][j]+beta*qdagC_0[j+1][i]+beta*qdagC_0[j+1][i];}
        }
        else if (j==(Nz-1)){
            for (i=0;i<=(Nz-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bCz[i]=gamma*qdagC_0[i][j]+beta*qdagC_0[j-1][i]+beta*qdagC_0[j-1][i];}
        }
        else {
            for (i=0;i<=(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bCz[i]=gamma*qdagC_0[i][j]+beta*qdagC_0[j+1][i]+beta*qdagC_0[j-1][i];}
        }
        DiagCzdz[j]=DiagCz[j];
        DiagCUzdz[j]=DiagCUz[j];
        DiagCLzdz[j]=DiagCLz[j];
        
        TDMA(Nr,DiagCLzdz,DiagCzdz,DiagCUzdz,bCz);
        
        for (i=0; i<=Nr-1; i++){
            qdagC[i][j][s]=bCz[i];}
    }
    
    for (i=0;i<=Nr-1;i++){
        bCr[i]=0;
    }
    
    for (i=0;i<=(Nr-1);i++){
        for (j=0; j<=Nz-1; j++){
            qdagC_0[i][j]=qdagC[i][j][s];}
    }
    
    destroy_2d_double_array(qdagC_0);
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









