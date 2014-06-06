

int C_Matrix_r(int ii){

    
    int i;
    //int j;
    //double alphaC,betaCL,betaCU;
    
    DiagCr=create_1d_double_array(Nr, "DiagCr");
    DiagCUr=create_1d_double_array(int(Nr-1), "DiagCUr");
    DiagCLr=create_1d_double_array(Nr-1, "DiagCLr");
    
    *DiagCr=0;
    *DiagCUr=0;
    *DiagCLr=0;
    
    for (i=0; i<=(Nr-1); i++){
        DiagCr[i]=1.0+(delt/(pow(delr,2))+((delt/2.0)*wC[i][ii]));}
    for (i=1; i<=(Nr-2);i++){
        DiagCUr[i]=-(delt/(2.0*pow(delr,2)))-(delt/(((delr*(i-1))+r_0)*4.0*delr));}
    for (i=0; i<=(Nr-2);i++){
        DiagCLr[i]=-(delt/(2.0*pow(delr,2)))+(delt/(((delr*(i+1-1))+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagCUr[0]=-(delt/(2.0*pow(delr,2)))-(delt/(((delr*(1-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow(delr,2)))+(delt/(((delr*(1-1))+r_0)*4.0*delr));
    
    DiagCLr[Nr-1]=-(delt/(2.0*pow(delr,2)))-(delt/(((delr*(Nr-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow(delr,2)))+(delt/(((delr*(Nr-1))+r_0)*4.0*delr));
    
    
    return C_Matrix_r(ii);
}

int C_Matrix_z(int ii){
    
    int i;
    //int j;
    //double alphaC,betaCL,betaCU;
    
    DiagCz=create_1d_double_array(Nr, "DiagCz");
    DiagCUz=create_1d_double_array(Nr-1, "DiagCUz");
    DiagCLz=create_1d_double_array(Nr-1, "DiagCLz");
    
    *DiagCz=0;
    *DiagCUz=0;
    *DiagCLz=0;
    
    for (i=0; i<=(Nz-1); i++){
        DiagCz[i]=1.0+delt/(pow(delz,2));}
    for (i=0; i<=(Nz-2); i++){
        DiagCUz[i]=-delt/(2.0*pow(delz,2));}
    for (i=0; i<=(Nz-2); i++){
        DiagCLz[i]=-delt/(2.0*pow(delz,2));}
    
    DiagCUz[0]=2.0*DiagCUz[1];
    DiagCLz[Nz-1]=2.0*DiagCLz[Nz-1];
    
    return C_Matrix_z(ii);
}


/**************************Finally time to define some propagators**********************/
int qC_forward(){
    //int errorflag;
    int s,i,j;
    //int errorflagC;
    double gamma,betaU,betaL,beta;
    
    qC=create_3d_double_array(Nr,Nz,Ns,"qC");
    qC_0=create_2d_double_array(Nr,Nz, "qC_0");
    bCr=create_1d_double_array(Nr, "bCr");
    bCz=create_1d_double_array(Nz, "bCz");
    DiagCr=create_1d_double_array(Nr, "DiagCr");
    DiagCUr=create_1d_double_array(int(Nr-1), "DiagCUr");
    DiagCLr=create_1d_double_array(Nr-1, "DiagCLr");
    DiagCz=create_1d_double_array(Nr, "DiagCz");
    DiagCUz=create_1d_double_array(Nr-1, "DiagCUz");
    DiagCLz=create_1d_double_array(Nr-1, "DiagCLz");
    DiagCrdr=create_1d_double_array(Nr, "DiagCrdr");
    DiagCUrdr=create_1d_double_array(Nr-1, "DiagCUrdr");
    DiagCLrdr=create_1d_double_array(Nr-1, "DiagCLrdr");
    DiagCzdz=create_1d_double_array(Nr, "DiagCzdz");
    DiagCUzdz=create_1d_double_array(Nr-1, "DiagCUzdz");
    DiagCLzdz=create_1d_double_array(Nr-1, "DiagCLzdz");
    
    
    **qC_0=0;
    ***qC=0;
    
    //Initialize the qs
    
    for (i=0;i<=(Nr-1);i++){
        for (j=1;j<=(Nz-1);j++){
            qC_0[i][j]=1.0;
            qC[i][j][0]=1.0;}}
    
    for (s=0;s<=(NC-1);s++){
        *bCr=0;
        /********************************scan over z***********************************************/
        for (i=0;i<=(Nr-1);i++){
            C_Matrix_z(i);
            if (i==0){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qC_0[i][j]+betaU*qC_0[i+1][j]+betaL*qC_0[i+1][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qC_0[i][j]+betaU*qC_0[i-1][j]+betaL*qC_0[i-1][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qC_0[i][j]+betaU*qC_0[i+1][j]+betaL*qC_0[i-1][j];}
            }
            DiagCrdr=DiagCr;
            DiagCUrdr=DiagCUr;
            DiagCLrdr=DiagCLr;
            
            TDMA(Nz,DiagCLrdr,DiagCrdr,DiagCUrdr,bCr);
            
            for (j=0; j<=Nz-1; j++){
                qC[i][j][s]=bCr[j];}
            
        }
        bCr=0;
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qC_0[i][j]=qC[i][j][s];}
        }
        *bCz=0;
        /***********************************scan over r***************************************************/
        for (j=0;j<=(Nz-1);j++){
            C_Matrix_r(j);
            if (j==0){
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow(delz,2)));
                    beta=(delt/(2.0*pow(delz,2)));
                    bCz[i]=gamma*qC_0[i][j]+beta*qC_0[j+1][i]+beta*qC_0[j+1][i];}
            }
            else if (j==(Nz-1)){
                for (i=0;i<=(Nz-1);i++){
                    gamma=1.0-(delt/(pow(delz,2)));
                    beta=(delt/(2.0*pow(delz,2)));
                    bCz[i]=gamma*qC_0[i][j]+beta*qC_0[j-1][i]+beta*qC_0[j-1][i];}
            }
            else {
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow(delz,2)));
                    beta=(delt/(2.0*pow(delz,2)));
                    bCr[j]=gamma*qC_0[i][j]+beta*qC_0[j+1][i]+beta*qC_0[j-1][i];}
            }
            DiagCzdz=DiagCz;
            DiagCUzdz=DiagCUz;
            DiagCLzdz=DiagCLz;
            
            TDMA(Nr,DiagCLzdz,DiagCzdz,DiagCUzdz,bCz);
            
            for (i=0; j<=Nr-1; i++){
                qC[i][j][s]=bCz[i];}
        }
        
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qC_0[i][j]=qC[i][j][s];}
        }
    }
    
    destroy_2d_double_array(qC_0);
    //destroy_3d_double_array(qC);
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
    return qC_forward();
}

/***********************************Define the complementary propagator***********************************/

int qdagC_forward(){
    //int errorflag;
    int s,i,j;
    //int errorflagA;
    double gamma,betaU,betaL,beta;
    
    qdagC=create_3d_double_array(Nr,Nz,Ns,"qC");
    qdagC_0=create_2d_double_array(Nr,Nz, "qC_0");
    bCr=create_1d_double_array(Nr, "bCr");
    bCz=create_1d_double_array(Nz, "bCz");
    DiagCr=create_1d_double_array(Nr, "DiagCr");
    DiagCUr=create_1d_double_array(int(Nr-1), "DiagCUr");
    DiagCLr=create_1d_double_array(Nr-1, "DiagCLr");
    DiagCz=create_1d_double_array(Nr, "DiagCz");
    DiagCUz=create_1d_double_array(Nr-1, "DiagCUz");
    DiagCLz=create_1d_double_array(Nr-1, "DiagCLz");
    DiagCrdr=create_1d_double_array(Nr, "DiagCrdr");
    DiagCUrdr=create_1d_double_array(Nr-1, "DiagCUrdr");
    DiagCLrdr=create_1d_double_array(Nr-1, "DiagCLrdr");
    DiagCzdz=create_1d_double_array(Nr, "DiagCzdz");
    DiagCUzdz=create_1d_double_array(Nr-1, "DiagCUzdz");
    DiagCLzdz=create_1d_double_array(Nr-1, "DiagCLzdz");
    
    **qdagC_0=0;
    ***qdagC=0;
    
    for (i=0;i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            qdagC[i][j][0]=qA[i][j][NA]*qB[i][j][NB];
            qdagC_0[i][j]=(qA[i][j][NA])*(qB[i][j][NB]);
        }
    }
    
    for (s=0;s<=NC-1;s++){
        *bCr=0;
        /********************************scan over z***********************************************/
        for (i=0;i<=(Nr-1);i++){
            C_Matrix_z(i);
            if (i==0){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qdagC_0[i][j]+betaU*qdagC_0[i+1][j]+betaL*qdagC_0[i+1][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qdagC_0[i][j]+betaU*qdagC_0[i-1][j]+betaL*qdagC_0[i-1][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wC[i][j]);
                    betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bCr[j]=gamma*qdagC_0[i][j]+betaU*qdagC_0[i+1][j]+betaL*qdagC_0[i-1][j];}
            }
            DiagCrdr=DiagCr;
            DiagCUrdr=DiagCUr;
            DiagCLrdr=DiagCLr;
            
            TDMA(Nz,DiagCLrdr,DiagCrdr,DiagCUrdr,bCr);
            
            for (j=0; j<=Nz-1; j++){
                qdagC[i][j][s]=bCr[j];}
            
        }
        *bCr=0;
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qdagC_0[i][j]=qdagC[i][j][s];}
        }
        
        
        
    }
    
    *bCz=0;
    /***********************************scan over r***************************************************/
    for (j=0;j<=(Nz-1);j++){
        C_Matrix_r(j);
        if (j==0){
            for (i=0;i<=(Nr-1);i++){
                gamma=1.0-(delt/(pow(delz,2)));
                beta=(delt/(2.0*pow(delz,2)));
                bCz[i]=gamma*qdagC_0[i][j]+beta*qdagC_0[j+1][i]+beta*qdagC_0[j+1][i];}
        }
        else if (j==(Nz-1)){
            for (i=0;i<=(Nz-1);i++){
                gamma=1.0-(delt/(pow(delz,2)));
                beta=(delt/(2.0*pow(delz,2)));
                bCz[i]=gamma*qdagC_0[i][j]+beta*qdagC_0[j-1][i]+beta*qdagC_0[j-1][i];}
        }
        else {
            for (i=0;i<=(Nr-1);i++){
                gamma=1.0-(delt/(pow(delz,2)));
                beta=(delt/(2.0*pow(delz,2)));
                bCr[j]=gamma*qdagC_0[i][j]+beta*qdagC_0[j+1][i]+beta*qdagC_0[j-1][i];}
        }
        DiagCzdz=DiagCz;
        DiagCUzdz=DiagCUz;
        DiagCLzdz=DiagCLz;
        
        TDMA(Nr,DiagCLzdz,DiagCzdz,DiagCUzdz,bCz);
        
        for (i=0; j<=Nr-1; i++){
            qdagC[i][j][s]=bCz[i];}
    }
    
    for (i=0;i<=(Nr-1);i++){
        for (j=0; j<=Nz-1; j++){
            qdagC_0[i][j]=qdagC[i][j][s];}
    }
    destroy_2d_double_array(qdagC_0);
    destroy_3d_double_array(qdagC);
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
    
    return qdagC_forward();
}













