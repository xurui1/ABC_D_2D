


void B_Matrix_r(int ii){
    
    
    int i;
    
    for (i=0;i<=Nr-1;i++){
        DiagBr[i]=0;}
    for (i=0;i<=Nr-2;i++){
        DiagBUr[i]=0;
        DiagBLr[i]=0;}
    
    
    for (i=0; i<=(Nr-1); i++){
        DiagBr[i]=1.0+(delt/(pow((double)delr,(int)2)))+((delt/2.0)*wB[i][ii]);}
    for (i=1; i<=(Nr-2);i++){
        DiagBUr[i]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(i-1))+r_0)*4.0*delr));}
    for (i=0; i<=(Nr-2);i++){
        DiagBLr[i]=-(delt/(2.0*pow((double) delr,(int)2)))+(delt/(((delr*(i+1-1))+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagBUr[0]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(1-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(1-1))+r_0)*4.0*delr));
    
    DiagBLr[Nr-1]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(Nr-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(Nr-1))+r_0)*4.0*delr));
    
}

void B_Matrix_z(int ii){
    
    int j;
    
    for (j=0;j<Nz;j++){
        DiagBz[j]=0;}
    for (j=0;j<(int)(Nz-1);j++){
        DiagBUz[j]=0;
        DiagBLz[j]=0;}
    
    
    for (j=0; j<=(Nz-1); j++){
        DiagBz[j]=1.0+delt/(pow((double)delz,(int)2));}
    for (j=0; j<=(Nz-2); j++){
        DiagBUz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    for (j=0; j<=(Nz-2); j++){
        DiagBLz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    
    DiagBUz[0]=2.0*DiagBUz[1];
    DiagBLz[Nz-1]=2.0*DiagBLz[Nz-1];
    
}


/**************************Finally time to define some propagators**********************/
double qB_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
    
    qB=create_3d_double_array(Nr,Nz,Ns,"qB");
    qB_0=create_2d_double_array(Nr,Nz, "qB_0");
    bBr=create_1d_double_array(Nr, "bBr");
    bBz=create_1d_double_array(Nz, "bBz");
    DiagBr=create_1d_double_array(Nr, "DiagBr");
    DiagBUr=create_1d_double_array(((int)Nr-1), "DiagBUr");
    DiagBLr=create_1d_double_array(((int)Nr-1), "DiagBLr");
    DiagBz=create_1d_double_array(Nr, "DiagBz");
    DiagBUz=create_1d_double_array(((int)Nr-1), "DiagBUz");
    DiagBLz=create_1d_double_array(((int)Nr-1), "DiagBLz");
    DiagBrdr=create_1d_double_array(Nr, "DiagBrdr");
    DiagBUrdr=create_1d_double_array(((int)Nr-1), "DiagBUrdr");
    DiagBLrdr=create_1d_double_array(((int)Nr-1), "DiagBLrdr");
    DiagBzdz=create_1d_double_array(Nr, "DiagBzdz");
    DiagBUzdz=create_1d_double_array(((int)Nr-1), "DiagBUzdz");
    DiagBLzdz=create_1d_double_array(((int)Nr-1), "DiagBLzdz");
    
    for (i=0; i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            qB_0[i][j]=0;
            for (s=0;s<=Ns-1;s++){
                qB[i][j][s]=0;}}}
    
    
    //Initialize the qs
    
    for (i=0;i<=(Nr-1);i++){
        for (j=1;j<=(Nz-1);j++){
            qB_0[i][j]=1.0;
            qB[i][j][0]=1.0;}}
    
    for (s=0;s<=(NB-1);s++){
        
        for (i=0;i<=Nr-1;i++){
            bBr[i]=0;}
        
        /********************************scan over z***********************************************/
        for (i=0;i<=(Nr-1);i++){
            B_Matrix_z(i);
            if (i==0){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qB_0[i][j]+betaU*qB_0[i+1][j]+betaL*qB_0[i+1][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qB_0[i][j]+betaU*qB_0[i-1][j]+betaL*qB_0[i-1][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qB_0[i][j]+betaU*qB_0[i+1][j]+betaL*qB_0[i-1][j];}
            }
            DiagBrdr[i]=DiagBr[i];
            DiagBUrdr[i]=DiagBUr[i];
            DiagBLrdr[i]=DiagBLr[i];
            
            TDMA(Nz,DiagBLrdr,DiagBrdr,DiagBUrdr,bBr);
            
            for (j=0; j<=Nz-1; j++){
                qB[i][j][s]=bBr[j];}
            
        }
        for (j=0;j<=Nr-1;j++){
            bBr[j]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qB_0[i][j]=qB[i][j][s];}
        }
        
        for (j=0;j<=Nz-1;j++){
            bBz[j]=0;
        }
        /***********************************scan over r***************************************************/
        for (j=0;j<=(Nz-1);j++){
            B_Matrix_r(j);
            if (j==0){
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bBz[i]=gamma*qB_0[i][j]+beta*qB_0[j+1][i]+beta*qB_0[j+1][i];}
            }
            else if (j==(Nz-1)){
                for (i=0;i<=(Nz-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bBz[i]=gamma*qB_0[i][j]+beta*qB_0[j-1][i]+beta*qB_0[j-1][i];}
            }
            else {
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bBz[i]=gamma*qB_0[i][j]+beta*qB_0[j+1][i]+beta*qB_0[j-1][i];}
            }
            DiagBzdz[j]=DiagBz[j];
            DiagBUzdz[j]=DiagBUz[j];
            DiagBLzdz[j]=DiagBLz[j];
            
            TDMA(Nr,DiagBLzdz,DiagBzdz,DiagBUzdz,bBz);
            
            for (i=0; i<=Nr-1; i++){
                qB[i][j][s]=bBz[i];}
        }
        
        for (i=0;i<=Nr-1;i++){
            bBz[i]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qB_0[i][j]=qB[i][j][s];}
        }
    }
    
    destroy_2d_double_array(qB_0);
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
    
    return ***qB;
}

/***********************************Define the complementary propagator***********************************/

double qdagB_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
    
    qdagB=create_3d_double_array(Nr,Nz,Ns,"qdagB");
    qdagB_0=create_2d_double_array(Nr,Nz, "qdagB_0");
    bBr=create_1d_double_array(Nr, "bBr");
    bBz=create_1d_double_array(Nz, "bBz");
    DiagBr=create_1d_double_array(Nr, "DiagBr");
    DiagBUr=create_1d_double_array(((int)Nr-1), "DiagBUr");
    DiagBLr=create_1d_double_array(((int)Nr-1), "DiagBLr");
    DiagBz=create_1d_double_array(Nr, "DiagBz");
    DiagBUz=create_1d_double_array(((int)Nr-1), "DiagBUz");
    DiagBLz=create_1d_double_array(((int)Nr-1), "DiagBLz");
    DiagBrdr=create_1d_double_array(Nr, "DiagBrdr");
    DiagBUrdr=create_1d_double_array(((int)Nr-1), "DiagBUrdr");
    DiagBLrdr=create_1d_double_array(((int)Nr-1), "DiagBLrdr");
    DiagBzdz=create_1d_double_array(Nr, "DiagBzdz");
    DiagBUzdz=create_1d_double_array(((int)Nr-1), "DiagBUzdz");
    DiagBLzdz=create_1d_double_array(((int)Nr-1), "DiagBLzdz");
    
    
    for (i=0; i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            qdagB_0[i][j]=0;
            for (s=0;s<=Ns-1;s++){
                qdagB[i][j][s]=0;
            }
        }
    }
    
    for (i=0;i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            qdagB[i][j][0]=qA[i][j][NA]*qC[i][j][NC];
            qdagB_0[i][j]=(qA[i][j][NA])*(qC[i][j][NC]);
        }
    }
    
    for (s=0;s<=NB-1;s++){
        
        for (i=0;i<=Nr-1;i++){
            bBr[i]=0;
        }
        
        /********************************scan over z***********************************************/
        for (i=0;i<=(Nr-1);i++){
            B_Matrix_z(i);
            if (i==0){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qdagB_0[i][j]+betaU*qdagB_0[i+1][j]+betaL*qdagB_0[i+1][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qdagB_0[i][j]+betaU*qdagB_0[i-1][j]+betaL*qdagB_0[i-1][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wB[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bBr[j]=gamma*qdagB_0[i][j]+betaU*qdagB_0[i+1][j]+betaL*qdagB_0[i-1][j];}
            }
            DiagBrdr[i]=DiagBr[i];
            DiagBUrdr[i]=DiagBUr[i];
            DiagBLrdr[i]=DiagBLr[i];
            
            TDMA(Nz,DiagBLrdr,DiagBrdr,DiagBUrdr,bBr);
            
            for (j=0; j<=Nz-1; j++){
                qdagB[i][j][s]=bBr[j];}
            
        }
        
        for (j=0;j<=Nz-1;j++){
            bBr[j]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qdagB_0[i][j]=qdagB[i][j][s];}
        }
        
        
    }
    
    for (j=0;j<=Nz-1;j++){
        bBz[j]=0;
    }
    /***********************************scan over r***************************************************/
    for (j=0;j<=(Nz-1);j++){
        B_Matrix_r(j);
        if (j==0){
            for (i=0;i<=(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bBz[i]=gamma*qdagB_0[i][j]+beta*qdagB_0[j+1][i]+beta*qdagB_0[j+1][i];}
        }
        else if (j==(Nz-1)){
            for (i=0;i<=(Nz-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bBz[i]=gamma*qdagB_0[i][j]+beta*qdagB_0[j-1][i]+beta*qdagB_0[j-1][i];}
        }
        else {
            for (i=0;i<=(Nr-1);i++){
                gamma=1.0-(delt/(pow((double)delz,(int)2)));
                beta=(delt/(2.0*pow((double)delz,(int)2)));
                bBz[i]=gamma*qdagB_0[i][j]+beta*qdagB_0[j+1][i]+beta*qdagB_0[j-1][i];}
        }
        DiagBzdz[j]=DiagBz[j];
        DiagBUzdz[j]=DiagBUz[j];
        DiagBLzdz[j]=DiagBLz[j];
        
        TDMA(Nr,DiagBLzdz,DiagBzdz,DiagBUzdz,bBz);
        
        for (i=0; i<=Nr-1; i++){
            qdagB[i][j][s]=bBz[i];}
    }
    
    for (i=0;i<=Nr-1;i++){
        bBz[i]=0;
    }
    
    for (i=0;i<=(Nr-1);i++){
        for (j=0; j<=Nz-1; j++){
            qdagB_0[i][j]=qdagB[i][j][s];}
    }
    
    destroy_2d_double_array(qdagB_0);
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
    
    return ***qdagB;
}









