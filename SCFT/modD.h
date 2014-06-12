

void D_Matrix_r(int ii){
    
    int i;
    
    for (i=0;i<=Nr-1;i++){
        DiagDr[i]=0;}
    for (i=0;i<=Nr-2;i++){
        DiagDUr[i]=0;
        DiagDLr[i]=0;}
    
    
    for (i=0; i<=(Nr-1); i++){
        DiagDr[i]=1.0+(delt/(pow((double)delr,(int)2)))+((delt/2.0)*wD[i][ii]);}
    for (i=1; i<=(Nr-2);i++){
        DiagDUr[i]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(i-1))+r_0)*4.0*delr));}
    for (i=0; i<=(Nr-2);i++){
        DiagDLr[i]=-(delt/(2.0*pow((double) delr,(int)2)))+(delt/(((delr*(i+1-1))+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagDUr[0]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(1-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(1-1))+r_0)*4.0*delr));
    
    DiagDLr[Nr-1]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(Nr-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(Nr-1))+r_0)*4.0*delr));
    
}

void D_Matrix_z(int ii){
    
    int j;
    
    for (j=0;j<Nz;j++){
        DiagDz[j]=0;}
    for (j=0;j<(int)(Nz-1);j++){
        DiagDUz[j]=0;
        DiagDLz[j]=0;}
    
    
    for (j=0; j<=(Nz-1); j++){
        DiagDz[j]=1.0+delt/(pow((double)delz,(int)2));}
    for (j=0; j<=(Nz-2); j++){
        DiagDUz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    for (j=0; j<=(Nz-2); j++){
        DiagDLz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    
    DiagDUz[0]=2.0*DiagDUz[1];
    DiagDLz[Nz-1]=2.0*DiagDLz[Nz-1];
    
}


/**************************Finally time to define some propagators**********************/
double qD_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
    
    qD=create_3d_double_array(Nr,Nz,Ns,"qD");
    qD_0=create_2d_double_array(Nr,Nz, "qD_0");
    bDr=create_1d_double_array(Nr, "bDr");
    bDz=create_1d_double_array(Nz, "bDz");
    DiagDr=create_1d_double_array(Nr, "DiagDr");
    DiagDUr=create_1d_double_array(((int)Nr-1), "DiagDUr");
    DiagDLr=create_1d_double_array(((int)Nr-1), "DiagDLr");
    DiagDz=create_1d_double_array(Nr, "DiagDz");
    DiagDUz=create_1d_double_array(((int)Nr-1), "DiagDUz");
    DiagDLz=create_1d_double_array(((int)Nr-1), "DiagDLz");
    DiagDrdr=create_1d_double_array(Nr, "DiagDrdr");
    DiagDUrdr=create_1d_double_array(((int)Nr-1), "DiagDUrdr");
    DiagDLrdr=create_1d_double_array(((int)Nr-1), "DiagDLrdr");
    DiagDzdz=create_1d_double_array(Nr, "DiagDzdz");
    DiagDUzdz=create_1d_double_array(((int)Nr-1), "DiagDUzdz");
    DiagDLzdz=create_1d_double_array(((int)Nr-1), "DiagDLzdz");
    
    for (i=0; i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            qD_0[i][j]=0;
            for (s=0;s<=Ns-1;s++){
                qD[i][j][s]=0;}}}
    
    //Initialize the qs
    
    for (i=0;i<=(Nr-1);i++){
        for (j=1;j<=(Nz-1);j++){
            qD_0[i][j]=1.0;
            qD[i][j][0]=1.0;}}
    
    for (s=0;s<=(ND-1);s++){
        
        for (i=0;i<=Nr-1;i++){
            bDr[i]=0;}
        
        /********************************scan over z***********************************************/
        for (i=0;i<=(Nr-1);i++){
            D_Matrix_z(i);
            if (i==0){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wD[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bDr[j]=gamma*qD_0[i][j]+betaU*qD_0[i+1][j]+betaL*qD_0[i+1][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wD[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bDr[j]=gamma*qD_0[i][j]+betaU*qD_0[i-1][j]+betaL*qD_0[i-1][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wD[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bDr[j]=gamma*qD_0[i][j]+betaU*qD_0[i+1][j]+betaL*qD_0[i-1][j];}
            }
            DiagDrdr[i]=DiagDr[i];
            DiagDUrdr[i]=DiagDUr[i];
            DiagDLrdr[i]=DiagDLr[i];
            
            TDMA(Nz,DiagDLrdr,DiagDrdr,DiagDUrdr,bDr);
            
            for (j=0; j<=Nz-1; j++){
                qD[i][j][s]=bDr[j];}
        }
        
        for (j=0;j<=Nz-1;j++){
            bDr[j]=0;}
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qD_0[i][j]=qD[i][j][s];}}
        
        for (j=0;j<=Nz-1;j++){
            bDz[j]=0;}
        
        /***********************************scan over r***************************************************/
        for (j=0;j<=(Nz-1);j++){
            D_Matrix_r(j);
            if (j==0){
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bDz[i]=gamma*qD_0[i][j]+beta*qD_0[j+1][i]+beta*qD_0[j+1][i];}
            }
            else if (j==(Nz-1)){
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bDz[i]=gamma*qD_0[i][j]+beta*qD_0[j-1][i]+beta*qD_0[j-1][i];}
            }
            else {
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bDz[i]=gamma*qD_0[i][j]+beta*qD_0[j+1][i]+beta*qD_0[j-1][i];}
            }
            DiagDzdz[j]=DiagDz[j];
            DiagDUzdz[j]=DiagDUz[j];
            DiagDLzdz[j]=DiagDLz[j];
            
            TDMA(Nr,DiagDLzdz,DiagDzdz,DiagDUzdz,bDz);
            
            for (i=0; i<=Nr-1; i++){
                qD[i][j][s]=bDz[i];}
        }
        
        for (i=0;i<=Nr-1;i++){
            bDz[i]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qD_0[i][j]=qD[i][j][s];}
        }
    }
    
    destroy_2d_double_array(qD_0);
    destroy_1d_double_array(DiagDLr);
    destroy_1d_double_array(DiagDUr);
    destroy_1d_double_array(DiagDr);
    destroy_1d_double_array(DiagDLrdr);
    destroy_1d_double_array(DiagDUrdr);
    destroy_1d_double_array(DiagDrdr);
    destroy_1d_double_array(DiagDLz);
    destroy_1d_double_array(DiagDUz);
    destroy_1d_double_array(DiagDz);
    destroy_1d_double_array(DiagDLzdz);
    destroy_1d_double_array(DiagDUzdz);
    destroy_1d_double_array(DiagDzdz);
    destroy_1d_double_array(bDr);
    destroy_1d_double_array(bDz);
    
    return ***qD;
    
}


