

void D_Matrix_r(int ii){
    
    int i;
    
    for (i=0;i<=int(Nr-1);i++){
        DiagDz[i]=0;}
    for (i=0;i<=int(Nr-2);i++){
        DiagDUz[i]=0;
        DiagDLz[i]=0;}
    
    
    for (i=0; i<=int(Nr-1); i++){
        DiagDz[i]=1.0+(delt/(pow((double)delr,(int)2)))+((delt/2.0)*wD[i][ii]);}
    for (i=1; i<=int(Nr-2);i++){
        DiagDUz[i]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(i-1))+r_0)*4.0*delr));}
    for (i=0; i<=int(Nr-3);i++){
        DiagDLz[i]=-(delt/(2.0*pow((double) delr,(int)2)))+(delt/(((delr*(i+1-1))+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagDUz[0]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(1-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(1-1))+r_0)*4.0*delr));
    
    DiagDLz[int(Nr-2)]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(int(Nr-1)))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(int(Nr-1)))+r_0)*4.0*delr));
    
}

void D_Matrix_z(int ii){
    
    int j;
    
    for (j=0;j<=int(Nz-1);j++){
        DiagDr[j]=0;}
    for (j=0;j<=int(Nz-2);j++){
        DiagDUr[j]=0;
        DiagDLr[j]=0;}
    
    
    for (j=0; j<=int(Nz-1); j++){
        DiagDr[j]=1.0+delt/(pow((double)delz,(int)2));}
    for (j=1; j<=int(Nz-2); j++){
        DiagDUr[j]=-delt/(2.0*pow((double)delz,(int)2));}
    for (j=0; j<=int(Nz-3); j++){
        DiagDLr[j]=-delt/(2.0*pow((double)delz,(int)2));}
    
    DiagDUr[0]=2.0*DiagDUr[0];
    DiagDLr[int(Nz-2)]=2.0*DiagDLr[int(Nz-2)];
    
}


/**************************Finally time to define some propagators**********************/
double qD_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
    
    qD=create_3d_double_array(Nr,Nz,Ns,"qD");
    qD_0=create_2d_double_array(Nr,Nz, "qD_0");
    bDr=create_1d_double_array(Nz, "bDr");
    bDz=create_1d_double_array(Nr, "bDz");
    DiagDr=create_1d_double_array(Nz, "DiagDr");
    DiagDUr=create_1d_double_array(((int)Nz-1), "DiagDUr");
    DiagDLr=create_1d_double_array(((int)Nz-1), "DiagDLr");
    DiagDz=create_1d_double_array(Nr, "DiagDz");
    DiagDUz=create_1d_double_array(((int)Nr-1), "DiagDUz");
    DiagDLz=create_1d_double_array(((int)Nr-1), "DiagDLz");
    DiagDrdr=create_1d_double_array(Nz, "DiagDrdr");
    DiagDUrdr=create_1d_double_array(((int)Nz-1), "DiagDUrdr");
    DiagDLrdr=create_1d_double_array(((int)Nz-1), "DiagDLrdr");
    DiagDzdz=create_1d_double_array(Nr, "DiagDzdz");
    DiagDUzdz=create_1d_double_array(((int)Nr-1), "DiagDUzdz");
    DiagDLzdz=create_1d_double_array(((int)Nr-1), "DiagDLzdz");
    
    for (i=0; i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            qD_0[i][j]=0;
            for (s=0;s<=int(Ns-1);s++){
                qD[i][j][s]=0;}}}
    
    //Initialize the qs
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=1;j<=int(Nz-1);j++){
            qD_0[i][j]=1.0;
            qD[i][j][0]=1.0;}}
    
    for (s=0;s<=int(ND-1);s++){
        
        for (j=0;j<=int(Nz-1);j++){
            bDr[j]=0;}
        
        /********************************scan over z***********************************************/
        for (i=0;i<=int(Nr-1);i++){
            D_Matrix_z(i);
            if (i==0){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wD[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bDr[j]=gamma*qD_0[i][j]+betaU*qD_0[int(i+1)][j]+betaL*qD_0[int(i+1)][j];}
            }
            else if (i==int(Nr-1)){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wD[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bDr[j]=gamma*qD_0[i][j]+betaU*qD_0[int(i-1)][j]+betaL*qD_0[int(i-1)][j];}
            }
            else {
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wD[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bDr[j]=gamma*qD_0[i][j]+betaU*qD_0[int(i+1)][j]+betaL*qD_0[int(i-1)][j];}
            }
            for (j=0;j<=int(Nz-1);j++){
                DiagDrdr[j]=DiagDr[j];}
            for (j=0;j<=int(Nz-2);j++){
                DiagDUrdr[j]=DiagDUr[j];
                DiagDLrdr[j]=DiagDLr[j];}
            
            TDMA(Nz,DiagDLrdr,DiagDrdr,DiagDUrdr,bDr);
            
            for (j=0; j<=int(Nz-1); j++){
                qD[i][j][s]=bDr[j];}
        }
        
        for (j=0;j<=int(Nz-1);j++){
            bDr[j]=0;}
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qD_0[i][j]=qD[i][j][s];}}
        
        for (i=0;i<=int(Nr-1);i++){
            bDz[i]=0;}
        
        /***********************************scan over r***************************************************/
        for (j=0;j<=int(Nz-1);j++){
            D_Matrix_r(j);
            if (j==0){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bDz[i]=gamma*qD_0[i][j]+beta*qD_0[i][int(j+1)]+beta*qD_0[i][int(j+1)];}
            }
            else if (j==int(Nz-1)){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bDz[i]=gamma*qD_0[i][j]+beta*qD_0[i][int(j-1)]+beta*qD_0[i][int(j-1)];}
            }
            else {
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bDz[i]=gamma*qD_0[i][j]+beta*qD_0[i][int(j+1)]+beta*qD_0[i][int(j-1)];}
            }
            for (i=0;i<=int(Nr-1);i++){
                DiagDzdz[i]=DiagDz[i];}
            for (i=0;i<=int(Nr-2);i++){
                DiagDUzdz[i]=DiagDUz[i];
                DiagDLzdz[i]=DiagDLz[i];}
            
            TDMA(Nr,DiagDLzdz,DiagDzdz,DiagDUzdz,bDz);
            
            for (i=0; i<=int(Nr-1); i++){
                qD[i][j][s]=bDz[i];}
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bDz[i]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qD_0[i][j]=qD[i][j][s];}
        }
    }
    
    return ***qD;
    
}


