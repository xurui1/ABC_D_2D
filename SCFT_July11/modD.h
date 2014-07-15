

void D_Matrix_z(int ii,double *DiagDr, double *DiagDLr, double *DiagDUr){
    
    int i;
    
    for (i=0;i<=int(Nr-1);i++){
        DiagDr[i]=0;}
    for (i=0;i<=int(Nr-2);i++){
        DiagDUr[i]=0;
        DiagDLr[i]=0;}
    
    
    for (i=0; i<=int(Nr-1); i++){
        DiagDr[i]=1.0+(delt/(pow((double)delr,(int)2)))+((delt/4.0)*wD[i][ii]);}
    for (i=1; i<=int(Nr-2);i++){
        DiagDUr[i]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*i)+r_0)*4.0*delr));}
    for (i=0; i<=int(Nr-3);i++){
        DiagDLr[i]=-(delt/(2.0*pow((double) delr,(int)2)))+(delt/(((delr*i)+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagDUr[0]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(0))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(0))+r_0)*4.0*delr));
    
    DiagDLr[int(Nr-2)]=-(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((delr*(int(Nr-2)))+r_0)*4.0*delr))
    -(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((delr*(int(Nr-2)))+r_0)*4.0*delr));
    
}

void D_Matrix_r(int ii,double *DiagDz,double *DiagDLz, double *DiagDUz){
    
    int j;
    
    for (j=0;j<=int(Nz-1);j++){
        DiagDz[j]=0;}
    for (j=0;j<=int(Nz-2);j++){
        DiagDUz[j]=0;
        DiagDLz[j]=0;}
    
    
    for (j=0; j<=int(Nz-1); j++){
        DiagDz[j]=1.0+delt/(pow((double)delz,(int)2));}
    for (j=1; j<=int(Nz-2); j++){
        DiagDUz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    for (j=0; j<=int(Nz-3); j++){
        DiagDLz[j]=-delt/(2.0*pow((double)delz,(int)2));}
    
    DiagDUz[0]=2.0*DiagDUz[0];
    DiagDLz[int(Nz-2)]=2.0*DiagDLz[int(Nz-2)];
    
}


/**************************Finally time to define some propagators**********************/
void qD_forward(){
    
    int s,i,j;
    double gamma,betaU,betaL,beta;
    double *DiagDr;
    double *DiagDUr;
    double *DiagDLr;
    double *DiagDrdr;
    double *DiagDUrdr;
    double *DiagDLrdr;
    double *DiagDz;
    double *DiagDUz;
    double *DiagDLz;
    double *DiagDzdz;
    double *DiagDUzdz;
    double *DiagDLzdz;
    
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
            for (s=0;s<=Ns;s++){
                qD[i][j][s]=0;}}}
    
    //Initialize the qs
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=1;j<=int(Nz-1);j++){
            qD_0[i][j]=1.0;
            qD[i][j][0]=1.0;}}
    
    for (s=0;s<=ND;s++){
        
        for (j=0;j<=int(Nz-1);j++){
            bDz[j]=0;}
        
        /********************************scan over z***********************************************/
        for (i=0;i<=int(Nr-1);i++){
            D_Matrix_r(i,DiagDz,DiagDLz,DiagDUz);
            if (i==0){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wD[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bDz[j]=gamma*qD_0[i][j]+betaU*qD_0[int(i+1)][j]+betaL*qD_0[int(i+1)][j];}
            }
            else if (i==int(Nr-1)){
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wD[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bDz[j]=gamma*qD_0[i][j]+betaU*qD_0[int(i-1)][j]+betaL*qD_0[int(i-1)][j];}
            }
            else {
                for (j=0;j<=int(Nz-1);j++){
                    gamma=1.0-(delt/(pow((double)delr,(int)2)))-((delt/2.0)*wD[i][j]);
                    betaL=(delt/(2.0*pow((double)delr,(int)2)))-(delt/(((i*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow((double)delr,(int)2)))+(delt/(((i*delr)+(r_0))*4.0*delr));
                    bDz[j]=gamma*qD_0[i][j]+betaU*qD_0[int(i+1)][j]+betaL*qD_0[int(i-1)][j];}
            }
            for (j=0;j<=int(Nz-1);j++){
                DiagDzdz[j]=DiagDz[j];}
            for (j=0;j<=int(Nz-2);j++){
                DiagDUzdz[j]=DiagDUz[j];
                DiagDLzdz[j]=DiagDLz[j];}
            
            TDMA(Nz,DiagDLzdz,DiagDzdz,DiagDUzdz,bDz);
            
            for (j=0; j<=int(Nz-1); j++){
                qD[i][j][s]=bDz[j];}
        }
        
        for (j=0;j<=int(Nz-1);j++){
            bDz[j]=0;}
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qD_0[i][j]=qD[i][j][s];}}
        
        for (i=0;i<=int(Nr-1);i++){
            bDr[i]=0;}
        
        /***********************************scan over r***************************************************/
        for (j=0;j<=int(Nz-1);j++){
            D_Matrix_z(j,DiagDr,DiagDLr,DiagDUr);
            if (j==0){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bDr[i]=gamma*qD_0[i][j]+beta*qD_0[i][int(j+1)]+beta*qD_0[i][int(j+1)];}
            }
            else if (j==int(Nz-1)){
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bDr[i]=gamma*qD_0[i][j]+beta*qD_0[i][int(j-1)]+beta*qD_0[i][int(j-1)];}
            }
            else {
                for (i=0;i<=int(Nr-1);i++){
                    gamma=1.0-(delt/(pow((double)delz,(int)2)));
                    beta=(delt/(2.0*pow((double)delz,(int)2)));
                    bDr[i]=gamma*qD_0[i][j]+beta*qD_0[i][int(j+1)]+beta*qD_0[i][int(j-1)];}
            }
            for (i=0;i<=int(Nr-1);i++){
                DiagDrdr[i]=DiagDr[i];}
            for (i=0;i<=int(Nr-2);i++){
                DiagDUrdr[i]=DiagDUr[i];
                DiagDLrdr[i]=DiagDLr[i];}
            
            TDMA(Nr,DiagDLrdr,DiagDrdr,DiagDUrdr,bDr);
            
            for (i=0; i<=int(Nr-1); i++){
                qD[i][j][s]=bDr[i];}
        }
        
        for (i=0;i<=int(Nr-1);i++){
            bDr[i]=0;
        }
        
        for (i=0;i<=int(Nr-1);i++){
            for (j=0; j<=int(Nz-1); j++){
                qD_0[i][j]=qD[i][j][s];}
        }
    }
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
    
    
}


