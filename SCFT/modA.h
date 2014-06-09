

int A_Matrix_r(int ii){


    int i;
    //double alphaA,betaAL,betaAU;
    
    DiagAr=create_1d_double_array(Nr, "DiagAr");
    DiagAUr=create_1d_double_array(int(Nr-1), "DiagAUr");
    DiagALr=create_1d_double_array(Nr-1, "DiagALr");
    
    for (i=0;i<=Nr-1;i++){
        DiagAr[i]=0;}
    for (i=0;i<=Nr-2;i++){
        DiagAUr[i]=0;
        DiagALr[i]=0;}
   
    
    for (i=0; i<=(Nr-1); i++){
        DiagAr[i]=1.0+(delt/(pow(delr,2))+((delt/2.0)*wA[i][ii]));}
    for (i=1; i<=(Nr-2);i++){
        DiagAUr[i]=-(delt/(2.0*pow(delr,2)))-(delt/(((delr*(i-1))+r_0)*4.0*delr));}
    for (i=0; i<=(Nr-2);i++){
        DiagALr[i]=-(delt/(2.0*pow(delr,2)))+(delt/(((delr*(i+1-1))+r_0)*4.0*delr));}
    
    //Corners are different due to the 0 order boundary condition used
    
    DiagAUr[0]=-(delt/(2.0*pow(delr,2)))-(delt/(((delr*(1-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow(delr,2)))+(delt/(((delr*(1-1))+r_0)*4.0*delr));
    
    DiagALr[Nr-1]=-(delt/(2.0*pow(delr,2)))-(delt/(((delr*(Nr-1))+r_0)*4.0*delr))
    -(delt/(2.0*pow(delr,2)))+(delt/(((delr*(Nr-1))+r_0)*4.0*delr));
    
    destroy_1d_double_array(DiagALr);
    destroy_1d_double_array(DiagAUr);
    destroy_1d_double_array(DiagAr);
    
    return A_Matrix_r(ii);
}

int A_Matrix_z(int ii){

    
    int j;
    //double alphaA,betaAL,betaAU;
    
  DiagAz=create_1d_double_array(Nr, "DiagAz");
  DiagAUz=create_1d_double_array(Nr-1, "DiagAUz");
  DiagALz=create_1d_double_array(Nr-1, "DiagALz");
    
    for (j=0;j<=Nz-1;j++){
       DiagAz[j]=0;}
    for (j=0;j<=Nz-2;j++){
        DiagAUz[j]=0;
        DiagALz[j]=0;}
   
    
    for (j=0; j<=(Nz-1); j++){
        DiagAz[j]=1.0+delt/(pow(delz,2));}
    for (j=0; j<=(Nz-2); j++){
        DiagAUz[j]=-delt/(2.0*pow(delz,2));}
    for (j=0; j<=(Nz-2); j++){
        DiagALz[j]=-delt/(2.0*pow(delz,2));}
    
    DiagAUz[0]=2.0*DiagAUz[1];
    DiagALz[Nz-1]=2.0*DiagALz[Nz-1];
    
    return A_Matrix_z(ii);
}


/**************************Finally time to define some propagators**********************/
int qA_forward(){
    //int errorflag;
    int s,i,j;
    //int errorflagA;
    double gamma,betaU,betaL,beta;
    
    qA=create_3d_double_array(Nr,Nz,Ns,"qA");
    qA_0=create_2d_double_array(Nr,Nz, "qA_0");
    bAr=create_1d_double_array(Nr, "bAr");
    bAz=create_1d_double_array(Nz, "bAz");
    DiagAr=create_1d_double_array(Nr, "DiagAr");
    DiagAUr=create_1d_double_array(int(Nr-1), "DiagAUr");
    DiagALr=create_1d_double_array(Nr-1, "DiagALr");
    DiagAz=create_1d_double_array(Nr, "DiagAz");
    DiagAUz=create_1d_double_array(Nr-1, "DiagAUz");
    DiagALz=create_1d_double_array(Nr-1, "DiagALz");
    DiagArdr=create_1d_double_array(Nr, "DiagArdr");
    DiagAUrdr=create_1d_double_array(Nr-1, "DiagAUrdr");
    DiagALrdr=create_1d_double_array(Nr-1, "DiagALrdr");
    DiagAzdz=create_1d_double_array(Nr, "DiagAzdz");
    DiagAUzdz=create_1d_double_array(Nr-1, "DiagAUzdz");
    DiagALzdz=create_1d_double_array(Nr-1, "DiagALzdz");
    
    for (i=0; i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            qA_0[i][j]=0;
            for (s=0;s<=Ns-1;s++){
                 qA[i][j][s]=0;}}}
   
    
//Initialize the qs
    
    for (i=0;i<=(Nr-1);i++){
        for (j=1;j<=(Nz-1);j++){
            qA_0[i][j]=1.0;
            qA[i][j][0]=1.0;}}
        
    for (s=0;s<=(NA-1);s++){
        
        for (i=0;i<=Nr-1;i++){
            bAr[i]=0;}
        
/********************************scan over z***********************************************/
        for (i=0;i<=(Nr-1);i++){
            A_Matrix_z(i);
            if (i==0){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wA[i][j]);
                    betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bAr[j]=gamma*qA_0[i][j]+betaU*qA_0[i+1][j]+betaL*qA_0[i+1][j];}
            }
            else if (i==(Nr-1)){
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wA[i][j]);
                    betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bAr[j]=gamma*qA_0[i][j]+betaU*qA_0[i-1][j]+betaL*qA_0[i-1][j];}
            }
            else {
                for (j=0;j<=(Nz-1);j++){
                    gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wA[i][j]);
                    betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                    bAr[j]=gamma*qA_0[i][j]+betaU*qA_0[i+1][j]+betaL*qA_0[i-1][j];}
            }
            DiagArdr=DiagAr;
            DiagAUrdr=DiagAUr;
            DiagALrdr=DiagALr;
            
            TDMA(Nz,DiagALrdr,DiagArdr,DiagAUrdr,bAr);
            
            for (j=0; j<=Nz-1; j++){
                qA[i][j][s]=bAr[j];}
            
            }
        for (i=0;i<=Nr-1;i++){
            bAr[i]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qA_0[i][j]=qA[i][j][s];}
            }
        
        for (j=0;j<=Nz-1;j++){
            bAz[j]=0;
        }
/***********************************scan over r***************************************************/
        for (j=0;j<=(Nz-1);j++){
            A_Matrix_r(j);
            if (j==0){
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow(delz,2)));
                    beta=(delt/(2.0*pow(delz,2)));
                    bAz[i]=gamma*qA_0[i][j]+beta*qA_0[j+1][i]+beta*qA_0[j+1][i];}
            }
            else if (j==(Nz-1)){
                for (i=0;i<=(Nz-1);i++){
                    gamma=1.0-(delt/(pow(delz,2)));
                    beta=(delt/(2.0*pow(delz,2)));
                    bAz[i]=gamma*qA_0[i][j]+beta*qA_0[j-1][i]+beta*qA_0[j-1][i];}
            }
            else {
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow(delz,2)));
                    beta=(delt/(2.0*pow(delz,2)));
                    bAr[j]=gamma*qA_0[i][j]+beta*qA_0[j+1][i]+beta*qA_0[j-1][i];}
            }
            DiagAzdz=DiagAz;
            DiagAUzdz=DiagAUz;
            DiagALzdz=DiagALz;
            
            TDMA(Nr,DiagALzdz,DiagAzdz,DiagAUzdz,bAz);
            
            for (i=0; j<=Nr-1; i++){
                qA[i][j][s]=bAz[i];}
        }
        
        for (i=0;i<=Nr-1;i++){
            bAr[i]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qA_0[i][j]=qA[i][j][s];}
        }
    }
    
    destroy_2d_double_array(qA_0);
    destroy_3d_double_array(qA);
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
    
    return qA_forward();
}
    
/***********************************Define the complementary propagator***********************************/

    int qdagA_forward(){
        //int errorflag;
        int s,i,j;
        //int errorflagA;
        double gamma,betaU,betaL,beta;
        
        qdagA=create_3d_double_array(Nr,Nz,Ns,"qdagA");
        qdagA_0=create_2d_double_array(Nr,Nz, "qdagA_0");
        bAr=create_1d_double_array(Nr, "bAr");
        bAz=create_1d_double_array(Nz, "bAz");
        DiagAr=create_1d_double_array(Nr, "DiagAr");
        DiagAUr=create_1d_double_array(int(Nr-1), "DiagAUr");
        DiagALr=create_1d_double_array(Nr-1, "DiagALr");
        DiagAz=create_1d_double_array(Nr, "DiagAz");
        DiagAUz=create_1d_double_array(Nr-1, "DiagAUz");
        DiagALz=create_1d_double_array(Nr-1, "DiagALz");
        DiagArdr=create_1d_double_array(Nr, "DiagArdr");
        DiagAUrdr=create_1d_double_array(Nr-1, "DiagAUrdr");
        DiagALrdr=create_1d_double_array(Nr-1, "DiagALrdr");
        DiagAzdz=create_1d_double_array(Nr, "DiagAzdz");
        DiagAUzdz=create_1d_double_array(Nr-1, "DiagAUzdz");
        DiagALzdz=create_1d_double_array(Nr-1, "DiagALzdz");

        
        for (i=0; i<=Nr-1;i++){
            for (j=0;j<=Nz-1;j++){
                qdagA_0[i][j]=0;
                for (s=0;s<=Ns-1;s++){
                    qdagA[i][j][s]=0;
                }
            }
        }
        
        for (i=0;i<=Nr-1;i++){
            for (j=0;j<=Nz-1;j++){
                qdagA[i][j][0]=qB[i][j][NB]*qC[i][j][NC];
                qdagA_0[i][j]=(qB[i][j][NB])*(qC[i][j][NC]);
            }
        }
        
        for (s=0;s<=NA-1;s++){
            
            for (i=0;i<=Nr-1;i++){
                bAr[i]=0;
            }
            
            /********************************scan over z***********************************************/
            for (i=0;i<=(Nr-1);i++){
                A_Matrix_z(i);
                if (i==0){
                    for (j=0;j<=(Nz-1);j++){
                        gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wA[i][j]);
                        betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        bAr[j]=gamma*qdagA_0[i][j]+betaU*qdagA_0[i+1][j]+betaL*qdagA_0[i+1][j];}
                }
                else if (i==(Nr-1)){
                    for (j=0;j<=(Nz-1);j++){
                        gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wA[i][j]);
                        betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        bAr[j]=gamma*qdagA_0[i][j]+betaU*qdagA_0[i-1][j]+betaL*qdagA_0[i-1][j];}
                }
                else {
                    for (j=0;j<=(Nz-1);j++){
                        gamma=1.0-(delt/(pow(delr,2)))-((delt/2.0)*wA[i][j]);
                        betaL=(delt/(2.0*pow(delr,2)))-(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        betaU=(delt/(2.0*pow(delr,2)))+(delt/((((i-1)*delr)+(r_0))*4.0*delr));
                        bAr[j]=gamma*qdagA_0[i][j]+betaU*qdagA_0[i+1][j]+betaL*qdagA_0[i-1][j];}
                }
                DiagArdr=DiagAr;
                DiagAUrdr=DiagAUr;
                DiagALrdr=DiagALr;
                
                TDMA(Nz,DiagALrdr,DiagArdr,DiagAUrdr,bAr);
                
                for (j=0; j<=Nz-1; j++){
                    qdagA[i][j][s]=bAr[j];}
                
            }
            
            for (i=0;i<=Nr-1;i++){
                bAr[i]=0;
            }
            
            for (i=0;i<=(Nr-1);i++){
                for (j=0; j<=Nz-1; j++){
                    qdagA_0[i][j]=qdagA[i][j][s];}
            }
            
            
            
        }
        
        for (j=0;j<=Nr-1;j++){
            bAz[j]=0;
        }
        /***********************************scan over r***************************************************/
        for (j=0;j<=(Nz-1);j++){
            A_Matrix_r(j);
            if (j==0){
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow(delz,2)));
                    beta=(delt/(2.0*pow(delz,2)));
                    bAz[i]=gamma*qdagA_0[i][j]+beta*qdagA_0[j+1][i]+beta*qdagA_0[j+1][i];}
            }
            else if (j==(Nz-1)){
                for (i=0;i<=(Nz-1);i++){
                    gamma=1.0-(delt/(pow(delz,2)));
                    beta=(delt/(2.0*pow(delz,2)));
                    bAz[i]=gamma*qdagA_0[i][j]+beta*qdagA_0[j-1][i]+beta*qdagA_0[j-1][i];}
            }
            else {
                for (i=0;i<=(Nr-1);i++){
                    gamma=1.0-(delt/(pow(delz,2)));
                    beta=(delt/(2.0*pow(delz,2)));
                    bAr[j]=gamma*qdagA_0[i][j]+beta*qdagA_0[j+1][i]+beta*qdagA_0[j-1][i];}
            }
            DiagAzdz=DiagAz;
            DiagAUzdz=DiagAUz;
            DiagALzdz=DiagALz;
            
            TDMA(Nr,DiagALzdz,DiagAzdz,DiagAUzdz,bAz);
            
            for (i=0; j<=Nr-1; i++){
                qdagA[i][j][s]=bAz[i];}
        }
        
        for (i=0;i<=Nr-1;i++){
            bAr[i]=0;
        }
        
        for (i=0;i<=(Nr-1);i++){
            for (j=0; j<=Nz-1; j++){
                qdagA_0[i][j]=qdagA[i][j][s];}
        }
        
        destroy_2d_double_array(qdagA_0);
        destroy_3d_double_array(qdagA);
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
        return qdagA_forward();
    }

        
        
        
        








