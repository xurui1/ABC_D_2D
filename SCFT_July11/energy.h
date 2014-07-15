void FreeEnergy(){
    
    int i,j,ii,jj;
    double F1,F2,F3,F4,F5,F6,F7,F8,F9;
    double FF1,FF2,FF3,FF4,FF5,FF6,FF7,FF8,FF9;
    double p_vect[4];
    double w_vect[4];
    
    fE=0.0;
    
    /*******************Bulk************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F1=0.0;
    FF1=0.0;
    
    for (i=1;i <= int(Nr-2);i++){
        for (j=1;j<= int(Nz-2); j++){
            p_vect[0]=pA[i][j];
            p_vect[1]=pB[i][j];
            p_vect[2]=pC[i][j];
            p_vect[3]=pD[i][j];
            
            w_vect[0]=wA[i][j];
            w_vect[1]=wB[i][j];
            w_vect[2]=wC[i][j];
            w_vect[3]=wD[i][j];
            
            for (ii=0; ii<=3; ii++){
                for (jj=0;jj<=3; jj++){
                    F1=F1+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*((i*delr)+(r_0))*delz;}
                FF1=FF1+p_vect[ii]*w_vect[ii]*delr*((i*delr)+(r_0))*delz;}
            }
        }

    /*******************Side 1*****************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F2=0.0;
    FF2=0.0;
    
    i=0;
        for (j=1;j<= int(Nz-2); j++){
            p_vect[0]=pA[i][j];
            p_vect[1]=pB[i][j];
            p_vect[2]=pC[i][j];
            p_vect[3]=pD[i][j];
            
            w_vect[0]=wA[i][j];
            w_vect[1]=wB[i][j];
            w_vect[2]=wC[i][j];
            w_vect[3]=wD[i][j];
            
            for (ii=0; ii<=3; ii++){
                for (jj=0;jj<=3; jj++){
                    F2=F2+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*((i*delr)+(r_0))*delz*0.5;}
                FF2=FF2+p_vect[ii]*w_vect[ii]*delr*((i*delr)+(r_0))*delz*0.5;}
        }
    
    /*******************Side 2*****************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F3=0.0;
    FF3=0.0;
    
    i=int(Nr-1);
        for (j=1;j<= int(Nz-2); j++){
            p_vect[0]=pA[i][j];
            p_vect[1]=pB[i][j];
            p_vect[2]=pC[i][j];
            p_vect[3]=pD[i][j];
            
            w_vect[0]=wA[i][j];
            w_vect[1]=wB[i][j];
            w_vect[2]=wC[i][j];
            w_vect[3]=wD[i][j];
            
            for (ii=0; ii<=3; ii++){
                for (jj=0;jj<=3; jj++){
                    F3=F3+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*((i*delr)+(r_0))*delz*0.5;}
                FF3=FF3+p_vect[ii]*w_vect[ii]*delr*((i*delr)+(r_0))*delz*0.5;}
        }
    
    /***************************************Side 3*****************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F4=0.0;
    FF4=0.0;
    
    for (i=1;i<=int(Nr-2);i++){
        j=0;
        p_vect[0]=pA[i][j];
        p_vect[1]=pB[i][j];
        p_vect[2]=pC[i][j];
        p_vect[3]=pD[i][j];
            
        w_vect[0]=wA[i][j];
        w_vect[1]=wB[i][j];
        w_vect[2]=wC[i][j];
        w_vect[3]=wD[i][j];
        
        for (ii=0; ii<=3; ii++){
            for (jj=0;jj<=3; jj++){
                F4=F4+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*((i*delr)+(r_0))*delz*0.5;}
            FF4=FF4+p_vect[ii]*w_vect[ii]*delr*((i*delr)+(r_0))*delz*0.5;}
    }
    
    /*******************************************Side 4*****************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F5=0.0;
    FF5=0.0;
    
    for (i=1;i<=int(Nr-2);i++){
        j=int(Nz-1);
        p_vect[0]=pA[i][j];
        p_vect[1]=pB[i][j];
        p_vect[2]=pC[i][j];
        p_vect[3]=pD[i][j];
        
        w_vect[0]=wA[i][j];
        w_vect[1]=wB[i][j];
        w_vect[2]=wC[i][j];
        w_vect[3]=wD[i][j];
        
        for (ii=0; ii<=3; ii++){
            for (jj=0;jj<=3; jj++){
                F5=F5+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*((i*delr)+(r_0))*delz*0.5;}
            FF5=FF5+p_vect[ii]*w_vect[ii]*delr*((i*delr)+(r_0))*delz*0.5;}
        
    }
    /********************************Corner 1******************************************/
    
    i=0;
    j=0;
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    F6=0.0;
    FF6=0.0;
    p_vect[0]=pA[i][j];
    p_vect[1]=pB[i][j];
    p_vect[2]=pC[i][j];
    p_vect[3]=pD[i][j];
    
    w_vect[0]=wA[i][j];
    w_vect[1]=wB[i][j];
    w_vect[2]=wC[i][j];
    w_vect[3]=wD[i][j];
    
    for (ii=0;ii<=3;ii++){
        for (jj=ii;jj<=3;jj++){
            F6=F6+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*((i*delr)+(r_0))*delz*0.25;}
        FF6=FF6+p_vect[ii]*w_vect[ii]*delr*((i*delr)+(r_0))*delz*0.25;
        }
    /********************************Corner 2******************************************/
    
    i=int(Nr-1);
    j=0;
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    F7=0.0;
    FF7=0.0;
    p_vect[0]=pA[i][j];
    p_vect[1]=pB[i][j];
    p_vect[2]=pC[i][j];
    p_vect[3]=pD[i][j];
    
    w_vect[0]=wA[i][j];
    w_vect[1]=wB[i][j];
    w_vect[2]=wC[i][j];
    w_vect[3]=wD[i][j];
    
    for (ii=0;ii<=3;ii++){
        for (jj=ii;jj<=3;jj++){
            F7=F7+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*((i*delr)+(r_0))*delz*0.25;}
        FF7=FF7+p_vect[ii]*w_vect[ii]*delr*((i*delr)+(r_0))*delz*0.25;
    }
    /********************************Corner 3******************************************/
    
    i=0;
    j=int(Nz-1);
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    F8=0.0;
    FF8=0.0;
    p_vect[0]=pA[i][j];
    p_vect[1]=pB[i][j];
    p_vect[2]=pC[i][j];
    p_vect[3]=pD[i][j];
    
    w_vect[0]=wA[i][j];
    w_vect[1]=wB[i][j];
    w_vect[2]=wC[i][j];
    w_vect[3]=wD[i][j];
    
    for (ii=0;ii<=3;ii++){
        for (jj=ii;jj<=3;jj++){
            F8=F8+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*((i*delr)+(r_0))*delz*0.25;}
        FF8=FF8+p_vect[ii]*w_vect[ii]*delr*((i*delr)+(r_0))*delz*0.25;
    }
    /********************************Corner 4******************************************/
    
    i=int(Nr-1);
    j=int(Nz-1);
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    F9=0.0;
    FF9=0.0;
    p_vect[0]=pA[i][j];
    p_vect[1]=pB[i][j];
    p_vect[2]=pC[i][j];
    p_vect[3]=pD[i][j];
    
    w_vect[0]=wA[i][j];
    w_vect[1]=wB[i][j];
    w_vect[2]=wC[i][j];
    w_vect[3]=wD[i][j];
    
    for (ii=0;ii<=3;ii++){
        for (jj=ii;jj<=3;jj++){
            F9=F9+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*((i*delr)+(r_0))*delz*0.25;}
        FF9=FF9+p_vect[ii]*w_vect[ii]*delr*((i*delr)+(r_0))*delz*0.25;
    }
    /******************Time to calculate the free energy!*******************************/
    
    fE=(F1+F2+F3+F4+F5+F6+F7+F8+F9)-(FF1+FF2+FF3+FF4+FF5+FF6+FF7+FF8+FF9);
    fE=fE*2.0*pi;
    fE=(fE/Vol)-(exp(muABC)*Q_ABC)-(exp(muD*kappaD)*Q_D/kappaD);

    dfffE=0.0;
    dfffE=abs(fE-fE_old);
    fE_old=fE;
    
    
}


void fE_homo(){
    
    int i,j;
    double pA_ave,pB_ave,pC_ave,pD_ave;
    double wA_ave,wB_ave,wC_ave,wD_ave;
    double dwA_ave,dwB_ave,dwC_ave,dwD_ave,dpp_ave;
    double eta_ave;
    double f_int, f_omeg;
    double p_vect[4],w_vect[4];
    
    f_int=0.0;
    f_omeg=0.0;
    
    dwA_ave=0.0;
    dwB_ave=0.0;
    dwC_ave=0.0;
    dwD_ave=0.0;
    
    eta_ave=0.0;
    
    pA_ave=0.002;
    pB_ave=pA_ave;
    pC_ave=pA_ave;
    pD_ave=1.0-(pA_ave+pB_ave+pC_ave);
    
    wA_ave=xAB*pB_ave+xAC*pC_ave+xAD*pD_ave+eta_ave;
    wB_ave=xAB*pA_ave+xBC*pC_ave+xBD*pD_ave+eta_ave;
    wC_ave=xAC*pA_ave+xBC*pB_ave+xCD*pD_ave+eta_ave;
    wD_ave=xAD*pA_ave+xBD*pB_ave+xCD*pC_ave+eta_ave;
    
    for (i=1; i<=10000000; i++){
        eta_ave=eta_ave-0.05*(1.0-(pA_ave+pB_ave+pC_ave+pD_ave));
        
        pA_ave=(exp(muABC-wA_ave*fracA-wB_ave*fracB-wC_ave*fracC))*(fracA);
        pB_ave=(exp(muABC-wA_ave*fracA-wB_ave*fracB-wC_ave*fracC))*(fracB);
        pC_ave=(exp(muABC-wA_ave*fracA-wB_ave*fracB-wC_ave*fracC))*(fracC);
        pD_ave=exp(kappaD*(muD-wD_ave));
        
        dwA_ave=(xAB*pB_ave+xAC*pC_ave+xAD*pD_ave+eta_ave)-wA_ave;
        dwB_ave=(xAB*pA_ave+xBC*pC_ave+xBD*pD_ave+eta_ave)-wB_ave;
        dwC_ave=(xAC*pA_ave+xBC*pB_ave+xCD*pD_ave+eta_ave)-wC_ave;
        dwD_ave=(xAD*pA_ave+xBD*pB_ave+xCD*pC_ave+eta_ave)-wD_ave;
        
        dpp_ave=1.0-(pA_ave+pB_ave+pC_ave+pD_ave);
        
        wA_ave=wA_ave+0.005*dwA_ave;
        wB_ave=wB_ave+0.005*dwB_ave;
        wC_ave=wC_ave+0.005*dwC_ave;
        wD_ave=wD_ave+0.005*dwD_ave;
    }
    
    phiABC_hom=pA_ave+pB_ave+pC_ave;
    phiD_hom=pD_ave;
    
    p_vect[0]=pA_ave;
    p_vect[1]=pB_ave;
    p_vect[2]=pC_ave;
    p_vect[3]=pD_ave;
    
    w_vect[0]=wA_ave;
    w_vect[1]=wB_ave;
    w_vect[2]=wC_ave;
    w_vect[3]=wD_ave;
    
    for (i=0; i<=3; i++){
        for (j=0; j<=3; j++){
            f_int=f_int+p_vect[i]*p_vect[j]*XM[i][j];
        }
        f_omeg=f_omeg+p_vect[i]*w_vect[i];
    }
    
    fE_hom=f_int-f_omeg-(exp(muABC-wA_ave*fracA-wB_ave*fracB-wC_ave*fracC))-(exp(kappaD*(muD-wD_ave))/kappaD);
    
}



