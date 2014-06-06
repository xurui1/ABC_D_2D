void FreeEnergy(){
    
    int i,j,ii,jj;
    double F1,F2,F3,F4,F5,F6,F7,F8,F9;
    double FF1,FF2,FF3,FF4,FF5,FF6,FF7,FF8,FF9;
    double p_vect[4]={};
    double w_vect[4]={};
    
    fE=0.0;
    
    /*******************Bulk************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F1=0.0;
    FF1=0.0;
    
    for (i=1;i <= (Nr-2);i++){
        for (j=1;j<= (Nz-2); j++){
            p_vect[0]=pA[i][j];
            p_vect[1]=pB[i][j];
            p_vect[2]=pC[i][j];
            p_vect[3]=pD[i][j];
            
            w_vect[0]=wA[i][j];
            w_vect[1]=wB[i][j];
            w_vect[2]=wC[i][j];
            w_vect[3]=wD[i][j];
            
            for (ii=0; i<=3; ii++){
                for (jj=0;i<=3; jj++){
                    F1=F1+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*(((i-1)*delr)+(r_0))*delz;}
                FF1=FF1+p_vect[ii]*w_vect[ii]*delr*(((i-1)*delr)+(r_0))*delz;}
            }
        }

    /*******************Side 1*****************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F2=0.0;
    FF2=0.0;
    
    for (i=0;;){
        for (j=1;j<= (Nz-2); j++){
            p_vect[0]=pA[i][j];
            p_vect[1]=pB[i][j];
            p_vect[2]=pC[i][j];
            p_vect[3]=pD[i][j];
            
            w_vect[0]=wA[i][j];
            w_vect[1]=wB[i][j];
            w_vect[2]=wC[i][j];
            w_vect[3]=wD[i][j];
            
            for (ii=0; i<=3; ii++){
                for (jj=0;i<=3; jj++){
                    F2=F2+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*(((i-1)*delr)+(r_0))*delz*0.5;}
                FF2=FF2+p_vect[ii]*w_vect[ii]*delr*(((i-1)*delr)+(r_0))*delz*0.5;}
        }
    }
    /*******************Side 2*****************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F3=0.0;
    FF3=0.0;
    
    for (i=(Nr-1);;){
        for (j=1;j<= (Nz-2); j++){
            p_vect[0]=pA[i][j];
            p_vect[1]=pB[i][j];
            p_vect[2]=pC[i][j];
            p_vect[3]=pD[i][j];
            
            w_vect[0]=wA[i][j];
            w_vect[1]=wB[i][j];
            w_vect[2]=wC[i][j];
            w_vect[3]=wD[i][j];
            
            for (ii=0; i<=3; ii++){
                for (jj=0;i<=3; jj++){
                    F3=F3+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*(((i-1)*delr)+(r_0))*delz*0.5;}
                FF3=FF3+p_vect[ii]*w_vect[ii]*delr*(((i-1)*delr)+(r_0))*delz*0.5;}
        }
    }
    /***************************************Side 3*****************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F4=0.0;
    FF4=0.0;
    
    for (i=1;i<=(Nr-2);i++){
        for (j=0;;){
            p_vect[0]=pA[i][j];
            p_vect[1]=pB[i][j];
            p_vect[2]=pC[i][j];
            p_vect[3]=pD[i][j];
            
            w_vect[0]=wA[i][j];
            w_vect[1]=wB[i][j];
            w_vect[2]=wC[i][j];
            w_vect[3]=wD[i][j];
            
            for (ii=0; i<=3; ii++){
                for (jj=0;i<=3; jj++){
                    F4=F4+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*(((i-1)*delr)+(r_0))*delz*0.5;}
                FF4=FF4+p_vect[ii]*w_vect[ii]*delr*(((i-1)*delr)+(r_0))*delz*0.5;}
        }
    }
    /*******************************************Side 4*****************************/
    memset(p_vect,0,4);
    memset(w_vect,0,4);
    
    F5=0.0;
    FF5=0.0;
    
    for (i=1;i<=(Nr-2);i++){
        for (j=(Nz-1);;){
            p_vect[0]=pA[i][j];
            p_vect[1]=pB[i][j];
            p_vect[2]=pC[i][j];
            p_vect[3]=pD[i][j];
            
            w_vect[0]=wA[i][j];
            w_vect[1]=wB[i][j];
            w_vect[2]=wC[i][j];
            w_vect[3]=wD[i][j];
            
            for (ii=0; i<=3; ii++){
                for (jj=0;i<=3; jj++){
                    F5=F5+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*(((i-1)*delr)+(r_0))*delz*0.5;}
                FF5=FF5+p_vect[ii]*w_vect[ii]*delr*(((i-1)*delr)+(r_0))*delz*0.5;}
        }
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
            F6=F6+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*(((i-1)*delr)+(r_0))*delz*0.25;}
        FF6=FF6+p_vect[ii]*w_vect[ii]*delr*(((i-1)*delr)+(r_0))*delz*0.25;
        }
    /********************************Corner 2******************************************/
    
    i=(Nr-1);
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
            F7=F7+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*(((i-1)*delr)+(r_0))*delz*0.25;}
        FF7=FF7+p_vect[ii]*w_vect[ii]*delr*(((i-1)*delr)+(r_0))*delz*0.25;
    }
    /********************************Corner 3******************************************/
    
    i=0;
    j=(Nz-1);
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
            F8=F8+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*(((i-1)*delr)+(r_0))*delz*0.25;}
        FF8=FF8+p_vect[ii]*w_vect[ii]*delr*(((i-1)*delr)+(r_0))*delz*0.25;
    }
    /********************************Corner 4******************************************/
    
    i=Nr-1;
    j=Nz-1;
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
            F9=F9+p_vect[ii]*p_vect[jj]*XM[ii][jj]*delr*(((i-1)*delr)+(r_0))*delz*0.25;}
        FF9=FF9+p_vect[ii]*w_vect[ii]*delr*(((i-1)*delr)+(r_0))*delz*0.25;
    }
    /******************Time to calculate the free energy!*******************************/
    
    fE=(F1+F2+F3+F4+F5+F6+F7+F8+F9)-(FF1+FF2+FF3+FF4+FF5+FF6+FF7+FF8+FF9);
    fE=fE*2.0*pi;
    fE=(fE/Vol)-(exp(muABC)*Q_ABC)-(exp(muD*kappaD)*Q_D/kappaD);
    
    //printing individual parts of fE
    //print*,-(exp(muAB)*Q_AB)-(exp(muC*kappaC)*Q_C/kappaC)-(exp(muED*kappaED)*Q_ED/kappaED), &
    //(F1+F2+F3+F4+F5+F6+F7+F8+F9)*2.0*pi/vol, &
    //(FF1+FF2+FF3+FF4+FF5+FF6+FF7+FF8+FF9)*2.0*pi/vol
    
    
    dfffE=0.0;
    dfffE=abs(fE-fE_old);
    fE_old=fE;
    
    destroy_2d_double_array(pA);
    destroy_2d_double_array(pB);
    destroy_2d_double_array(pC);
    destroy_2d_double_array(pD);
    destroy_2d_double_array(wA);
    destroy_2d_double_array(wB);
    destroy_2d_double_array(wC);
    destroy_2d_double_array(wD);
    
    

    return FreeEnergy();
}




