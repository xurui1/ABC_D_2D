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
    
    for (i=1; i<=10000; i++){
        eta_ave=eta_ave-0.05*(1.0-(pA_ave+pB_ave+pC_ave+pD_ave));
        
        pA_ave=(exp(muABC-wA_ave*fracA-wB_ave*fracB-wC_ave*fracC))*(1.0-fracB-fracC);
        pB_ave=(exp(muABC-wA_ave*fracA-wB_ave*fracB-wC_ave*fracC))*(1.0-fracA-fracC);
        pC_ave=(exp(muABC-wA_ave*fracA-wB_ave*fracB-wC_ave*fracC))*(1.0-fracA-fracB);
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