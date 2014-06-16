double totalphi(){
    
    int i,j;
    
    phiA=0.0;
    phiB=0.0;
    phiC=0.0;
    phiD=0.0;
    ptot_in=0.0;
    phistar=0.0;
    
    for (i=1; i<=(Nr-2); i++){
        for (j=1; j<=(Nz-2); j++){
            phiA=phiA+(pA[i][j])*delr*(float(i-1)*delr+r_0)*delz;
            phiB=phiB+(pB[i][j])*delr*(float(i-1)*delr+r_0)*delz;
            phiC=phiC+(pC[i][j])*delr*(float(i-1)*delr+r_0)*delz;
            phiD=phiD+(pD[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        }}
    
    i=0;
    for (j=1; j<=(Nz-2); j++){
        phiA=phiA+0.5*(pA[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiB=phiB+0.5*(pB[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiC=phiC+0.5*(pC[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiD=phiD+0.5*(pD[i][j])*delr*(float(i-1)*delr+r_0)*delz;
    }
    
    i=(Nr-1);
    for (j=1; j<=(Nz-2); j++){
        phiA=phiA+0.5*(pA[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiB=phiB+0.5*(pB[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiC=phiC+0.5*(pC[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiD=phiD+0.5*(pD[i][j])*delr*(float(i-1)*delr+r_0)*delz;
    }
    
    j=0;
    for (i=1; i<=(Nr-2); i++){
        phiA=phiA+0.5*(pA[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiB=phiB+0.5*(pB[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiC=phiC+0.5*(pC[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiD=phiD+0.5*(pD[i][j])*delr*(float(i-1)*delr+r_0)*delz;
    }
    
    j=(Nz-1);
    for (i=1; i<=(Nr-2); i++){
        phiA=phiA+0.5*(pA[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiB=phiB+0.5*(pB[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiC=phiC+0.5*(pC[i][j])*delr*(float(i-1)*delr+r_0)*delz;
        phiD=phiD+0.5*(pD[i][j])*delr*(float(i-1)*delr+r_0)*delz;
    }
    
    phiA=phiA+0.25*delr*delz*((pA[0][0]*(float(i-1)*delr+r_0))+(pA[0][int(Nz-1)]*(float(i-1)*delr+r_0))+
                              (pA[int(Nr-1)][0]*(float(Nr-1)*delr+r_0))+(pA[int(Nr-1)][int(Nz-1)]*(float(Nr-1)*delr+r_0)));
    
    phiB=phiB+0.25*delr*delz*((pB[0][0]*(float(i-1)*delr+r_0))+(pB[0][int(Nz-1)]*(float(i-1)*delr+r_0))+
                              (pB[int(Nr-1)][0]*(float(Nr-1)*delr+r_0))+(pB[int(Nr-1)][int(Nz-1)]*(float(Nr-1)*delr+r_0)));
    
    phiC=phiC+0.25*delr*delz*((pC[0][0]*(float(i-1)*delr+r_0))+(pC[0][int(Nz-1)]*(float(i-1)*delr+r_0))+
                              (pC[int(Nr-1)][0]*(float(Nr-1)*delr+r_0))+(pC[int(Nr-1)][int(Nz-1)]*(float(Nr-1)*delr+r_0)));
  
    phiD=phiD+0.25*delr*delz*((pD[0][0]*(float(i-1)*delr+r_0))+(pD[0][int(Nz-1)]*(float(1-1)*delr+r_0))+
                              (pD[int(Nr-1)][0]*(float(Nr-1)*delr+r_0))+(pD[int(Nr-1)][int(Nz-1)]*(float(Nr-1)*delr+r_0)));
    
    phiA=2.0*pi*phiA/Vol;
    phiB=2.0*pi*phiB/Vol;
    phiC=2.0*pi*phiC/Vol;
    phiD=2.0*pi*phiD/Vol;
    
    ptot_in=phiA+phiB+phiC+phiD;
    
    return phiD;
}
