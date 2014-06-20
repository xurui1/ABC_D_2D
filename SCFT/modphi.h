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
            phiA=phiA+(pA[i][j])*delr*(float(i-1)*delr+r_0)*delz; //why does the code require float(i-1)?
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



/**********************************Function for calculating concentrations***********/
void phi(){
    int i,j,s;
    
    for (i=0;i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            pA[i][j]=0;
            pB[i][j]=0;
            pC[i][j]=0;
            pD[i][j]=0;
        }
    }
    /******************************Calculating phiA***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=int(NA-1);s++){
                if (s==0 or s==int(NA-1)){
                    pA[i][j]+= ((qA[i][j][s])*(qdagA[i][j][int(NA-s)]));}
                else if ((s%2) !=0 and s!=int(NA-1)){
                    pA[i][j]+= 4*((qA[i][j][s])*(qdagA[i][j][int(NA-s)]));}
                else if ((s%2) ==0 and s!=int(NA-1) and s!=0){
                    pA[i][j]+= 2*((qA[i][j][s])*(qdagA[i][j][int(NA-s)]));}
                else {cout<< "Error calculating phiA"<< endl;break;}
            }
            pA[i][j]=pA[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiB***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=int(NB-1);s++){
                if (s==0 or s==int(NB-1)){
                    pB[i][j]+= ((qB[i][j][s])*(qdagB[i][j][int(NB-s)]));}
                else if (s%2 !=0 and s!=int(NB-1)){
                    pB[i][j]+= 4*((qB[i][j][s])*(qdagB[i][j][int(NB-s)]));}
                else if (s%2 ==0 and s!=int(NB-1) and s!=0){
                    pB[i][j]+= 2*((qB[i][j][s])*(qdagB[i][j][int(NB-s)]));}
                else {cout<< "Error calculating phiB"<< endl;break;}
            }
            pB[i][j]=pB[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiC***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=int(NC-1);s++){
                if (s==0 or s==int(NC-1)){
                    pC[i][j]+= ((qC[i][j][s])*(qdagC[i][j][int(NC-s)]));}
                else if (s%2 !=0 and s!=int(NC-1)){
                    pC[i][j]+= 4*((qC[i][j][s])*(qdagC[i][j][int(NC-s)]));}
                else if (s%2 ==0 and s!=int(NC-1) and s!=0){
                    pC[i][j]+= 2*((qC[i][j][s])*(qdagC[i][j][int(NC-s)]));}
                else {cout<< "Error calculating phiC"<< endl;break;}
            }
            pC[i][j]=pC[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiD***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=int(ND-1);s++){
                if (s==0 or s==int(ND-1)){
                    pD[i][j]+= ((qD[i][j][s])*(qD[i][j][int(ND-s)]));}
                else if (s%2 !=0 and s!=int(ND-1)){
                    pD[i][j]+= 4*((qD[i][j][s])*(qD[i][j][int(ND-s)]));}
                else if (s%2 ==0 and s!=int(ND-1) and s!=0){
                    pD[i][j]+= 2*((qD[i][j][s])*(qD[i][j][int(ND-s)]));}
                else {cout<< "Error calculating phiD"<< endl;break;}
            }
            pD[i][j]=pD[i][j]*delt/3.0;
        }
    }
    for (i=0;i<=int(Nr-1);i++){
        for (j=0;j<=int(Nz-1);j++){
            pA[i][j]=exp(muABC)*(pA[i][j]);
            pB[i][j]=exp(muABC)*(pB[i][j]);
            pC[i][j]=exp(muABC)*(pC[i][j]);
            pD[i][j]=exp(kappaD*muD)*(pD[i][j])/kappaD;
        }
    }
}
