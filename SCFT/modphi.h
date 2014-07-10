void totalphi(){
    
    int i,j;
    
    phiA=0.0;
    phiB=0.0;
    phiC=0.0;
    phiD=0.0;
    ptot_in=0.0;

    for (i=1; i<=int(Nr-2); i++){
        for (j=1; j<=int(Nz-2); j++){
            phiA+=(pA[i][j])*delr*(i*delr+r_0)*delz;
            phiB+=(pB[i][j])*delr*(i*delr+r_0)*delz;
            phiC+=(pC[i][j])*delr*(i*delr+r_0)*delz;
            phiD+=(pD[i][j])*delr*(i*delr+r_0)*delz;
        }}
    
    i=0;
    for (j=1; j<=int(Nz-2); j++){
        phiA+=0.5*(pA[i][j])*delr*(i*delr+r_0)*delz;
        phiB+=0.5*(pB[i][j])*delr*(i*delr+r_0)*delz;
        phiC+=0.5*(pC[i][j])*delr*(i*delr+r_0)*delz;
        phiD+=0.5*(pD[i][j])*delr*(i*delr+r_0)*delz;
    }
    
    i=int(Nr-1);
    for (j=1; j<=int(Nz-2); j++){
        phiA+=0.5*(pA[i][j])*delr*(i*delr+r_0)*delz;
        phiB+=0.5*(pB[i][j])*delr*(i*delr+r_0)*delz;
        phiC+=0.5*(pC[i][j])*delr*(i*delr+r_0)*delz;
        phiD+=0.5*(pD[i][j])*delr*(i*delr+r_0)*delz;
    }
    
    j=0;
    for (i=1; i<=int(Nr-2); i++){
        phiA+=0.5*(pA[i][j])*delr*(i*delr+r_0)*delz;
        phiB+=0.5*(pB[i][j])*delr*(i*delr+r_0)*delz;
        phiC+=0.5*(pC[i][j])*delr*(i*delr+r_0)*delz;
        phiD+=0.5*(pD[i][j])*delr*(i*delr+r_0)*delz;
    }
    
    j=int(Nz-1);
    for (i=1; i<=int(Nr-2); i++){
        phiA+=0.5*(pA[i][j])*delr*(i*delr+r_0)*delz;
        phiB+=0.5*(pB[i][j])*delr*(i*delr+r_0)*delz;
        phiC+=0.5*(pC[i][j])*delr*(i*delr+r_0)*delz;
        phiD+=0.5*(pD[i][j])*delr*(i*delr+r_0)*delz;
    }
    
    phiA+=0.25*delr*delz*((pA[0][0]*(double(0)*delr+r_0))+
                              (pA[0][int(Nz-1)]*(double(0)*delr+r_0))+
                              (pA[int(Nr-1)][0]*(double(Nr-1)*delr+r_0))+
                              (pA[int(Nr-1)][int(Nz-1)]*(double(Nr-1)*delr+r_0)));
    
    phiB+=0.25*delr*delz*((pB[0][0]*(double(0)*delr+r_0))+
                              (pB[0][int(Nz-1)]*(double(0)*delr+r_0))+
                              (pB[int(Nr-1)][0]*(double(Nr-1)*delr+r_0))+
                              (pB[int(Nr-1)][int(Nz-1)]*(double(Nr-1)*delr+r_0)));
    
    phiC+=0.25*delr*delz*((pC[0][0]*(double(0)*delr+r_0))+
                            (pC[0][int(Nz-1)]*(double(0)*delr+r_0))+
                              (pC[int(Nr-1)][0]*(double(Nr-1)*delr+r_0))+
                              (pC[int(Nr-1)][int(Nz-1)]*(double(Nr-1)*delr+r_0)));
  
    phiD+=0.25*delr*delz*((pD[0][0]*(double(0)*delr+r_0))+
                              (pD[0][int(Nz-1)]*(double(0)*delr+r_0))+
                              (pD[int(Nr-1)][0]*(double(Nr-1)*delr+r_0))+
                              (pD[int(Nr-1)][int(Nz-1)]*(double(Nr-1)*delr+r_0)));
    

    
    phiA=(2.0*pi*phiA)/Vol;
    phiB=(2.0*pi*phiB)/Vol;
    phiC=(2.0*pi*phiC)/Vol;
    phiD=(2.0*pi*phiD)/Vol;
    
    ptot_in=phiA+phiB+phiC+phiD;
    //cout<<ptot_in<<" "<<phiA<<" "<<phiB<<" "<<phiC<<" "<<phiD<<endl;

    
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
            for (s=0;s<=NA;s++){
                if (s==0 or s==NA){
                    pA[i][j]+= ((qA[i][j][s])*(qdagA[i][j][int(NA-s)]));}
                else if ((s%2) !=0 and s!=NA){
                    pA[i][j]+= 4*((qA[i][j][s])*(qdagA[i][j][int(NA-s)]));}
                else if ((s%2) ==0 and s!=NA and s!=0){
                    pA[i][j]+= 2*((qA[i][j][s])*(qdagA[i][j][int(NA-s)]));}
                else {cout<< "Error calculating phiA"<< endl;break;}
            }
            pA[i][j]=pA[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiB***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=NB;s++){
                if (s==0 or s==NB){
                    pB[i][j]+= ((qB[i][j][s])*(qdagB[i][j][int(NB-s)]));}
                else if (s%2 !=0 and s!=NB){
                    pB[i][j]+= 4*((qB[i][j][s])*(qdagB[i][j][int(NB-s)]));}
                else if (s%2 ==0 and s!=NB and s!=0){
                    pB[i][j]+= 2*((qB[i][j][s])*(qdagB[i][j][int(NB-s)]));}
                else {cout<< "Error calculating phiB"<< endl;break;}
            }
            pB[i][j]=pB[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiC***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=NC;s++){
                if (s==0 or s==NC){
                    pC[i][j]+= ((qC[i][j][s])*(qdagC[i][j][int(NC-s)]));}
                else if (s%2 !=0 and s!=NC){
                    pC[i][j]+= 4*((qC[i][j][s])*(qdagC[i][j][int(NC-s)]));}
                else if (s%2 ==0 and s!=NC and s!=0){
                    pC[i][j]+= 2*((qC[i][j][s])*(qdagC[i][j][int(NC-s)]));}
                else {cout<< "Error calculating phiC"<< endl;break;}
            }
            pC[i][j]=pC[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiD***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=ND;s++){
                if (s==0 or s==ND){
                    pD[i][j]+= ((qD[i][j][s])*(qD[i][j][int(ND-s)]));}
                else if (s%2 !=0 and s!=ND){
                    pD[i][j]+= 4*((qD[i][j][s])*(qD[i][j][int(ND-s)]));}
                else if (s%2 ==0 and s!=ND and s!=0){
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
