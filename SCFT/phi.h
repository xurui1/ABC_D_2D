

/**********************************Function for calculating concentrations***********/
int phi(){
    int i,j,s;
    **pA=0;**pB=0;**pC=0;**pD=0;
    
    /******************************Calculating phiA***********************************/
    for (i=1;i<=Nr;i++) {
        for (j=1;j<=Nz;j++){
            for (s=0;s<=NA;s++){
                if (s==0 or s==NA){
                    pA[i][j]=pA[i][j]+((qA[i][j][s])*(qdagA[i][j][NA-s]));}
                else if ((s%2) !=0 and s!=NA){
                    pA[i][j]=pA[i][j]+4*((qA[i][j][s])*(qdagA[i][j][NA-s]));}
                else if ((s%2) ==0 and s!=NA and s!=0){
                    pA[i][j]=pA[i][j]+2*((qA[i][j][s])*(qdagA[i][j][NA-s]));}
                else {cout<< "Error calculating phiA"<< endl;}
                return 0;
            }
            **pA=**pA*delt/3.0;
        }
    }
    /******************************Calculating phiB***********************************/
    for (i=1;i<=Nr;i++) {
        for (j=1;j<=Nz;j++){
            for (s=0;s<=NB;s++){
                if (s==0 or s==NB){
                    pB[i][j]=pB[i][j]+((qB[i][j][s])*(qdagB[i][j][NB-s]));}
                else if (s%2 !=0 and s!=NB){
                    pB[i][j]=pB[i][j]+4*((qB[i][j][s])*(qdagB[i][j][NB-s]));}
                else if (s%2 ==0 and s!=NB and s!=0){
                    pB[i][j]=pB[i][j]+2*((qB[i][j][s])*(qdagB[i][j][NB-s]));}
                else {cout<< "Error calculating phiB"<< endl;}
                return 0;
            }
            pB[i][j]=pB[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiC***********************************/
    for (i=1;i<=Nr;i++) {
        for (j=1;j<=Nz;j++){
            for (s=0;s<=NC;s++){
                if (s==0 or s==NC){
                    pC[i][j]=pC[i][j]+((qC[i][j][s])*(qdagC[i][j][NC-s]));}
                else if (s%2 !=0 and s!=NC){
                    pC[i][j]=pC[i][j]+4*((qC[i][j][s])*(qdagC[i][j][NC-s]));}
                else if (s%2 ==0 and s!=NC and s!=0){
                    pC[i][j]=pC[i][j]+2*((qC[i][j][s])*(qdagC[i][j][NC-s]));}
                else {cout<< "Error calculating phiB"<< endl;}
                return 0;
            }
            **pC=**pC*delt/3.0;
        }
    }
    /******************************Calculating phiD***********************************/
    for (i=1;i<=Nr;i++) {
        for (j=1;j<=Nz;j++){
            for (s=0;s<=ND;s++){
                if (s==0 or s==ND){
                    pD[i][j]=pD[i][j]+((qD[i][j][s])*(qD[i][j][ND-s]));}
                else if (s%2 !=0 and s!=ND){
                    pD[i][j]=pD[i][j]+4*((qD[i][j][s])*(qD[i][j][ND-s]));}
                else if (s%2 ==0 and s!=ND and s!=0){
                    pD[i][j]=pD[i][j]+2*((qD[i][j][s])*(qD[i][j][ND-s]));}
                else {cout<< "Error calculating phiB"<< endl;}
                return 0;
            }
            **pD=**pD*delt/3.0;
        }
    }
    
    **pA=exp(muABC)*(**pA);
    **pB=exp(muABC)*(**pB);
    **pC=exp(muABC)*(**pC);
    **pD=exp(kappaD*muD)*(**pD)/kappaD;
    
    destroy_2d_double_array(pA);
    destroy_2d_double_array(pB);
    destroy_2d_double_array(pC);
    destroy_2d_double_array(pD);
    destroy_3d_double_array(qA);
    destroy_3d_double_array(qB);
    destroy_3d_double_array(qC);
    destroy_3d_double_array(qD);
    destroy_3d_double_array(qdagA);
    destroy_3d_double_array(qdagB);
    destroy_3d_double_array(qdagC);
    
    
    
    return phi();
}