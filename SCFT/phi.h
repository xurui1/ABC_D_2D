

/**********************************Function for calculating concentrations***********/
double phi(){
    int i,j,s;
    
    for (i=0;i<=Nr;i++){
        for (j=0;j<=Nz;j++){
            pA[i][j]=0;
            pB[i][j]=0;
            pC[i][j]=0;
            pD[i][j]=0;
        }
    }
    
    
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
            pA[i][j]=pA[i][j]*delt/3.0;
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
                else {cout<< "Error calculating phiC"<< endl;}
                return 0;
            }
            pC[i][j]=pC[i][j]*delt/3.0;
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
                else {cout<< "Error calculating phiD"<< endl;}
                return 0;
            }
            pD[i][j]=pD[i][j]*delt/3.0;
        }
    }
    for (i=0;i<=Nr;i++){
        for (j=0;j<=Nz;j++){
            pA[i][j]=exp(muABC)*(pA[i][j]);
            pB[i][j]=exp(muABC)*(pB[i][j]);
            pC[i][j]=exp(muABC)*(pC[i][j]);
            pD[i][j]=exp(kappaD*muD)*(pD[i][j])/kappaD;
        }
    }
   
    return **pA, **pB, **pC, **pD;
}