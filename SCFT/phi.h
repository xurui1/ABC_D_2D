

/**********************************Function for calculating concentrations***********/
double phi(){
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
                    pA[i][j]=pA[i][j]+((qA[i][j][s])*(qdagA[i][j][int(NA-s)]));}
                else if ((s%2) !=0 and s!=int(NA-1)){
                    pA[i][j]=pA[i][j]+4*((qA[i][j][s])*(qdagA[i][j][int(NA-s)]));}
                else if ((s%2) ==0 and s!=int(NA-1) and s!=0){
                    pA[i][j]=pA[i][j]+2*((qA[i][j][s])*(qdagA[i][j][int(NA-s)]));}
                else {cout<< "Error calculating phiA"<< endl;return 0;}
            }
            pA[i][j]=pA[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiB***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=int(NB-1);s++){
                if (s==0 or s==int(NB-1)){
                    pB[i][j]=pB[i][j]+((qB[i][j][s])*(qdagB[i][j][int(NB-s)]));}
                else if (s%2 !=0 and s!=int(NB-1)){
                    pB[i][j]=pB[i][j]+4*((qB[i][j][s])*(qdagB[i][j][int(NB-s)]));}
                else if (s%2 ==0 and s!=int(NB-1) and s!=0){
                    pB[i][j]=pB[i][j]+2*((qB[i][j][s])*(qdagB[i][j][int(NB-s)]));}
                else {cout<< "Error calculating phiB"<< endl;return 0;}
            }
            pB[i][j]=pB[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiC***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=int(NC-1);s++){
                if (s==0 or s==int(NC-1)){
                    pC[i][j]=pC[i][j]+((qC[i][j][s])*(qdagC[i][j][int(NC-s)]));}
                else if (s%2 !=0 and s!=int(NC-1)){
                    pC[i][j]=pC[i][j]+4*((qC[i][j][s])*(qdagC[i][j][int(NC-s)]));}
                else if (s%2 ==0 and s!=int(NC-1) and s!=0){
                    pC[i][j]=pC[i][j]+2*((qC[i][j][s])*(qdagC[i][j][int(NC-s)]));}
                else {cout<< "Error calculating phiC"<< endl;return 0;}
            }
            pC[i][j]=pC[i][j]*delt/3.0;
        }
    }
    /******************************Calculating phiD***********************************/
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            for (s=0;s<=int(ND-1);s++){
                if (s==0 or s==int(ND-1)){
                    pD[i][j]=pD[i][j]+((qD[i][j][s])*(qD[i][j][int(ND-s)]));}
                else if (s%2 !=0 and s!=int(ND-1)){
                    pD[i][j]=pD[i][j]+4*((qD[i][j][s])*(qD[i][j][int(ND-s)]));}
                else if (s%2 ==0 and s!=int(ND-1) and s!=0){
                    pD[i][j]=pD[i][j]+2*((qD[i][j][s])*(qD[i][j][int(ND-s)]));}
                else {cout<< "Error calculating phiD"<< endl;return 0;}
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
   
    return **pD;
}