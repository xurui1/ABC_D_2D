/************************calculating eta (incompressibility condition)***********************************/

int pressure(){
    int i,j;
    for (i=1;i<=Nr;i++) {
        for (j=1;j<=Nz;j++){
            eta[i][j]=eta[i][j]-(1.0-(pA[i][j]+pB[i][j]+pC[i][j]+pD[i][j]));
        }
    }
    return pressure();
}
/***********************calculating eta2 (pinning condition)***********************************/
int pressure2(){
    int i,j;
    for (i=Ntip;;){
        for (j=Mtip;;){
            eta2[i][j]=eta2[i][j]-10*((pA[i][j]-pD[i][j])-(pB[i][j]+pC[i][j]));
        }
    }
    
    return pressure2();
}