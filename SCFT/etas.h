/************************calculating eta (incompressibility condition)***********************************/

void pressure(){
    int i,j;
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            eta[i][j]=eta[i][j]-(1.0-(pA[i][j]+pB[i][j]+pC[i][j]+pD[i][j]));
        }
    }
}
/***********************calculating eta2 (pinning condition)***********************************/
void pressure2(){
    cout<< "pressure2"<<endl;
    int i,j;
    for (i=Ntip;;){
        for (j=Mtip;;){
            eta2[i][j]=eta2[i][j]-10*((pA[i][j]+pD[i][j])-(pB[i][j]+pC[i][j]));
        }
    }
}