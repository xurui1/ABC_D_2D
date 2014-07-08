double unifRand()
{
    return double(rand()) / double(RAND_MAX);
}

void rand_field(int initial){
    
#include <ctime>      //Call system time libraries to define the integer seed for random numbers

    int i,j,dummy;
    srand((unsigned)time(NULL));
    cout<<unifRand()<<endl;
    cout<<unifRand()<<endl;

    
    //This is for reading the initial omega fields from a files.
    //The M and N should match the size of the data file
    if (initial==1){ //Data for bilayer
        if (bilayer==1){
            if(iter==1){
                fstream myfile;
                myfile.open("results/shapes/bilayer_M50_N50.dat");
                for (i=0;i<=int(Nr-1);i++) {
                    for (j=0;j<=int(Nz-1);j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        wA[i][j]=wA[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*unifRand()-1.0)*1.0;
                    }
                    
                }
                myfile.close();
            }
            else if(iter != 1){
                fstream myfile;
                myfile.open("results/omega.bilayer");
                for (i=0;i<=int(Nr-1);i++) {
                    for (j=0;j<=int(Nz-1);j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        wA[i][j]=wA[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*unifRand()-1.0)*1.0;                    }
                }
                myfile.close();
            }
        }
        if (disk==1){
            if(iter==1){
                fstream myfile;
                myfile.open("results/shapes/disk_M50_N50.dat");
                for (i=0;i<=int(Nr-1);i++) {
                    for (j=1;j<=int(Nz-1);j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        wA[i][j]=wA[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*unifRand()-1.0)*1.0;
                    }
                }
                myfile.close();
            }
            else if(iter != 1){
                fstream myfile;
                myfile.open("results/omega.bilayer");
                for (i=0;i<=int(Nr-1);i++) {
                    for (j=1;j<=int(Nz-1);j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        wA[i][j]=wA[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*unifRand()-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*unifRand()-1.0)*1.0;
                    }
                }
                myfile.close();
            }
        }
    }
    else if (initial ==0){
        for (i=0;i<=int(Nr-1);i++) {
            for (j=0;j<=int(Nz-1);j++){
                wA[i][j]=0.5*unifRand();
                wB[i][j]=0.5*unifRand();
                wC[i][j]=0.5*unifRand();
                wD[i][j]=0.5*unifRand();}}
    }
    else if (initial ==2) { //for special cases
        for (i=0;i<=int(Nr-1);i++) {
            for (j=0;j<=int(Nz-1);j++){
                wA[i][j]=1;
                wB[i][j]=1;
                wC[i][j]=-1;
                wD[i][j]=1;}}
        
        for (i=0;i <= int(Nr-1);i++) {
            for (j=1;j <= 3;j++){
                wB[i][j]=-1.0;
                wB[i][j]=1.0;
            }
        }
        for (i=0;i<=int(Nr-1);i++) {
            for (j=0;j<=int(Nz-1);j++){
                wD=wB;}}
    }
}







/**********************************Calculate new omega fields***************************************/

void new_fields (){
    int i,j;
    Conv_w=0.0;
    Conv_p=0.0;
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            dwA[i][j]= xAB*(pB[i][j])+xAC*(pC[i][j])+xAD*(pD[i][j])-(eta2[i][j])+(eta[i][j])-(wA[i][j]);
            dwB[i][j]= xAB*(pB[i][j])+xBC*(pC[i][j])+xBD*(pD[i][j])+(eta2[i][j])+(eta[i][j])-(wB[i][j]);
            dwC[i][j]= xAC*(pA[i][j])+xBC*(pB[i][j])+xCD*(pD[i][j])+(eta2[i][j])+(eta[i][j])-(wC[i][j]);
            dwD[i][j]= xAD*(pA[i][j])+xBD*(pB[i][j])+xCD*(pC[i][j])-(eta2[i][j])+(eta[i][j])-(wD[i][j]);
            dpp[i][j]=1.0-(pA[i][j]+pB[i][j]+pC[i][j]+pD[i][j]);
            
            //changed (-) to (+) for eta2, according to Kyles paper (for hydrophobic components).
            
            Conv_w=Conv_w+abs((dwA[i][j])+(dwB[i][j])+(dwC[i][j])+(dwD[i][j]));
            Conv_p=Conv_p+abs(dpp[i][j]);
            
            // updating the omega condition
            wA[i][j]=(wA[i][j])+sig*(dwA[i][j])-sig2*(dpp[i][j]);
            wB[i][j]=(wB[i][j])+sig*(dwB[i][j])-sig2*(dpp[i][j]);
            wC[i][j]=(wC[i][j])+sig*(dwC[i][j])-sig2*(dpp[i][j]);
            wD[i][j]=(wD[i][j])+sig*(dwD[i][j])-sig2*(dpp[i][j]);
            
        }
    }
    Conv_w=Conv_w/(Nr*Nz);
    Conv_p=Conv_p/(Nr*Nz);
}
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
    i=int(Ntip-1);
    j=int(Mtip-1);
    eta2[i][j]=eta2[i][j]-10*((pA[i][j]+pD[i][j])-(pB[i][j]+pC[i][j]));
}