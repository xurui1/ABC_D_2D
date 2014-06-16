void rand_field(int initial){
    
#include <ctime>      //Call system time libraries to define the integer seed for random numbers

    int i,j,dummy;
    srand (static_cast <unsigned> (time(0)));
    double random=(double)rand()/(double)RAND_MAX;
    
    cout<<random<<endl;
    cout<<random<<endl;
    //This is for reading the initial omega fields from a files.
    //The M and N should match the size of the data file
    if (initial==1){ //Data for bilayer
        if (bilayer==1){
            if(iter==1){
                fstream myfile;
                myfile.open("./results/shapes/bilayer_M50_N50.dat");
                for (i=0;i<=int(Nr-1);i++) {
                    for (j=0;j<=int(Nz-1);j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        wA[i][j]=wA[i][j]+(2.0*random-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*random-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*random-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*random-1.0)*1.0;
                    }
                    
                }
                myfile.close();
            }
            else if(iter != 1){
                fstream myfile;
                myfile.open("./results/omega.bilayer");
                for (i=0;i<=int(Nr-1);i++) {
                    for (j=0;j<=int(Nz-1);j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        wA[i][j]=wA[i][j]+(2.0*random-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*random-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*random-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*random-1.0)*1.0;                    }
                }
                myfile.close();
            }
        }
        if (disk==1){
            if(iter==1){
                fstream myfile;
                myfile.open("./results/shapes/bilayer_M50_N50.dat");
                for (i=0;i<=int(Nr-1);i++) {
                    for (j=1;j<=int(Nz-1);j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        wA[i][j]=wA[i][j]+(2.0*random-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*random-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*random-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*random-1.0)*1.0;
                    }
                    
                    
                }
            }
            else if(iter != 1){
                fstream myfile;
                myfile.open("./results/omega.bilayer");
                for (i=0;i<=int(Nr-1);i++) {
                    for (j=1;j<=int(Nz-1);j++){
                        myfile << dummy<<**wA<<**wB<<**wC<<**wD;
                        wA[i][j]=wA[i][j]+(2.0*random-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*random-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*random-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*random-1.0)*1.0;
                    }
                }
                myfile.close();
            }
        }
    }
    else if (initial ==0){
        for (i=0;i<=int(Nr-1);i++) {
            for (j=0;j<=int(Nz-1);j++){
                wA[i][j]=0.5*random;
                wB[i][j]=0.5*random;
                wC[i][j]=0.5*random;
                wD[i][j]=0.5*random;
            }
        }
    }
    else if (initial ==2) { //for special cases
        for (i=0;i<=int(Nr-1);i++) {
            for (j=0;j<=int(Nz-1);j++){
                wA[i][j]=1;
                wB[i][j]=1;
                wC[i][j]=-1;
                wD[i][j]=1;
            }
        }
        
        
        for (i=0;i <= int(Nr-1);i++) {
            for (j=1;j <= 3;j++){
                wB[i][j]=-1.0;
                wB[i][j]=1.0;
            }
        }
        wD=wB;
    }
    
    
}