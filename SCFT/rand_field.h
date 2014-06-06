int rand_field(int initial){
    
     //Defines the type of initial cond
    int i,j,dummy; //dummy is a garbage index
    //int IOstatus1;
    double rand;// dum_rad; //used for random number
    
    
    //This is for reading the initial omega fields from a files.
    //The M and N should match the size of the data file
    if (initial==1){ //Data for bilayer
        if (bilayer==1){
            if(iter==1){
                fstream myfile;
                myfile.open("./results/shapes/bilayer_M50_N50.dat");
                for (i=1;i<=Nr;i++) {
                    for (j=1;j<=Nz;j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        wA[i][j]=wA[i][j]+(2.0*rand-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*rand-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*rand-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*rand-1.0)*1.0;
                    }
                    
                }
                myfile.close();
            }
            else if(iter != 1){
                fstream myfile;
                myfile.open("./results/omega.bilayer");
                for (i=1;i<=Nr;i++) {
                    for (j=1;j<=Nz;j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        wA[i][j]=wA[i][j]+(2.0*rand-1.0)*1.0;
                        wB[i][j]=wB[i][j]+(2.0*rand-1.0)*1.0;
                        wC[i][j]=wC[i][j]+(2.0*rand-1.0)*1.0;
                        wD[i][j]=wB[i][j]+(2.0*rand-1.0)*1.0;                    }
                }
                myfile.close();
            }
        }
        if (disk==1){
            if(iter==1){
                fstream myfile;
                myfile.open("./results/shapes/bilayer_M50_N50.dat");
                for (i=1;i<=Nr;i++) {
                    for (j=1;j<=Nz;j++){
                        myfile << dummy<<wA[i][j]<<wB[i][j]<<wC[i][j]<<wD[i][j];
                        **wA=**wA+(2.0*rand-1.0)*1.0;
                        **wB=**wB+(2.0*rand-1.0)*1.0;
                        **wC=**wC+(2.0*rand-1.0)*1.0;
                        **wD=**wB+(2.0*rand-1.0)*1.0;
                    }
                    
                    
                }
            }
            else if(iter != 1){
                fstream myfile;
                myfile.open("./results/omega.bilayer");
                for (i=1;i<=Nr;i++) {
                    for (j=1;j<=Nz;j++){
                        myfile << dummy<<**wA<<**wB<<**wC<<**wD;
                        **wA=**wA+(2.0*rand-1.0)*1.0;
                        **wB=**wB+(2.0*rand-1.0)*1.0;
                        **wC=**wC+(2.0*rand-1.0)*1.0;
                        **wD=**wB+(2.0*rand-1.0)*1.0;
                    }
                }
                myfile.close();
            }
        }
    }
    else if (initial ==0){
        for (i=1;i<=Nr;i++) {
            for (j=1;j<=Nz;j++){
                **wA=0.5*rand;
                **wB=0.5*rand;
                **wC=0.5*rand;
                **wD=0.5*rand;
            }
        }
    }
    else if (initial ==2) { //for special cases
        **wA=1;
        **wB=1;
        **wC=-1;
        **wD=1;
        
        
        for (i=1;i <= Nr;i++) {
            for (j=1;j <= 3;j++){
                **wB=-1.0;
                **wB=1.0;
            }
        }
        wD=wB;
    }
    
    return rand_field(initial);
}