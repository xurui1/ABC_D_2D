void profile(int msg){
    int  i,j,ii;
    float cr,cz,ct,epsct;
    float px,py,pz;
    string outfile,orderpar;
    
    if (msg==1){
        epsct=(2.0*pi)/Nr;
        fstream myfile;
        myfile.open("./results/profile1.dat");
        for (i=0;i<=(Nr-1);i++){
            for (j=0;j<=(Nz-1);j++){
                myfile<<(i*delr)<<(j*delz)<<pA[i][j]<<pB[i][j]<<pC[i][j]<<pD[i][j]<<endl;
            }
        }
        myfile.close();
    }
    else if (msg==2){
        epsct=(2.0*pi)/Nr;
        fstream myfile;
        myfile.open("./results/outphi/profile.dat");
        for (i=0;i<=(Nr-1);i++){
            for (j=0;j<=(Nz-1);j++){
                myfile<<(i*delr)<<(j*delz)<<pA[i][j]<<pB[i][j]<<pC[i][j]<<pD[i][j]<<endl;
            }
        }
        myfile.close();
        myfile.open("./results/xyz.dat");
        ct=0.0;
        for (i=0;i<=Nr-1;i++){
            px=(((i-1)*delr+r_0))*cos(ct);
            py=(((i-1)*delr+r_0))*sin(ct);
            pz=(i-1)*delz;
            myfile<<ct<<(i-1)*delr<<(i-1)*delz;
        }
        myfile.close();
        myfile.open("./results/ABCD.dat");
        for (ii=0;ii<=(Nr-1);ii++){
            for (i=0;i<=(Nr-1);i++){
                for(j=0;j<=(Nr-1);j++){
                    myfile<<pA[i][j]<<pB[i][j]<<pC[i][j]<<pD[i][j]<<endl;
                }
            }
        }
        myfile.close();
        
    
    }
    
    
    
}