

void cleanme(){
    ofstream myfile;
    myfile.open("results/profile1.dat", ofstream::trunc);
    if (myfile.is_open()){cout<<"results/profile1.dat good"<<endl;}
    myfile.close();
    myfile.open("results/xyz.dat", ofstream::trunc);
    if (myfile.is_open()){cout<<"results/xyz.dat good"<<endl;}
    myfile.close();
    myfile.open("results/ABCDE.dat", ofstream::trunc);
    if (myfile.is_open()){cout<<"results/ABCDE.dat good"<<endl;}
    myfile.close();
    myfile.open("results/profile4.dat",ofstream::trunc);
    if (myfile.is_open()){cout<<"results/profile4.dat good"<<endl;}
    myfile.close();
    myfile.open("results/phi1Dr.dat",ofstream::trunc);
    if (myfile.is_open()){cout<<"results/phi1Dr.dat good"<<endl;}
    myfile.close();
    myfile.open("results/phi1Dz.dat",ofstream::trunc);
    if (myfile.is_open()){cout<<"results/phi1Dz.dat good"<<endl;}
    myfile.close();
    myfile.open("results/phi.dat",ofstream::trunc);
    if (myfile.is_open()){cout<<"results/phi.dat good"<<endl;}
    myfile.close();
    myfile.open("results/omega.dat",ofstream::trunc);
    if (myfile.is_open()){cout<<"results/omega.dat good"<<endl;}
    myfile.close();
    myfile.open("results/a.dat",ofstream::trunc);
    if (myfile.is_open()){cout<<"results/a.dat good"<<endl;}
    myfile.close();
    myfile.open("results/b.dat",ofstream::trunc);
    if (myfile.is_open()){cout<<"results/b.dat good"<<endl;}
    myfile.close();
    
}


void clean_data(){
    if (disk==1) {
        ofstream myfile;
        myfile.open("results/fE_disk.dat", ofstream::trunc);
        myfile.close();
        myfile.open("results/omega.disk",ofstream::trunc);
        myfile.close();
    }
    if (bilayer==1) {
        ofstream myfile;
        myfile.open("results/fE_bilayer.dat", ofstream::trunc);
        myfile.close();
        myfile.open("results/omega.bilayer",ofstream::trunc);
        myfile.close();
    }
}

void show_data(double Area){
    cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
    
    if (bilayer==1){
        cout<<"Shape chosen is: Bilayer"<< endl;
        //cout<<"order_parameter"<<OP<< endl;
        cout<<"muABC="<<muABC<<"   "<<"muD="<< muD<< endl;
        cout<<"fE="<<((fE-fE_hom)*Vol)/Area<<endl;
    }
    
    else if (disk==1){
        cout<<"Shape chosen is: Disk"<< endl;
        //cout<<"order_parameter"<<OP<< endl;
        cout<<"muABC="<<muABC<<"   "<<"muD="<< muD<< endl;
        cout<<"fE="<<((fE-fE_hom)*Vol)/Area<<endl;
    }
    cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
}

void save_data(double Area, double phiA, double phiB, double phiC, double phiD, double Tip_R){
    if (disk==1){
        fstream myfile;
        myfile.open("results/fE_disk.dat");
        myfile<<ND<<"   "<<((fE-fE_hom)*Vol)/Area<<"    "<<((phiA+phiB+phiD)-phiABC_hom)*Vol/Area
        <<" "<<(r_0+Tip_R)<<"   "<<(phiA+phiB+phiC)<<"  "<<(phiD)<<"    "<<Area<<"  "<<muABC<<" "<<muD<< endl;
        myfile.close();
    }
    
    if (bilayer==1){
        fstream myfile;
        myfile.open("results/fE_bilayer.dat");
        myfile<<ND<<"   "<<((fE-fE_hom)*Vol)/Area<<"    "<<((phiA+phiB+phiD)-phiABC_hom)*Vol/Area
        <<" "<< (r_0+Tip_R)<<"  "<<(phiA+phiB+phiC)<<"  "<<(phiD)<<"    "<<Area<<"  "<<muABC<<" "<<muD<< endl;
        myfile.close();
    }
}

void write_data(){
    int i,j;
    fstream myfile;
    
    myfile.open("results/phi.dat");
    for (i=0;i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            myfile <<(i*delr)<<" "<<(j*delz)<<" "<<pA[i][j]<<" "<<pB[i][j]<<" "<<pC[i][j]<<" "<<pD[i][j]<<endl;}}
    myfile.close();
    
    myfile.open("results/phi1Dr.dat");
    for (i=0;i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            myfile <<(i*delr)<<" "<<(j*delz)<<" "<<pA[i][Nz/2]<<" "<<pB[i][Nz/2]<<" "<<pC[i][Nz/2]<<" "<<pD[i][Nz/2]<<endl;}}
    myfile.close();
    
    myfile.open("results/phi1Dz.dat");
    for (i=0;i<=Nr-1;i++){
        for (j=0;j<=Nz-1;j++){
            myfile <<(i*delr)<<" "<<(j*delz)<<" "<<pA[Nr/2][j]<<" "<<pB[Nr/2][j]<<" "<<pC[Nr/2][j]<<" "<<pD[Nr/2][j]<<endl;}}
    myfile.close();
    
    if (disk==1){
        myfile.open("results/omega.disk");
        for (i=0;i<=Nr-1;i++){
            for (j=0;j<=Nz-1;j++){
                myfile<<i*delr<<" "<<j*delz<<" "<<wA[i][j]<<" "<<wB[i][j]<<" "<<wC[i][j]<<" "<<wD[i][j]<<endl;}}
        myfile.close();}
    
    if (bilayer==1){
        myfile.open("results/omega.bilayer");
        for (i=0;i<=Nr-1;i++){
            for (j=0;j<=Nz-1;j++){
                myfile<<i*delr<<" "<<j*delz<<" "<<wA[i][j]<<" "<<wB[i][j]<<" "<<wC[i][j]<<" "<<wD[i][j]<<endl;}}
        myfile.close();}
    
}

void profile(int msg){
    int  i,j,ii;
    double ct,epsct;
    double px,py,pz;
    string outfile,orderpar;
    
    if (msg==1){
        epsct=(2.0*pi)/Nr;
        fstream myfile;
        myfile.open("results/profile1.dat");
        for (i=0;i<=(Nr-1);i++){
            for (j=0;j<=(Nz-1);j++){
                myfile<<(i*delr)<<" "<<(j*delz)<<" "<<pA[i][j]<<" "<<pB[i][j]<<" "<<pC[i][j]<<" "<<pD[i][j]<<endl;
            }
        }
        myfile.close();
    }
    else if (msg==2){
        epsct=(2.0*pi)/Nr;
        fstream myfile;
        myfile.open("results/outphi/profile.dat");
        for (i=0;i<=(Nr-1);i++){
            for (j=0;j<=(Nz-1);j++){
                myfile<<(i*delr)<<" "<<(j*delz)<<" "<<pA[i][j]<<" "<<pB[i][j]<<" "<<pC[i][j]<<" "<<pD[i][j]<<endl;
            }
        }
        myfile.close();
        myfile.open("results/xyz.dat");
        ct=0.0;
        for (i=0;i<=Nr-1;i++){
            px=(((i-1)*delr+r_0))*cos(ct);
            py=(((i-1)*delr+r_0))*sin(ct);
            pz=(i-1)*delz;
            myfile<<ct<<(i-1)*delr<<" "<<(i-1)*delz;
        }
        myfile.close();
        myfile.open("results/ABCD.dat");
        for (ii=0;ii<=(Nr-1);ii++){
            for (i=0;i<=(Nr-1);i++){
                for(j=0;j<=(Nr-1);j++){
                    myfile<<pA[i][j]<<" "<<pB[i][j]<<" "<<pC[i][j]<<" "<<pD[i][j]<<endl;
                }
            }
        }
        myfile.close();
    }
}



