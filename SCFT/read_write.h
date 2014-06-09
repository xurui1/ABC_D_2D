int read(){
    // Reading in the parameters
    int n;
    n=0;
    fstream myfile;
    myfile.open("data.txt");
    if (myfile.is_open()){
        n=n+1;
        while (!myfile.eof()){
            myfile >> D_r>> D_z;
            myfile >> Nr>> Nz;
            myfile >> NA>>NB>>NC;
            myfile >> ND;
            myfile >> sig>>sig2;
            myfile >> xAB>>xAC>>xAD;
            myfile >> xBC>>xBD;
            myfile >> xCD;
            myfile >> muABC>>muD;
        }
        myfile.close();
        cout<<D_r<< D_z<< Nr<< Nz<< NA<< NB << NC<< ND << sig<< sig2<< xAB<< xAC<< xAD<< xBC<< xBD<< xCD<< muABC<< muD<< endl;
    }
    
    else {
        cout<< "Unable to read file data.dat "<<n<< endl;
        exit(1);
    }
    
        return read();
}

void cleanme(){
    ofstream myfile;
    myfile.open("./results/profile1.dat", ofstream::trunc);
    myfile.close();
    myfile.open("./results/xyz.dat", ofstream::trunc);
    myfile.close();
    myfile.open("./results/ABCDE.dat", ofstream::trunc);
    myfile.close();
    myfile.open("./results/profile4.dat",ofstream::trunc);
    myfile.close();
    myfile.open("./results/phi1Dr.dat",ofstream::trunc);
    myfile.close();
    myfile.open("./results/phi1Dz.dat",ofstream::trunc);
    myfile.close();
    myfile.open("./results/phi.dat",ofstream::trunc);
    myfile.close();
    myfile.open("./results/omega.dat",ofstream::trunc);
    myfile.close();
    myfile.open("./results/a.dat",ofstream::trunc);
    myfile.close();
    myfile.open("./results/b.dat",ofstream::trunc);
    myfile.close();
    
}


void clean_data(){
    if (disk==1) {
        ofstream myfile;
        myfile.open("./results/fE_disk.dat", std::ofstream::trunc);
        myfile.close();
    }
    if (bilayer==1) {
        ofstream myfile;
        myfile.open("./results/fE_bilayer.dat", std::ofstream::trunc);
        myfile.close();
    }
}

void show_data(double Area){
    cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
    
    if (bilayer==1){
        cout<<"Shape chosen is: Bilayer"<< endl;
        //cout<<"order_parameter"<<OP<< endl;
        cout<<"muABC="<<muABC<<"   "<<"muD="<< muD<< endl;
        cout<<"fE="<<((fE-fE_hom)*Vol)/Area;
    }
    
    else if (disk==1){
        cout<<"Shape chosen is: Disk"<< endl;
        //cout<<"order_parameter"<<OP<< endl;
        cout<<"muABC="<<muABC<<"   "<<"muD="<< muD<< endl;
        cout<<"fE="<<((fE-fE_hom)*Vol)/Area;
    }
    cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
}

void save_data(double Area, double phiA, double phiB, double phiC, double phiD, double Tip_R){
    if (disk==1){
        fstream myfile;
        myfile.open("./results/fE_disk.dat");
        myfile<<ND<<((fE-fE_hom)*Vol)/Area<<((phiA+phiB+phiD)-phiABC_hom)*Vol/Area<< (r_0+Tip_R)<<
        (phiA+phiB+phiC)<<(phiD)<<Area<<muABC<<muD<< endl;
        myfile.close();
    }
    
    if (bilayer==1){
        fstream myfile;
        myfile.open("./results/fE_bilayer.dat");
        myfile<<ND<<((fE-fE_hom)*Vol)/Area<<((phiA+phiB+phiD)-phiABC_hom)*Vol/Area<< (r_0+Tip_R)<<
        (phiA+phiB+phiC)<<(phiD)<<Area<<muABC<<muD<< endl;
        myfile.close();
    }
}



