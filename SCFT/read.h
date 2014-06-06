int read(){
    // Reading in the parameters
    int n;
    n=0;
    fstream myfile;
    myfile.open("./data.dat");
    if (myfile.is_open()){
        n=n+1;
        while (!myfile.eof()){
            myfile >> D_r>>D_z;
            myfile >> Nr>>Nz;
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
    
    else cout<< "Unable to read file, attempt "<<n<< endl;
    
        return read();
}
