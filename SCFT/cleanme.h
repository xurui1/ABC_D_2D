
/**************************************Clean files****************************************************/
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
