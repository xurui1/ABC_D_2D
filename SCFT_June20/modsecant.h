/******************************************************************************************/
/**************************************The Secant method***********************************/

int secant(){
    
    
    int s,s2,ii,msg,msg1,i,j;
    double mu1,mu2,mu3;
    double fE1,fE2,fE3;
    int initial;
    
    
    
    mu1=muD;
    mu2=muD+0.01;
    
    bilayer=1;
    msg1=1;
    
    rand_field(initial);
    ii=1;
    msg=0;
    
    for (;;){
        for (i=0;i<=int(Nr-1);i++){
            for (j=0;j<=int(Nz-1);j++){
                eta2[i][j]=0;
                eta[i][j]=0;
            }
        }
        
        
        fE_homo( );
        Vol=2.0*pi*(0.5*(Nz-1)*delz*(pow(((Nr-1)*delr+r_0),2)-pow((r_0),2)));
        
        s2=0;
        fE_old=0.0;
        cleanme();
        
        for (s=1;s<=10000;s++){
            qA_forward();
            qB_forward();
            qC_forward();
            qD_forward();
            
            qdagA_forward();
            qdagB_forward();
            qdagC_forward();
            Q_partition( );
            phi( );
            pressure( );
            FreeEnergy( );
            new_fields( );
            
            
            cout<< "Field convergence: "<<Conv_w<<"Free deltaE: "<< fE-fE_hom<< muD << ii << endl;
            write_data();
            
            if (Conv_p<1.0e-4 and Conv_w<1.0e-4 and dfffE<1.0e-4){break;}
        }
        
        
        if (ii==1){
            fE1=fE-fE_hom;
            muD=mu2;
        }
        
        if (ii==2) {
            fE2=fE-fE_hom;
            mu3=mu2-(fE2*((mu2-mu1)/(fE2-fE1)));
            muD=mu3;
        }
        
        if (ii==3){
            fE3=fE-fE_hom;
            if (abs(fE3)<1.0e-5){
                break;}
            else{
                mu1=mu2;
                mu2=mu3;
                muD=mu1;
                ii=1;
                msg=1;
            }
        }
        cout<<muD<<endl;
        
        if (msg==0){
            ii=ii+1;
        }
        msg=0;
        
        if (ii>3){cout<<"something is wrong in the secant method mod"<<endl;}
        
        break;
        
    }
    
    if ((msg1==1) and(disk==0)){
        disk=1;
        bilayer=0;
    }
    
    return secant();
}