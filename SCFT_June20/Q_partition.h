/******************************Calculate chain partition function**********************************/
double Q_partition(){
    int i,j;
    Q_ABC=0.0;
    Q_D=0.0;
    
    for (i=1; i<=int(Nr-2); i++){
        for (j=1; j<=int(Nz-2); j++){
            Q_ABC=Q_ABC+(qdagA[i][j][int(NA-1)])*((int(i-1)*delr)+(r_0))*delr*delz;
            Q_D=Q_D+(qD[i][j][0])*((int(i-1)*delr)+(r_0))*delr*delz;
        }
    }
    i=0;
    for (j=1; j<=int(Nz-2); j++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][int(NA-1)])*((int(i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*((int(i-1)*delr)+(r_0))*delr*delz;
    }
    i=Nr-1;
    for (j=1; j<=int(Nz-2); j++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][int(NA-1)])*((int(i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*((int(i-1)*delr)+(r_0))*delr*delz;
    }
    j=0;
    for (i=1; i<=int(Nr-2); i++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][int(NA-1)])*((int(i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*((int(i-1)*delr)+(r_0))*delr*delz;
    }
    j=Nz-1;
    for (i=1; i<=int(Nr-2); i++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][int(NA-1)])*((int(i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*((int(i-1)*delr)+(r_0))*delr*delz;
    }
    
    Q_ABC=Q_ABC+0.25*delr*delz*(qdagA[0][0][int(NA-1)]*(float(1-1)*delr+r_0)+
                                qdagA[0][int(Nz-1)][int(NA-1)]*(float(1-1)*delr+r_0)+
                                qdagA[int(Nr-1)][0][int(NA-1)]*(float(Nr-1)*delr+(r_0))+
                                qdagA[int(Nr-1)][int(Nz-1)][int(NA-1)]*(float(Nr-1)*delr+(r_0)));
    
    Q_D=Q_D+0.25*delr*delz*(qD[0][0][int(ND-1)]*(float(1-1)*delr+r_0)+
                            qD[0][int(Nz-1)][int(ND-1)]*(float(1-1)*delr+r_0)+
                            qD[int(Nr-1)][0][int(ND-1)]*(float(Nr-1)*delr+(r_0))+
                            qD[int(Nr-1)][int(Nz-1)][int(ND-1)]*(float(Nr-1)*delr+(r_0)));

    
    Q_ABC=Q_ABC*2.0*pi;
    Q_D=Q_D*2.0*pi;
    Q_ABC=Q_ABC*2.0/Vol;
    Q_D=Q_D*2.0/Vol;
    
    cout<<"Q_ABC: "<<Q_ABC<<" Q_D:"<<Q_D<<endl;
    
    return Q_D;
}