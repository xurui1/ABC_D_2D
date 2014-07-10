/******************************Calculate chain partition function**********************************/
double Q_partition(){
    int i,j;
    Q_ABC=0.0;
    Q_D=0.0;
    
    for (i=1; i<=int(Nr-2); i++){
        for (j=1; j<=int(Nz-2); j++){
            Q_ABC=Q_ABC+(qdagB[i][j][NB])*((i*delr)+(r_0))*delr*delz;
            Q_D=Q_D+(qD[i][j][0])*((i*delr)+(r_0))*delr*delz;
        }
    }
    i=0;
    for (j=1; j<=int(Nz-2); j++){
        Q_ABC=Q_ABC+0.5*(qdagB[i][j][NB])*((i*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*((i*delr)+(r_0))*delr*delz;
    }
    i=int(Nr-1);
    for (j=1; j<=int(Nz-2); j++){
        Q_ABC=Q_ABC+0.5*(qdagB[i][j][NB])*((i*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*((i*delr)+(r_0))*delr*delz;
    }
    j=0;
    for (i=1; i<=int(Nr-2); i++){
        Q_ABC=Q_ABC+0.5*(qdagB[i][j][NB])*((i*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*((i*delr)+(r_0))*delr*delz;
    }
    j=int(Nz-1);
    for (i=1; i<=int(Nr-2); i++){
        Q_ABC=Q_ABC+0.5*(qdagB[i][j][NB])*((i*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*((i*delr)+(r_0))*delr*delz;
    }
    
    Q_ABC=Q_ABC+0.25*delr*delz*(qdagB[0][0][NB]*(double(0.0)*delr+r_0)+
                                qdagB[0][int(Nz-1)][NB]*(double(0.0)*delr+r_0)+
                                qdagB[int(Nr-1)][0][NB]*(double(Nr-1)*delr+(r_0))+
                                qdagB[int(Nr-1)][int(Nz-1)][NB]*(double(Nr-1)*delr+(r_0)));
    
    Q_D=Q_D+0.25*delr*delz*(qD[0][0][ND]*(double(0.0)*delr+r_0)+
                            qD[0][int(Nz-1)][ND]*(double(0.0)*delr+r_0)+
                            qD[int(Nr-1)][0][ND]*(double(Nr-1)*delr+(r_0))+
                            qD[int(Nr-1)][int(Nz-1)][ND]*(double(Nr-1)*delr+(r_0)));

    
    Q_ABC=Q_ABC*2.0*pi;
    Q_D=Q_D*2.0*pi;
    Q_ABC=Q_ABC*2.0/Vol;
    Q_D=Q_D*2.0/Vol;
    
    //cout<<"Q_ABC: "<<Q_ABC<<" Q_D:"<<Q_D<<endl;
    
    return Q_D;
}