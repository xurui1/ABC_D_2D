/******************************Calculate chain partition function**********************************/
void Q_partition(){
    int i,j;
    Q_ABC=0.0;
    Q_D=0.0;
    
    for (i=1; i<=Nr-2; i++){
        for (j=1; j<=Nz-2; j++){
            Q_ABC=Q_ABC+(qdagA[i][j][NA-1])*(((i-1)*delr)+(r_0))*delr*delz;
            Q_D=Q_D+(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
        }
    }
    i=0;
    for (j=1; j<=Nz-2; j++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][NA-1])*(((i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
    }
    i=Nr-1;
    for (j=1; j<=Nz-2; j++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][NA-1])*(((i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
    }
    j=0;
    for (i=1; i<=Nr-2; i++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][NA-1])*(((i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
    }
    j=Nz-1;
    for (i=1; i<=Nr-2; i++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][NA-1])*(((i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
    }
    
    Q_ABC=Q_ABC*2.0*pi;
    Q_D=Q_D*2.0*pi;
    Q_ABC=Q_ABC*2.0/Vol;
    Q_D=Q_D*2.0/Vol;
    
    
}