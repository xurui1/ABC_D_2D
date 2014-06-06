/******************************Calculate chain partition function**********************************/
int Q_partition(){
    int i,j;
    Q_ABC=0.0;
    Q_D=0.0;
    
    for (i=2; i<=Nr-1; i++){
        for (j=2; j<=Nz-1; j++){
            Q_ABC=Q_ABC+(qdagA[i][j][NA])*(((i-1)*delr)+(r_0))*delr*delz;
            Q_D=Q_D+(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
        }
    }
    i=1;
    for (j=2; j<=Nz-1; j++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][NA])*(((i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
    }
    i=Nr;
    for (j=2; j<=Nz-1; j++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][NA])*(((i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
    }
    j=1;
    for (i=2; i<=Nr-1; i++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][NA])*(((i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
    }
    j=Nz;
    for (i=2; i<=Nr-1; i++){
        Q_ABC=Q_ABC+0.5*(qdagA[i][j][NA])*(((i-1)*delr)+(r_0))*delr*delz;
        Q_D=Q_D+0.5*(qD[i][j][0])*(((i-1)*delr)+(r_0))*delr*delz;
    }
    
    Q_ABC=Q_ABC*2.0*pi;
    Q_D=Q_D*2.0*pi;
    Q_ABC=Q_ABC*2.0/Vol;
    Q_D=Q_D*2.0/Vol;
    
    //clear memory
    destroy_3d_double_array(qdagA);
    destroy_3d_double_array(qD);
    
    return Q_partition();
}