double TDMA(const size_t N_in, const double DiagL_in[], const double Diag_in[], const double DiagU_in[],double b_in[]) {
    size_t i;
    
    /* Allocate scratch space. */
    double* c = (double*)malloc(sizeof(double) * N_in);
    double* d = (double*)malloc(sizeof(double) * N_in);
    //double* c = create_1d_double_array(N_in, "C");
    //double* d = create_1d_double_array(N_in, "D");
    
    
    if (!c) {
        cout<< "TDMA error"<< endl;
    }
    
    c[0] = DiagU_in[0] / Diag_in[0];
    d[0] = b_in[0] / Diag_in[0];
    
    /* loop from 1 to N - 1 inclusive */
    
    
    for (i = 1; i <= int(N_in-2); i++) {
        c[i]=DiagU_in[i]/(Diag_in[i]-(c[i-1]*DiagL_in[i-1]));
        b_in[i] = (b_in[i] - DiagL_in[i] * b_in[i-1]);
    }
    for (i=1; i<= int(N_in-1);i++){
        d[i]=(b_in[i]-(d[i-1]*DiagL_in[i-1]))/(Diag_in[i]-(c[i-1]*DiagL_in[i-1]));
    }
    
    b_in[int(N_in-1)]=d[int(N_in-1)];
    
    /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
    for (i = int(N_in - 1); i --> 0;){
        b_in[i] =d[i] - c[i] * b_in[i + 1];
    }
    /* free scratch space */
    free(c);free(d);
    
    
    return b_in[N_in];
}