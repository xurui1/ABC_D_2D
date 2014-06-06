int TDMA(const size_t N_in, const double DiagL_in[], const double Diag_in[], const double DiagU_in[],double b_in[]) {
    size_t i;
    
    /* Allocate scratch space. */
    double* c = (double*)malloc(sizeof(double) * N_in);
    double* d = (double*)malloc(sizeof(double) * N_in);
    
    if (!c) {
        cout<< "TDMA error"<< endl;
    }
    
    c[0] = DiagU_in[0] / Diag_in[0];
    d[0] = b_in[0] / Diag_in[0];
    
    /* loop from 1 to N - 1 inclusive */
    for (i = 1; i < N_in; i++) {
        double m = 1.0 / (Diag_in[i] - (DiagL_in[i] * c[i - 1]));
        c[i] = c[i] * m;
        b_in[i] = (b_in[i] - DiagL_in[i] * b_in[i - 1]) * m;
    }
    
    /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
    for (i = N_in - 1; i-- > 0; )
        b_in[i] = b_in[i] - c[i] * b_in[i + 1];
    
    /* free scratch space */
    free(c);
    free(d);
    
    return TDMA(N_in, DiagL_in,Diag_in,DiagU_in,b_in);
    
}