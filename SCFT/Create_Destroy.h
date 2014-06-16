void Destroy(){
    destroy_2d_double_array(eta);
    destroy_2d_double_array(dpp);
    destroy_2d_double_array(eta2);
    
    destroy_2d_double_array(wA);
    destroy_2d_double_array(dwA);
    destroy_2d_double_array(pA);
    
    destroy_2d_double_array(wB);
    destroy_2d_double_array(dwB);
    destroy_2d_double_array(pB);
    
    destroy_2d_double_array(wC);
    destroy_2d_double_array(dwC);
    destroy_2d_double_array(pC);
    
    destroy_2d_double_array(wD);
    destroy_2d_double_array(dwD);
    destroy_2d_double_array(pD);
    
    destroy_3d_double_array(qA);
    destroy_3d_double_array(qB);
    destroy_3d_double_array(qC);
    destroy_3d_double_array(qD);
    
    destroy_2d_double_array(qA_0);
    destroy_2d_double_array(qB_0);
    destroy_2d_double_array(qC_0);
    destroy_2d_double_array(qD_0);
    
    destroy_3d_double_array(qdagA);
    destroy_3d_double_array(qdagB);
    destroy_3d_double_array(qdagC);
    
    destroy_2d_double_array(qdagA_0);
    destroy_2d_double_array(qdagB_0);
    destroy_2d_double_array(qdagC_0);
    
    destroy_1d_double_array(DiagALr);
    destroy_1d_double_array(DiagAUr);
    destroy_1d_double_array(DiagAr);
    destroy_1d_double_array(DiagALrdr);
    destroy_1d_double_array(DiagAUrdr);
    destroy_1d_double_array(DiagArdr);
    destroy_1d_double_array(DiagALz);
    destroy_1d_double_array(DiagAUz);
    destroy_1d_double_array(DiagAz);
    destroy_1d_double_array(DiagALzdz);
    destroy_1d_double_array(DiagAUzdz);
    destroy_1d_double_array(DiagAzdz);
    destroy_1d_double_array(bAr);
    destroy_1d_double_array(bAz);
    
    destroy_1d_double_array(DiagBLr);
    destroy_1d_double_array(DiagBUr);
    destroy_1d_double_array(DiagBr);
    destroy_1d_double_array(DiagBLrdr);
    destroy_1d_double_array(DiagBUrdr);
    destroy_1d_double_array(DiagBrdr);
    destroy_1d_double_array(DiagBLz);
    destroy_1d_double_array(DiagBUz);
    destroy_1d_double_array(DiagBz);
    destroy_1d_double_array(DiagBLzdz);
    destroy_1d_double_array(DiagBUzdz);
    destroy_1d_double_array(DiagBzdz);
    destroy_1d_double_array(bBr);
    destroy_1d_double_array(bBz);
    
    destroy_1d_double_array(DiagCLr);
    destroy_1d_double_array(DiagCUr);
    destroy_1d_double_array(DiagCr);
    destroy_1d_double_array(DiagCLrdr);
    destroy_1d_double_array(DiagCUrdr);
    destroy_1d_double_array(DiagCrdr);
    destroy_1d_double_array(DiagCLz);
    destroy_1d_double_array(DiagCUz);
    destroy_1d_double_array(DiagCz);
    destroy_1d_double_array(DiagCLzdz);
    destroy_1d_double_array(DiagCUzdz);
    destroy_1d_double_array(DiagCzdz);
    destroy_1d_double_array(bCr);
    destroy_1d_double_array(bCz);
    
    destroy_1d_double_array(DiagDLr);
    destroy_1d_double_array(DiagDUr);
    destroy_1d_double_array(DiagDr);
    destroy_1d_double_array(DiagDLrdr);
    destroy_1d_double_array(DiagDUrdr);
    destroy_1d_double_array(DiagDrdr);
    destroy_1d_double_array(DiagDLz);
    destroy_1d_double_array(DiagDUz);
    destroy_1d_double_array(DiagDz);
    destroy_1d_double_array(DiagDLzdz);
    destroy_1d_double_array(DiagDUzdz);
    destroy_1d_double_array(DiagDzdz);
    destroy_1d_double_array(bDr);
    destroy_1d_double_array(bDz);
}

void create(){
    eta =create_2d_double_array(Nr,Nz, "eta");             //incompressibility condition
    dpp =create_2d_double_array(Nr,Nz, "dpp");             //not sure what dpp is used for (update w fields?)
    eta2 =create_2d_double_array(Nr,Nz, "eta2");           //pinning condition
    
    wA=create_2d_double_array(Nr,Nz, "wA");                //interaction potential of chain A
    dwA=create_2d_double_array(Nr,Nz, "dwA");              //differential of wA
    pA=create_2d_double_array(Nr,Nz, "pA");                //probability amplitude of chain A
    
    wB=create_2d_double_array(Nr,Nz, "wB");                //interaction potential of chain B
    dwB=create_2d_double_array(Nr,Nz, "dwB");              //differential of wB
    pB=create_2d_double_array(Nr,Nz, "pB");                //probability amplitude of chain B
    
    wC=create_2d_double_array(Nr,Nz, "wC");                //interaction potential of chain C
    dwC=create_2d_double_array(Nr,Nz, "dwC");              //differential of wC
    pC=create_2d_double_array(Nr,Nz, "pC");                //probability amplitude of chain C
    
    wD=create_2d_double_array(Nr,Nz, "wD");                //interaction potential of chain D
    dwD=create_2d_double_array(Nr,Nz, "dwD");              //differential of wD
    pD=create_2d_double_array(Nr,Nz, "pD");                //probability amplitude of chain D
    
    
    qA=create_3d_double_array(Nr,Nz,Ns,"qA");
    qA_0=create_2d_double_array(Nr,Nz, "qA_0");
    qdagA=create_3d_double_array(Nr,Nz,Ns,"qdagA");
    qdagA_0=create_2d_double_array(Nr,Nz, "qdagA_0");
    bAr=create_1d_double_array(Nz, "bAr");
    bAz=create_1d_double_array(Nr, "bAz");
    DiagAr=create_1d_double_array(Nz, "DiagAr");
    DiagAUr=create_1d_double_array(((int)Nz-1), "DiagAUr");
    DiagALr=create_1d_double_array(((int)Nz-1), "DiagALr");
    DiagAz=create_1d_double_array(Nr, "DiagAz");
    DiagAUz=create_1d_double_array(((int)Nr-1), "DiagAUz");
    DiagALz=create_1d_double_array(((int)Nr-1), "DiagALz");
    DiagArdr=create_1d_double_array(Nz, "DiagArdr");
    DiagAUrdr=create_1d_double_array(((int)Nz-1), "DiagAUrdr");
    DiagALrdr=create_1d_double_array(((int)Nz-1), "DiagALrdr");
    DiagAzdz=create_1d_double_array(Nr, "DiagAzdz");
    DiagAUzdz=create_1d_double_array(((int)Nr-1), "DiagAUzdz");
    DiagALzdz=create_1d_double_array(((int)Nr-1), "DiagALzdz");
    
    
    qB=create_3d_double_array(Nr,Nz,Ns,"qB");
    qB_0=create_2d_double_array(Nr,Nz, "qB_0");
    qdagB=create_3d_double_array(Nr,Nz,Ns,"qdagB");
    qdagB_0=create_2d_double_array(Nr,Nz, "qdagB_0");
    bBr=create_1d_double_array(Nz, "bBr");
    bBz=create_1d_double_array(Nr, "bBz");
    DiagBr=create_1d_double_array(Nz, "DiagBr");
    DiagBUr=create_1d_double_array(((int)Nz-1), "DiagBUr");
    DiagBLr=create_1d_double_array(((int)Nz-1), "DiagBLr");
    DiagBz=create_1d_double_array(Nr, "DiagBz");
    DiagBUz=create_1d_double_array(((int)Nr-1), "DiagBUz");
    DiagBLz=create_1d_double_array(((int)Nr-1), "DiagBLz");
    DiagBrdr=create_1d_double_array(Nz, "DiagBrdr");
    DiagBUrdr=create_1d_double_array(((int)Nz-1), "DiagBUrdr");
    DiagBLrdr=create_1d_double_array(((int)Nz-1), "DiagBLrdr");
    DiagBzdz=create_1d_double_array(Nr, "DiagBzdz");
    DiagBUzdz=create_1d_double_array(((int)Nr-1), "DiagBUzdz");
    DiagBLzdz=create_1d_double_array(((int)Nr-1), "DiagBLzdz");
    
    
    qC=create_3d_double_array(Nr,Nz,Ns,"qC");
    qC_0=create_2d_double_array(Nr,Nz, "qC_0");
    qdagC=create_3d_double_array(Nr,Nz,Ns,"qdagC");
    qdagC_0=create_2d_double_array(Nr,Nz, "qdagC_0");
    bCr=create_1d_double_array(Nz, "bCr");
    bCz=create_1d_double_array(Nr, "bCz");
    DiagCr=create_1d_double_array(Nz, "DiagCr");
    DiagCUr=create_1d_double_array(((int)Nz-1), "DiagCUr");
    DiagCLr=create_1d_double_array(((int)Nz-1), "DiagCLr");
    DiagCz=create_1d_double_array(Nr, "DiagCz");
    DiagCUz=create_1d_double_array(((int)Nr-1), "DiagCUz");
    DiagCLz=create_1d_double_array(((int)Nr-1), "DiagCLz");
    DiagCrdr=create_1d_double_array(Nz, "DiagCrdr");
    DiagCUrdr=create_1d_double_array(((int)Nz-1), "DiagCUrdr");
    DiagCLrdr=create_1d_double_array(((int)Nz-1), "DiagCLrdr");
    DiagCzdz=create_1d_double_array(Nr, "DiagCzdz");
    DiagCUzdz=create_1d_double_array(((int)Nr-1), "DiagCUzdz");
    DiagCLzdz=create_1d_double_array(((int)Nr-1), "DiagCLzdz");
    
    
    qD=create_3d_double_array(Nr,Nz,Ns,"qD");
    qD_0=create_2d_double_array(Nr,Nz, "qD_0");
    bDr=create_1d_double_array(Nz, "bDr");
    bDz=create_1d_double_array(Nr, "bDz");
    DiagDr=create_1d_double_array(Nz, "DiagDr");
    DiagDUr=create_1d_double_array(((int)Nz-1), "DiagDUr");
    DiagDLr=create_1d_double_array(((int)Nz-1), "DiagDLr");
    DiagDz=create_1d_double_array(Nr, "DiagDz");
    DiagDUz=create_1d_double_array(((int)Nr-1), "DiagDUz");
    DiagDLz=create_1d_double_array(((int)Nr-1), "DiagDLz");
    DiagDrdr=create_1d_double_array(Nz, "DiagDrdr");
    DiagDUrdr=create_1d_double_array(((int)Nz-1), "DiagDUrdr");
    DiagDLrdr=create_1d_double_array(((int)Nz-1), "DiagDLrdr");
    DiagDzdz=create_1d_double_array(Nr, "DiagDzdz");
    DiagDUzdz=create_1d_double_array(((int)Nr-1), "DiagDUzdz");
    DiagDLzdz=create_1d_double_array(((int)Nr-1), "DiagDLzdz");
    
}