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
    
    qA=create_3d_double_array(Nr,Nz,int(Ns+1),"qA");
    qA_0=create_2d_double_array(Nr,Nz, "qA_0");
    qdagA=create_3d_double_array(Nr,Nz,int(Ns+1),"qdagA");
    qdagA_0=create_2d_double_array(Nr,Nz, "qdagA_0");
    
    qB=create_3d_double_array(Nr,Nz,int(Ns+1),"qB");
    qB_0=create_2d_double_array(Nr,Nz, "qB_0");
    qdagB=create_3d_double_array(Nr,Nz,int(Ns+1),"qdagB");
    qdagB_0=create_2d_double_array(Nr,Nz, "qdagB_0");
    
    qC=create_3d_double_array(Nr,Nz,int(Ns+1),"qC");
    qC_0=create_2d_double_array(Nr,Nz, "qC_0");
    qdagC=create_3d_double_array(Nr,Nz,int(Ns+1),"qdagC");
    qdagC_0=create_2d_double_array(Nr,Nz, "qdagC_0");
    
    qD=create_3d_double_array(Nr,Nz,int(Ns+1),"qD");
    qD_0=create_2d_double_array(Nr,Nz, "qD_0");
}