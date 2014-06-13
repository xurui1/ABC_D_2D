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