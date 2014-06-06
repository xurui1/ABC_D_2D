/**********************************Calculate new omega fields***************************************/

int new_fields (){
    int i,j;
    Conv_w=0.0;
    Conv_p=0.0;
    for (i=1;i<=Nr;i++) {
        for (j=1;j<=Nz;j++){
            **dwA= xAB*(**pB)+xAC*(**pC)+xAD*(**pD)-(**eta2)+(**eta)-(**wA);
            **dwB= xAB*(**pB)+xBC*(**pC)+xBD*(**pD)-(**eta2)+(**eta)-(**wB);
            **dwC= xAC*(**pA)+xBC*(**pB)+xCD*(**pD)-(**eta2)+(**eta)-(**wC);
            **dwD= xAD*(**pA)+xBD*(**pB)+xCD*(**pC)-(**eta2)+(**eta)-(**wD);
            **dpp=1.0-(**pA+**pB+**pC+**pD);
            
            Conv_w=Conv_w+abs((**dwA)+(**dwB)+(**dwC)+(**dwD));
            Conv_p=Conv_p+abs(**dpp);
            
            // updating the omega condition
            **wA=(**wA)+sig*(**dwA)-sig2*(**dpp);
            **wB=(**wB)+sig*(**dwB)-sig2*(**dpp);
            **wC=(**wC)+sig*(**dwC)-sig2*(**dpp);
            **wD=(**wD)+sig*(**dwD)-sig2*(**dpp);
            
        }
    }
    Conv_w=Conv_w/(Nr*Nz);
    Conv_p=Conv_p/(Nr*Nz);
    return new_fields();
}