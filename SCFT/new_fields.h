/**********************************Calculate new omega fields***************************************/

void new_fields (){
    int i,j;
    Conv_w=0.0;
    Conv_p=0.0;
    for (i=0;i<=int(Nr-1);i++) {
        for (j=0;j<=int(Nz-1);j++){
            dwA[i][j]= xAB*(pB[i][j])+xAC*(pC[i][j])+xAD*(pD[i][j])-(eta2[i][j])+(eta[i][j])-(wA[i][j]);
            dwB[i][j]= xAB*(pB[i][j])+xBC*(pC[i][j])+xBD*(pD[i][j])-(eta2[i][j])+(eta[i][j])-(wB[i][j]);
            dwC[i][j]= xAC*(pA[i][j])+xBC*(pB[i][j])+xCD*(pD[i][j])-(eta2[i][j])+(eta[i][j])-(wC[i][j]);
            dwD[i][j]= xAD*(pA[i][j])+xBD*(pB[i][j])+xCD*(pC[i][j])-(eta2[i][j])+(eta[i][j])-(wD[i][j]);
            dpp[i][j]=1.0-(pA[i][j]+pB[i][j]+pC[i][j]+pD[i][j]);
            
            Conv_w=Conv_w+abs((dwA[i][j])+(dwB[i][j])+(dwC[i][j])+(dwD[i][j]));
            Conv_p=Conv_p+abs(dpp[i][j]);
            
            // updating the omega condition
            wA[i][j]=(wA[i][j])+sig*(dwA[i][j])-sig2*(dpp[i][j]);
            wB[i][j]=(wB[i][j])+sig*(dwB[i][j])-sig2*(dpp[i][j]);
            wC[i][j]=(wC[i][j])+sig*(dwC[i][j])-sig2*(dpp[i][j]);
            wD[i][j]=(wD[i][j])+sig*(dwD[i][j])-sig2*(dpp[i][j]);
            
        }
    }
    Conv_w=Conv_w/(Nr*Nz);
    Conv_p=Conv_p/(Nr*Nz);
    }