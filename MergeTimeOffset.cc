void MergeTimeOffset(){

  std::string::ROOTFILE_COSMIC = std::getenv("ROOTFILE_COSMIC");
  std::ifsstream ifsPi0Offet(Form("%s/Pi0Peak.dat",ROOTFILE_COSMIC.c_str()));
  if( !ifsPi0Offet.is_open() ){
    std::cout<< "Pi0Peak.dat is not existed" << std::endl; 
    return;
  }
  Double_t TimeOffsetTotal[2716]={0};
  Double_t TimeOffsetTotalSigma[2716] = {0xFFFF};
  Double_t TimeOffset[2716]={0};
  Double_t TimeOffsetSigma[2716] ={0xFFFF};
  Doubel_t TimeOffsetCrystalPosition[2716]={0};
  Double_t TimeOffsetCosmic[2716] = {0};
  Double_t TimeOffsetCosmicSigma[2716] = {0xFFFF};

  Double_t ResolutionT[2716]={0xFFFF};
  Int_t tID;
  Double_t tOffset;
  Double_t tOffsetSigma;
  while( ifsPi0Offet >> tID >> tOffset >> tOffsetSigma ){
    TimeOffsetSigma[tID] = tOffsetSigma;
    if( TimeOffsetSigma > 0 && TimeOffsetSigma < 0xFFFF){
      TimeOffset[tID] = tOffset;
    }else{
      TimeOffset[tID] = 0;
    }
  }
  for( int i = 0; i< 2716; i++){
    double x,y;
    handler->GetMetricPosition( i,x,y);
    TimeOffsetCrystalPosition[i] = (TMath::Sqrt(2624*2624+x*x+y*y)-2624)/299.7;
  }
  
  std::ifstream ifsCosmicOffset("CosmicOutTimeDeltaResolution_19.dat");
  while( ifsCosmicOffset >> tID >> tOffset >> tOffsetSigma ){
    if( tOffsetSigma < 0xFFFF ){
      TimeOffsetCosmic[tID] = tOffset;
      TimeOffsetCosmicSigma[tID] = tOffsetSigma;
    }else{
      TimeOffsetCosmic[tID] = 0;
      TimeOffsetCosmicSigma[tID] = 0xFFFF;
    }
  }
}
