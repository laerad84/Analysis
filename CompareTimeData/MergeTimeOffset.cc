void MergeTimeOffset(){

  std::string ANALYSISLIB = std::getenv("ANALYSISLIB");
  gSystem->Load(Form("%s/lib/libAnalysisLib.so",ANALYSISLIB.c_str()));
  gSystem->Load("/home/jwlee/local/Analysis/DrawLib/lib/libDrawLib.so");
  CsIPoly* csi = new CsIPoly("DeltaT","DeltaT");

  IDHandler* handler = new IDHandler();

  std::ifstream ifsPi0Offet(Form("%s/Data/TimeOffset/Pi0Peak.dat",ANALYSISLIB.c_str()));
  if( !ifsPi0Offet.is_open() ){
    std::cout<< "Pi0Peak.dat is not existed" << std::endl; 
    return;
  }

  Double_t TimeOffsetTotal[2716]={0};
  Double_t TimeOffsetTotalSigma[2716] = {0xFFFF};

  Double_t TimeOffsetPi0[2716]={0};
  Double_t TimeOffsetPi0Sigma[2716] ={0xFFFF};
  Double_t TimeOffsetCrystalPosition[2716]={0};
  Double_t TimeOffsetCosmic[2716] = {0};
  Double_t TimeOffsetCosmicSigma[2716] = {0xFFFF};
  
  Double_t ResolutionT[2716]={0xFFFF};
  Int_t tID;
  Double_t tOffset;
  Double_t tOffsetSigma;

  while( ifsPi0Offet >> tID >> tOffset >> tOffsetSigma ){
    TimeOffsetPi0Sigma[tID] = tOffsetSigma;
    if( TimeOffsetPi0Sigma[tID] > 0 && TimeOffsetPi0Sigma[tID] < 0xFFFF){
      TimeOffsetPi0[tID] = tOffset;
    }else{
      TimeOffsetPi0[tID] = 0;
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


  for( int i = 0 ;i< 2716; i++){
    if( TimeOffsetPi0Sigma[i]    < 0xFFFF && TimeOffsetPi0Sigma[i] >= 0  &&
	TimeOffsetCosmicSigma[i] < 0xFFFF && TimeOffsetCosmicSigma[i] >= 0 ){
      TimeOffsetTotal[i]       = TimeOffsetCosmic[i] + TimeOffsetPi0[i] ;//- TimeOffsetCrystalPosition[i];
      TimeOffsetTotalSigma[i]  = TimeOffsetCosmicSigma[i];
    }else{
      TimeOffsetTotal[i]       = 0; 
      TimeOffsetTotalSigma[i]  = 0xFFFF;
    }    
  }

  std::ofstream ofs("TimeOffset_with_cosmic.dat");
  for( int i = 0; i< 2716; i++){
    ofs << i << "\t" << TimeOffsetTotal[i]  << "\t" << TimeOffsetTotalSigma[i] << "\n";
  }

  ofs.close();
}
