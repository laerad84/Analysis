void MergeResult(){

  std::ifstream ifsPi0("Pi0Peak.dat");
  std::ifstream ifsCsm("CosmicOut_V4_TimeDeltaResolution.dat");
  
  if( !ifsPi0.is_open() || !ifsCsm.is_open() ){
    std::cerr <<" File Open Error " << std::endl;
    return ;
  }
  int tmpID; 
  Double_t tmpTime;
  Double_t tmpResolution;
  
  Double_t Pi0Time[2716] = {0};
  Double_t CSMTime[2716] = {0};
  Double_t Pi0Res[2716] = {0};
  Double_t CSMRes[2716] = {0};

  while( ifsPi0 >> tmpID >> tmpTime >> tmpResolution ){
    Pi0Time[tmpID] = tmpTime;
    Pi0Res[tmpID] = tmpResolution;
  }
  while( ifsCsm >> tmpID >> tmpTime >> tmpResolution ){
    CSMTime[tmpID] = tmpTime; 
    CSMRes[tmpID] = tmpResolution;
  }

  std::ofstream ofs("ResultTimeDelta.dat");
  for( int i = 0; i < 2716; i++){
    if( CSMRes[i] < 0xFFFF ){
      ofs << i                       << "\t" 
	  << Pi0Time[i] + CSMTime[i] << "\t"
	  << CSMRes[i]               << "\n";
    }else{
      ofs << i                       << "\t"
	  <<  0                      << "\t"
          << 0xFFFF                  << "\n";
    }
  }
  ofs.close();
}



  
