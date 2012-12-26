#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <fstream>
int main(){
  std::string TimeOffsetFile = "/home/jwlee/local/Analysis/CosmicAnalysisPrograms/CosmicRayAnaSignal/testNewWORKCompileOffset.txt";

  std::ifstream* ifs(TimeOffsetFile.c_str());
  
  Double_t TimeOffset[2716];
  
  Int_t ID;
  Double_t Offset;
  Double_t RMS;

  if( !ifs.is_open() ){
    std::cerr << "File is not exist." << std::endl; 
    return ; 
  }
  while ( ifs >>  ID >> Offset >> RMS ){
    TimeOffset[ID] = Offset;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  std::string wavDir = std::getenv("ROOTFILE_WAV");
  std::cout<< wavDir << std::endl;
}

  
