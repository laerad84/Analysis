#ifndef __E14ReadSumFile__h__
#include "E14ReadSumFile.h"

E14ReadSumFile::E14ReadSumFile(Int_t Flag = 0){
  ch = new TChain("T");
  this->SetBranchAddress();
  if( Flag == 1 ){
    this->ReadTrack();
  }
}

E14ReadSumFile::~E14ReadSumFile(){
  delete ch;
}

bool E14ReadSumFile::SetBranchAddress(){
  
  ch->SetBranchAddress("CsiNumber",&CsiNumber);
  ch->SetBranchAddress("CsiModID",CsiModID);//CsiNumber
  ch->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
  ch->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  ch->SetBranchAddress("CsiIntegratedADC",CsiIntegratedADC);//CsiNumber
  
  ch->SetBranchAddress("CC03Number",&CC03Number);
  ch->SetBranchAddress("CC03ModID",CC03ModID);//CC03Number
  ch->SetBranchAddress("CC03Ene",CC03Ene);//CC03Number
  ch->SetBranchAddress("CC03Time",CC03Time);//CC03Number
  
  ch->SetBranchAddress("CVNumber",&CVNumber);
  ch->SetBranchAddress("CVModID",CVModID);//CVNumber
  ch->SetBranchAddress("CVEne",CVEne);//CVNumber
  ch->SetBranchAddress("CVTime",CVTime);//CVNumber
  
  ch->SetBranchAddress("LaserNumber",&LaserNumber);
  ch->SetBranchAddress("LaserModID",LaserModID);//LaserNumber
  ch->SetBranchAddress("LaserEne",LaserEne);//LaserNumber
  ch->SetBranchAddress("LaserTime",LaserTime);//LaserNumber
  
  ch->SetBranchAddress("CosmicNumber",&CosmicNumber);
  ch->SetBranchAddress("CosmicModID",CosmicModID);//CosmicNumber
  ch->SetBranchAddress("CosmicEne",CosmicEne);//CosmicNumber
  ch->SetBranchAddress("CosmicTime",CosmicTime);//CosmicNumber
  
  ch->SetBranchAddress("EtcNumber",&EtcNumber);
  ch->SetBranchAddress("EtcModID",EtcModID);//EtcNumber
  ch->SetBranchAddress("EtcEne",EtcEne);//EtcNumber
  ch->SetBranchAddress("EtcTime",EtcTime);//EtcNumber
  
  ch->SetBranchAddress("SciNumber",&SciNumber);
  ch->SetBranchAddress("SciModID",SciModID);//SciNumber
  ch->SetBranchAddress("SciEne",SciEne);//SciNumber
  ch->SetBranchAddress("SciTime",SciTime);//SciNumber

  return true;
}

bool E14ReadSumFile::ReadTrack(){
  
  ch->SetBranchAddress("nTrack",&nTrack);
  ch->SetBranchAddress("track",track);//nTrack
  ch->SetBranchAddress("mother",mother);//nTrack
  ch->SetBranchAddress("pid",pid);//nTrack
  ch->SetBranchAddress("mass",mass);//nTrack
  ch->SetBranchAddress("ek",ek);//nTrack
  ch->SetBranchAddress("end_ek",end_ek);//nTrack
  ch->SetBranchAddress("p",p);//nTrack
  ch->SetBranchAddress("v",v);//nTrack
  ch->SetBranchAddress("end_p",end_p);//nTrack
  ch->SetBranchAddress("end_v",end_v);//nTrack
  
  return true;
}

int E14ReadSumFile::GetEntry(long eventNumber){
  return ch->GetEntry(eventNumber);
}

long E14ReadSumFile::GetEntries(){
  return ch->GetEntries();
}

int E14ReadSumFile::Add(const char* filename){
  return ch->Add(filename);
}

#endif
