#ifndef __E14ReadSumFile__h__
#include "E14ReadSumFile.h"
#endif

E14ReadSumFile::E14ReadSumFile(){
  ch = new TChain("T");
  this->SetBranchAddress();
}

E14ReadSumFile::~E14ReadSumFile(){
  delete ch;
}

bool E14ReadSumFile::SetBranchAddress(){
  
  ch->SetBranchAddress("CsiNumber",&CsiNumber);
  ch->SetBranchAddress("CsiModID",CsiModID);
  ch->SetBranchAddress("CsiEne",CsiEne);
  ch->SetBranchAddress("CsiTime",CsiTime);

  ch->SetBranchAddress("CC03Number",&CC03Number);
  ch->SetBranchAddress("CC03ModID",CC03ModID);
  ch->SetBranchAddress("CC03Ene",CC03Ene);
  ch->SetBranchAddress("CC03Time",CC03Time);

  ch->SetBranchAddress("CVNumber",&CVNumber);
  ch->SetBranchAddress("CVModID",CVModID);
  ch->SetBranchAddress("CVEne",CVEne);
  ch->SetBranchAddress("CVTime",CVTime);

  ch->SetBranchAddress("LaserNumber",&LaserNumber);
  ch->SetBranchAddress("LaserModID",LaserModID);
  ch->SetBranchAddress("LaserEne",LaserEne);
  ch->SetBranchAddress("LaserTime",LaserTime);

  ch->SetBranchAddress("CosmicNumber",&CosmicNumber);
  ch->SetBranchAddress("CosmicModID",CosmicModID);
  ch->SetBranchAddress("CosmicEne",CosmicEne);
  ch->SetBranchAddress("CosmicTime",CosmicTime);

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


