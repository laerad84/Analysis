#ifndef __E14ReadSumFile__h__
#include "E14_CLUSTER_BUILDER/E14ReadSumFile.h"

E14ReadSumFile::E14ReadSumFile(){
  ch = new TChain("T");
  this->SetBranchAddress();
  for(int  iLaser  = 0;  iLaser < N_TOTAL_Laser; iLaser++){
    LaserThreshold[iLaser] = 500;
  }
  for( int iCosmic = 0; iCosmic < N_TOTAL_Cosmic; iCosmic++){
    CosmicThreshold[iCosmic]  = 50000;
  }
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

bool E14ReadSumFile::SetOutputFile(const char* filename){
  tfOut = new TFile(filename,"recreate");
  trOut = new TTree("Tree","Output Of E14ReadSumFile");
  if( tfOut != NULL && trOut != NULL ){ return true;}
  else{ return false;}
}


bool E14ReadSumFile::BranchtoTree(TTree* tr){
  if( tr == NULL ){ return false;}
  trOut = tr;

  trOut->Branch("EventNumber",&EventNumber,"EventNumber/I");

  trOut->Branch("CsiTotalE",&CsiTotalE,"CsiTotalE/D");
  trOut->Branch("CsiNumber",&CsiNumber,"CsiNumber/I");
  trOut->Branch("CsiModID" ,CsiModID  ,"CsiModID[CsiNumber]/I");//CsiNumber
  trOut->Branch("CsiEne"   ,CsiEne    ,"CsiEne[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiTime"  ,CsiTime   ,"CsiTime[CsiNumber]/D");//CsiNumber

  trOut->Branch("CC03TotalE",&CC03TotalE,"CC03TotalE/D");
  trOut->Branch("CC03HitNumber",&CC03HitNumber,"CC03HitNumber/I");
  trOut->Branch("CC03Number",&CC03Number,"CC03Number/I");
  trOut->Branch("CC03ModID" ,CC03ModID  ,"CC03ModID[CC03Number]/I");//CC03Number
  trOut->Branch("CC03Ene"   ,CC03Ene    ,"CC03Ene[CC03Number]/D");//CC03Number
  trOut->Branch("CC03Time"  ,CC03Time   ,"CC03Time[CC03Number]/D");//CC03Number

  trOut->Branch("CVTotalE",&CVTotalE,"CVTotalE/D");
  trOut->Branch("CVHitNumber",&CVHitNumber,"CVHitNumber/I");
  trOut->Branch("CVNumber",&CVNumber,"CVNumber/I");
  trOut->Branch("CVModID" ,CVModID  ,"CVModID[CVNumber]/I");//CVNumber
  trOut->Branch("CVEne"   ,CVEne    ,"CVEne[CVNumber]/D");//CVNumber
  trOut->Branch("CVTime"  ,CVTime   ,"CVTime[CVNumber]/D");//CVNumber

  trOut->Branch("OEVTotalE",&OEVTotalE,"OEVTotalE/D");
  trOut->Branch("OEVHitNumber",&OEVHitNumber,"OEVHitNumber/I");
  trOut->Branch("OEVNumber",&OEVNumber,"OEVNumber/I");
  trOut->Branch("OEVModID" ,OEVModID  ,"OEVModID[OEVNumber]/I");//OEVNumber
  trOut->Branch("OEVEne"   ,OEVEne    ,"OEVEne[OEVNumber]/D");//OEVNumber
  trOut->Branch("OEVTime"  ,OEVTime   ,"OEVTime[OEVNumber]/D");//OEVNumber

  trOut->Branch("CrateTotalE",&CrateTotalE,"CrateTotalE/D");
  trOut->Branch("CrateHitNumber",&CrateHitNumber,"CrateHitNumber/I");
  trOut->Branch("CrateNumber",&CrateNumber,"CrateNumber/I");
  trOut->Branch("CrateModID" ,CrateModID  ,"CrateModID[CrateNumber]/I");//CrateNumber
  trOut->Branch("CrateEne"   ,CrateEne    ,"CrateEne[CrateNumber]/D");//CrateNumber
  trOut->Branch("CrateTime"  ,CrateTime   ,"CrateTime[CrateNumber]/D");//CrateNumber

  trOut->Branch("CosmicTotalE",&CosmicTotalE,"CosmicTotalE/D");
  trOut->Branch("CosmicNumber",&CosmicNumber,"CosmicNumber/I");
  trOut->Branch("CosmicModID" ,CosmicModID  ,"CosmicModID[CosmicNumber]/I");//CosmicNumber
  trOut->Branch("CosmicEne"   ,CosmicEne    ,"CosmicEne[CosmicNumber]/D");//CosmicNumber
  trOut->Branch("CosmicTime"  ,CosmicTime   ,"CosmicTime[CosmicNumber]/D");//CosmicNumber

  trOut->Branch("LaserTotalE",&LaserTotalE,"LaserTotalE/D");
  trOut->Branch("LaserNumber",&LaserNumber,"LaserNumber/I");
  trOut->Branch("LaserModID" ,LaserModID  ,"LaserModID[LaserNumber]/I");//LaserNumber
  trOut->Branch("LaserEne"   ,LaserEne    ,"LaserEne[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserTime"  ,LaserTime   ,"LaserTime[LaserNumber]/D");//LaserNumber

  trOut->Branch("LaserBit",&LaserBit,"LaserBit/I");
  trOut->Branch("CosmicBitUp",&CosmicBitUp,"CosmicBitUp/I");
  trOut->Branch("CosmicBitDn",&CosmicBitDn,"CosmicBitDn/I");

  return true;
}
 
bool E14ReadSumFile::SumData(){
  if( trOut ==NULL ){ return false;}
  CsiTotalE=0.;
  for( int iCsi = 0; iCsi < CsiNumber; iCsi++){
    CsiTotalE += CsiEne[iCsi];
  }  
  CC03TotalE=0.;
  for( int iCC03 = 0; iCC03 < CC03Number; iCC03++){
    CC03TotalE += CC03Ene[iCC03];
  }  
  CVTotalE=0.;
  for( int iCV = 0; iCV < CVNumber; iCV++){
    CVTotalE += CVEne[iCV];
  }
  OEVTotalE=0.;
  for( int iOEV = 0; iOEV < OEVNumber; iOEV++){
    OEVTotalE += OEVEne[iOEV];
  }
  
  CrateTotalE=0.;
  for( int iCrate = 0; iCrate < CrateNumber; iCrate++){
    CrateTotalE += CrateEne[iCrate];
  }
  
  CosmicTotalE=0.;
  for( int iCosmic = 0; iCosmic < CosmicNumber; iCosmic++){
    CosmicTotalE += CosmicEne[iCosmic];
  }
  
  LaserTotalE=0.;
  for( int iLaser = 0; iLaser < LaserNumber; iLaser++){
    LaserTotalE += LaserEne[iLaser];
  }
  return true;
}


bool E14ReadSumFile::SumCrate(){
  return true;
}

bool E14ReadSumFile::Write(){

  trOut->Write();
  tfOut->Close();
  return true;
}
bool E14ReadSumFile::Fill(){
  trOut->Fill();
  return true;
}

bool E14ReadSumFile::SetTriggerBit(){
  LaserBit =0;
  CosmicBitUp = 0; 
  CosmicBitDn = 0; 
  if(LaserEne[0] > LaserThreshold[0]){
    LaserBit = 1;
  }

  for( int iCosmic = 0; iCosmic < 5 ; iCosmic++){
    if( CosmicEne[iCosmic]    > CosmicThreshold[iCosmic]   ||
	CosmicEne[iCosmic+10] > CosmicThreshold[iCosmic+10]){
      CosmicBitUp |= 1 << iCosmic;
    }
    if( CosmicEne[iCosmic+5]  > CosmicThreshold[iCosmic+5]  ||
	CosmicEne[iCosmic+15] > CosmicThreshold[iCosmic+15] ){
      CosmicBitDn |= 1 << iCosmic;
     }
  }
  if( LaserBit == 0 && (CosmicBitUp & CosmicBitDn) ==0 ){
    return true;
  }else{
    return false;
  }
}
  
#endif
