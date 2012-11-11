#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <list>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "TFile.h"
#include "TTree.h"

Int_t main( int argc , char** argv ){
  
  if( argc != 2 ){ 
    std::cerr <<" Input File type " << std::endl;
    return -1; 
  }

  int FileType = atoi( argv[1] ); 

  if( FileType == 0 ){
    std::cout<< "Read Simulation File" << std::endl;
  }else if( FileType == 1 ){
    std::cout<< "Read Data( Wave Analysis ) File" << std::endl;
  }else if( FileType == 2 ){
    std::cout<< "Read Data( Sumup ) File" << std::endl; 
  }

  std::string ROOTFILE_SUMUP = std::getenv("ROOTFILE_SUMUP");
  std::string ROOTFILE_WAV   = std::getenv("ROOTFILE_WAV");

  std::string ROOTFILE_3PI0CALIBRATIONWAV = std::getenv("ROOTFILE_3PI0CALIBRATIONWAV");
  std::string ROOTFILE_3PI0CALIBRATIONSUM = std::getenv("ROOTFILE_3PI0CALIBRATIONSUM");
  std::string ROOTFILE_3PI0CALIBRATIONSIM = std::getenv("ROOTFILE_3PI0CALIBRATIONSIM");
      
  TChain* ch = new TChain("trCalibration");
  if( FileType == 0){
    for( int i = 0; i< 400; i++){
      ch->Add(Form("%s/Calibration_with_4e9/Calibration_%03d0_15.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
    }
  }else if( FileType == 1 ){
    std::string HOMEDIR = std::getenv("HOME");
    std::ifstream ifsRunNumber(Form("%s/local/Analysis/RunList/KLRunList_2.txt",HOMEDIR.c_str()));
    if( !ifsRunNumber.is_open() ){ 
      std::cerr << "File dosen't exist" << std::endl;
      return -1; 
    }
    int tmpRunNumber; 
    while( ifsRunNumber >> tmpRunNumber ){
      ch->Add(Form("%s/CalibrationADV_%d_15.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  }else if( FileType ==2 ){
    std::string HOMEDIR  =std::getenv("HOME");
    std::ifstream ifsRunNumber(Form("%s/local/Analysis/RunList/KLRunList_2.txt",HOMEDIR.c_str()));
    if( !ifsRunNumber.is_open() ){
      std::cerr << "File dosen't exist" << std::endl;
      return -1;
    }
    int tmpRunNumber;
    while( ifsRunNumber >> tmpRunNumber ){
      ch->Add(Form("%s/CalibrationADV_%d_15.root",ROOTFILE_3PI0CALIBRATIONSUM.c_str(),tmpRunNumber));
    }
  }
  std::cout<< "Total Event Number : " << ch->GetEntries() << std::endl;  
  E14GNAnaDataContainer data;
  data.setBranchAddress( ch );

  //// Set Output File //// 
  char *RunName[3] = {"SIM","WAV","SUM"};
  TFile* tfout = new TFile(Form("Kl_Total_%s.root",RunName[FileType]),"RECREATE");
  TTree* trKL = new TTree("trKL","Klong Tree");
  Double_t KLMass;
  Double_t KLChisq;
  Double_t KLE;
  Double_t KLPos[3];
  
  Double_t GammaE[6];
  Double_t GammaPos[6][3];
  trKL->Branch("KLMass",&KLMass,"KLMass/D");
  trKL->Branch("KLChisq",&KLChisq,"KLChisq/D");
  trKL->Branch("KLE",&KLE,"KLE/D");
  trKL->Branch("KLPos",KLPos,"KLPos[3]/D");
  trKL->Branch("GammaE",GammaE,"GammaE[6]/D");
  trKL->Branch("GammaPos",GammaPos,"GammaPos[6][3]/D");

  for( int ievent = 0; ievent < ch->GetEntries() ; ievent++){
    ch->GetEntry( ievent );
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    data.getData( clist );
    data.getData( glist );
    data.getData( klVec );
    
    if( klVec.size() == 0 ){ continue; }
    if( clist.size() != 0 ){ continue; }
    int GammaID = 0;
    for( std::list<Gamma>::iterator itGamma = glist.begin(); 
	 itGamma != glist.end();
	 itGamma++,GammaID++){
      GammaE[GammaID] = (*itGamma).e();
      GammaPos[GammaID][0] = (*itGamma).x();
      GammaPos[GammaID][1] = (*itGamma).y();
      GammaPos[GammaID][2] = (*itGamma).z();
      
      if( GammaID >=  6 ){ break; }
    }
    
    KLMass = klVec[0].m();
    KLE    = klVec[0].e();
    KLPos[0] = klVec[0].vx();
    KLPos[1] = klVec[0].vy();
    KLPos[2] = klVec[0].vz();

    trKL->Fill();
  }
  trKL->Write();
  tfout->Close();
  return 0; 
}

