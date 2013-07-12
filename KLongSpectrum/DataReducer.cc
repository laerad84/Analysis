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
  }else if( FileType == 3 ){
    std::cout<< "Read Data( SIMFast ) File" << std::endl; 
  }

  std::string ROOTFILE_SUMUP = std::getenv("ROOTFILE_SUMUP");
  std::string ROOTFILE_WAV   = std::getenv("ROOTFILE_WAV");

  //std::string ROOTFILE_3PI0CALIBRATIONWAV = std::getenv("ROOTFILE_3PI0CALIBRATIONWAV");
  std::string ROOTFILE_3PI0CALIBRATIONWAV = std::getenv("ROOTFILE_WAV");
  std::string ROOTFILE_3PI0CALIBRATIONSIM = "/group/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/3pi0Run/SIM3PI0";  
  TChain* ch;
  if( FileType == 0 ){ 
    ch = new TChain("T");
  }else{
    ch = new TChain("T");
  }

  if( FileType == 0){
    //for( int i = 0; i< 120; i++){
    for( int i = 0; i< 60; i++){
      //ch->Add(Form("%s/Calibration_with_4e9/Calibration_%03d0_15.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
      //ch->Add(Form("%s/out_KL3pi0.mac_1000000_%d_FEB_CL_KL.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
      //ch->Add(Form("%s/Sim3pi0_wav_KL_RES_LY_%d.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
      ch->Add(Form("%s/Sim3pi0_wav_KL_RES_LY_pe_%d.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
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
      if( tmpRunNumber < 4249 ){ continue; }
      if( tmpRunNumber > 4624 ){ continue; }
      //ch->Add(Form("%s/CalibrationADV_%d_15.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
      ch->Add(Form("%s/run_wav_%d_Cal_FNL_COS_newTimeOffset.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  }else if(FileType == 3){ //Fsim event
    for( int i = 0; i< 20; i++){
      ch->Add(Form("%s/Sim3pi0_wav_fast_KL_RES_LY_pe_%d.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
    }
  }
  std::cout<< "Total Event Number : " << ch->GetEntries() << std::endl;  

  int CsiL1nTrig;
  double CsiL1TrigCount[20];
  ch->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  ch->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);

  E14GNAnaDataContainer data;
  data.setBranchAddress( ch );
  
  //// Set Output File //// 
  char *RunName[4] = {"SIM","WAV","SUM","SIMFAST"};
  TFile* tfout = new TFile(Form("Kl_Total_%s.root",RunName[FileType]),"RECREATE");
  TTree* trKL = new TTree("trKL","Klong Tree");
  Double_t KLMass;
  Double_t KLChisq;
  Double_t KLSecChisq;
  Double_t KLE;
  Double_t KLPos[3];
  Double_t KLMom[3];
  Double_t GammaE[6];
  Double_t GammaPos[6][3];
  Double_t GammaTime[6];
  trKL->Branch("CsiL1nTrig",&CsiL1nTrig,"CsiL1nTrig/I");
  trKL->Branch("CsiL1TrigCount",CsiL1TrigCount,"CsiL1TrigCount[20]/D");
  E14GNAnaDataContainer dataCopy;
  dataCopy.branchOfKlong(trKL);
  /*
  trKL->Branch("KLMass"  ,&KLMass  ,"KLMass/D");
  trKL->Branch("KLChisq" ,&KLChisq ,"KLChisq/D");
  trKL->Branch("KLSecChisq",&KLSecChisq,"KLSecChisq/D");
  trKL->Branch("KLE"     ,&KLE     ,"KLE/D");
  trKL->Branch("KLPos"   ,KLPos    ,"KLPos[3]/D");
  trKL->Branch("GammaE"  ,GammaE   ,"GammaE[6]/D");
  trKL->Branch("GammaPos",GammaPos ,"GammaPos[6][3]/D");
  trKL->Branch("GammaTime",GammaTime,"GammaTime[6]/D");
  trKL->Branch("KLMom"   ,KLMom    ,"KLMom[3]/D");
  */

  std::cout<< ch->GetEntries() << std::endl;

  for( int ievent = 0; ievent < ch->GetEntries() ; ievent++){
    ch->GetEntry( ievent );
    data.reset();
    //if(ievent > 1000){ break; }
    
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    data.getData( clist );
    data.getData( glist );
    data.getData( klVec );
    dataCopy.setData( clist );
    dataCopy.setData( glist );
    dataCopy.setData( klVec );
    
    //std::cout<< klVec.size() << "\t" << clist.size() << "\t" << glist.size() << std::endl;
    if( klVec.size() == 0 ){ continue; }
    //if( clist.size() == 0 ){ continue; }
    if( glist.size() <  6 ){ continue; }
    //if( glist.size() != 6 ){ continue; }
    //std::cout<< clist.size() << std::endl;
    /*
    int GammaID = 0;

    for( std::list<Gamma>::iterator itGamma = glist.begin();itGamma != glist.end();itGamma++,GammaID++){
      //std::cout<<"G:"<< GammaID << std::endl;
      if( GammaID >=  6 ){ break; }
      GammaE[GammaID]      = (*itGamma).e();
      GammaPos[GammaID][0] = (*itGamma).x();
      GammaPos[GammaID][1] = (*itGamma).y();
      GammaPos[GammaID][2] = (*itGamma).z();
      GammaTime[GammaID]   = (*itGamma).t();
    }
    KLMass   = klVec[0].m();
    KLE      = klVec[0].e();
    KLPos[0] = klVec[0].vx();
    KLPos[1] = klVec[0].vy();
    KLPos[2] = klVec[0].vz();
    KLMom[0] = (klVec[0].p3()).x();
    KLMom[1] = (klVec[0].p3()).y();
    KLMom[2] = (klVec[0].p3()).z();
    KLChisq  = (klVec[0]).chisqZ();
    if( klVec.size()==2){
      KLSecChisq = (klVec[1]).chisqZ();
    }else{
      KLSecChisq = 0xFFFF;
    }
    */
    trKL->Fill();
  }
  
  trKL->Write();
  tfout->Close();
  return 0; 
}

