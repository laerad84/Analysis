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

#include "User_Functions.h"
Int_t main( int argc , char** argv ){
  
  if( argc != 2 ){ 
    std::cerr <<" Input File type " << std::endl;
    return -1; 
  }
  const double KLMassdata = 497.614;

  const int nFileType = 18;
  int FileType = atoi( argv[1] ); 
  char *RunName[nFileType] = {"SIM","WAV","SUM",
			      "SIMFAST","WAVNOCV","WAVNEWCOMPNONCAL",
			      "3pi0_OldComp","3pi0_LaserComp","3pi0_3pi0Comp",
			      "3pi0_noComp_wopi0","3pi0_OldComp_wopi0","3pi0_NoCompNoCal",
			      "3pi0_OldComp_NOCV","SIMFULL_NEW_7G","DATA_NONTIMECAL1",
			      "DATA_NONTIMECALNOCV","test","WAVNOCV_TIME_5nsCut"};
  std::cout<< RunName[FileType] << std::endl;

  if( FileType == 0 ){
    std::cout<< "Read Simulation File" << std::endl;
  }else if( FileType == 1 ){
    std::cout<< "Read Data( Wave Analysis ) File" << std::endl;
  }else if( FileType == 2 ){
    std::cout<< "Read Data( Sumup ) File" << std::endl; 
  }else if( FileType == 3 ){
    std::cout<< "Read Data( SIMFast ) File" << std::endl; 
  }else if( FileType == 4 ){
    std::cout<< "Read Data( Wave Analysis:W/O CV) File" << std::endl;
  }else if( FileType == 5 ){
    std::cout<< "Read Data( Wave Analysis:Old Compensated) File" << std::endl;
  }else if( FileType == 6 ){
    std::cout<< "ReadData( Wave Analysis:Old Compensation) File" << std::endl;
  }else if( FileType == 7 ){
    std::cout<< "ReadData( Wave Analysis:Laser Compensation) File" << std::endl;
  }else if( FileType == 8 ){
    std::cout<< "ReadData( Wave Analysis:3pi0 Compensation) File" << std::endl; 
  }else if( FileType == 9 ){
    std::cout<< "ReadData( Wave Analysis:Non Compensation) File" << std::endl;
  }else if( FileType == 10 ){
    std::cout<< "ReadData( Wave Analysis:Old Comp&No pi0 ) File" << std::endl;
  }else if( FileType == 11 ){
    std::cout<< "ReadData( Wave Analysis: No Comp No Cal) File" << std::endl;
  }else if( FileType == 12 ){
    std::cout<< "ReadData( Wave Analysis: W/O CV OldComp) File" << std::endl;
  }else if( FileType == 13 ){
    std::cout<< "ReadData( SIMulation FULL ) File" << std::endl;
  }else if( FileType == 14 ){
    std::cout << "ReadData(NonTime Cal(CV)) File" << std::endl;
  }else if( FileType == 15 ){
    std::cout << "ReadData(NonTime Cal(NOCV)) File" << std::endl;
  }else if( FileType == 16 ){
    std::cout << "ReadData(Time Cal(NOCV) File" << std::endl;
  }else if( FileType == 17 ){
    std::cout << "ReadData(NOCV_FULL) File" << std::endl;
  }else{
    return -1;
  }

  std::string ROOTFILE_SUMUP = std::getenv("ROOTFILE_SUMUP");
  std::string ROOTFILE_WAV   = std::getenv("ROOTFILE_WAV");

  //std::string ROOTFILE_3PI0CALIBRATIONWAV = std::getenv("ROOTFILE_3PI0CALIBRATIONWAV");
  std::string ROOTFILE_3PI0CALIBRATIONWAV = std::getenv("ROOTFILE_WAV");
  std::string ROOTFILE_3PI0CALIBRATIONSIM = "/group/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/3pi0Run/SIM3PI0";  
  TChain* ch;
  switch ( FileType ){
  case 0:
    ch = new TChain("T");
    break;
  default:
    ch = new TChain("T");
    break;
  }

  std::cout<< "Read File" << std::endl;
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
  }else if(FileType == 4){
    std::string HOMEDIR = std::getenv("HOME");
    std::ifstream ifsRunNumber(Form("%s/local/Analysis/RunList/3pi0WOCVRunList.txt",HOMEDIR.c_str()));
    if( !ifsRunNumber.is_open() ){ 
      std::cerr << "File dosen't exist" << std::endl;
      return -1; 
    }
    int tmpRunNumber; 
    while( ifsRunNumber >> tmpRunNumber ){
      //if( tmpRunNumber < 4249 ){ continue; }
      //if( tmpRunNumber > 4624 ){ continue; }
      //ch->Add(Form("%s/CalibrationADV_%d_15.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
      ch->Add(Form("%s/run_wav_%d_Cal_FNL_COS_newTimeOffset.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }    
  }else if(FileType == 5){
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
      ch->Add(Form("%s/run_wav_%d_Cal_FNL_COS.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  }else if(FileType == 6){
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
      ch->Add(Form("%s/run_wav_%d_3pi0_OldComp.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  }else if(FileType == 7){
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
      ch->Add(Form("%s/run_wav_%d_3pi0_LaserComp.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  }else if(FileType == 8){
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
      ch->Add(Form("%s/run_wav_%d_3pi0_3pi0Comp.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  }else if(FileType == 9){
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
      ch->Add(Form("%s/run_wav_%d_3pi0_noCal.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  }else if(FileType == 10){
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
      ch->Add(Form("%s/run_wav_%d_3pi0_OldComp_wopi0.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  } else if(FileType == 11){
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
      ch->Add(Form("%s/run_wav_%d_3pi0_NoCompNoCal.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  }else if(FileType == 12){
    std::string HOMEDIR = std::getenv("HOME");
    std::ifstream ifsRunNumber(Form("%s/local/Analysis/RunList/RunList_3pi0_wo_CV.csv",HOMEDIR.c_str()));
    if( !ifsRunNumber.is_open() ){ 
      std::cerr << "File dosen't exist" << std::endl;
      return -1; 
    }
    int tmpRunNumber; 
    while( ifsRunNumber >> tmpRunNumber ){
      ch->Add(Form("%s/run_wav_%d_Cal_FNL_COS_newTimeOffset.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
    }
  }else if( FileType == 13){
    //for( int i = 0; i< 120; i++){
    for( int i = 0; i< 8; i++){
      //ch->Add(Form("%s/Calibration_with_4e9/Calibration_%03d0_15.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
      //ch->Add(Form("%s/out_KL3pi0.mac_1000000_%d_FEB_CL_KL.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
      //ch->Add(Form("%s/Sim3pi0_wav_KL_RES_LY_%d.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
      //ch->Add(Form("%s/Sim3pi0_wav_ALCV_KL_RES_LY_pe_5E8_KL_%d.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
      ch->Add(Form("%s/Sim3pi0_wav_ALCV_KL_RES_LY_pe_5E8_KL_%d_7G.root",ROOTFILE_3PI0CALIBRATIONSIM.c_str(),i));
    }
  }else if( FileType == 14){
    std::string HOMEDIR = std::getenv("HOME");
    //std::ifstream ifsRunNumber(Form("%s/local/Analysis/RunList/KLRunList_3_W_CV.txt",HOMEDIR.c_str()));
    std::ifstream ifsRunNumber(Form("%s/local/Analysis/RunList/RunList_3pi0_5Crate150.csv",HOMEDIR.c_str()));
    int tmpRunNumber;
    if( !ifsRunNumber.is_open()){ std::cout<< "No RunList file" << std::endl;return -1;}
    while( ifsRunNumber >> tmpRunNumber ){
      ch->Add(Form("%s/run_wav_%d_2.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
    }
  }else if( FileType == 15){
    std::string HOMEDIR = std::getenv("HOME");
    std::ifstream ifsRunNumber(Form("%s/local/Analysis/RunList/RunList_3pi0_wo_CV.csv",HOMEDIR.c_str()));
    int tmpRunNumber;
    if( !ifsRunNumber.is_open()){ std::cout<< "No RunList file" << std::endl;return -1;}
    while( ifsRunNumber >> tmpRunNumber ){
      ch->Add(Form("%s/run_wav_%d_2.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
    }
  }else if( FileType == 16){
    std::string HOMEDIR = std::getenv("HOME");
    std::ifstream ifsRunNumber(Form("%s/local/Analysis/RunList/RunList_3pi0_wo_CV.csv",HOMEDIR.c_str()));
    int tmpRunNumber;
    if( !ifsRunNumber.is_open()){ std::cout<< "No RunList file" << std::endl;return -1;}
    while( ifsRunNumber >> tmpRunNumber ){
      ch->Add(Form("%s/run_wav_%d_2.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
    }
  }else if(FileType == 17){
    std::string HOMEDIR = std::getenv("HOME");
    std::ifstream ifsRunNumber(Form("%s/local/Analysis/RunList/3pi0WOCVRunList.txt",HOMEDIR.c_str()));
    if( !ifsRunNumber.is_open() ){ 
      std::cerr << "File dosen't exist" << std::endl;
      return -1; 
    }
    int tmpRunNumber; 
    while( ifsRunNumber >> tmpRunNumber ){
      //if( tmpRunNumber < 4249 ){ continue; }
      //if( tmpRunNumber > 4624 ){ continue; }
      //ch->Add(Form("%s/CalibrationADV_%d_15.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
      ch->Add(Form("%s/run_wav_%d_GammaTime_5nsCut.root",ROOTFILE_3PI0CALIBRATIONWAV.c_str(),tmpRunNumber));
      //ch->Add(Form("/media/3TB_1/DataAll/DataAll/Data/run_wav_%d_GammaTime_TCut.root",tmpRunNumber));
    }    
  }
  
  CsiCut* csiCut = new CsiCut();
  GammaCut* gammaCut = new GammaCut();

  std::cout<< "Total Event Number : " << ch->GetEntries() << std::endl;  
  int CsiNumber;
  int CsiModID[3000];
  double CsiEne[3000];
  double CsiTime[3000];
  double CsiSignal[3000];

  int cCsiNumber;
  int cCsiModID[3000];
  double cCsiEne[3000];
  double cCsiTime[3000];
  double cCsiSignal[3000];

  int s_arrSize = 120;
  Int_t    GamClusNumbers;
  Int_t    GamClusSizes[120];
  Double_t GamClusCsiSignal[120][120];
  Double_t GamClusCsiChisq[120][120];
  Int_t    GamClusCsiL1[120][120];
  Int_t    GamClusCsiCrate[120][120];
  int RunNumber;
  int EventNumber;
  int CsiL1nTrig;
  double CsiL1TrigCount[20];
  E14GNAnaDataContainer data;
  data.setBranchAddress( ch );
  csiCut->SetBranchAddress( ch );
  
  ch->SetBranchAddress("RunNumber",&RunNumber);
  ch->SetBranchAddress("EventNumber",&EventNumber);
  /*
  ch->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  ch->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  ch->SetBranchAddress("CsiNumber",&CsiNumber);
  ch->SetBranchAddress("CsiModID",CsiModID);
  ch->SetBranchAddress("CsiEne",CsiEne);
  ch->SetBranchAddress("CsiTime",CsiTime);
  ch->SetBranchAddress("CsiSignal",CsiSignal);
  */
  /*
  ch->SetBranchAddress("GamClusNumbers",&GamClusNumbers);
  ch->SetBranchAddress("GamClusSizes",GamClusSizes);//GamClusNumbers
  ch->SetBranchAddress("GamClusCsiSignal",GamClusCsiSignal);//GamClusNumbers
  ch->SetBranchAddress("GamClusCsiChisq",GamClusCsiChisq);//GamClusNumbers
  ch->SetBranchAddress("GamClusCsiL1",GamClusCsiL1);//GamClusNumbers
  ch->SetBranchAddress("GamClusCsiCrate",GamClusCsiCrate);//GamClusNumbers
  */

  //// Set Output File //// 
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

  //trKL->Branch("CsiL1nTrig",&CsiL1nTrig,"CsiL1nTrig/I");
  //trKL->Branch("CsiL1TrigCount",CsiL1TrigCount,"CsiL1TrigCount[20]/D");

  E14GNAnaDataContainer dataCopy;
  dataCopy.branchOfClusterList(trKL);
  //dataCopy.branchOfGammaList(trKL);
  //dataCopy.branchOfPi0List(trKL);
  dataCopy.branchOfKlong(trKL);
  csiCut->Branch(trKL);
  gammaCut->Branch(trKL);
  int DataCutCondition;
  trKL->Branch("DataCutCondition",&DataCutCondition,"DataCutCondition/I");
  trKL->Branch("RunNumber",&RunNumber,"RunNumber/I");
  trKL->Branch("EventNumber",&EventNumber,"EventNumber/I");

  /*
  trKL->Branch("CsiNumber",&cCsiNumber,"CsiNumber/I");
  trKL->Branch("CsiModID",cCsiModID,"CsiModID[CsiNumber]/I");  
  trKL->Branch("CsiEne",cCsiEne,"CsiEne[CsiNumber]/D");//CsiNumber
  trKL->Branch("CsiTime",cCsiTime,"CsiTime[CsiNumber]/D");//CsiNumber
  trKL->Branch("CsiSignal",cCsiSignal,"CsiSignal[CsiNumber]/D");//CsiNumber  

  trKL->Branch("CsiL1nTrig",&CsiL1nTrig,"CsiL1nTrig/I");
  trKL->Branch("CsiL1TrigCount",CsiL1TrigCount,"CsiL1TrigCount[20]/D");
  */
  
  /*
  trKL->Branch("GamClusNumbers",&GamClusNumbers,"GamClusNumbers/I");
  trKL->Branch("GamClusSizes",GamClusSizes,"GamClusSizes[GamClusNumbers]/I");//GamClusNumbers 
  trKL->Branch("GamClusCsiSignal",GamClusCsiSignal,Form("GamClusCsiSignal[GamClusNumbers][%d]/D",s_arrSize));//GamClusNumbers 
  trKL->Branch("GamClusCsiChisq",GamClusCsiChisq,Form("GamClusCsiChisq[GamClusNumbers][%d]/D",s_arrSize));//GamClusNumbers 
  trKL->Branch("GamClusCsiL1",GamClusCsiL1,Form("GamClusCsiL1[GamClusNumbers][%d]/I",s_arrSize));//GamClusNumbers 
  trKL->Branch("GamClusCsiCrate",GamClusCsiCrate,Form("GamClusCsiCrate[GamClusNumbers][%d]/I",s_arrSize));//GamClusNumbers 
  */
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
    //std::cout<< ievent << std::endl;
    //data.reset();
    dataCopy.reset();
    //if(ievent > 1000){ break; }
    /*
    cCsiNumber = 0; 
    for( int i = 0; i< 3000; i++){
      cCsiModID[i] = 0;
      cCsiEne[i]   = 0; 
      cCsiTime[i]  = 0;
      cCsiSignal[i] = 0; 
    }

    for( int i = 0; i< csiCut->CsiNumber; i++){
      cCsiModID[i] = CsiCut->CsiID[i];
      cCsiEne[i] = CsiCut->CsiEne[i];
      cCsiTime[i] = CsiCut->CsiTime[i];
      cCsiSignal[i] = CsiCut->CsiSignal[i];
      cCsiNumber++;
    }
    */
    

    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    data.getData( clist );
    data.getData( glist );
    data.getData( klVec );
   

    //dataCopy.setData( clist );
    //dataCopy.setData( glist );
    if( klVec.size() == 0 ){ continue; }
    dataCopy.setData( clist );
    dataCopy.setData( klVec );
    gammaCut->Decision( klVec );
    //std::cout<< klVec.size() << "\t" << clist.size() << "\t" << glist.size() << std::endl;
    //if( clist.size() == 0 ){ continue; }
    if( glist.size() <  6 ){ continue; }
    /*
    bool bGPosition  = false;
    bool bGEne       = false;
    bool bPi0pt      = false;

    bool bKlongMass  = false;
    bool bKlongChisq = false;
    bool bKlongChisqD= false;
    bool bKlongVtxZ  = false;
    bool bKlongE     = false;

    double pi0MaxPt = 0;

    for( int i = 0; i< klVec[0].pi0().size();i++){
      double x(0);
      double y(0);
      x = klVec[0].pi0()[i].g1().x();
      y = klVec[0].pi0()[i].g1().y();
      if( TMath::Abs(x) < 150 && TMath::Abs(y)< 150 ){
	bGPosition = true;
      }   
      if(TMath::Sqrt( x*x +y*y) > 850 ){continue; }
      if(TMath::Abs(y) > 550 ){ continue; }
      if(klVec[0].pi0()[i].g1().e() < 200 ){
	bGEne = true; 
      }
    }
    for( int i = 0; i< klVec[0].pi0().size();i++){
      double x(0);
      double y(0);
      x = klVec[0].pi0()[i].g2().x();
      y = klVec[0].pi0()[i].g2().y();
      if( TMath::Abs(x) < 150 && TMath::Abs(y)< 150 ){
	bGPosition = true;
      }   
      if(TMath::Sqrt( x*x +y*y) > 850 ){continue; }
      if(TMath::Abs(y) > 550 ){ continue; }
      if( klVec[0].pi0()[i].g2().e() < 200 ){
	bGEne = true; 
      }
    }


    for( int i = 0; i< klVec[0].pi0().size(); i++){
      double px = klVec[0].pi0()[i].p3()[0];
      double py = klVec[0].pi0()[i].p3()[1];
      double pt = TMath::Sqrt( px*px+py*py );
      if( pt > pi0MaxPt ){ pi0MaxPt = pt ;}
    }    
    if( pi0MaxPt > 150 ){ bPi0pt  = true; }


    if( TMath::Abs( klVec[0].m() - KLMassdata ) > 20 ){
      bKlongMass = true;
    }

    if( klVec[0].chisqZ() > 50 ){ bKlongChisq = true; }
    if( klVec[1].chisqZ() - klVec[0].chisqZ() > 50 ){ bKlongChisqD = true; }
    if( klVec[0].e() > 5000 ){ bKlongE = true; }
    if( klVec[0].vz() > 5000 || klVec[0].vz() < 1500 ){ bKlongVtxZ = true; }

    DataCutCondition =0;
    if( bGPosition ){
      DataCutCondition |= 1<< 1;
    }
    if( bGEne ){
      DataCutCondition |= 1<< 2;
    }
    if( bPi0pt ){
      DataCutCondition |= 1<< 3;      
    }
    if( bKlongChisq ){
      DataCutCondition |= 1<< 4;
    }
    if( bKlongChisqD ){
      DataCutCondition |= 1<< 5;
    }
    if( bKlongE ){
      DataCutCondition |= 1 << 6;
    }
    if( bKlongVtxZ ){
      DataCutCondition |= 1 << 7;
    }
    if( bKlongMass ){
      DataCutCondition |= 1 << 8;
    }
    */
    trKL->Fill();
  }
  
  trKL->Write();
  tfout->Close();
  return 0; 
}
