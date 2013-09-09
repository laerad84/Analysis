#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <string>
#include <list>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
//#include "cluster/ClusterFinder.h"
#include "ClusterFinder_EDIT.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "User_Function.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "IDHandler.h"

//#include "E14WavReader.h"
#include "E14WavReader_V1.h"
//#include "E14WaveReader_V2.h"
#include "EnergyConverter.h"

#include "L1TrigCounter.h"

int
main( int argc ,char ** argv ){
  
  int RunNumber = atoi( argv[1]);
  int TypeIndex    = atoi( argv[2]);
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string HOME         = std::getenv("HOME");
  std::string iFileForm="%s/run_wav_%d.root";
  std::string oFileForm="%s/run_wav_%d_%s.root";
  std::string TCalFile = Form("%s/Data/TimeOffset/TimeOffset_with_cosmic.dat",ANALYSISLIB.c_str());  
  std::string ECalFile = Form("%s/local/Analysis/K3pi0Producer/Data/CalibrationFactorADV_15.dat",HOME.c_str());
  std::string TempCalibrationFilename = Form("%s/Data/Temperature_Factor/TemperatureCorrectionFactor.root",ANALYSISLIB.c_str());  

  const int nFileType = 6; 
  char* FileTypes[nFileType]={"3pi0_OldComp"  ,
			      "3pi0_LaserComp",
			      "3pi0_3pi0Comp",
			      "3pi0_noCal",//Missed noComp//
			      "3pi0_OldComp_wopi0",
			      "3pi0_NoCompNoCal"};
  int EnergyConvInt;
  switch( TypeIndex ){
  case 0:
    EnergyConvInt = 0; 
    break;
  case 1:
    EnergyConvInt = 1;
    break;
  case 2:
    EnergyConvInt = 3;
    break;
  case 3:
    EnergyConvInt = 2;
    break;
  case 4:
    EnergyConvInt = 0;
    break;
  case 5:
    EnergyConvInt = 2;
    break;
  default:
    return -1;
    break;
  }

  TFile* tfTempCorr  =new TFile(TempCalibrationFilename.c_str());
  TTree* trTempCorr  =(TTree*)tfTempCorr->Get("TemperatureCorrectionCsI");
  Double_t TempCorFactor=0;
  trTempCorr->SetBranchAddress("CorrectionFactor",&TempCorFactor);
  trTempCorr->GetEntry(RunNumber);
  if( TempCorFactor == 0){
    TempCorFactor = 1;
  }
  std::cout<< TempCorFactor << std::endl;
  EnergyConverter* Converter = new EnergyConverter(EnergyConvInt);// 1=new version // 2=non Height convertion // 3 = withCalibrawtion Result
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",
					  ANALYSISLIB.c_str()));

  TChain* trin = new TChain("Tree"); 
  trin->Add(Form(iFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber));
  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber,FileTypes[TypeIndex]),"recreate");
  TTree* trout = new TTree("T", "Output from Time zero" );  
  
  int EventNumber;
  int nCSIDigi = 0;
  int CSIDigiID[2716]={-1};
  double CSIDigiE[2716] = {0.};
  double CSIDigiTime[2716] ;
  double CSIDigiHHTime[2716];
  double CSIDigiSignal[2716]={0};
  double CSIL1TrigCount[20];
  int    CSIL1nTrig;

  trout->Branch("RunNumber"  ,&RunNumber   ,"RunNumber/I");
  trout->Branch("EventNumber",&EventNumber ,"EventNumber/I");
  trout->Branch("CsiNumber"  ,&nCSIDigi    ,"CsiNumber/I");
  trout->Branch("CsiModID"   ,CSIDigiID    ,"CSIModID[CsiNumber]/I");//nCSIDigi
  trout->Branch("CsiEne"     ,CSIDigiE     ,"CsiEne[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiTime"    ,CSIDigiTime  ,"CsiTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiHHTime"  ,CSIDigiHHTime,"CsiHHTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiSignal"  ,CSIDigiSignal,"CsiSignal[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiL1nTrig" ,&CSIL1nTrig  ,"CsiL1nTrig/I");
  trout->Branch("CsiL1TrigCount",CSIL1TrigCount,"CsiL1TrigCount[20]/D");
  /*
  trout->Branch("nCSIDigi",&nCSIDigi,"nCSIDigi/I");
  trout->Branch("CSIDigiID",CSIDigiID,"CSIDigiID[nCSIDigi]/I");//nCSIDigi
  trout->Branch("CSIDigiE" ,CSIDigiE,"CSIDigiE[nCSIDigi]/D");//nCSIDigi
  trout->Branch("CSIDigiTime",CSIDigiTime,"CSIDigiTime[nCSIDigi]/D");//nCSIDigi
  trout->Branch("CSIDigiHHTime",CSIDigiHHTime,"CSIDigiHHTime[nCSIDigi]/D");//nCSIDigi
  */
  double CSIL1TrigCountThreshold[20] = {1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
					1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};

  E14GNAnaDataContainer data; 

  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  data.branchOfKlong( trout );

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  double TimeDelta[2716];
  double TimeDeltaSig[2716];
  double CalibrationFactor[2716] = {1};
  double TimeDeltaLength[2716]={0};
  int tmpID;
  double tmpDelta;
  double tmpDeltaSig;
  double tmpCalFactor; 
  std::string ANAFILEDIR = std::getenv("HOME");
  //std::ifstream ifs(Form("%s/local/Analysis/K3pi0Producer/Data/Pi0Peak.dat",ANAFILEDIR.c_str()));
  std::ifstream ifsTCal(Form(TCalFile.c_str(),ANALYSISLIB.c_str()));
  if( !ifsTCal.is_open() ) { std::cerr <<"File does not exist."<< Form(TCalFile.c_str(),ANALYSISLIB.c_str())  << std::endl; return -1;}
  while( ifsTCal >> tmpID >> tmpDelta >> tmpDeltaSig ){
    TimeDelta[ tmpID ]    = tmpDelta;
    TimeDeltaSig[ tmpID ] = tmpDeltaSig; 
  }
  std::ifstream ifsECal(Form(ECalFile.c_str(),ANAFILEDIR.c_str()));
  if( !ifsECal.is_open() ){ std::cerr << "File does not exist." << Form(ECalFile.c_str(),ANAFILEDIR.c_str()) << std::endl; return -1; }
  while( ifsECal >> tmpID >> tmpCalFactor ){
    CalibrationFactor[ tmpID ] = tmpCalFactor;
  }
  
  IDHandler* handler = new IDHandler();
  double x,y;
  double AlzPosition = 2624;//mm
  double sol = 299.792458;//mm/ns
  for( int i = 0 ; i< 2716; i++){
    handler->GetMetricPosition( i, x, y );
    double length = TMath::Sqrt( AlzPosition*AlzPosition + x*x + y*y );
    TimeDeltaLength[i] = length/sol;
  }
  Double_t Pi0PeakCorFactor = 0.9937;
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////


  E14WavReader_V1* reader = new E14WavReader_V1(trin);
  GammaFinder   gFinder;
  ClusterFinder_EDIT clusterFinder;
  L1TrigCounter* l1 = new L1TrigCounter();
  l1->ReadMapFile();
  l1->SetThreshold(1800);
  l1->Reset();


  int nCsI               = 0;
  int CsIID[2716]        = {-1};
  double CsIEnergy[2716] = {0.};
  double CsISignal[2716] = {0.};
  double CsITime[2716]   = {-1};
  double CsIHHTime[2716] = {-1};
  
  std::cout<< __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;
  Long_t entries =  reader->fChain->GetEntries();
  for( int ievent  = 0; ievent < entries ; ievent++){
    nCsI = 0;
    nCSIDigi = 0;
    CSIL1nTrig =0;
    l1->Reset();
    for( int ich = 0; ich < 2716; ich++){
      CsIID[ich]     = -1; 
      CsIEnergy[ich] = 0.;
      CsITime[ich]   = -1;
      CsIHHTime[ich] = -1;
      CsISignal[ich] = 0.;

      CSIDigiID[ich] = 0;
      CSIDigiE[ich] = 0;
      CSIDigiTime[ich]   = 0;
      CSIDigiHHTime[ich] = 0;
      CSIDigiSignal[ich] = 0.;
    }    
    if( (ievent%100) ==0 && ievent ){ std::cout<< ievent << std::endl;}
    reader->GetEntry( ievent  );
    data.reset();
    EventNumber = reader->EventNo;
    if(( reader->TrigFlag & 3 ) != 0 ){ continue; }
    if( reader->CsinTimeCluster == 0 ){ continue; }
    if( reader->CsinTimeCluster  > 2  ){ continue ;}    
    if( reader->CsinTimeCluster == 2 && reader->CsiTimeClusterHead[0] < 50 ){ continue; }

    //// assumption :: All REAL csi Event if TimeCluster #0 
    //    std::cout<< "Analysis" << std::endl;

    nCsI = 0;
    for( int ich  = 0; ich < reader->CsiNumber; ich++){      
      int CsiID        = reader->CsiID[ich];

      double CsiTime   = reader->CsiTime[ich];
      double CsiSignal = reader->CsiSignal[ich]; 
      //double CsiEnergy = reader->CsiEne[ich]*CalibrationFactor[ reader->CsiID[ich]]/TempCorFactor;      
      //double CsiEnergy = reader->CsiEne[ich]; oldinary
      Double_t CsiEnergy=0;
      if( TypeIndex < 3 ){
	CsiEnergy = Converter->ConvertToEnergy( CsiID , CsiSignal)*CalibrationFactor[CsiID]/TempCorFactor*Pi0PeakCorFactor;
      }else if( TypeIndex == 3 ){
	CsiEnergy = Converter->ConvertToEnergy( CsiID , CsiSignal)*CalibrationFactor[CsiID]/TempCorFactor;	
      }else if( TypeIndex == 4 ){
	CsiEnergy = Converter->ConvertToEnergy( CsiID , CsiSignal)*CalibrationFactor[CsiID]/TempCorFactor;	
      }else if( TypeIndex ==5 ){
	CsiEnergy = Converter->ConvertToEnergy( CsiID , CsiSignal)/TempCorFactor;
      }

      double CsiHHTime = reader->CsiHHTime[ich];
      int CsiTimeClusterID = reader->CsiTimeClusterID[ich];
      l1->Fill(CsiID,CsiSignal);
      if( CsiTimeClusterID == 0){
	CsIID[nCsI]     =  CsiID;
	CsISignal[nCsI] =  CsiSignal;
	CsIEnergy[nCsI] =  CsiEnergy;
	CsITime[nCsI]   =  CsiTime;
	CsIHHTime[nCsI] =  CsiHHTime;
	nCsI++;
      }
    }
    CSIL1nTrig =0;
    std::vector<double> vecCount    = l1->GetCount(); 
    if( vecCount.size() != 20 ){ std::cout << "Vector Size Error" << std::endl;}
    for( int i = 0; i< 20; i++){
      CSIL1TrigCount[i] = vecCount.at(i);
      if(CSIL1TrigCount[i] > CSIL1TrigCountThreshold[i]){
	CSIL1nTrig++;
      }
    }
    /// Adjustment ///

    nCSIDigi=0;
    for( int idigi = 0; idigi< nCsI; idigi++){
      if( CsISignal[idigi] > 5 && CsIEnergy[idigi]>3){
	CSIDigiID[ nCSIDigi ]     = CsIID[idigi];
	CSIDigiE[ nCSIDigi ]      = CsIEnergy[idigi];
	CSIDigiTime[ nCSIDigi ]   = CsITime[idigi]-TimeDelta[CsIID[idigi]];
	CSIDigiHHTime[ nCSIDigi ] = CsIHHTime[idigi];
	CSIDigiSignal[nCSIDigi]   = CsISignal[idigi];
	nCSIDigi++;
      }
    }
    
    if( nCSIDigi<5){ continue;}

    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    
    /*
    std::cout<< "Cluster" << std::endl;
    std::cout << nCSIDigi << std::endl;
    for( int i = 0; i< nCSIDigi; i++ ){

      std::cout<< i            << " / "
	       << nCSIDigi     << " : " 
	       << CSIDigiID[i] << " : " 
	       << CSIDigiE[i]  << " : "
	       << CSIDigiTime[i] << std::endl;
    }
    */
    
    clist = clusterFinder.findCluster( nCSIDigi,CSIDigiID,CSIDigiE,CSIDigiTime);

    gFinder.findGamma( clist, glist );
    if( clist.size() < 6 ){ continue; }
    if( glist.size() ==6 ){ 
      if( user_rec(glist,klVec)){
	data.setData( clist );
	data.setData( glist );
	user_cut( data, klVec );
	data.setData(klVec);    
	trout->Fill();
      }
    }
  }
  std::cout<< "End" << std::endl;
  trout->Write();
  std::cout<< "Write" << std::endl;
  tfout->Close();
  return 0;
}
