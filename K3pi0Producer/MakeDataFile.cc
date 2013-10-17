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
#include "cluster/ClusterFinder.h"
//#include "ClusterFinder_EDIT.h"
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
#include "L1TrigCounter.h"


int
main( int argc ,char ** argv ){
  int FileLevel=-1;
  if( argc < 2 ){ return -1; }
  FileLevel = std::atoi(argv[2]); 
  
  switch( FileLevel){
  case 0 :
    std::cout << "RawData Conversion, NO Calibration" << std::endl;
    break;
  case 1 :
    std::cout << "3pi0 Calibration Applied" << std::endl;
    break;
  case 2 :
    std::cout << "All calibration Applied" << std::endl;
    break;
  default :
    std::cout << "Invaelid Number" << std::endl;
    return -1;
  }

  int RunNumber = atoi( argv[1]);
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string HOME         = std::getenv("HOME");

  std::string iFileForm="%s/run_wav_%d.root";//Run ID
  std::string oFileForm="%s/run_wav_%d_%d.root";// Calibration Type , Run ID

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// READ CALIBRATION CONSTANTS
  Double_t Pi0PeakCorFactor = 0.9937;
  std::cout<< Pi0PeakCorFactor << std::endl;
  std::cout<< "READ CALIBRATION FILE" << std::endl;
  //Old Calibration File.
  //std::string TCalFile = Form("%s/Data/TimeOffset/testNewWORKCompileOffset.txt",ANALYSISLIB.c_str());  
  //std::string ECalFile = Form("%s/local/Analysis/K3pi0Producer/Data/CalibrationFactorADV_15.dat",HOME.c_str());
  //std::string TCalFile = Form("%s/Data/TimeOffset/TimeOffset_with_cosmic.dat",ANALYSISLIB.c_str());  
  std::string TCalFile = Form("%s/Data/CalibrationFile/TimeOffset_Shower_10.dat",ANALYSISLIB.c_str());
  std::string ECalFile = Form("%s/Data/CalibrationFile/CalibrationFactorWithoutNOCV.dat",ANALYSISLIB.c_str());

  std::cout<< "Temperature Calibration" << std::endl;
  std::string TempCalibrationFilename = Form("%s/Data/Temperature_Factor/TemperatureCorrectionFactor.root",ANALYSISLIB.c_str());  
  TFile* tfTempCorr  =new TFile(TempCalibrationFilename.c_str());
  TTree* trTempCorr  =(TTree*)tfTempCorr->Get("TemperatureCorrectionCsI");
  Double_t TempCorFactor=0;
  trTempCorr->SetBranchAddress("CorrectionFactor",&TempCorFactor);
  trTempCorr->GetEntry(RunNumber);
  if( TempCorFactor == 0){
    TempCorFactor = 1;
  }
  std::cout<< TempCorFactor << std::endl;
  if( tfTempCorr == NULL ){
    std::cout<< tfTempCorr->GetName() << "is not opened" << std::endl;
    return -1;
  }

  //Set Calibration constant 
  double TimeDelta[2716]={0};
  double CalibrationFactor[2716] = {1};
  for( int i = 0; i< 2716; i++){
    CalibrationFactor[i] = 1;
    TimeDelta[i]         = 0;
  }

  double TimeDeltaLength[2716]={0};
  int    tmpID;
  double tmpDelta;
  double tmpSig;
  double tmpCalFactor; 
  std::string ANAFILEDIR = std::getenv("HOME");
  //std::ifstream ifs(Form("%s/local/Analysis/K3pi0Producer/Data/Pi0Peak.dat",ANAFILEDIR.c_str()));
  std::cout<< "Read TimeOffset" << std::endl;
  std::ifstream ifsTCal(Form(TCalFile.c_str(),ANALYSISLIB.c_str()));
  if( !(ifsTCal.is_open())){ 
    std::cout << TCalFile << "is not opened" << std::endl;
    return -1;
  }
  while( ifsTCal >> tmpID >> tmpDelta >> tmpSig ){
    TimeDelta[ tmpID ]    = tmpDelta;
  }
  std::cout<< "Read Energy Calibration File" << std::endl;
  std::ifstream ifsECal(Form(ECalFile.c_str(),ANAFILEDIR.c_str()));
  while( ifsECal >> tmpID >> tmpCalFactor ){
    CalibrationFactor[ tmpID ] = tmpCalFactor;
  }
  if( !(ifsECal.is_open())){
    std::cout<< ECalFile << "is not opened" << std::endl;
    return -1; 
  }
  
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Set L1 Trigger Scheme
  std::cout<< "Set Trigger Scheme" << std::endl;
  double CSIL1TrigCountThreshold[20] = {1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
					1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};

  L1TrigCounter* l1 = new L1TrigCounter();
  l1->ReadMapFile();
  l1->SetThreshold(1800);
  l1->Reset();

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  TChain* trin = new TChain("Tree"); 
  trin->Add(Form(iFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber));
  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber,FileLevel),"recreate");
  TTree* trout = new TTree("T",Form("Output from Time zero;%s;%s",ECalFile.c_str(),TCalFile.c_str()) );  
  
  int EventNumber;
  int nCSIDigi = 0;
  int CSIDigiID[2716]={-1};
  double CSIDigiE[2716] = {0.};
  double CSIDigiTime[2716] ;
  double CSIDigiHHTime[2716];
  double CSIDigiSignal[2716]={0};
  double CSIL1TrigCount[20];
  int    CSIL1nTrig;

  int nCC03Digi = 0;
  int CC03DigiID[32]={-1};
  double CC03DigiE[32] = {0.};
  double CC03DigiTime[32] ;
  double CC03DigiHHTime[32];
  double CC03IntegratedADC[32];
  double CC03Signal[32];

  int nCVDigi = 0;
  int CVDigiID[10]={-1};
  double CVDigiE[10] = {0.};
  double CVDigiTime[10];
  double CVDigiHHTime[10];
  double CVIntegratedADC[10];
  double CVSignal[10];

  int nOEVDigi = 0;
  int OEVDigiID[44]={-1};
  double OEVDigiE[44] = {0.};
  double OEVDigiTime[44] ;
  double OEVDigiHHTime[44];
  double OEVIntegratedADC[44];
  double OEVSignal[44];

  int nCosmicDigi = 0;
  int CosmicDigiID[20]={-1};
  double CosmicDigiE[20] = {0.};
  double CosmicDigiTime[20] ;
  double CosmicDigiHHTime[20];
  double CosmicIntegratedADC[20];
  double CosmicSignal[20];

  int nLaserDigi = 0;
  int LaserDigiID[5]={-1};
  double LaserDigiE[5] = {0.};
  double LaserDigiTime[5] ;
  double LaserDigiHHTime[5];
  double LaserIntegratedADC[5];
  double LaserSignal[5];

  int nEtcDigi = 0;
  int EtcDigiID[16]={-1};
  double EtcDigiE[16] = {0.};
  double EtcDigiTime[16] ;
  double EtcDigiHHTime[16];
  double EtcIntegratedADC[16];
  double EtcSignal[16];

  trout->Branch("RunNumber"  ,&RunNumber   ,"RunNumber/I");
  trout->Branch("EventNumber",&EventNumber ,"EventNumber/I");
  trout->Branch("CsiNumber"  ,&nCSIDigi    ,"CsiNumber/I");
  trout->Branch("CsiModID"   ,CSIDigiID    ,"CsiModID[CsiNumber]/I");//nCSIDigi
  trout->Branch("CsiEne"     ,CSIDigiE     ,"CsiEne[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiTime"    ,CSIDigiTime  ,"CsiTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiHHTime"  ,CSIDigiHHTime,"CsiHHTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiSignal"  ,CSIDigiSignal,"CsiSignal[CsiNumber]/D");//nCSIDigi

  trout->Branch("CsiL1nTrig",&CSIL1nTrig,"CsiL1nTrig/I");
  trout->Branch("CsiL1TrigCount",CSIL1TrigCount,"CsiL1TrigCount[20]/D");

  trout->Branch("CC03Number"       ,&nCC03Digi,"CC03Number/I");
  trout->Branch("CC03ModID"        ,CC03DigiID,"CC03ModID[CC03Number]/I");//nCC03Digi
  trout->Branch("CC03Ene"          ,CC03DigiE ,"CC03Ene[CC03Number]/D");//nCC03Digi
  trout->Branch("CC03Time"         ,CC03DigiTime  ,"CC03Time[CC03Number]/D");//nCC03Digi
  trout->Branch("CC03HHTime"       ,CC03DigiHHTime,"CC03HHTime[CC03Number]/D");//nCC03Digi
  trout->Branch("CC03IntegratedADC",CC03IntegratedADC,"CC03IntegratedADC[CC03Number]/D");//nCC03Digi
  trout->Branch("CC03Signal"       ,CC03Signal,       "CC03Signal[CC03Number]/D");//nCC03Digi

  trout->Branch("CVNumber"       ,&nCVDigi,"CVNumber/I");
  trout->Branch("CVModID"        ,CVDigiID,"CVModID[CVNumber]/I");//nCVDigi
  trout->Branch("CVEne"          ,CVDigiE ,"CVEne[CVNumber]/D");//nCVDigi
  trout->Branch("CVTime"         ,CVDigiTime  ,"CVTime[CVNumber]/D");//nCVDigi
  trout->Branch("CVHHTime"       ,CVDigiHHTime,"CVHHTime[CVNumber]/D");//nCVDigi
  trout->Branch("CVIntegratedADC",CVIntegratedADC,"CVIntegratedADC[CVNumber]/D");//nCVDigi
  trout->Branch("CVSignal"       ,CVSignal       ,"CVSignal[CVNumber]/D");//nCVDigi

  trout->Branch("CosmicNumber"       ,&nCosmicDigi,"CosmicNumber/I");
  trout->Branch("CosmicModID"        ,CosmicDigiID,"CosmicModID[CosmicNumber]/I");//nCosmicDigi
  trout->Branch("CosmicEne"          ,CosmicDigiE ,"CosmicEne[CosmicNumber]/D");//nCosmicDigi
  trout->Branch("CosmicTime"         ,CosmicDigiTime  ,"CosmicTime[CosmicNumber]/D");//nCosmicDigi
  trout->Branch("CosmicHHTime"       ,CosmicDigiHHTime,"CosmicHHTime[CosmicNumber]/D");//nCosmicDigi
  trout->Branch("CosmicIntegratedADC",CosmicIntegratedADC,"CosmicIntegratedADC[CosmicNumber]/D");//nCosmicDigi
  trout->Branch("CosmicSignal"       ,CosmicSignal       ,"CosmicSignal[CosmicNumber]/D");//nCosmicDigi

  trout->Branch("EtcNumber"       ,&nEtcDigi,"EtcNumber/I");
  trout->Branch("EtcModID"        ,EtcDigiID,"EtcModID[EtcNumber]/I");//nEtcDigi
  trout->Branch("EtcEne"          ,EtcDigiE ,"EtcEne[EtcNumber]/D");//nEtcDigi
  trout->Branch("EtcTime"         ,EtcDigiTime  ,"EtcTime[EtcNumber]/D");//nEtcDigi
  trout->Branch("EtcHHTime"       ,EtcDigiHHTime,"EtcHHTime[EtcNumber]/D");//nEtcDigi
  trout->Branch("EtcIntegratedADC",EtcIntegratedADC,"EtcIntegratedADC[EtcNumber]/D");//nEtcDigi
  trout->Branch("EtcSignal"       ,EtcSignal       ,"EtcSignal[EtcNumber]/D");//nEtcDigi

  trout->Branch("OEVNumber"       ,&nOEVDigi,"OEVNumber/I");
  trout->Branch("OEVModID"        ,OEVDigiID,"OEVModID[OEVNumber]/I");//nOEVDigi
  trout->Branch("OEVEne"          ,OEVDigiE ,"OEVEne[OEVNumber]/D");//nOEVDigi
  trout->Branch("OEVTime"         ,OEVDigiTime  ,"OEVTime[OEVNumber]/D");//nOEVDigi
  trout->Branch("OEVHHTime"       ,OEVDigiHHTime,"OEVHHTime[OEVNumber]/D");//nOEVDigi
  trout->Branch("OEVIntegratedADC",OEVIntegratedADC,"OEVIntegratedADC[OEVNumber]/D");//nOEVDigi
  trout->Branch("OEVSignal"       ,OEVSignal       ,"OEVSignal[OEVNumber]/D");//nOEVDigi

  trout->Branch("LaserNumber"       ,&nLaserDigi,"LaserNumber/I");
  trout->Branch("LaserModID"        ,LaserDigiID,"LaserModID[LaserNumber]/I");//nLaserDigi
  trout->Branch("LaserEne"          ,LaserDigiE ,"LaserEne[LaserNumber]/D");//nLaserDigi
  trout->Branch("LaserTime"         ,LaserDigiTime  ,"LaserTime[LaserNumber]/D");//nLaserDigi
  trout->Branch("LaserHHTime"       ,LaserDigiHHTime,"LaserHHTime[LaserNumber]/D");//nLaserDigi
  trout->Branch("LaserIntegratedADC",LaserIntegratedADC,"LaserIntegratedADC[LaserNumber]/D");//nLaserDigi
  trout->Branch("LaserSignal"       ,LaserSignal       ,"LaserSignal[LaserNumber]/D");//nLaserDigi

  E14GNAnaDataContainer data;

  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  data.branchOfKlong( trout );

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  IDHandler* handler = new IDHandler();
  double x,y;
  double AlzPosition = 2624;//mm
  double sol = 299.792458;//mm/ns

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////


  E14WavReader_V1* reader = new E14WavReader_V1(trin);
  GammaFinder   gFinder;
  ClusterFinder clusterFinder;

  int nCsI               = 0;
  int CsIID[2716]        = {-1};
  double CsIEnergy[2716] = {0.};
  double CsISignal[2716] = {0.};
  double CsITime[2716]   = {-1};
  double CsIHHTime[2716] = {-1};
  
  //std::cout<< __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;
  Long_t entries =  reader->fChain->GetEntries();
  for( int ievent  = 0; ievent < entries ; ievent++){
    /// Initialize /// 
    nCsI = 0;
    nCSIDigi = 0;
    CSIL1nTrig=0;
    nOEVDigi = 0;
    nCC03Digi = 0;
    nCVDigi   = 0;
    nCosmicDigi = 0;
    nLaserDigi  = 0;
    nEtcDigi    = 0;

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
    for( int ich = 0; ich < 44; ich++){
      OEVDigiID[ich] = 0;
      OEVDigiE[ich] = 0.;
      OEVDigiTime[ich] = 0.;
      OEVDigiHHTime[ich] = 0.;
      OEVIntegratedADC[ich] = 0;
      OEVSignal[ich]= 0;
    }      
    for( int ich = 0; ich < 10; ich++){
      CVDigiID[ich] = 0;
      CVDigiE[ich] = 0.;
      CVDigiTime[ich] = 0.;
      CVDigiHHTime[ich] = 0.;
      CVIntegratedADC[ich] = 0;
      CVSignal[ich]= 0;
    }
    for( int ich =0; ich < 32; ich++){
      CC03DigiID[ich] = 0;
      CC03DigiE[ich] = 0.;
      CC03DigiTime[ich] = 0.;
      CC03DigiHHTime[ich] = 0.;
      CC03IntegratedADC[ich] = 0;
      CC03Signal[ich]= 0;
    }
    for( int ich = 0; ich < 5; ich++){
      LaserDigiID[ich] = 0;
      LaserDigiE[ich] = 0.;
      LaserDigiTime[ich] = 0.;
      LaserDigiHHTime[ich] = 0.;
      LaserIntegratedADC[ich] = 0;
      LaserSignal[ich]= 0;
    }
    for( int ich = 0; ich < 20; ich++){      
      CosmicDigiID[ich] = 0;
      CosmicDigiE[ich] = 0.;
      CosmicDigiTime[ich] = 0.;
      CosmicDigiHHTime[ich] = 0.;     
      CosmicIntegratedADC[ich] = 0;
      CosmicSignal[ich]= 0;
    }
    for( int ich = 0; ich < 16; ich++){
      EtcDigiID[ich] = 0;
      EtcDigiE[ich] = 0.;
      EtcDigiTime[ich] = 0.;
      EtcDigiHHTime[ich] = 0.;     
      EtcIntegratedADC[ich] = 0;
      EtcSignal[ich]= 0;
    }
    // End initialize // 
    /////////////////////////////////////////////////////////////

    if( (ievent%100) ==0 && ievent ){ std::cout<< ievent << std::endl;}
    reader->GetEntry( ievent  );
    data.reset();
    EventNumber = reader->EventNo;

    //if(( reader->TrigFlag & 3 ) != 0 ){ continue; }
    //if( reader->CsinTimeCluster == 0 ){ continue; }
    //if( reader->CsinTimeCluster  > 2  ){ continue ;}    
    if( reader->CsinTimeCluster == 2 && reader->CsiTimeClusterHead[0] < 50 ){ continue; }
    //// assumption :: All REAL csi Event if TimeCluster #0 
    //   std::cout<< "Analysis" << std::endl;

    for (int idigi = 0; idigi<reader->CVNumber; idigi++){
      CVSignal[idigi] = reader->CVSignal[idigi];
      CVDigiID[idigi] = reader->CVID[idigi];
      CVDigiE[idigi]  = reader->CVEne[idigi];
      CVDigiTime[idigi] = reader->CVTime[idigi];
      nCVDigi++;
    }
    //std::cout<< __PRETTY_FUNCTION__<< " : " << __LINE__ << std::endl;
    for (int idigi = 0; idigi<reader->EtcNumber; idigi++){
      EtcSignal[idigi] = reader->EtcSignal[idigi];
      EtcDigiID[idigi] = reader->EtcID[idigi];
      EtcDigiE[idigi]  = reader->EtcEne[idigi];
      EtcDigiTime[idigi] = reader->EtcTime[idigi];
      nEtcDigi++;
    }
    //std::cout<< __PRETTY_FUNCTION__<< " : " << __LINE__ << std::endl;
    for (int idigi = 0; idigi<reader->CC03Number; idigi++){
      CC03Signal[idigi] = reader->CC03Signal[idigi];
      CC03DigiID[idigi] = reader->CC03ID[idigi];
      CC03DigiE[idigi]  = reader->CC03Ene[idigi];
      CC03DigiTime[idigi] = reader->CC03Time[idigi];
      nCC03Digi++;
    }
    //std::cout<< __PRETTY_FUNCTION__<< " : " << __LINE__ << std::endl;
    for (int idigi = 0; idigi<reader->LaserNumber; idigi++){
      LaserSignal[idigi] = reader->LaserSignal[idigi];
      LaserDigiID[idigi] = reader->LaserID[idigi];
      LaserDigiE[idigi]  = reader->LaserEne[idigi];
      LaserDigiTime[idigi] = reader->LaserTime[idigi];
      nLaserDigi++;
    }
    for (int idigi = 0; idigi< reader->OEVNumber; idigi++){
      OEVSignal[idigi] = reader->OEVSignal[idigi];
      OEVDigiID[idigi] = reader->OEVID[idigi];
      OEVDigiE[idigi]  = reader->OEVEne[idigi];
      OEVDigiTime[idigi] = reader->OEVTime[idigi];
      nOEVDigi++;
    }
    for (int idigi = 0; idigi<reader->CosmicNumber; idigi++){
      CosmicSignal[idigi] = reader->CosmicSignal[idigi];
      CosmicDigiID[idigi] = reader->CosmicID[idigi];
      CosmicDigiE[idigi]  = reader->CosmicEne[idigi];
      CosmicDigiTime[idigi] = reader->CosmicTime[idigi];
      nCosmicDigi++;
    }
    
    for( int ich  = 0; ich < reader->CsiNumber; ich++){      
      int CsiID        = reader->CsiID[ich];
      double CsiTime   = reader->CsiTime[ich];
      double CsiSignal = reader->CsiSignal[ich]; 
      double CsiEnergy =0;
      switch( FileLevel ){
      case 0:
	CsiEnergy = reader->CsiEne[ich]/TempCorFactor;      
	break;
      case 1:
	CsiEnergy = reader->CsiEne[ich]*CalibrationFactor[reader->CsiID[ich]]/TempCorFactor;      
	break;
      case 2:
	CsiEnergy = reader->CsiEne[ich]*CalibrationFactor[reader->CsiID[ich]]/TempCorFactor*Pi0PeakCorFactor;
	break;
      default:
	break;
      }
      double CsiHHTime = reader->CsiHHTime[ich];
      int CsiTimeClusterID = reader->CsiTimeClusterID[ich];
      l1->Fill(CsiID, CsiSignal );
      if( CsiTimeClusterID == 0){
	CsIID[nCsI]     =  CsiID;
	CsISignal[nCsI] =  CsiSignal;
	CsIEnergy[nCsI] =  CsiEnergy;
	CsITime[nCsI]   =  CsiTime;
	CsIHHTime[nCsI] =  CsiHHTime;
	nCsI++;
      }
    }
    CSIL1nTrig = 0;
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
      if( CsISignal[idigi] > 5 && CsIEnergy[idigi]>0.5){
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
    //if( glist.size() ==6 ){ 
    if( user_rec(glist,klVec)){
      data.setData( clist );
      data.setData( glist );
      user_cut( data, klVec );
      data.setData(klVec);    
      trout->Fill();
    }
  }

  std::cout<< "End" << std::endl;
  trout->Write();
  std::cout<< "Write" << std::endl;
  tfout->Close();
  return 0;
}
