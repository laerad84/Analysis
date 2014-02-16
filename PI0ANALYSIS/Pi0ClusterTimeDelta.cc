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
#include "EnergyConverter.h"


int
main( int argc ,char ** argv ){
  
  int RunNumber = atoi(argv[1]);
  int Types     = atoi(argv[2]);
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string HOME         = std::getenv("HOME");

  std::string iFileForm="%s/run_wav_%d.root";
  std::string oFileForm;
  switch( Types ){
  case 0:
    oFileForm ="%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0_nocal_nopi0peak_notempcorr.root";
    break;
  case 1:
    oFileForm ="%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0_nopi0peak_notempcorr.root";
    break;
  case 2:
    oFileForm ="%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0.root";
    break;
  case 3:
    oFileForm ="%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0_nocal.root";
    break;
  case 4:
    oFileForm ="%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0_nopi0peak.root";
    break;
  case 5:
    oFileForm ="%s/run_wav_%d_Cal_FNL_COS_newCompensate_pi0.root";
    break;
  default : 
    return -1;
  }
  //std::string TCalFile = Form("%s/Data/TimeOffset/TimeOffset_with_cosmic.dat",ANALYSISLIB.c_str());  
  //std::string TCalFile = Form("%s/Data/TimeOffset/testNewWORKCompileOffset.txt",ANALYSISLIB.c_str());  
  std::string TCalFile = Form("~/local/Analysis/KLongSpectrum/Data/TimeOffset_Shower_10.dat");  
  std::string ECalFile = Form("%s/local/Analysis/K3pi0Producer/Data/CalibrationFactorADV_15.dat",HOME.c_str());
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
  Double_t Pi0PeakCorFactor = 0.9937;

  EnergyConverter* Converter = new EnergyConverter(1);
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",ANALYSISLIB.c_str()));
  TChain* trin = new TChain("Tree"); 
  trin->Add(Form(iFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber));
  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber),"recreate");
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
  trout->Branch("CsiModID"   ,CSIDigiID    ,"CsiModID[CsiNumber]/I");//CsiNumber
  trout->Branch("CsiEne"     ,CSIDigiE     ,"CsiEne[CsiNumber]/D");//CsiNumber
  trout->Branch("CsiTime"    ,CSIDigiTime  ,"CsiTime[CsiNumber]/D");//CsiNumber
  trout->Branch("CsiHHTime"  ,CSIDigiHHTime,"CsiHHTime[CsiNumber]/D");//CsiNumber
  trout->Branch("CsiSignal"  ,CSIDigiSignal,"CsiSignal[CsiNumber]/D");//CsiNumber

  trout->Branch("CsiL1nTrig",&CSIL1nTrig,"CsiL1nTrig/I");
  trout->Branch("CsiL1TrigCount",CSIL1TrigCount,"CsiL1TrigCount[20]/D");

  int     CVNumber      = 0;
  Short_t CVModID[100]  = {-1};
  double  CVEne[100] = {0};
  double  CVTime[100]   = {0};
  int    SciNumber;
  double SciEne[1];
  double SciTime[1];

  trout->Branch("SciNumber",&SciNumber,"SciNumber/I");
  trout->Branch("SciEne"   ,&SciEne  ,"SciEne/D");//SciNumber
  trout->Branch("SciTime"  ,&SciTime ,"SciTime/D");//SciNumber
  trout->Branch("CVNumber" ,&CVNumber,"CVNumber/I");
  trout->Branch("CVModID"  ,CVModID  ,"CVModID[CVNumber]/S");//CVNumber
  trout->Branch("CVEne"    ,CVEne    ,"CVEne[CVNumber]/D");//CVNumber
  trout->Branch("CVTime"   ,CVTime   ,"CVTime[CVNumber]/D");//CVNumber


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
  for( int i = 0; i< 2716; i++){
    CalibrationFactor[i] = 1;
  }
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
  
  double AlzPosition = 3526;//mm

  L1TrigCounter* l1 = new L1TrigCounter();
  l1->ReadMapFile();
  l1->SetThreshold(1800);
  l1->Reset();

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
  
  std::cout<< __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;
  Long_t entries =  reader->fChain->GetEntries();
  for( int ievent  = 0; ievent < entries ; ievent++){
    nCsI = 0;
    nCSIDigi = 0;
    CSIL1nTrig=0;
    l1->Reset();
    for( int ich = 0; ich < 2716; ich++){
      CsIID[ich]     = -1; 
      CsIEnergy[ich] = 0.;
      CsITime[ich]   = 0;
      CsIHHTime[ich] = 0;
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
    //if(( reader->TrigFlag & 3 ) != 0 ){ continue; }
    //if( reader->CsinTimeCluster == 0 ){ continue; }
    //if( reader->CsinTimeCluster  > 2  ){ continue ;}    
    //if( reader->CsinTimeCluster == 2 && reader->CsiTimeClusterHead[0] < 50 ){ continue; }
    //// assumption :: All REAL csi Event if TimeCluster #0 
    //    std::cout<< "Analysis" << std::endl;

    nCsI = 0;
    for( int ich  = 0; ich < reader->CsiNumber; ich++){      
      int CsiID        = reader->CsiID[ich];
      double CsiTime   = reader->CsiTime[ich];
      double CsiSignal = reader->CsiSignal[ich]; 
      double CsiEnergy =0;
      switch( Types ){
      case 0:
	CsiEnergy = reader->CsiEne[ich];
	break;
      case 1:
	CsiEnergy = reader->CsiEne[ich]*CalibrationFactor[reader->CsiID[ich]];//TempCorFactor*Pi0PeakorFactor;      
	break;
      case 2:
	CsiEnergy = reader->CsiEne[ich]*CalibrationFactor[reader->CsiID[ich]]/TempCorFactor*Pi0PeakCorFactor;      
	break;
      case 3:
	CsiEnergy = reader->CsiEne[ich]/TempCorFactor*Pi0PeakCorFactor;
	break;
      case 4:
	CsiEnergy = reader->CsiEne[ich]*CalibrationFactor[reader->CsiID[ich]]/TempCorFactor;
	break;
      case 5:
	CsiEnergy =  Converter->ConvertToEnergy( CsiID, CsiSignal)/TempCorFactor*Pi0PeakCorFactor*CalibrationFactor[reader->CsiID[ich]];
	//std::cout << CsiID << "\t" << CsiSignal << "\t " << CsiEnergy << std::endl;
	break;
      default :
	return -1;
      }

      double CsiHHTime = reader->CsiHHTime[ich];
      int CsiTimeClusterID = reader->CsiTimeClusterID[ich];
      l1->Fill(CsiID, CsiSignal );
      //if( CsiTimeClusterID == 0){
      CsIID[nCsI]     =  CsiID;
      CsISignal[nCsI] =  CsiSignal;
      CsIEnergy[nCsI] =  CsiEnergy;
      CsITime[nCsI]   =  CsiTime;
      CsIHHTime[nCsI] =  CsiHHTime;
      nCsI++;
      //}
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
	//CSIDigiTime[ nCSIDigi ]   = CsITime[idigi]+TimeDelta[CsIID[idigi]];
	CSIDigiTime[ nCSIDigi ]   = CsITime[idigi]-TimeDelta[CsIID[idigi]];
	CSIDigiHHTime[ nCSIDigi ] = CsIHHTime[idigi];
	CSIDigiSignal[nCSIDigi]   = CsISignal[idigi];
	nCSIDigi++;
	//std::cout<< CSIDigiID[nCSIDigi-1]  << "\t" << CSIDigiE[nCSIDigi-1] << std::endl;
      }
    }
    
    
    if( nCSIDigi<5){ continue;}

    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::list<Pi0>     plist;
    clist = clusterFinder.findCluster( nCSIDigi,CSIDigiID,CSIDigiE,CSIDigiTime);
    gFinder.findGamma( clist, glist );
    if( clist.size() < 2 ){ continue; }
    if( glist.size() != 2 ){ continue; }



    /////////////////////////////////////////////////////////////
    // Copy CV & Sci Data
    /////////////////////////////////////////////////////////////
    for( int icv = 0; icv < 100; icv++){
      CVModID[icv] = -1;
      CVEne[icv] = 0; 
      CVTime[icv] = 0;
    }
    CVNumber  = 0;
    SciNumber = 0;
    SciEne[0] = 0;
    SciTime[0] = 0;
    for( int icv  =0; icv < reader->CVNumber; icv++){
      CVTime[icv]   = reader->CVTime[icv];
      CVEne[icv]    = reader->CVSignal[icv];
      CVModID[icv]  = reader->CVID[icv];
      CVNumber++;
    }
    for( int ir = 0; ir < reader->EtcNumber; ir++){
      if( reader->EtcID[ir]==1 ){
	SciNumber = 1;
	SciEne[0] = reader->EtcSignal[ir];
	SciTime[0] = reader->EtcTime[ir];
      }
      break;
    }
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    if( user_rec(glist,plist)){
      std::list<Pi0>::iterator it = plist.begin();
      (*it).setVtx(0,0,AlzPosition);
      (*it).updateVars();
      user_cut(data,plist);
      data.setData(plist);      
      trout->Fill();
      data.eventID++;
    }
  }

  std::cout<< "End" << std::endl;
  trout->Write();
  std::cout<< "Write" << std::endl;
  tfout->Close();
  return 0;
}
