//Standard Libs
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

//E14 Libs
#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"
#include "csimap/CsiMap.h"

//ROOT & CLHEP Libs
#include "CLHEP/Vector/ThreeVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"

//Extra Libs
#include "IDHandler.h"
#include "E14WavReader_V1.h"
//#include "L1TrigCounter.h"
//#include "CrateIDHandler.h"
#include "CosmicTriggerTree.h"
#include "User_Function.h"
#include "User_Functions.h"
#include "GeneralFunctions.h"
#include "EnergyConverter.h"
int
main( int argc ,char ** argv ){
  
  int RunNumber;// = atoi( argv[1]);
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string HOME         = std::getenv("HOME");

  std::string iFileForm="%s/run_wav_%d.root";
  std::string oFileForm="%s/run_wav_All_GammaTime_2G_NewComp.root";
  std::string Pi0RunList = Form("%s/local/Analysis/RunList/Pi0RunList.txt",HOME.c_str());
  std::ifstream ifs( Pi0RunList.c_str());
  if( !ifs.is_open() ){ return -1;}
  std::vector<int> runList;
  int tmpRun;
  while( ifs >> tmpRun ){
    runList.push_back( tmpRun );
  }


  std::string TCalFile = Form("%s/Data/TimeOffset/testNewWORKCompileOffset.txt",ANALYSISLIB.c_str());  
  std::string ECalFile = Form("%s/local/Analysis/K3pi0Producer/Data/CalibrationFactorADV_15.dat",HOME.c_str());


  CsiMap*            map        = CsiMap::getCsiMap();
  CosmicTriggerTree* cosmicTrig = new CosmicTriggerTree();
  ClusterFinder      cFinder;
  GammaFinder        gFinder;
  EnergyConverter*    Converter = new EnergyConverter(3);
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",
					  ANALYSISLIB.c_str()));
  ///////////////////////////////////////////////////////
  GammaCut* gammaCut = new GammaCut();
  CsiCut*   csiCut   = new CsiCut();

  std::string TempCalibrationFilename = Form("%s/Data/Temperature_Factor/TemperatureCorrectionFactor.root",ANALYSISLIB.c_str());  
  TFile* tfTempCorr  =new TFile(TempCalibrationFilename.c_str());
  TTree* trTempCorr  =(TTree*)tfTempCorr->Get("TemperatureCorrectionCsI");
  Double_t TempCorFactor=0;
  trTempCorr->SetBranchAddress("CorrectionFactor",&TempCorFactor);
  trTempCorr->GetEntry(RunNumber);
  if( TempCorFactor == 0){
    TempCorFactor = 1;
  }

  TF1* THCorFunc = new TF1("THCorFunc",THCorrectionFunction,0,25000,3);
  THCorFunc->SetParameters(0, 1.672,0.0319651);
  
  std::cout<< TempCorFactor << std::endl;
  //Double_t Pi0PeakCorFactor = 0.9937;//Original
  //Double_t Pi0PeakCorFactor = 1;//20140316 RAW
  Double_t Pi0PeakCorFactor = 1.001011;//20140316 Fit

  TChain* trin = new TChain("Tree"); 
  //trin->Add(Form(iFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber));
  for( int i = 0; i< runList.size() ; i++){
    trin->Add(Form(iFileForm.c_str(),ROOTFILE_WAV.c_str(),runList.at(i)));
  }


  //TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber),"recreate");
  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_WAV.c_str()),"recreate");
  TTree* trout = new TTree("T", "Output from Time zero" );  
  
  int EventNumber;
  
  int CsiNumber = 0;
  int CSIDigiID[2716]={-1};
  double CSIDigiE[2716] = {0.};
  double CSIDigiTime[2716] ;
  double CSIDigiHHTime[2716];
  double CSIDigiSignal[2716]={0};
  double CSIL1TrigCount[20];
  int    CSIL1nTrig;
  double CsiChisq[2716];
  short  CsiID[2716];
  short  CsiNDF[2716];
  short  CsiCrate[2716];
  short  CsiL1[2716];
  short  CsiGB[2716];
  short  CsiPosID[2716];
  
  Int_t LaserNumber;
  Double_t LaserSignal[10];
  Double_t LaserHHTime[10];
  Double_t LaserTime[10];
  Double_t LaserChisq[10];
  Short_t  LaserID[10];
  Short_t  LaserNDF[10];
  
  Int_t OEVNumber;
  Double_t OEVSignal[50];
  Double_t OEVTime[50];
  Double_t OEVChisq[50];
  Short_t  OEVID[50];
  Short_t  OEVNDF[50];

  Int_t CC03Number;
  Double_t CC03Signal[50];
  Double_t CC03Time[50];
  Double_t CC03Chisq[50];
  Short_t  CC03ID[50];
  Short_t  CC03NDF[50];

  Int_t CVNumber;
  Double_t CVSignal[50];
  Double_t CVTime[50];
  Double_t CVEne[50];
  Double_t CVChisq[50];
  Short_t  CVID[50];
  Short_t  CVNDF[50];

  Int_t CosmicNumber;
  Double_t CosmicSignal[20];
  Double_t CosmicTime[20];
  Double_t CosmicChisq[20];
  Short_t  CosmicID[20];
  Short_t  CosmicNDF[20];
  
  Int_t EtcNumber;
  Double_t EtcSignal[10];
  Double_t EtcTime[10];
  Double_t EtcHHTime[10];
  Double_t EtcChisq[10];
  Short_t  EtcID[10];
  Short_t  EtcNDF[10];

  Int_t SciNumber;
  Double_t SciSignal[10];
  Double_t SciEne[10];
  Double_t SciTime[10];
  Double_t SciHHTime[10];
  Double_t SciChisq[10];
  Short_t  SciID[10];
  Short_t  SciNDF[10];

  int s_arrSize = 120;
  Int_t    GamClusNumbers;
  Int_t    GamClusSizes[120];
  Double_t GamClusCsiSignal[120][120];
  Double_t GamClusCsiChisq[120][120];
  Int_t    GamClusCsiL1[120][120];
  Int_t    GamClusCsiCrate[120][120];

  ;;
  //Branch
  {
    trout->Branch("RunNumber"  ,&RunNumber   ,"RunNumber/I");
    trout->Branch("EventNumber",&EventNumber ,"EventNumber/I");
    csiCut->Branch(trout);
    gammaCut->Branch(trout);
    /*
    trout->Branch("CsiNumber"  ,&CsiNumber    ,"CsiNumber/I");
    trout->Branch("CsiModID"   ,CSIDigiID    ,"CsiModID[CsiNumber]/I");//CsiNumber
    trout->Branch("CsiEne"     ,CSIDigiE     ,"CsiEne[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiTime"    ,CSIDigiTime  ,"CsiTime[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiHHTime"  ,CSIDigiHHTime,"CsiHHTime[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiSignal"  ,CSIDigiSignal,"CsiSignal[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiChisq"   ,CsiChisq     ,"CsiChisq[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiNDF"     ,CsiNDF       ,"CsiNDF[CsiNumber]/S");//CsiNumber
    trout->Branch("CsiID"      ,CsiID        ,"CsiID[CsiNumber]/S");//CsiNumber
    trout->Branch("CsiPosID"   ,CsiPosID     ,"CsiPosID[CsiNumber]/S");//CsiNumber
    trout->Branch("CsiGB"      ,CsiGB        ,"CsiGB[CsiNumber]/S");//CsiNumber
    trout->Branch("CsiCrate"   ,CsiCrate     ,"CsiCrate[CsiNumber]/S");//CsiNumber
    
    trout->Branch("CsiL1nTrig" ,&CSIL1nTrig  ,"CsiL1nTrig/I");
    trout->Branch("CsiL1TrigCount",CSIL1TrigCount,"CsiL1TrigCount[20]/D");
    */
    trout->Branch("GamClusNumbers",&GamClusNumbers,"GamClusNumbers/I");
    trout->Branch("GamClusSizes",GamClusSizes,"GamClusSizes[GamClusNumbers]/I");  
    trout->Branch("GamClusCsiSignal",GamClusCsiSignal,Form("GamClusCsiSignal[GamClusNumbers][%d]/D",s_arrSize));
    trout->Branch("GamClusCsiChisq",GamClusCsiChisq,Form("GamClusCsiChisq[GamClusNumbers][%d]/D",s_arrSize));
    trout->Branch("GamClusCsiL1",GamClusCsiL1,Form("GamClusCsiL1[GamClusNumbers][%d]/I",s_arrSize));
    trout->Branch("GamClusCsiCrate",GamClusCsiCrate,Form("GamClusCsiCrate[GamClusNumbers][%d]/I",s_arrSize));
    
    trout->Branch("OEVNumber",&OEVNumber,"OEVNumber/I");
    trout->Branch("OEVSignal",OEVSignal,"OEVSignal[OEVNumber]/D");//OEVNumber
    trout->Branch("OEVTime",OEVTime,"OEVTime[OEVNumber]/D");//OEVNumber
    trout->Branch("OEVChisq",OEVChisq,"OEVChisq[OEVNumber]/D");//OEVNumber
    trout->Branch("OEVNDF",OEVNDF,"OEVNDF[OEVNumber]/S");//OEVNumber
    trout->Branch("OEVID",OEVID,"OEVID[OEVNumber]/S");//OEVNumber
    
    trout->Branch("CC03Number",&CC03Number,"CC03Number/I");
    trout->Branch("CC03Signal",CC03Signal,"CC03Signal[CC03Number]/D");//CC03Number
    trout->Branch("CC03Time",CC03Time,"CC03Time[CC03Number]/D");//CC03Number
    trout->Branch("CC03Chisq",CC03Chisq,"CC03Chisq[CC03Number]/D");//CC03Number
    trout->Branch("CC03NDF",CC03NDF,"CC03NDF[CC03Number]/S");//CC03Number
    trout->Branch("CC03ID",CC03ID,"CC03ID[CC03Number]/S");//CC03Number
    
    trout->Branch("CVNumber",&CVNumber,"CVNumber/I");
    trout->Branch("CVSignal",CVSignal,"CVSignal[CVNumber]/D");//CVNumber
    trout->Branch("CVTime",CVTime,"CVTime[CVNumber]/D");//CVNumber
    trout->Branch("CVChisq",CVChisq,"CVChisq[CVNumber]/D");//CVNumber
    trout->Branch("CVNDF",CVNDF,"CVNDF[CVNumber]/S");//CVNumber
    trout->Branch("CVID",CVID,"CVID[CVNumber]/S");//CVNumber
    
    trout->Branch("EtcNumber",&EtcNumber,"EtcNumber/I");
    trout->Branch("EtcSignal",EtcSignal,"EtcSignal[EtcNumber]/D");//EtcNumber
    trout->Branch("EtcTime",EtcTime,"EtcTime[EtcNumber]/D");//EtcNumber
    trout->Branch("EtcHHTime",EtcHHTime,"EtcHHTime[EtcNumber]/D");//EtcNumber
    trout->Branch("EtcChisq",EtcChisq,"EtcChisq[EtcNumber]/D");//EtcNumber
    trout->Branch("EtcID",EtcID,"EtcID[EtcNumber]/S");//EtcNumber
    trout->Branch("EtcNDF",EtcNDF,"EtcNDF[EtcNumber]/S");//EtcNumber  

    trout->Branch("SciNumber",&SciNumber,"SciNumber/I");
    trout->Branch("SciSignal",SciSignal,"SciSignal[SciNumber]/D");//SciNumber
    trout->Branch("SciTime",SciTime,"SciTime[SciNumber]/D");//SciNumber
    trout->Branch("SciHHTime",SciHHTime,"SciHHTime[SciNumber]/D");//SciNumber
    trout->Branch("SciChisq",SciChisq,"SciChisq[SciNumber]/D");//SciNumber
    trout->Branch("SciID",SciID,"SciID[SciNumber]/S");//SciNumber
    trout->Branch("SciNDF",SciNDF,"SciNDF[SciNumber]/S");//SciNumber  
    
    trout->Branch("LaserNumber",&LaserNumber,"LaserNumber/I");
    trout->Branch("LaserSignal",LaserSignal,"LaserSignal[LaserNumber]/D");//LaserNumber
    trout->Branch("LaserTime",LaserTime,"LaserTime[LaserNumber]/D");//LaserNumber
    trout->Branch("LaserHHTime",LaserHHTime,"LaserHHTime[LaserNumber]/D");//LaserNumber
    trout->Branch("LaserChisq",LaserChisq,"LaserChisq[LaserNumber]/D");//LaserNumber
    trout->Branch("LaserID",LaserID,"LaserID[LaserNumber]/S");//LaserNumber
    trout->Branch("LaserNDF",LaserNDF,"LaserNDF[LaserNumber]/S");//LaserNumber  
    
  }
  cosmicTrig->Branch(trout);
  double CSIL1TrigCountThreshold[20] = {1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
					1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};

  E14GNAnaDataContainer data; 

  data.branchOfClusterList( trout );
  data.branchOfPi0List( trout );
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
  
  IDHandler* handler = new IDHandler();
  double x,y;
  //double AlzPosition = 2624;//mm
  double AlzPosition = 3526;//mm
  double sol = 299.792458;//mm/ns
  for( int i = 0 ; i< 2716; i++){
    handler->GetMetricPosition( i, x, y );
    double length = TMath::Sqrt( AlzPosition*AlzPosition + x*x + y*y );
    TimeDeltaLength[i] = length/sol;
  }




  /*
  L1TrigCounter* l1 = new L1TrigCounter();
  l1->ReadMapFile();
  l1->SetThreshold(1800);
  l1->Reset();
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  E14WavReader_V1*   reader = new E14WavReader_V1(trin);

  int nCsI               = 0;
  int CsIID[2716]        = {-1};
  double CsIEnergy[2716] = {0.};
  double CsISignal[2716] = {0.};
  double CsITime[2716]   = {-1};
  double CsIHHTime[2716] = {-1};
  
  std::cout<< __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;
  Long_t entries =  reader->fChain->GetEntries();
  for( int ievent  = 0; ievent < entries ; ievent++){
    RunNumber   = reader->RunNo;
    EventNumber = reader->EventNo;
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////    



    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

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

    ///////////////////////////////////////////////
    // init 
    ///////////////////////////////////////////////
    {

      nCsI = 0;
      CsiNumber = 0;
      CSIL1nTrig=0;
      for( int ich = 0; ich< 20; ich++){
	CSIL1TrigCount[ich] = 0;
      }
      //l1->Reset();
      for( int ich = 0; ich < 2716; ich++){
	CsIID[ich]     = -1; 
	CsIEnergy[ich] = 0.;
	CsITime[ich]   = 0;
	CsIHHTime[ich] = 0;
	CsISignal[ich] = 0.;
	
	CSIDigiID[ich]     = 0;
	CSIDigiE[ich]      = 0;
	CSIDigiTime[ich]   = 0;
	CSIDigiHHTime[ich] = 0;
	CSIDigiSignal[ich] = 0.;
	
	CsiChisq[ich]=0;
	CsiID[ich]=-1;
	CsiNDF[ich]=0;
	CsiCrate[ich]=0;
	CsiL1[ich]=-1;
	CsiGB[ich]=0;
	CsiPosID[ich]=-1;	
      }

      LaserNumber = 0;
      for( int ich  =0; ich < 10; ich++){
	LaserSignal[ich] = 0;
	LaserHHTime[ich] = 0;
	LaserTime[ich]   = 0;
	LaserChisq[ich]  = 0;
	LaserID[ich]     = -1;
	LaserNDF[ich]    = 0;
      }
      OEVNumber =0;
      for( int ich = 0; ich < 50; ich++){
	OEVSignal[ich] = 0;
	OEVTime[ich]   = 0;
	OEVID[ich]     = -1;
	OEVChisq[ich]  = 0;
	OEVNDF[ich]    = 0;
      }
      CC03Number =0;
      for( int ich = 0; ich < 50; ich++){
	CC03Signal[ich] = 0;
	CC03Time[ich]   = 0;
	CC03ID[ich]     = -1;
	CC03Chisq[ich]  = 0;
	CC03NDF[ich]    = 0;
      }
      CVNumber =0;
      for( int ich = 0; ich < 50; ich++){
	CVSignal[ich] = 0;
	CVTime[ich]   = 0;
	CVID[ich]     = -1;
	CVChisq[ich]  = 0;
	CVNDF[ich]    = 0;
      }
      /*
      CosmicNumber =0;
      for( int ich = 0; ich < 20; ich++){
	CosmicSignal[ich] = 0;
	CosmicTime[ich]   = 0;
	CosmicID[ich]     = -1;
	CosmicChisq[ich]  = 0;
	CosmicNDF[ich]    = 0;
      }
      */
      EtcNumber =0;
      for( int ich = 0; ich < 10; ich++){
	EtcSignal[ich] = 0;
	EtcTime[ich]   = 0;
	EtcID[ich]     = -1;
	EtcChisq[ich]  = 0;
	EtcNDF[ich]    = 0;
      }
      SciNumber =0;
      for( int ich = 0; ich < 10; ich++){
	SciSignal[ich] = 0;
	SciTime[ich]   = 0;
	SciID[ich]     = -1;
	SciChisq[ich]  = 0;
	SciNDF[ich]    = 0;
      }
    }
    ///////////////////////////////////////////////


    for( int ich = 0; ich < reader->CVNumber; ich++){
      CVSignal[ich] = reader->CVSignal[ich];
      CVID[ich]     = reader->CVID[ich];
      CVTime[ich]   = reader->CVTime[ich];
      CVChisq[ich]  = reader->CVChisq[ich];
      CVNDF[ich]    = reader->CVNDF[ich];
      CVNumber++;
    }
    for( int ich = 0; ich < reader->LaserNumber; ich++){
      LaserSignal[ich] = reader->LaserSignal[ich];
      LaserID[ich]     = reader->LaserID[ich];
      LaserTime[ich]   = reader->LaserTime[ich];
      LaserChisq[ich]  = reader->LaserChisq[ich];
      LaserNDF[ich]    = reader->LaserNDF[ich];
      LaserNumber++;
    }
    for( int ich = 0; ich < reader->OEVNumber; ich++){
      OEVSignal[ich] = reader->OEVSignal[ich];
      OEVID[ich]     = reader->OEVID[ich];
      OEVTime[ich]   = reader->OEVTime[ich];
      OEVChisq[ich]  = reader->OEVChisq[ich];
      OEVNDF[ich]    = reader->OEVNDF[ich];
      OEVNumber++;
    }
    for( int ich = 0; ich < reader->CC03Number; ich++){
      CC03Signal[ich] = reader->CC03Signal[ich];
      CC03ID[ich]     = reader->CC03ID[ich];
      CC03Time[ich]   = reader->CC03Time[ich];
      CC03Chisq[ich]  = reader->CC03Chisq[ich];
      CC03NDF[ich]    = reader->CC03NDF[ich];
      CC03Number++;
    }
    
    cosmicTrig->InitVar();
    for( int ich = 0; ich < reader->CosmicNumber; ich++){
      /*
      CosmicSignal[ich] = reader->CosmicSignal[ich];
      CosmicID[ich]     = reader->CosmicID[ich];
      CosmicTime[ich]   = reader->CosmicTime[ich];
      CosmicChisq[ich]  = reader->CosmicChisq[ich];
      CosmicNDF[ich]    = reader->CosmicNDF[ich];
      CosmicNumber++;
      */
      cosmicTrig->CosmicSignal[ich]=reader->CosmicSignal[ich];
      cosmicTrig->CosmicTime[ich]  =reader->CosmicTime[ich];
      cosmicTrig->CosmicID[ich]    =reader->CosmicID[ich];
    }
    cosmicTrig->CosmicNumber = reader->CosmicNumber;
    cosmicTrig->TriggerDecision();

    for( int ich = 0; ich < reader->EtcNumber; ich++){
      EtcSignal[ich] = reader->EtcSignal[ich];
      EtcID[ich]     = reader->EtcID[ich];
      EtcTime[ich]   = reader->EtcTime[ich];
      EtcChisq[ich]  = reader->EtcChisq[ich];
      EtcNDF[ich]    = reader->EtcNDF[ich];
      EtcNumber++;
    }
    for( int ich = 0; ich < reader->EtcNumber; ich++){
      if( reader->EtcID[ich] == 1){
	SciNumber =1;
	SciEne[0] = reader->EtcSignal[ich];
	SciTime[0] = reader->EtcTime[ich];
	SciSignal[0] = reader->EtcSignal[ich];	  
	break;
      }
    }
    GamClusNumbers = 0;
    for( int i = 0; i< 120; i++){
      GamClusSizes[i] = 0;
      for( int j = 0; j< 120; j++){
	GamClusCsiSignal[i][j] = 0;
	GamClusCsiChisq[i][j]  = 0;
	GamClusCsiL1[i][j]     = 0;
	GamClusCsiCrate[i][j]  = -1;
      }
    }

    if( LaserSignal[0] > 50 ){ continue; }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Convert
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    nCsI = 0;
    for( int ich  = 0; ich < reader->CsiNumber; ich++){      
      if( reader->CsiSignal[ich] > 5 && reader->CsiEne[ich]>0.5){
	CSIDigiID[CsiNumber]     = reader->CsiID[ich];
	//CSIDigiE[CsiNumber]    = reader->CsiEne[ich]*Pi0PeakCorFactor/TempCorFactor*CalibrationFactor[reader->CsiID[ich]]
	CSIDigiE[CsiNumber]      = Converter->ConvertToEnergy(reader->CsiID[ich], reader->CsiSignal[ich])*Pi0PeakCorFactor/TempCorFactor*CalibrationFactor[reader->CsiID[ich]];//ForNewComp
	CSIDigiTime[CsiNumber]   = reader->CsiTime[ich]-TimeDelta[reader->CsiID[ich]]-THCorFunc->Eval(reader->CsiSignal[ich]);
	CSIDigiHHTime[CsiNumber] = reader->CsiHHTime[ich];
	CSIDigiSignal[CsiNumber] = reader->CsiSignal[ich];
	CsiChisq[CsiNumber]      = reader->CsiChisq[ich];
	CsiID[CsiNumber]         = reader->CsiID[ich];
	CsiNDF[CsiNumber]        = reader->CsiNDF[ich];
	CsiNumber++;
      }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( CsiNumber<4){ continue;}

    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::list<Gamma>   glistTCut;
    std::list<Pi0>     plist;

    csiCut->DecisionForPi0Run( CsiNumber, CSIDigiID, CSIDigiE,CSIDigiTime, CSIDigiSignal, CsiChisq,CsiNDF);
    //clist = cFinder.findCluster( CsiNumber,CSIDigiID,CSIDigiE,CSIDigiTime);
    clist = cFinder.findCluster( csiCut->CsiNumber,csiCut->CsiID,csiCut->CsiEne,csiCut->CsiTime);
    gFinder.findGamma( clist, glist );
    gammaCut->Decision( glist );
    SetGammaTime( glist );


    //////////////////////    //////////////////////
    //////////////////////    //////////////////////
    //////////////////////    //////////////////////
    //////////////////////    //////////////////////



    if( clist.size() < 2 ){ continue; }
    if( glist.size() < 2 ){ continue; }
    
    //std::cout << "MaxDeltaT:" << gammaCut->GMaxDeltaT << std::endl;
    //std::cout << "MinDist  :" << gammaCut->GMinDist   << std::endl;
    std::list<Gamma>::iterator git = glist.begin();
    for( ; git != glist.end(); git++){
      SetGammaTime( (*git));
    }
    //GammaTimeDeltaCut( glist, glistTCut,2);
    GammaTimeDeltaCutEventTime( glist, glistTCut,csiCut->CsiEventTime,3);
    data.setData( clist );
    data.setData( glistTCut );
    std::list<Gamma>::iterator gitT = glistTCut.begin();
    //std::cout<<glistTCut.size() << "\t" <<  (*gitT).t() << std::endl;
    /*
    if( glist.size() < 10 ){ 
      if (User_RecG6( glist, klVec ) ){
	data.setData( klVec );
      }
    }
    */
    if( glistTCut.size() ==2 ){
      if( User_RecG2(glist,plist)){
	std::list<Pi0>::iterator it = plist.begin();
	(*it).setVtx( 0,0,AlzPosition );
	(*it).updateVars();
	data.setData( clist );
	data.setData( glist );
	//user_cut( data, klVec );
	data.setData(plist);    
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
