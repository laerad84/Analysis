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

//Extra Libs
#include "IDHandler.h"
#include "E14WavReader_V1.h"
#include "L1TrigCounter.h"
#include "CrateIDHandler.h"
#include "CosmicTriggerTree.h"
#include "User_Function.h"

int
main( int argc ,char ** argv ){
  
  int RunNumber = atoi( argv[1]);
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string HOME         = std::getenv("HOME");

  std::string iFileForm="%s/run_wav_%d.root";
  std::string oFileForm="%s/run_wav_%d_Cal_FNL_COS_newTimeOffset.root";

  //std::string TCalFile = Form("%s/Data/TimeOffset/TimeOffset_with_cosmic.dat",ANALYSISLIB.c_str());  
  std::string TCalFile = Form("%s/Data/TimeOffset/testNewWORKCompileOffset.txt",ANALYSISLIB.c_str());  
  std::string ECalFile = Form("%s/local/Analysis/K3pi0Producer/Data/CalibrationFactorADV_15.dat",HOME.c_str());

  CsiMap*            map        = CsiMap::getCsiMap();
  CrateIDHandler*    CIDHandler = new CrateIDHandler();
  CosmicTriggerTree* cosmicTrig = new CosmicTriggerTree();
  ClusterFinder      cFinder;
  GammaFinder        gFinder;

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

  TChain* trin = new TChain("Tree"); 
  trin->Add(Form(iFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber));
  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber),"recreate");
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
    
    trout->Branch("GamClusNumbers",&GamClusNumbers,"GamClusNumbers/I");
    trout->Branch("GamClusSizes",&GamClusSizes,"GamClusSizes[GamClusNumbers]/I");  
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
  
  IDHandler* handler = new IDHandler();
  double x,y;
  double AlzPosition = 2624;//mm
  double sol = 299.792458;//mm/ns
  for( int i = 0 ; i< 2716; i++){
    handler->GetMetricPosition( i, x, y );
    double length = TMath::Sqrt( AlzPosition*AlzPosition + x*x + y*y );
    TimeDeltaLength[i] = length/sol;
  }

  L1TrigCounter* l1 = new L1TrigCounter();
  l1->ReadMapFile();
  l1->SetThreshold(1800);
  l1->Reset();

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
      l1->Reset();
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
    


    nCsI = 0;
    for( int ich  = 0; ich < reader->CsiNumber; ich++){      
      l1->Fill(reader->CsiID[ich], reader->CsiSignal[ich] );
      if( reader->CsiSignal[ich] > 5 && reader->CsiEne[ich]>0.5){
	CSIDigiID[CsiNumber]     = reader->CsiID[ich];
	CSIDigiE[CsiNumber]      = reader->CsiEne[ich]*CalibrationFactor[reader->CsiID[ich]]/TempCorFactor*Pi0PeakCorFactor;
	CSIDigiTime[CsiNumber]   = reader->CsiTime[ich]-TimeDelta[reader->CsiID[ich]];
	CSIDigiHHTime[CsiNumber] = reader->CsiHHTime[ich];
	CSIDigiSignal[CsiNumber] = reader->CsiSignal[ich];
	CsiChisq[CsiNumber]      = reader->CsiChisq[ich];
	CsiID[CsiNumber]         = reader->CsiID[ich];
	CsiNDF[CsiNumber]        = reader->CsiNDF[ich];
	double x(0), y(0);
	x = map->getX( reader->CsiID[ich] );
	y = map->getY( reader->CsiID[ich] );
	if( x > 0 ){
	  if( y > 0 ){ CsiPosID[CsiNumber]  = 0; }
	  else{ CsiPosID[CsiNumber] = 1; }
	}else{
	  if( y > 0 ){ CsiPosID[CsiNumber] = 2; }
	  else{ CsiPosID[CsiNumber]  =3; }
	}
	if( y > 0 ){
	  if( x < -100 ){
	    CsiGB[CsiNumber] = -1;
	  }else{
	    CsiGB[CsiNumber] = 1;
	  }
	}else{
	  if( x < 75 ){
	    CsiGB[CsiNumber] = -1;
	  }else{
	    CsiGB[CsiNumber] = 1;
	  }
	}

	CsiCrate[CsiNumber] = CIDHandler->GetCrate(reader->CsiID[ich]);
	CsiL1[CsiNumber]    = CIDHandler->GetL1(reader->CsiID[ich]);
	CsiNumber++;
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
    if( CsiNumber<5){ continue;}

    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    
    clist = cFinder.findCluster( CsiNumber,CSIDigiID,CSIDigiE,CSIDigiTime);
    gFinder.findGamma( clist, glist );
    if( clist.size() < 6 ){ continue; }
    if( glist.size() ==6 ){ 
      if( user_rec(glist,klVec)){
	data.setData( clist );
	data.setData( glist );
	user_cut( data, klVec );
	data.setData(klVec);    
	std::list<Gamma>::iterator git = glist.begin();
	int clNumber = 0;
	GamClusNumbers = glist.size();
	for(; git != glist.end(); git++){
	  GamClusSizes[clNumber] = (*git).clusterIdVec().size();
	  for( int i = 0; i < (*git).clusterIdVec().size(); i++){
	    for( int in = 0; in < CsiNumber; in++){
	      if( CsiID[in] == (*git).clusterIdVec()[i]){
		GamClusCsiSignal[clNumber][i] = CSIDigiSignal[in];
		GamClusCsiChisq[clNumber][i]  = CsiChisq[in];
		GamClusCsiCrate[clNumber][i]  = CsiCrate[in];
		GamClusCsiL1[clNumber][i]     = CsiL1[in];
		break;
	      }
	    }
	  }
	  clNumber++;
	}
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
