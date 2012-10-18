#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "E14WavReader.h"
#include <cstring>
#include <string>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "User_Function.h"
#include "EnergyConverter.h"

#include "ClusterFinder_EDIT.h"

int
main( int argc ,char ** argv ){
  
  //int RunNumber = atoi( argv[1]);
  std::string WAVFILE   = std::getenv("ROOTFILE_WAV");
  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  EnergyConverter *Converter = new EnergyConverter();  
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",
					  ANALIBDIR.c_str()));
  //TFile* tfin = new TFile(Form("%s/TEMPLATE_FIT_RESULT_1_%d.root",WAVFILE.c_str(),RunNumber));
  /*
  TTree* trin = (TTree*)tfin->Get("WFTree");
  TFile* tfout = new TFile( Form("%s/TEMPLATE_FIT_RESULT_1_%d_TIME.root",WAVFILE.c_str(), RunNumber),"RECREATE");
  */
  TChain* trin = new TChain("WFTree");
  for( int i = 0; i< 22; i++){
    trin->Add(Form("%s/TEMPLATE_FIT_RESULT_1_%d.root",WAVFILE.c_str(),4503+i));
  }
  TFile* tfout = new TFile("Pi0Out.root","RECREATE");
  TH2D* hisTimeDelta = new TH2D("hisTimeDelta","hisTimeDelta",2716, 0, 2716,
				400, -100, 100 );  
  TH2D* hisEnergyTimeDelta[2716];
  TH1D* hisTimeDeltaCH[2716];
  for( int i = 0; i< 2716; i++){
    hisTimeDeltaCH[i] = new TH1D(Form("hisTimeDelta%d",i),
				 Form("hisTimeDelta%d",i),
				 400,-100,100);
    hisEnergyTimeDelta[i] = new TH2D(Form("hisEnergyTimeDelta%d",i),
				     Form("hisEnergyTimeDelta%d",i),
				     40,0,16000,
				     400,-100,100);
  }

  E14WavReader* reader = new E14WavReader(trin);
  Long_t entries =  reader->fChain->GetEntries();
  
  TTree* trout = new TTree("trOut","");

  const int nCSI = 2716;
  Int_t    RunNumber;
  Int_t    EventNumber;
  Double_t ScintiSignal = 0;
  Double_t ScintiHHTime = -500.;
  Double_t ScintiTime   =-500.;
  Int_t    nCSIDigi     = 0;
  Double_t CSIDigiE[nCSI];
  Double_t CSIDigiTime[nCSI];
  Double_t CSIDigiHHTime[nCSI];
  Int_t    CSIDigiID[nCSI];
  Double_t CSIDigiSignal[nCSI];

  trout->Branch( "RunNumber"     , &RunNumber     , "RunNumber/I");
  trout->Branch( "EventNumber"   , &EventNumber   , "EventNumber/I");
  trout->Branch( "ScintiSignal"  , &ScintiSignal  , "ScintiSignal/D");
  trout->Branch( "ScintiHHTimne" , &ScintiHHTime  , "ScintiHHTime/D");
  trout->Branch( "ScintiTime"    , &ScintiTime    , "ScintiTime/D");
  trout->Branch( "nCSIDigi"      , &nCSIDigi     , "nCSIDigi/I" );
  trout->Branch( "CSIDigiE"      , CSIDigiE      , "CSIDidgiE[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CSIDigiTime"   , CSIDigiTime   , "CSIDigiTime[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CSIDigiHHTime" , CSIDigiHHTime , "CSIDigiHHTime[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CSIDigiID"     , CSIDigiID     , "CSIDigiID[nCSIDigi]/I");//nCSIDigi
  trout->Branch( "CSIDigiSignal" , CSIDigiSignal , "CSIDigiSignal[nCSIDigi]/D");//nCSIDigi

  E14GNAnaDataContainer data;
  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  data.branchOfPi0List( trout );
  
  ClusterFinder_EDIT clusterFinder;
  GammaFinder gFinder;


  TH1D* stepHist = new TH1D("hisStep","Step;Step;Survived Event",20,0,20);
  for( int ievent  = 0; ievent < entries ; ievent++){
    //std::cout << ievent << "/" << entries << std::endl;
    reader->GetEntry( ievent  );
    ScintiTime = -500;
    ScintiHHTime = -500.;
    nCSIDigi = 0;

    for( Int_t iCSI = 0; iCSI < nCSI; iCSI++ ){
      CSIDigiE[iCSI]      = 0;
      CSIDigiTime[iCSI]   = 0;
      CSIDigiHHTime[iCSI] = 0;
      CSIDigiID[iCSI]     = -1;
    }

    stepHist->Fill(0);   
    //std::cout<< "Scinti" << std::endl;
    RunNumber = reader->RunNo;
    EventNumber = reader->EventNo;
    Int_t    ScintiID;
    Bool_t   ScintiOn = false;
    for( int i = 0; i< reader->EtcNumber ; i++){    
      if( reader->EtcID[i] == 1 ){
	ScintiID     = reader->EtcID[i]; 
	ScintiHHTime = reader->EtcHHTime[i];
	ScintiTime   = reader->EtcTime[i];
	ScintiSignal = reader->EtcSignal[i];
	ScintiOn     = true;
	break;
      }
    }
    if( !ScintiOn || reader->EtcSignal[ScintiID] < 50 ){ continue; }
    stepHist->Fill(1);
    //std::cout << "CsI" << std::endl;

    for( int ich  = 0; ich < reader->CsiNumber; ich++){
      
      if( reader->CsiSignal[ich] < 5 ){ 
	continue;
      }

      //std::cout << "Converter " << std::endl;
      if( Converter->ConvertToEnergy( reader->CsiID[ich], reader->CsiSignal[ich] ) > 1.5){
	CSIDigiID[nCSIDigi]     = reader->CsiID[ich];
	CSIDigiTime[nCSIDigi]   = reader->CsiTime[ich];
	CSIDigiHHTime[nCSIDigi] = reader->CsiHHTime[ich];
	CSIDigiSignal[nCSIDigi]     = reader->CsiSignal[ich];
	CSIDigiE[nCSIDigi]      = Converter->ConvertToEnergy( reader->CsiID[ich] ,reader->CsiSignal[ich] );
	hisTimeDelta->Fill( CSIDigiID[nCSIDigi], CSIDigiHHTime[nCSIDigi] - ScintiHHTime );
	hisEnergyTimeDelta[CSIDigiID[nCSIDigi]]->Fill( CSIDigiSignal[nCSIDigi], CSIDigiHHTime[nCSIDigi]- ScintiHHTime );
	nCSIDigi++;
      }
    }
    
    //std::cout<< "Clustering" << std::endl;
    //std::cout << nCSIDigi << std::endl;
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    clist = clusterFinder.findCluster( nCSIDigi, CSIDigiID, CSIDigiE, CSIDigiHHTime);
    gFinder.findGamma( clist , glist );

    if( glist.size() != 2 ) continue;
    stepHist->Fill(2);
    //std::cout <<"Pi0 Reconstruction" << std::endl;
    std::list<Pi0> piList; 
    double mass =0;
    double position = 3526;
    if( !user_rec(glist,piList,mass,position )) continue; 
    stepHist->Fill(3);
    Gamma g1 = glist.front();
    Gamma g2 = glist.back();
    /*
    g1.clusterIdVec();
    g1.clusterEVec();
    g1.clusterTimeVec();
    */
    for( std::list<Gamma>::iterator it = glist.begin(); it != glist.end(); it++ ){ 
      for( int icluster = 0; icluster < (*it).clusterIdVec().size(); icluster++ ){
	if( (*it).clusterEVec()[icluster] > 30 ){
	  hisTimeDeltaCH[ (*it).clusterIdVec()[icluster] ]->Fill( (*it).clusterTimeVec()[icluster] - ScintiHHTime );
	}
      }
    }



    //std::cout<< mass <<std::endl; 
    data.setData( piList );
    trout->Fill();
  }

  for( int i = 0; i< 2716; i++){
    hisEnergyTimeDelta[i]->Write();
    hisTimeDeltaCH[i]->Write();
  }

  trout->Write();
  stepHist->Write();
  hisTimeDelta->Write();
  tfout->Close();
  return 0;
}
