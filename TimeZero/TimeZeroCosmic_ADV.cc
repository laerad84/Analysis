// TimeZeroTestWithCosmic_V1
// Fit TimeData of Cosmic event and calcultate Delta of Fitted - Data 
// 

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
#include "TF1.h"
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
#include "IDHandler.h"
#include "CsIImage.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "HoughCsI.h"
#include "Chisq_cosmic.h"
#include "E14CosmicAnalyzer.h"

#include "TChain.h"
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// TimeZeroCosmic [RunNumberList] [IterationNumber]
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

int
main( int argc ,char ** argv ){

  gStyle->SetOptFit(111111111);

  std::string WAVFILE   = std::getenv("ROOTFILE_WAV");
  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  std::string ANALYSISDIR = std::getenv("ANALYSISDIR");
  std::string DIRNAME   = "TimeZero";
  std::string ROOTFILE_COSMIC= std::getenv("ROOTFILE_COSMIC");
  EnergyConverter *Converter = new EnergyConverter();  
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",
					  ANALIBDIR.c_str()));
   
  TChain* trin = new TChain("Tree");
  int fRunNumber = atoi(argv[1]);
  int IterationNumber = atoi( argv[2]);  
  std::cout << "////////////////////////////////////////////////////\n";
  std::cout << "RunNumber       : " << fRunNumber << "\n";
  std::cout << "IterationNumber : " << IterationNumber << "\n";
  std::cout << "////////////////////////////////////////////////////\n";    
  trin->Add(Form("%s/run_wav_%d.root",WAVFILE.c_str(),fRunNumber));  
  
  IDHandler* handler        = new IDHandler();
  TH2D* TriggerMap          = new TH2D("Trigger","Trigger",5,0,5,5,0,5);
  CsIImage* TimeMap         = new CsIImage( handler );
  CsIImage* EnergyMap       = new CsIImage( handler );
  E14CosmicAnalyzer* cosmicAnalyzer = new E14CosmicAnalyzer();
  TGraphErrors* grHeightTimeNoADJ = new TGraphErrors();
  TGraphErrors* grHeightTimePi0   = new TGraphErrors();
  TGraphErrors* grHeightTimeADJ   = new TGraphErrors();
  TGraphErrors* grTrack           = new TGraphErrors();
  grHeightTimeNoADJ->SetMarkerStyle( 6 );
  grHeightTimePi0->SetMarkerStyle( 6 );
  grHeightTimeADJ->SetMarkerStyle( 6 );



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Time Offset 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Double_t TimeOffsetTotal[2716] = {0};
  Double_t TimeOffsetTotalSigma[2716]= {0xFFFF};

  Double_t TimeOffset[2716]={500};
  Double_t TimeOffsetSigma[ 2716 ] = {0xFFFF};
  Double_t TimeOffsetCrystalPosition[2716] = {0};
  Double_t DeltaT[2716]={0};
  Double_t ResolutionT[2716]={0xFFFF};

  ////////////////////////////
  /* Time offset of Pi0Peak */ 
  ////////////////////////////
  std::ifstream ifs(Form("%s/Pi0Peak.dat",ROOTFILE_COSMIC.c_str()));
  if( !ifs.is_open()){ std::cout<< "Pi0Peak.dat is not existed" << std::endl; return -1;}
  Int_t tID;
  Double_t tOffset;
  Double_t tOffsetSigma;
  while( ifs >> tID >> tOffset >> tOffsetSigma ){
    TimeOffsetSigma[ tID ] = tOffsetSigma;
    if( tOffsetSigma > 0 ){
      TimeOffset[ tID ] = tOffset;      
    }else{
      TimeOffset[ tID ] = -384 ;
    }
  }

  ///////////////////////////////////////
  /* Time offset from crystal position */
  ///////////////////////////////////////
  for( int i = 0; i< 2716; i++){
    double x,y; 
    handler->GetMetricPosition( i, x, y );
    TimeOffsetCrystalPosition[i] = (TMath::Sqrt( 2624*2624 + x*x +y*y ) - 2624 )/ 299.7 ; // ns 
  }

  //////////////////////////////////////
  /*Time Offset from iteration*/
  //////////////////////////////////////

  if( IterationNumber >0 ){
    std::cout<< "IterationNumber is bigger than 0" << std::endl;
    std::string offsetFilename=Form("%s/CosmicOutTimeDeltaResolution_%d.dat",ROOTFILE_COSMIC.c_str(),IterationNumber-1);
    std::cout<< offsetFilename << std::endl;
    std::ifstream ifst(Form("%s/CosmicOutTimeDeltaResolution_%d.dat",ROOTFILE_COSMIC.c_str(),IterationNumber-1));
    if( ifst.is_open() ){
      std::cout << Form("CosmicOutTimeDeltaResolution_%d.dat is opened.",IterationNumber-1) << std::endl;
    }else{
      std::cout << "File isn't opened" << std::endl;
      return -1;
    }

    int tmpID;
    double tmpDelta;
    double tmpResolution;
    while( ifst >> tmpID >> tmpDelta >> tmpResolution ){
      DeltaT[tmpID] = tmpDelta;
      ResolutionT[tmpID] = tmpResolution; 
      std::cout<< tmpID << " : " << tmpDelta << " : " << tmpResolution << std::endl;
    }
  }

  ////////////////////////////////////
  /* Sum up All Offsets */ 
  ///////////////////////////////////

  for( int idIndex  = 0; idIndex  < 2716; idIndex++){
    TimeOffsetTotal[idIndex] += TimeOffset[idIndex];
    TimeOffsetTotal[idIndex] -= TimeOffsetCrystalPosition[idIndex];
    if( IterationNumber >0 ){
      if( ResolutionT[ idIndex ] > 0 && ResolutionT[ idIndex ] < 0xFFFF ){
	TimeOffsetTotal[idIndex] += DeltaT[idIndex];
	TimeOffsetTotalSigma[idIndex] = ResolutionT[ idIndex ];
      }
    }
  }
  for( int idIndex = 0; idIndex < 2716; idIndex++){
    std::cout<< idIndex << " : " << TimeOffset[idIndex] << " : " <<  TimeOffsetTotal[idIndex] << std::endl;
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  /*
  for( int i = 0; i< 2716; i++){
    std::cout<< i<< " : " << DeltaT[i] << " : " << ResolutionT[i] << std::endl;
  }
  */
  
  TFile* tfout = new TFile(Form("%s/CosmicOut_TimeCalibration_ADV_%d_%d.root",ROOTFILE_COSMIC.c_str(),fRunNumber,IterationNumber),"RECREATE");

  E14WavReader* reader = new E14WavReader(trin);
  Long_t entries =  reader->fChain->GetEntries();
  
  TTree* trout = new TTree("trOut","");
  
  const int nCSI = 2716;
  Int_t    RunNumber;
  Int_t    EventNumber;
  Double_t ScintiSignal;
  Double_t ScintiHHTime;
  Double_t ScintiTime;
  Int_t    nCSIDigi;
  Double_t CSIDigiE[nCSI];//nCSIDigi
  Double_t CSIDigiTime[nCSI];//nCSIDigi
  Double_t CSIDigiHHTime[nCSI];//nCSIDigi
  Int_t    CSIDigiID[nCSI];//nCSIDigi
  Double_t CSIDigiSignal[nCSI];//nCSIDigi
  Double_t FitP0[nCSI];//nCSIDigi
  Double_t FitP1[nCSI];//nCSIDigi
  Double_t FitChisq[nCSI];//nCSIDigi
  Double_t CSIDigiDeltaT0[nCSI];//nCSIDigi
  Double_t CSIDigiDeltaT1[nCSI];//nCSIDigi
  Int_t    CosmicTrigUp;
  Int_t    CosmicTrigDn;
  Double_t Roh;
  Double_t Theta;

  trout->Branch( "RunNumber"     , &RunNumber      , "RunNumber/I");
  trout->Branch( "EventNumber"   , &EventNumber    , "EventNumber/I");
  trout->Branch( "ScintiSignal"  , &ScintiSignal   , "ScintiSignal/D");
  trout->Branch( "ScintiHHTimne" , &ScintiHHTime   , "ScintiHHTime/D");
  trout->Branch( "ScintiTime"    , &ScintiTime     , "ScintiTime/D");
  trout->Branch( "nCSIDigi"      , &nCSIDigi       , "nCSIDigi/I" );
  trout->Branch( "CSIDigiE"      , CSIDigiE        , "CSIDigiE[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CSIDigiTime"   , CSIDigiTime     , "CSIDigiTime[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CSIDigiHHTime" , CSIDigiHHTime   , "CSIDigiHHTime[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CSIDigiID"     , CSIDigiID       , "CSIDigiID[nCSIDigi]/I");//nCSIDigi
  trout->Branch( "CSIDigiSignal" , CSIDigiSignal   , "CSIDigiSignal[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CSIDigiDeltaT0", CSIDigiDeltaT0  , "CSIDigiDeltaT0[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CSIDigiDeltaT1", CSIDigiDeltaT1  , "CSIDigiDeltaT1[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "FitP0"         , FitP0           , "FitP0[nCSIDigi]/D" );//nCSIDigi
  trout->Branch( "FitP1"         , FitP1           , "FitP1[nCSIDigi]/D" );//nCSIDigi
  trout->Branch( "FitChisq"      , FitChisq        , "FitChisq[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CosmicTrigUp"  , &CosmicTrigUp   , "CosmicTrigUp/I");
  trout->Branch( "CosmicTrigDn"  , &CosmicTrigDn   , "CosmicTrigDn/I");
  trout->Branch( "Roh"           , &Roh            , "Roh/D");
  trout->Branch( "Theta"         , &Theta          , "Theta/D");

  TCanvas* can = new TCanvas("can","can",1200,800);
  can->Divide( 3,2 );
  for( int i = 0; i< 6; i++){
    can->cd( i +1 ); 
    gPad->SetGridx();
    gPad->SetGridy();
  }


  TH1D* stepHist = new TH1D("hisStep","Step;Step;Survived Event",20,0,20);
  for( int ievent  = 0; ievent < entries ; ievent++ ){
    if( ievent % 1000 == 0){std::cout << ievent << "/" << entries << std::endl;}
    //if( ievent > 10000 ){ break; }
    reader->GetEntry( ievent  );
    TimeMap->Reset();
    EnergyMap->Reset();
    TriggerMap->Reset();

    grHeightTimeNoADJ->Set(0);
    grHeightTimePi0->Set(0);
    grHeightTimeADJ->Set(0);
    grTrack->Set(0);

    cosmicAnalyzer->Reset();

    ScintiTime = -500;
    ScintiHHTime = -500.;
    nCSIDigi = 0;
    CosmicTrigUp = -1;
    CosmicTrigDn = -1;

    for( Int_t iCSI = 0; iCSI < nCSI; iCSI++ ){
      CSIDigiE[iCSI]      = 0;
      CSIDigiTime[iCSI]   = 0;
      CSIDigiHHTime[iCSI] = 0;
      CSIDigiID[iCSI]     = -1;
    }
    for( int i = 0; i< 2; i++){
      FitP0[i] = 0xFFFF;
      FitP1[i] = 0xFFFF;
      FitChisq[i] = 0xFFFF;
    }    

    stepHist->Fill(0);   
    //std::cout<< "Scinti" << std::endl;
    RunNumber = reader->RunNo;
    EventNumber = reader->EventNo;
    if( reader->CosmicTrigFlagUp == 0 && 
	reader->CosmicTrigFlagDn == 0 ){
      continue;
    }
    if( reader->TrigFlag != 2 ){
      continue;
    }
    int nTrigUp = 0;
    int nTrigDn = 0;

    for( int i = 0; i< 5; i++){
      Int_t TrigFlag = 1 << i ;
      if( (TrigFlag & reader->CosmicTrigFlagUp ) != 0){
	nTrigUp++;
	TriggerMap->Fill( i, 4 );
	CosmicTrigUp = i;
      }
      if( (TrigFlag & reader->CosmicTrigFlagDn ) != 0){
	nTrigDn++;
	TriggerMap->Fill( i, 0 );
	CosmicTrigDn = i;
      }
    }
    if( nTrigUp != 1  || nTrigDn != 1 ){      
      continue;
    }
	
      

    stepHist->Fill(1);
    //std::cout << "CsI" << std::endl;

    Double_t CosmicEventEnergyThreshold= 3.;
    Double_t CosmicEnergyThreshold = 5.;

    //std::cout<< "Channel Loop" << std::endl;
    for( int ich  = 0; ich < reader->CsiNumber; ich++){
      if( reader->CsiSignal[ich] < 5 ){ 
	continue;
      }

      Double_t x,y;
      handler->GetMetricPosition( reader->CsiID[ich] , x, y );
      
      //std::cout << "Converter " << std::endl;
      double Energy  = Converter->ConvertToEnergy( reader->CsiID[ich], reader->CsiSignal[ich] );         

      if( Energy > CosmicEventEnergyThreshold ){//&& TimeOffsetSigma[reader->CsiID[ich]] > 0 ){
	//std::cout << ich << " : " <<   Energy << std::endl; 	
	CSIDigiID[nCSIDigi]     = reader->CsiID[ich];
	CSIDigiTime[nCSIDigi]   = reader->CsiTime[ich];
	CSIDigiHHTime[nCSIDigi] = reader->CsiHHTime[ich];
	CSIDigiSignal[nCSIDigi] = reader->CsiSignal[ich];
	CSIDigiE[nCSIDigi]      = Converter->ConvertToEnergy( reader->CsiID[ich] ,reader->CsiSignal[ich] );
	nCSIDigi++;
	grHeightTimeNoADJ->SetPoint( grHeightTimeNoADJ->GetN(),
				     y,
				     reader->CsiTime[ ich ] );
	grHeightTimePi0  ->SetPoint( grHeightTimePi0->GetN(),
				     y,
				     reader->CsiTime[ ich ] - TimeOffsetTotal[ reader->CsiID[ ich ] ]);
	grHeightTimeADJ  ->SetPoint( grHeightTimeADJ->GetN(), 
				     y,
				     reader->CsiTime[ ich ] - TimeOffsetTotal[ reader->CsiID[ ich ] ]);
	
	if( IterationNumber == 0 ){
	  grHeightTimeNoADJ->SetPointError( grHeightTimeNoADJ->GetN()-1,
					    0,
					    1/TMath::Sqrt(Energy/14));
	  grHeightTimePi0  ->SetPointError( grHeightTimePi0->GetN()-1,
					    0,
					    1/TMath::Sqrt(Energy/14));
	  
	  grHeightTimeADJ  ->SetPointError( grHeightTimeADJ->GetN()-1,
					    0,
					    1/TMath::Sqrt(Energy/14));
	}else{ 
	  grHeightTimeNoADJ->SetPointError( grHeightTimeNoADJ->GetN()-1,
					    0,
					    ResolutionT[reader->CsiID[ich]]/TMath::Sqrt(Energy/14));
	  grHeightTimePi0  ->SetPointError( grHeightTimePi0->GetN()-1,
					    0,
					    ResolutionT[reader->CsiID[ich]]/TMath::Sqrt(Energy/14));  
	  grHeightTimeADJ  ->SetPointError( grHeightTimeADJ->GetN()-1,
					    0,
					    ResolutionT[reader->CsiID[ich]]/TMath::Sqrt(Energy/14));
	}
	grTrack->SetPoint( grTrack->GetN(),x,y);
	
	TimeMap->Fill( reader->CsiID[ ich ]   , reader->CsiTime[ ich ] - TimeOffset[ reader->CsiID[ ich ] ] );
	EnergyMap->Fill( reader->CsiID[ ich ] , Converter->ConvertToEnergy( reader->CsiID[ ich ] , reader->CsiSignal[ ich ] ));
	
      }
    }
    
    if( nCSIDigi < 8 ){ 
      //std::cout<< "nCSIDigi:" << nCSIDigi <<  std::endl;
      continue; 
    }
    
    cosmicAnalyzer->GetResult( grTrack, Roh, Theta );    
    grHeightTimePi0->Fit( "pol1", "Q","", -1000, 1000);
    TF1* funcTime0 = grHeightTimePi0->GetFunction("pol1");    
    grHeightTimeADJ->Fit( "pol1", "Q","", -1000,1000);
    TF1* funcTime1 = grHeightTimeADJ->GetFunction("pol1");
    FitP0[0] = funcTime0->GetParameter(0);
    FitP0[1] = funcTime1->GetParameter(0);
    FitP1[0] = funcTime0->GetParameter(1);
    FitP1[1] = funcTime1->GetParameter(1);
    FitChisq[0] = funcTime0->GetChisquare() / funcTime0->GetNDF();
    FitChisq[1] = funcTime1->GetChisquare() / funcTime1->GetNDF();
    for( int idigi = 0; idigi < nCSIDigi; idigi++){
      Double_t x,y;
      handler->GetMetricPosition( CSIDigiID[ idigi ] , x, y);
      CSIDigiDeltaT0[idigi] = CSIDigiTime[idigi] - TimeOffsetTotal[CSIDigiID[idigi]] - funcTime0->Eval(y);
      CSIDigiDeltaT1[idigi] = CSIDigiTime[idigi] - TimeOffsetTotal[CSIDigiID[idigi]] - funcTime1->Eval(y);
    }
    grHeightTimePi0->GetListOfFunctions()->Delete();
    grHeightTimeADJ->GetListOfFunctions()->Delete();

    //std::cout<<FitP0[0] << std::endl;
    /*
    can->cd(1);
    TimeMap->Draw("colz");
    can->cd(2);
    EnergyMap->Draw("colz");
    can->cd(3);
    TriggerMap->Draw("col");
    can->cd(4);
    grHeightTimeNoADJ->Draw("AP");
    grHeightTimeNoADJ->Fit("pol1","","",-600,600);
    can->cd(5);
    grHeightTimePi0->Draw("AP");
    grHeightTimePi0->Fit("pol1","","",-600,600);
    can->cd(6);
    grHeightTimeADJ->Draw("AP");
    grHeightTimeADJ->Fit("pol1","","",-600,600);
    can->Update();
    can->Modified();

    getchar();
    */

    trout->Fill();
  }

  trout->Write();
  stepHist->Write();
  tfout->Close();
  //app->Run();
  return 0;
}
