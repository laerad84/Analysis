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
#include "TVector.h"
#include "TCanvas.h"
#include "TApplication.h"

#include "CsIPoly.h"
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// TimeZeroCosmic [RunNumberList] [IterationNumber]
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
double linearFunc( double *x , double *par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  return p0+x0*p1;
}

Double_t DFromLine( double x, double y, double roh, double theta){
  TVector2 vec(x,y);
  TVector2 vecp = vec.Rotate( -1*theta );
  double DistFromLine = TMath::Abs( vecp.X() - roh );  
  return DistFromLine;
}
Double_t HFromLine( double x, double y, double roh, double theta ){
  TVector2 vec(x,y);
  TVector2 vecp = vec.Rotate( -1*theta );
  double HeightFromLine = vecp.Y();
  return HeightFromLine;
}
void DHFromLine( double x, double y, double roh, double theta, double& H, double& D){
  TVector2 vec(x,y);
  TVector2 vecp = vec.Rotate( -1*theta );
  H = vecp.Y();
  D = TMath::Abs( vecp.X() - roh );  
}

#define LINE() std::cout<< __PRETTY_FUNCTION__ << std::endl;std::cout<<__LINE__ << std::endl;
#define STATUS( msg ) LINE(); std::cout<< msg << std::endl;

int
main( int argc ,char ** argv ){
  gStyle->SetOptFit(111111111); 
 std::string WAVFILE   = std::getenv("ROOTFILE_WAV");
  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  //std::string ANALYSISDIR = std::getenv("ANALYSISDIR");
  std::string DIRNAME   = "TimeZero";
  std::string ROOTFILE_COSMIC= std::getenv("ROOTFILE_COSMIC");
  EnergyConverter *Converter = new EnergyConverter();  
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",
					  ANALIBDIR.c_str()));
  STATUS("Check_TFile");
  TChain* trin = new TChain("Tree");
  int fRunNumber = atoi(argv[1]);
  int IterationNumber = atoi( argv[2]);
  std::cout << "////////////////////////////////////////////////////\n";
  std::cout << "RunNumber       : " << fRunNumber << "\n";
  std::cout << "IterationNumber : " << IterationNumber << "\n";
  std::cout << "////////////////////////////////////////////////////\n";    
  trin->Add(Form("%s/run_wav_%d.root",WAVFILE.c_str(),fRunNumber));  
  
  IDHandler* handler       = new IDHandler();
  TH2D* TriggerMap         = new TH2D("Trigger","Trigger",5,0,5,5,0,5);
  CsIPoly* TimeMap         = new CsIPoly("TimeMap","TimeMap");
  CsIPoly* EnergyMap       = new CsIPoly("EnergyMap","EnergyMap");
  CsIPoly* hisTrack        = new CsIPoly("hisTrack","hisTrack");
  CsIPoly* hisTrackD       = new CsIPoly("hisTrackD","hisTrackD");
  CsIPoly* hisTrackL       = new CsIPoly("hisTrackL","hisTrackH");

  E14CosmicAnalyzer* cosmicAnalyzer = new E14CosmicAnalyzer();
  TGraphErrors* grHeightTimeNoADJ = new TGraphErrors();
  TGraphErrors* grHeightTimePi0   = new TGraphErrors();
  TGraphErrors* grHeightTimeADJ   = new TGraphErrors();
  TGraphErrors* grTrack           = new TGraphErrors();
  
  grTrack->SetMarkerStyle(21);
  grHeightTimeNoADJ->SetMarkerStyle(21);
  grHeightTimePi0->SetMarkerStyle(21);
  grHeightTimeADJ->SetMarkerStyle(21);
  std::cout << __PRETTY_FUNCTION__ << std::endl;


  TApplication* app = new TApplication("app",&argc, argv);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Time Offset 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Crystal Position Dep. is vanished 20130305 of clearance
  /// 

  Double_t TimeOffsetTotal[2716] = {0};
  Double_t TimeOffsetTotalSigma[2716]= {0xFFFF};

  Double_t TimeOffset[2716]={500};
  Double_t TimeOffsetSigma[ 2716 ] = {0xFFFF};
  Double_t TimeOffsetCrystalPosition[2716] = {0};
  Double_t DeltaT[2716]={0};
  Double_t ResolutionT[2716]={0xFFFF};

  ////////////////////////////
  // Time offset of Pi0Peak // 
  ////////////////////////////
  /*
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
      TimeOffset[ tID ] = 0;
    }
  }
  */
  
  for( int i = 0; i< 2716; i++){
    if( i < 2240 ){
      TimeOffset[i] = 0.;
      TimeOffsetSigma[i] = 1;
    }else{
      TimeOffset[i] = 16;
      TimeOffsetSigma[i] = 1;
    }
  }

  ///////////////////////////////////////
  // Time offset from crystal position //
  ///////////////////////////////////////
  for( int i = 0; i< 2716; i++){
    double x,y; 
    handler->GetMetricPosition( i, x, y );
    TimeOffsetCrystalPosition[i] = (TMath::Sqrt( 2624*2624 + x*x +y*y ) - 2624 )/ 299.7 ; // ns 
  }

  //////////////////////////////////////
  //Time Offset from iteration//
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
  // Sum up All Offsets // 
  ///////////////////////////////////

  for( int idIndex  = 0; idIndex  < 2716; idIndex++){
    TimeOffsetTotal[idIndex] += TimeOffset[idIndex];
    //TimeOffsetTotal[idIndex] -= TimeOffsetCrystalPosition[idIndex];
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
  
  //  TFile* tfout = new TFile(Form("%s/CosmicOut_TimeCalibration_FNL_FIXV_%d_%d.root",ROOTFILE_COSMIC.c_str(),fRunNumber,IterationNumber),"RECREATE");

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
  Double_t CSIDigiDeltaT[nCSI];//nCSIDigi
  Double_t DistFromLine[nCSI];//nCSIDigi;
  Double_t HeightFromLine[nCSI];//CSIDigi;

  Int_t    CosmicTrigUp;
  Int_t    CosmicTrigDn;
  Double_t Roh;
  Double_t Theta;
  
  {
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
  trout->Branch( "CSIDigiDeltaT" , CSIDigiDeltaT   , "CSIDigiDeltaT[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "FitP0"         , FitP0           , "FitP0[nCSIDigi]/D" );//nCSIDigi
  trout->Branch( "FitP1"         , FitP1           , "FitP1[nCSIDigi]/D" );//nCSIDigi
  trout->Branch( "FitChisq"      , FitChisq        , "FitChisq[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "DistFromLine"  , DistFromLine    , "DistFromLine[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "HeightFromLine", HeightFromLine  , "HeightFromLine[nCSIDigi]/D");//nCSIDigi
  trout->Branch( "CosmicTrigUp"  , &CosmicTrigUp   , "CosmicTrigUp/I");
  trout->Branch( "CosmicTrigDn"  , &CosmicTrigDn   , "CosmicTrigDn/I");
  trout->Branch( "Roh"           , &Roh            , "Roh/D");
  trout->Branch( "Theta"         , &Theta          , "Theta/D");
  }
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
    hisTrack->Reset();
    cosmicAnalyzer->Reset();

    ScintiTime   = -500;
    ScintiHHTime = -500.;
    nCSIDigi     = 0;
    CosmicTrigUp = -1;
    CosmicTrigDn = -1;

    for( Int_t iCSI = 0; iCSI < nCSI; iCSI++ ){
      CSIDigiE[iCSI]      = 0;
      CSIDigiTime[iCSI]   = 0;
      CSIDigiHHTime[iCSI] = 0;
      CSIDigiID[iCSI]     = -1;
    }
    for( int i = 0; i< 2; i++){
      FitP0[i]    = 0xFFFF;
      FitP1[i]    = 0xFFFF;
      FitChisq[i] = 0xFFFF;
    }    

    stepHist->Fill(0);   
    //std::cout<< "Scinti" << std::endl;
    RunNumber   = reader->RunNo;
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

    bool SameScintiTriggered = false;
    for( int i = 0; i< 5; i++){
      Int_t TrigFlag = 1 << i; 
      if(( TrigFlag & reader->CosmicTrigFlagUp ) !=0 &&
	 ( TrigFlag & reader->CosmicTrigFlagDn ) !=0 ){
	SameScintiTriggered = true; 
      }						     
    }

    /*
    if( nTrigUp >= 1  || nTrigDn >= 1 ){      
      continue;
    }
    */

    if( !SameScintiTriggered ){ continue; }

    stepHist->Fill(1);
    //std::cout << "CsI" << std::endl;
    ////////////////////////////////////////////////////////////////////////////////////////
    //std::cout<< "Channel Loop" << std::endl;
    // For Analyzing Track
    // Set Data to track which signal height is higher than 5 cnt and Energy Bigger than 3.
    ////////////////////////////////////////////////////////////////////////////////////////
    const Double_t CosmicEventEnergyThreshold = 3.;
    const Double_t CosmicEventHeightThreshold = 5.;
    ////////////////////////////////////////////////////////////////////////////////////////
    for( int ich = 0; ich < reader->CsiNumber; ich++){
      if( reader->CsiSignal[ich] < CosmicEventHeightThreshold ){ continue; }
      Double_t x,y;
      double Energy = Converter->ConvertToEnergy( reader->CsiID[ich], reader->CsiSignal[ich] );
      if( Energy < CosmicEventEnergyThreshold ) { continue; }
      handler->GetMetricPosition( reader->CsiID[ich] , x, y );
      grTrack->SetPoint( grTrack->GetN(),x,y);

    }
    cosmicAnalyzer->GetResult( grTrack, Roh, Theta );    
    //std::cout<< Roh << "\t" << Theta << std::endl;
    Double_t c  = TMath::Cos(Theta/180*TMath::Pi());
    Double_t s  = TMath::Sin(Theta/180*TMath::Pi());    
    Double_t y0 = -1000;
    Double_t x0 = ( Roh*c*c - s*(y0-Roh*s))/c;
    Double_t y1 = 1000;
    Double_t x1 = ( Roh*c*c - s*(y1-Roh*s))/c;

    ////////////////////////////////////////////////////////////////////////////////////////
    // For Analyzing Time-Track relation 
    ////////////////////////////////////////////////////////////////////////////////////////
    for( int ich  = 0; ich < reader->CsiNumber; ich++){
      if( reader->CsiSignal[ich] < 5 ){ 
	continue;
      }
      Double_t x,y;
      Double_t D,H;
      handler->GetMetricPosition( reader->CsiID[ich] , x, y );      
      DHFromLine( x,y,Roh,Theta/180*TMath::Pi(),H,D);

      double Energy  = Converter->ConvertToEnergy( reader->CsiID[ich], reader->CsiSignal[ich] );         
      //Distance cut //////////////////
      if( reader->CsiID[ich] < 2240 ){
	if( TMath::Abs(D) > 50 ){ continue; }
      }else{
	if( TMath::Abs(D) > 100 ){ continue; }
      }
      //std::cout<<reader->CsiID[ich] << "\t" << x << "\t" << y << "\t" << D << "\t"<< H << std::endl;

      /////////////////////////////////

      if( Energy > CosmicEventEnergyThreshold ){//&& TimeOffsetSigma[reader->CsiID[ich]] > 0 ){

	CSIDigiID[nCSIDigi]     = reader->CsiID[ich];
	CSIDigiTime[nCSIDigi]   = reader->CsiTime[ich];
	CSIDigiHHTime[nCSIDigi] = reader->CsiHHTime[ich];
	CSIDigiSignal[nCSIDigi] = reader->CsiSignal[ich];
	CSIDigiE[nCSIDigi]      = Converter->ConvertToEnergy( reader->CsiID[ich] ,reader->CsiSignal[ich] );
	DistFromLine[nCSIDigi]  = D;
	HeightFromLine[nCSIDigi]= H;
	nCSIDigi++;

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set Time Data 
	// grheightTimeNoADJ : No correction
	// grHeightTimeADJ   : Correct time with Full Adjust [ pi0 + crystal position ]
	////////////////////////////////////////////////////////////////////////////////////////////////////
	hisTrack->Fill(reader->CsiID[ich], CSIDigiE[nCSIDigi-1]);
	grHeightTimeNoADJ->SetPoint( grHeightTimeNoADJ->GetN(),
				     H,
				     reader->CsiTime[ ich ] );
	grHeightTimeADJ  ->SetPoint( grHeightTimeADJ->GetN(), 
				     H,
				     reader->CsiTime[ ich ] - TimeOffsetTotal[ reader->CsiID[ ich ] ]);

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set Error 
	// If Iteration Number == 0 , Assume time resolution of each channel as 1/sqrt(e).
	// If Iteration Number >  1 , Time resolutions of each channels are set to [acquired value]/sqrt(e);
	////////////////////////////////////////////////////////////////////////////////////////////////////
	if( IterationNumber == 0 ){
	  grHeightTimeNoADJ->SetPointError( grHeightTimeNoADJ->GetN()-1,
					    0,
					    1/TMath::Sqrt(Energy/14));
	  grHeightTimeADJ  ->SetPointError( grHeightTimeADJ->GetN()-1,
					    0,
					    1/TMath::Sqrt(Energy/14));
	}else{ 
	  grHeightTimeNoADJ->SetPointError( grHeightTimeNoADJ->GetN()-1,
					    0,
					    ResolutionT[reader->CsiID[ich]]/TMath::Sqrt(Energy/14));
	  grHeightTimeADJ  ->SetPointError( grHeightTimeADJ->GetN()-1,
					    0,
					    ResolutionT[reader->CsiID[ich]]/TMath::Sqrt(Energy/14));
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////	
	TimeMap  ->Fill( reader->CsiID[ich],
			 reader->CsiTime[ich] - TimeOffset[ reader->CsiID[ich]] );
	EnergyMap->Fill( reader->CsiID[ich],
			 Converter->ConvertToEnergy( reader->CsiID[ich] , reader->CsiSignal[ich] ));	
      }
    }
    
    if( nCSIDigi < 8 ){ 
      //std::cout<< "nCSIDigi:" << nCSIDigi <<  std::endl;
      continue; 
    }

    ///////////////////////////////////////////////////////////////////
    // Fitting with linear data 
    ///////////////////////////////////////////////////////////////////    
    for( Int_t ipoint = 0; ipoint < grHeightTimeADJ->GetN(); ipoint++){
      TF1* funcTime = new TF1("linearFunc",linearFunc,-1000,1000,2);
      funcTime->SetParameter(1,-1./300);
      funcTime->SetParLimits(1,-1./300,-1./300);
      Double_t tmpErr = grHeightTimeADJ->GetEY()[ipoint];
      // For get off effect of point //
      grHeightTimeADJ->SetPointError(ipoint, 0, 0xFFFF);
      grHeightTimeADJ->Fit(funcTime,"Q","",-1000,1000);

      FitP0[ipoint]    = funcTime->GetParameter(0);
      FitP1[ipoint]    = funcTime->GetParameter(1);
      FitChisq[ipoint] = funcTime->GetChisquare() / funcTime->GetNDF();     
      // restore error point // 
      grHeightTimeADJ->SetPointError( ipoint, 0, tmpErr);
    }

    TLine* line  = new TLine(x0,y0,x1,y1);
    can->cd(1);
    TimeMap->Draw("col");
    can->cd(2);
    EnergyMap->Draw("col");
    can->cd(3);    
    can->cd(4);
    
    can->cd(5);
    gPad->SetLogz();
    hisTrack->Draw("col");

    grTrack->Draw("L");
    line->SetLineWidth(2);
    line->SetLineColor(5);
    line->Draw("same");

    can->cd(6);
    grHeightTimeADJ->Draw("AP");
    grHeightTimeADJ->GetYaxis()->SetRangeUser(0,200);
    can->Update();
    can->Modified();
    can->Update();
    getchar();
    line->Delete();


    ///////////////////////////////////////////////////////////////////
    // Get Time Offset
    ///////////////////////////////////////////////////////////////////
    
    for( int idigi = 0; idigi < nCSIDigi; idigi++){
      Double_t x,y;
      Double_t D,H;
      handler->GetMetricPosition( CSIDigiID[ idigi ],x,y);      
      DHFromLine(x,y,Roh,Theta*TMath::Pi()/180,H,D);      
      Double_t FuncFitTime = FitP0[idigi] + FitP1[idigi]*HeightFromLine[idigi];
      CSIDigiDeltaT[idigi] = CSIDigiTime[idigi] - TimeOffsetTotal[CSIDigiID[idigi]] - FuncFitTime; 

    }
    

    grHeightTimePi0->GetListOfFunctions()->Delete();
    grHeightTimeADJ->GetListOfFunctions()->Delete();
    trout->Fill();
  }
  //trout->Write();
  //stepHist->Write();
  //tfout->Close();
  app->Run();
  return 0;
}
