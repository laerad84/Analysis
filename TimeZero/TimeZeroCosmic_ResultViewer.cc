#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include "TFile.h"
#include "TGraph.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TApplication.h"


int
main( int argc, char** argv ){

  gStyle->SetOptFit(11111111);
  gStyle->SetOptStat(11111111);
  
  std::string RunList = argv[1];
  int IterationNumber = atoi( argv[2] );


  //TApplication* app = new TApplication("app",&argc, argv);
  TChain* trin = new TChain("trOut");

  std::string ROOTFILE_WAV=std::getenv("ROOTFILE_WAV");
  std::string DIRNAME = "TimeZero";
  std::string ANALYSISDIR = std::getenv("ANALYSISDIR");
  std::string ROOTFILE_COSMIC=std::getenv("ROOTFILE_COSMIC");
  std::ifstream ifsRunList(RunList.c_str());
  int tmpRunNumber;
  while( ifsRunList >> tmpRunNumber ){
    //trin->Add(Form("%s/CosmicOut_TimeCalibration_FNL_FIXV_%d_%d.root",ROOTFILE_COSMIC.c_str(),tmpRunNumber,IterationNumber));
    trin->Add(Form("%s/CosmicOut_TimeCalibration_FNL_FIXV_%d_%d.root",ROOTFILE_COSMIC.c_str(),tmpRunNumber,IterationNumber));
  }  

  int    ID[2716];
  double Delta[2716];
  double Resolution[2716];
  for( int i = 0; i< 2716; i++){
    Resolution[i] = 0xFFFF;
    Delta[i]      = 0;
  }

  const int nCSI = 2716;
  Int_t    RunNumber;
  Int_t    EventNumber;
  Double_t ScintiSignal = 0;
  Double_t ScintiHHTime = -500.;
  Double_t ScintiTime   =-500.;
  Int_t    nCSIDigi     = 0;
  Double_t CSIDigiE[nCSI];//nCSIDigi
  Double_t CSIDigiTime[nCSI];//nCSIDigi
  Double_t CSIDigiHHTime[nCSI];//nCSIDigi
  Int_t    CSIDigiID[nCSI];//nCSIDigi
  Double_t CSIDigiSignal[nCSI];//nCSIDigi
  Double_t FitP0[nCSI];
  Double_t FitP1[nCSI];
  Double_t FitChisq[nCSI];
  Double_t CSIDigiDeltaT[nCSI];//nCSIDigi
  Int_t    CosmicTrigUp;
  Int_t    CosmicTrigDn;
  Double_t Roh;
  Double_t Theta;
  Double_t DistFromLine[nCSI];//nCSIDigi
  Double_t HeightFromLine[nCSI];//nCSIDigi
  


  trin->SetBranchStatus( "*", 0 );
  trin->SetBranchStatus( "RunNumber"      );
  trin->SetBranchStatus( "EventNumber"    );
  trin->SetBranchStatus( "ScintiSignal"   );
  trin->SetBranchStatus( "ScintiHHTimne"  );
  trin->SetBranchStatus( "ScintiTime"     );
  trin->SetBranchStatus( "nCSIDigi"       );
  trin->SetBranchStatus( "CSIDigiE"       );
  trin->SetBranchStatus( "CSIDigiTime"    );
  trin->SetBranchStatus( "CSIDigiHHTime"  );
  trin->SetBranchStatus( "CSIDigiID"      );
  trin->SetBranchStatus( "CSIDigiSignal"  );
  trin->SetBranchStatus( "CSIDigiDeltaT"  );
  trin->SetBranchStatus( "FitP0"          );
  trin->SetBranchStatus( "FitP1"          );
  trin->SetBranchStatus( "FitChisq"       );
  trin->SetBranchStatus( "CosmicTrigUp"   );
  trin->SetBranchStatus( "CosmicTrigDn"   );
  trin->SetBranchStatus( "Roh"            );
  trin->SetBranchStatus( "Theta"          );
  trin->SetBranchStatus( "DistFromLine"   );
  trin->SetBranchStatus( "HeightFromLine" );
  

  trin->SetBranchAddress( "RunNumber"      , &RunNumber      );
  trin->SetBranchAddress( "EventNumber"    , &EventNumber    );
  trin->SetBranchAddress( "ScintiSignal"   , &ScintiSignal   );
  trin->SetBranchAddress( "ScintiHHTimne"  , &ScintiHHTime   );
  trin->SetBranchAddress( "ScintiTime"     , &ScintiTime     );
  trin->SetBranchAddress( "nCSIDigi"       , &nCSIDigi       );
  trin->SetBranchAddress( "CSIDigiE"       , CSIDigiE        );
  trin->SetBranchAddress( "CSIDigiTime"    , CSIDigiTime     );
  trin->SetBranchAddress( "CSIDigiHHTime"  , CSIDigiHHTime   );
  trin->SetBranchAddress( "CSIDigiID"      , CSIDigiID       );
  trin->SetBranchAddress( "CSIDigiSignal"  , CSIDigiSignal   );
  trin->SetBranchAddress( "CSIDigiDeltaT"  , CSIDigiDeltaT   );
  trin->SetBranchAddress( "FitP0"          , FitP0           );
  trin->SetBranchAddress( "FitP1"          , FitP1           );
  trin->SetBranchAddress( "FitChisq"       , FitChisq        );
  trin->SetBranchAddress( "CosmicTrigUp"   , &CosmicTrigUp   );
  trin->SetBranchAddress( "CosmicTrigDn"   , &CosmicTrigDn   );
  trin->SetBranchAddress( "Roh"            , &Roh            );
  trin->SetBranchAddress( "Theta"          , &Theta          );
  trin->SetBranchAddress( "DistFromLine"   , DistFromLine    );
  trin->SetBranchAddress( "HeightFromLine" , HeightFromLine  );

  TFile* tfout = new TFile(Form("%s/CosmicOuthist_%d.root",ROOTFILE_COSMIC.c_str(),IterationNumber), "recreate");

  TH1D* hisResult = new TH1D("hisResult","hisResult;DeltaTime[ns];nEnetries[0.1ns]",200,-10,10);
  TH2D* hisDeltaChannel = new TH2D("hisDeltaChannel","hisDeltaChannel",2716,0,2716,400,-40,40);
  TH1D* hisDelta[2716];
  TH2D* hisDeltaAll = new TH2D("hisDeltaAll","hisDeltaAll;ChannelNo;DeltaT",2716,0,2716,400,-40,40);
  TH1D* hisDeltaNoCut[2716];

  TGraphErrors* grDelta = new TGraphErrors();
  TGraphErrors* grRES   = new TGraphErrors();   

  TPostScript* ps  = new TPostScript(Form("%s/CosmicOuthist_%d.ps",ROOTFILE_COSMIC.c_str(),IterationNumber),111);
  TCanvas *can = new TCanvas("can","",1000*TMath::Sqrt(0.5),1000);
  can->Divide(2,3);
  for( int i  = 0; i< 6; i++){
    can->cd(i+1);
    gPad->SetGridx();
    gPad->SetGridy();
  }

  //ps->NewPage();

  for( int i = 0; i< 2716; i++){
    hisDelta[i] = new TH1D(Form("hisDelta%d",i ),Form("hisDelta%d",i),400,-40,40);
    hisDeltaNoCut[i] = new TH1D(Form("hisDeltaNoCut%d",i),Form("hisDeltaNoCut%d",i),400,-40,40);
  }

  for( int ievent = 0; ievent < trin->GetEntries(); ievent++){

    trin->GetEntry(ievent);
    if( FitChisq[1] > 3 ){ continue; }
    for( int idigi = 0; idigi < nCSIDigi ; idigi++){
      hisDeltaNoCut[ CSIDigiID[ idigi ] ]->Fill( CSIDigiDeltaT[ idigi ] );
      hisDeltaAll->Fill( CSIDigiID[ idigi  ] , CSIDigiDeltaT[ idigi ] );
    
      if( CSIDigiID[idigi] < 2240 ){ // Case of Small Crystal 
	if( CSIDigiE[idigi] > 11 && CSIDigiE[idigi] < 17 ){
	  hisDelta[ CSIDigiID[ idigi ] ]->Fill( CSIDigiDeltaT[ idigi ] );
	  hisDeltaChannel->Fill( CSIDigiID[ idigi  ] , CSIDigiDeltaT[ idigi ] );
	}
      }else{ // Case of Large Crystal ... 
	if( CSIDigiE[idigi] > 22 && CSIDigiE[idigi] < 34 ){
	  hisDelta[ CSIDigiID[ idigi ] ]->Fill( CSIDigiDeltaT[ idigi ] );
	  hisDeltaChannel->Fill( CSIDigiID[ idigi  ] , CSIDigiDeltaT[ idigi ] );
	}
      }
      //std::cout  << CSIDigiID[ idigi ] << std::endl;
    }
  }

  Int_t iCan = 1; 
  for( int i = 0; i< 2716; i++){
    //std::cout << hisDelta[i]->GetEntries() << std::endl;
    if( hisDelta[i]->GetEntries() > 10){
      int rst = hisDelta[i]->Fit("gaus","Q","",
				 hisDelta[i]->GetBinCenter( hisDelta[i]->GetMaximumBin() ) - 6, 
				 hisDelta[i]->GetBinCenter( hisDelta[i]->GetMaximumBin() ) + 6);
      TF1* func = NULL;
      func = hisDelta[i]->GetFunction("gaus");
      if( func != NULL ){
	double mean = func->GetParameter(1);
	double sigma = func->GetParameter(2);
	int rst_1 = hisDelta[i]->Fit("gaus","Q","",mean-3*sigma, mean+3*sigma);
	if( rst_1 != 0 ){ std::cout<< i << std::endl; }
	func = hisDelta[i]->GetFunction("gaus");
	Delta[i] = func->GetParameter(1);
	Resolution[i] = func->GetParameter(2);
	if( TMath::Abs(Delta[i] - mean ) > 3*sigma || rst_1 != 0 ){ 
	  rst_1 = hisDelta[i]->Fit("gaus","Q","",mean-4*sigma, mean+4*sigma);
	  if( rst_1 !=0 ){ std::cout<< i << std::endl; }
	  func = hisDelta[i]->GetFunction("gaus");
	  Delta[i] = func->GetParameter(1);
	  Resolution[i] = func->GetParameter(2);
	}
	grDelta->SetPoint( grDelta->GetN(), i, func->GetParameter(1));
	grDelta->SetPointError( grDelta->GetN()-1, 0, func->GetParError(2));
	grRES->SetPoint( grRES->GetN() , i , func->GetParameter(2));
	hisResult->Fill(Delta[i]);
      }
    }
 
    can->cd(iCan);
    hisDelta[i]->DrawCopy(); 
    if( iCan == 6 ){
      can->Update();
      can->Modified();
      ps->NewPage();
    }
    if( iCan == 6 ){ iCan = 0;}
    iCan++;
    hisDelta[i]     ->Write(); 
    hisDeltaNoCut[i]->Write();
  }

  int tmpID;
  double tmpDelta;
  double tmpResolution;
  if( IterationNumber > 0){
    std::ifstream ifs(Form("%s/CosmicOutTimeDeltaResolution_%d.dat",ROOTFILE_COSMIC.c_str(),IterationNumber-1));
    if( !ifs.is_open()){ std::cout<< "Error:" << Form("CosmicOutTimeDeltaResolution_%d.dat",IterationNumber-1) << " is not exist"<< std::endl; return -1;}
    while ( ifs >> tmpID >> tmpDelta >> tmpResolution ){
      Delta[ tmpID ] += tmpDelta;
      if( Resolution[tmpID] == 0xFFFF ){
	Resolution[ tmpID ] = tmpResolution;
      }
    }
    ifs.close();
  }

  std::ofstream ofs(Form("%s/CosmicOutTimeDeltaResolution_%d.dat",ROOTFILE_COSMIC.c_str(),IterationNumber));
  /*
  for( int i = 0; i< grRES->GetN(); i++){
    hisResult->Fill(grDelta->GetY()[i]);
    Delta[(int)(grDelta->GetX()[i])] += grDelta->GetY()[i];
    Resolution[(int)(grRES->GetX()[i])] = grRES->GetY()[i] ;
  }
  */
  for( int i = 0; i< 2716; i++){
    ofs << i             << "\t" 
	<< Delta[i]      << "\t"
	<< Resolution[i] << "\n";
  }
  
  grDelta->SetNameTitle("grDelta","grDelta");
  grRES->SetNameTitle("grRES","grRES");
  grDelta->Write();
  grRES->Write();
  hisDeltaAll->Write();
  hisDeltaChannel->Write();
  hisResult->Write();
  tfout->Close();
  ps->Close();
  ofs.close();

  //app->Run();
}
