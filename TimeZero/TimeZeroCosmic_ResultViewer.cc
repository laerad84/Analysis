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

int
main( int argc, char** argv ){

  gStyle->SetOptFit(11111111);
  gStyle->SetOptStat(11111111);
  
  std::string RunList = argv[1];
  int IterationNumber = atoi( argv[2] );

  TChain* trin = new TChain("trOut");

  std::string ROOTFILE_WAV=std::getenv("ROOTFILE_WAV");
  std::ifstream ifsRunList(RunList.c_str());
  int tmpRunNumber;
  while( ifsRunList >> tmpRunNumber ){
    trin->Add(Form("%s/CosmicOut_TimeCalibration_%d_%d.root",ROOTFILE_WAV.c_str(),tmpRunNumber,IterationNumber));
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
  Double_t FitP0[2];
  Double_t FitP1[2];
  Double_t FitChisq[2];
  Double_t CSIDigiDeltaT0[nCSI];//nCSIDigi
  Double_t CSIDigiDeltaT1[nCSI];//nCSIDigi
  Int_t    CosmicTrigUp;
  Int_t    CosmicTrigDn;
  Double_t Roh;
  Double_t Theta;


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
  trin->SetBranchStatus( "CSIDigiDeltaT0" );
  trin->SetBranchStatus( "CSIDigiDeltaT1" );
  trin->SetBranchStatus( "FitP0"          );
  trin->SetBranchStatus( "FitP1"          );
  trin->SetBranchStatus( "FitChisq"       );
  trin->SetBranchStatus( "CosmicTrigUp"   );
  trin->SetBranchStatus( "CosmicTrigDn"   );
  trin->SetBranchStatus( "Roh"            );
  trin->SetBranchStatus( "Theta"          );


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
  trin->SetBranchAddress( "CSIDigiDeltaT0" , CSIDigiDeltaT0  );
  trin->SetBranchAddress( "CSIDigiDeltaT1" , CSIDigiDeltaT1  );
  trin->SetBranchAddress( "FitP0"          , FitP0           );
  trin->SetBranchAddress( "FitP1"          , FitP1           );
  trin->SetBranchAddress( "FitChisq"       , FitChisq        );
  trin->SetBranchAddress( "CosmicTrigUp"   , &CosmicTrigUp   );
  trin->SetBranchAddress( "CosmicTrigDn"   , &CosmicTrigDn   );
  trin->SetBranchAddress( "Roh"            , &Roh            );
  trin->SetBranchAddress( "Theta"          , &Theta          );

  TFile* tfout = new TFile(Form("CosmicOuthist_%d.root",IterationNumber), "recreate");

  TH2D* hisDeltaChannel = new TH2D("hisDeltaChannel","hisDeltaChannel",2716,0,2716,100,-10,10);
  TH1D* hisDelta[2716];
  TH2D* hisDeltaAll = new TH2D("hisDeltaAll","hisDeltaAll;ChannelNo;DeltaT",2716,0,2716,100,-10,10);
  TH1D* hisDeltaNoCut[2716];

  TGraphErrors* grDelta = new TGraphErrors();
  TGraphErrors* grRES   = new TGraphErrors();   

  TPostScript* ps  = new TPostScript(Form("CosmicOuthist_%d.ps",IterationNumber),111);
  TCanvas *can = new TCanvas("can","",1000*TMath::Sqrt(0.5),1000);
  can->Divide(2,3);
  for( int i  = 0; i< 6; i++){
    can->cd(i+1);
    gPad->SetGridx();
    gPad->SetGridy();
  }
  ps->NewPage();

  for( int i = 0; i< 2716; i++){
    hisDelta[i] = new TH1D(Form("hisDelta%d",i ),Form("hisDelta%d",i),100,-10,10);
    hisDeltaNoCut[i] = new TH1D(Form("hisDeltaNoCut%d",i),Form("hisDeltaNoCut%d",i),100,-10,10);
  }

  for( int ievent = 0; ievent < trin->GetEntries(); ievent++){

    trin->GetEntry(ievent);
    //if( FitChisq[1] > 7 ){ continue; }
    for( int idigi = 0; idigi < nCSIDigi ; idigi++){
      hisDeltaNoCut[ CSIDigiID[ idigi ] ]->Fill( CSIDigiDeltaT1[ idigi ] );
      hisDeltaAll->Fill( CSIDigiID[ idigi  ] , CSIDigiDeltaT1[ idigi ] );
    
      if( CSIDigiID[idigi] < 2240 ){ // Case of Small Crystal 
	if( CSIDigiE[idigi] > 11 && CSIDigiE[idigi] < 17 ){
	  hisDelta[ CSIDigiID[ idigi ] ]->Fill( CSIDigiDeltaT1[ idigi ] );
	  hisDeltaChannel->Fill( CSIDigiID[ idigi  ] , CSIDigiDeltaT1[ idigi ] );
	}
      }else{ // Case of Large Crystal ... 
	if( CSIDigiE[idigi] > 22 && CSIDigiE[idigi] < 34 ){
	  hisDelta[ CSIDigiID[ idigi ] ]->Fill( CSIDigiDeltaT1[ idigi ] );
	  hisDeltaChannel->Fill( CSIDigiID[ idigi  ] , CSIDigiDeltaT1[ idigi ] );
	}
      }
      //std::cout  << CSIDigiID[ idigi ] << std::endl;
    }
  }

  for( int i = 0; i< 2716; i++){
    //std::cout << hisDelta[i]->GetEntries() << std::endl;
    if( hisDelta[i]->GetEntries() > 10){
      int rst = hisDelta[i]->Fit("gaus","Q","",hisDelta[i]->GetBinCenter( hisDelta[i]->GetMaximumBin() ) - 2, hisDelta[i]->GetBinCenter( hisDelta[i]->GetMaximumBin() ) + 2);
      TF1* func = NULL;
      func = hisDelta[i]->GetFunction("gaus");
      if( func != NULL ){
	double mean = func->GetParameter(1);
	double sigma = func->GetParameter(2);
	int rst_1 = hisDelta[i]->Fit("gaus","Q","",mean-1.5*sigma, mean+1.5*sigma);
	func = hisDelta[i]->GetFunction("gaus");
	grDelta->SetPoint( grDelta->GetN(), i, func->GetParameter(1));
	grDelta->SetPointError( grDelta->GetN()-1, 0, func->GetParError(2));
	grRES->SetPoint( grRES->GetN() , i , func->GetParameter(2));
      }
    }
 
    can->cd( i%6 + 1);
    hisDelta[i]->Draw();
    if( (i+1)%6 == 0){
      can->Modified();
      can->Update();
      ps->NewPage();
    }

    hisDelta[ i ]  ->Write(); 
    hisDeltaNoCut[ i]->Write();
  }

  int    ID[2716];
  double Delta[2716];
  double Resolution[2716];
  for( int i = 0; i< 2716; i++){
    Resolution[i] = 0xFFFF;
    Delta[i]      = 0;
  }


  int tmpID;
  double tmpDelta;
  double tmpResolution;
  if( IterationNumber > 0){
    std::ifstream ifs(Form("CosmicTimeDeltaResolution_%d.dat",IterationNumber));
    while ( ifs >> tmpID >> tmpDelta >> tmpResolution ){
      Delta[ tmpID ] = tmpDelta;
      Resolution[ tmpID ] = tmpResolution;
    }
    ifs.close();
  }

  std::ofstream ofs(Form("CosmicOutTimeDeltaResolution_%d.dat",IterationNumber));
  for( int i = 0; i< grRES->GetN(); i++){
    Delta[(int)(grDelta->GetX()[i])] += grDelta->GetY()[i];
    Resolution[(int)(grRES->GetX()[i])] = grRES->GetY()[i] ;
  }
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
  tfout->Close();
  ps->Close();
  ofs.close();

}
