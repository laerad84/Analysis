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
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat(111111111);
  
  TFile* tf = new TFile("CosmicOut_V4.root");
  TTree* trin = (TTree*)tf->Get("trOut");

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

  TFile* tfout = new TFile("CosmicOut_V4_hist.root", "recreate");

  TH2D* hisDeltaChannel = new TH2D("hisDeltaChannel","hisDeltaChannel",2716,0,2716,100,-10,10);
  TH1D* hisDelta[2716];
  TGraphErrors* grDelta = new TGraphErrors();
  TGraphErrors* grRES   = new TGraphErrors(); 
  

  TPostScript* ps  = new TPostScript("CosmicOut_V4_hist.ps",111);
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
  }

  for( int ievent = 0; ievent < trin->GetEntries(); ievent++){

    trin->GetEntry(ievent);
    
    for( int idigi = 0; idigi < nCSIDigi ; idigi++){
      hisDelta[ CSIDigiID[ idigi ] ]->Fill( CSIDigiDeltaT1[ idigi ] );
      hisDeltaChannel->Fill( CSIDigiID[ idigi  ] , CSIDigiDeltaT1[ idigi ] );
      //std::cout  << CSIDigiID[ idigi ] << std::endl;
    }
  }

  for( int i = 0; i< 2716; i++){
    //std::cout << hisDelta[i]->GetEntries() << std::endl;
    if( hisDelta[i]->GetEntries() > 10){
      int rst = hisDelta[i]->Fit("gaus","Q","",
				 hisDelta[i]->GetBinCenter( hisDelta[i]->GetMaximumBin() )-hisDelta[i]->GetRMS(), 
				 hisDelta[i]->GetBinCenter( hisDelta[i]->GetMaximumBin() )+hisDelta[i]->GetRMS());
      TF1* func = NULL;
      func = hisDelta[i]->GetFunction("gaus");
      if( func != NULL ){
	double mean  = func->GetParameter(1);
	double sigma = func->GetParameter(2);
	int rst_1 = hisDelta[i]->Fit("gaus","Q","",mean-1.5*sigma, mean+1.5*sigma );
	func = hisDelta[i]->GetFunction("gaus");
	grDelta->SetPoint( grDelta->GetN(), i, func->GetParameter(1));
	grDelta->SetPointError( grDelta->GetN()-1, 0, func->GetParError(2));
	grRES->SetPoint( grRES->GetN() , i , func->GetParameter(2));
      }
    }
 
    can->cd( (i%6) + 1);
    hisDelta[i]->Draw();
    if( ((i+1)%6) == 0){
      can->Modified();
      can->Update();
      ps->NewPage();
    }

    hisDelta[ i ]  ->Write(); 
    
  }

  std::ifstream ifs("CosmicOut_V3_TimeDeltaResolution.dat");
  std::ofstream ofs("CosmicOut_V4_TimeDeltaResolution.dat");
  int    ID[2716];
  double Delta[2716];
  double Resolution[2716];
  for( int i = 0; i< 2716; i++){
    Resolution[i] = 0xFFFF;
    Delta[i]      = 0xFFFF;
  }
  int tmpID;
  double tmpDelta;
  double tmpResolution;
  while ( ifs >> tmpID >> tmpDelta >> tmpResolution ){
    Delta[ tmpID ] = tmpDelta;
    Resolution[ tmpID ] = tmpResolution;
  }
  ifs.close();
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
  hisDeltaChannel->Write();
  tfout->Close();
  ps->Close();
  ofs.close();

}
