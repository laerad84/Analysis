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

  TFile* tfout = new TFile("CosmicOut_Resolution.root", "recreate");

  TH2D* hisDeltaChannel = new TH2D("hisDeltaChannel","hisDeltaChannel",2716,0,2716,100,-10,10);

  const int nRegion = 4;

  TH1D* hisDelta[2716][nRegion];
  Double_t Mean[2716][nRegion]  = {{0.}}; 
  Double_t Sigma[2716][nRegion] = {{0XFFFF}};

  TGraphErrors* grSigmaChannel[2716];
  TGraphErrors* grDelta[nRegion];
  TGraphErrors* grRes[nRegion];  

  for( int i= 0 ; i< nRegion; i++){
    grDelta[i] = new TGraphErrors();
    grRes[i]   = new TGraphErrors(); 
    grDelta[i]->SetNameTitle(Form("TimeDelta%d",i),Form("TimeDelta%d",i));
    grRes[i]->SetNameTitle(Form("TimeResolution%d",i),Form("TimeResolution%d",i));
  }
  Double_t E[nRegion+1] ={4,8,12,16,20};
  Double_t LE[nRegion+1]={8,16,24,32,40};
  
  for( int i = 0; i< 2716; i++){
      grSigmaChannel[i]= new TGraphErrors();
      grSigmaChannel[i]->SetNameTitle(Form("grSigmaChannel%d",i),Form("SigmaChannel%d",i));
      for( int iE = 0 ; iE < nRegion; iE++){	
      hisDelta[i][iE] = new TH1D(Form("hisDelta_%d_%d",i,iE),Form("hisDelta_%d_%d",i,iE),
				 100,-10,10);
    }
  }
  
  std::ifstream ifs("CosmicOut_V4_TimeDeltaResolution.dat");
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
  
  TPostScript* ps       = new TPostScript("CosmicOut_Resolution.ps",111);
  TCanvas *can          = new TCanvas("can","",700*TMath::Sqrt(0.5),700);
  can->Divide(2,3);  
  for( int i  = 0; i< 6; i++){
     can->cd(i+1);
     gPad->SetGridx();
     gPad->SetGridy();
   }
   ps->NewPage();

   for( int ievent = 0; ievent < trin->GetEntries(); ievent++){
     trin->GetEntry(ievent);
     for( int idigi = 0; idigi < nCSIDigi ; idigi++){
       if( CSIDigiID[ idigi ] < 2240 ){
	 for( int iE = 0; iE < nRegion; iE++){
	   if( CSIDigiE[ idigi ] <= E[ iE ] ) break; 
	   if( CSIDigiE[ idigi ] >  E[ iE ] &&
	       CSIDigiE[ idigi ] <= E[ iE+1] ){
	     hisDelta[ CSIDigiID[ idigi ] ][iE]->Fill( CSIDigiDeltaT1[ idigi ] );
	   }else{ continue; }
	 }
       }else{
	 for( int iE = 0; iE < nRegion; iE++){
	   if( CSIDigiE[ idigi ] <= LE[ iE ] ) break; 
	   if( CSIDigiE[ idigi ] >  LE[ iE ] &&
	       CSIDigiE[ idigi ] <= LE[ iE+1 ] ){
	     hisDelta[ CSIDigiID[ idigi ] ][iE]->Fill( CSIDigiDeltaT1[ idigi ] );
	   }else{ continue; }
	 }
       }
       hisDeltaChannel->Fill( CSIDigiID[ idigi  ] , CSIDigiDeltaT1[ idigi ] );
     }
   }
   
   for( int i = 0; i< 2716; i++){
     //std::cout << hisDelta[i]->GetEntries() << std::endl;
     for( int iE = 0; iE < nRegion; iE++ ){
       	 can->cd( iE + 1  );
       if( hisDelta[i][iE]->GetEntries() > 10){
	 int rst = hisDelta[i][iE]->Fit("gaus","Q","",
					hisDelta[i][iE]->GetBinCenter( hisDelta[i][iE]->GetMaximumBin() )-2*hisDelta[i][iE]->GetRMS(), 
					hisDelta[i][iE]->GetBinCenter( hisDelta[i][iE]->GetMaximumBin() )+2*hisDelta[i][iE]->GetRMS());
	 TF1* func = NULL;
	 func = hisDelta[i][iE]->GetFunction("gaus");
	 if( func != NULL ){
	   double mean  = func->GetParameter(1);
	   double sigma = func->GetParameter(2);
	   int rst_1 = hisDelta[i][iE]->Fit("gaus","Q","",mean-2*sigma, mean+2*sigma );
	   func = hisDelta[i][iE]->GetFunction("gaus");
	   grDelta[iE]->SetPoint( grDelta[iE]->GetN(), i, func->GetParameter(1));
	   grDelta[iE]->SetPointError( grDelta[iE]->GetN()-1, 0, func->GetParError(2));
	   grRes[iE]->SetPoint( grRes[iE]->GetN() , i , func->GetParameter(2));
	   Mean[i][iE] = func->GetParameter(1);
	   Sigma[i][iE] = func->GetParameter(2);
	   if( i <2240){
	     grSigmaChannel[i]->SetPoint(grSigmaChannel[i]->GetN(),(E[iE]+E[iE+1])/2,func->GetParameter(2));
	     grSigmaChannel[i]->SetPointError(grSigmaChannel[i]->GetN()-1,(E[iE+1]-E[iE])/2,func->GetParError(2));
	   }else{
	     grSigmaChannel[i]->SetPoint(grSigmaChannel[i]->GetN(),(LE[iE]+LE[iE+1])/2,func->GetParameter(2));
	     grSigmaChannel[i]->SetPointError(grSigmaChannel[i]->GetN()-1,(LE[iE+1]-LE[iE])/2,func->GetParError(2));
	   }
	 }
	 hisDelta[i][iE]->Draw();
       }
     }
     can->cd( nRegion + 1 );
     grSigmaChannel[i]->Draw("AP");
     
     
     can->Modified();
     can->Update();
     ps->NewPage();
    
  }

   std::ofstream ofs("Delta_Resolution_by_Energy.dat");
   for( int i = 0; i< 2716; i++){
     ofs << i ;
     for( int iE = 0; iE < nRegion; iE++){       
       ofs << "\t" << Mean[i][iE];
     } 
     for( int iE = 0; iE < nRegion; iE++){
       ofs << "\t" << Sigma[i][iE];
     }
     ofs << "\n";
   }

  for( int i = 0; i< 2716; i++){
    grSigmaChannel[i]->Write();
    for( int iE = 0; iE < nRegion; iE++ ){
      hisDelta[i][iE] ->Write();      
    }
  }
  
  for( int iE = 0; iE < nRegion; iE++){
    grDelta[iE]->Write();
    grRes[iE]->Write();
  }

  hisDeltaChannel->Write();
  tfout->Close();
  ps->Close();
  ofs.close();

}
