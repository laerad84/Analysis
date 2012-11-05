#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "EventTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TApplication.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TProfile.h"

void RotationTheta( Double_t  Theta, Double_t  x, Double_t  y, Double_t & nx,Double_t & ny){
  nx = x*TMath::Cos( Theta ) - y*TMath::Sin( Theta );
  ny = x*TMath::Sin( Theta ) + y*TMath::Cos( Theta );
} 
void ConvertIndex( Int_t RTNIndex, Int_t& iBinsX, Int_t& iBinsY, Int_t NxBin, Int_t NyBin){
  if( RTNIndex <= 0 ){
    iBinsX = -1;
    iBinsY = -1; 
    return;
  }
  iBinsY = (int)(RTNIndex / (NxBin+2));
  iBinsX = (int)(RTNIndex % (NyBin+2));  
}

int main( int argc , char** argv ){

  TApplication* app = new TApplication("app",&argc, argv);

  TChain* chain = new TChain("eventTree00");
  chain->Add(Form("/Volume0/gamma/template_gamma_210MeV_10deg-1E5-0.root"));
  EventTree* trin = new EventTree( chain );

  trin->Show(0);
  std::cout<<trin->fChain->GetEntries() << std::endl;

  TCanvas* can = new TCanvas("can","can",800,800);
  can->Divide(2,2);
  
  TH2D* hisEdep = new TH2D("hisEdep","hisEdep",14,0-7*25,14*25-7*25, 14,0-7*25,14*25-7*25);
  TH2D* hisEdepZ = new TH2D("hisEdepZ","hisEdepZ",14,0-7*25,14*25-7*25, 500,0,500);
  Double_t DepEnergy[14][14]={{0}};
  Double_t XArr[14] = {0};
  Double_t YArr[14] = {0};
  for( int i = 0; i< 14; i++){
    XArr[i] = -7*25 + 25*i + 12.5;
    YArr[i] = -7*25 + 25*i + 12.5; 
  }
    
  Double_t MeanX=0;
  Double_t MeanY=0;
  
  for( int ievent  = 0; ievent < trin->fChain->GetEntries(); ievent++ ){
    trin->GetEntry(ievent);
    hisEdep->Reset();
    hisEdepZ->Reset();

    double TotalEnergy = 0;
    for( int ihit = 0; ihit < trin->CSI_hits_; ihit++){
      for( int iclusterx = 0; iclusterx < 14; iclusterx++){
	for( int iclustery  =0; iclustery < 14; iclustery++){
	  DepEnergy[iclusterx][iclustery] = 0;
	}
      }
      std::cout<< trin->CSI_hits_r_fX[ihit] <<  " : " 
	       << trin->CSI_hits_r_fY[ihit] << std::endl;
      Int_t RTNIndex = hisEdep->Fill( trin->CSI_hits_r_fX[ihit],
				      trin->CSI_hits_r_fY[ihit], 
				      trin->CSI_hits_edep[ihit]);
      hisEdepZ->Fill( trin->CSI_hits_r_fX[ihit],
		      trin->CSI_hits_r_fZ[ihit],
		      trin->CSI_hits_edep[ihit]);
      Int_t iBinsX;
      Int_t iBinsY;
      ConvertIndex( RTNIndex, iBinsX, iBinsY, 14,14);
      TotalEnergy += trin->CSI_hits_edep[ihit];
      DepEnergy[iBinsX][iBinsY] += trin->CSI_hits_edep[ihit];
      MeanX += trin->CSI_hits_edep[ihit]*XArr[iBinsX-1];
      MeanY += trin->CSI_hits_edep[ihit]*YArr[iBinsY-1];

    }
    MeanX = MeanX/TotalEnergy;
    MeanY = MeanY/TotalEnergy;

    can->cd(1);
    gPad->SetLogz();
    hisEdep->Draw("colz");
    TMarker* mk = new TMarker( MeanX, MeanY,6);
    mk->Draw();
    can->cd(2);
    gPad->SetLogz();
    hisEdepZ->Draw("colz");
    hisEdepZ->ProfileX()->Draw("same");
    can->Update();
    can->Modified();
    getchar();

  }
  app->Run();
  app->Terminate();
  return 0; 
}
