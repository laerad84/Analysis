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
#include "CsI_Module.h"
#include "CsIPoly.h"
#include "IDHandler.h"


void RotationTheta( Double_t  Theta, Double_t  x, Double_t  y, Double_t & nx,Double_t & ny){
  nx = x*TMath::Cos( Theta ) - y*TMath::Sin( Theta );
  ny = x*TMath::Sin( Theta ) + y*TMath::Cos( Theta );
} 

void ConvertPosition( Double_t Radius, Double_t Theta, Double_t x, Double_t y, Double_t& nx, Double_t& ny ){
  nx = (x+Radius)*TMath::Cos( Theta ) - y*TMath::Sin( Theta );
  ny = (x+Radius)*TMath::Sin( Theta ) - y*TMath::Cos( Theta );
}

int main( int argc , char** argv ){

  TApplication* app = new TApplication("app",&argc, argv);

  CsIPoly* CsIEne = new CsIPoly("CsIEnergy","CsIEnergy");
  CsIPoly* CsIPos = new CsIPoly("CsIPosition","CsIPosition");
  CsIPoly* CsIEneT = new CsIPoly("CsIEnergyT","CsIEnergyT");

  TChain* chain = new TChain("eventTree00");
  chain->Add(Form("/Volume0/gamma/template_gamma_210MeV_10deg-1E5-0.root"));
  EventTree* trin = new EventTree( chain );

  trin->Show(0);
  std::cout<<trin->fChain->GetEntries() << std::endl;

  TCanvas* can = new TCanvas("can","can",800,800);
  can->Divide(2,2);
  
  TH2D* hisEdep = new TH2D("hisEdep","hisEdep",14,0-7*25,14*25-7*25, 14,0-7*25,14*25-7*25);
  TH2D* hisEdepZ = new TH2D("hisEdepZ","hisEdepZ",14,0-7*25,14*25-7*25, 500,0,500);
  
  Double_t MeanX=0;
  Double_t MeanY=0;
  
  for( int ievent  = 0; ievent < trin->fChain->GetEntries(); ievent++ ){
    trin->GetEntry(ievent);
    hisEdep->Reset();
    hisEdepZ->Reset();
    CsIEne->Reset();
    CsIPos->Reset();
    CsIEneT->Reset();
    double TotalEnergy = 0;
    std::vector<int> CrystalIDVec;
    std::vector<double> CrystalEnergy;
    for( int ihit = 0; ihit < trin->CSI_hits_; ihit++){

      std::cout<< trin->CSI_hits_r_fX[ihit] <<  " : " 
	       << trin->CSI_hits_r_fY[ihit] << std::endl;
      Int_t RTNIndex = hisEdep->Fill( trin->CSI_hits_r_fX[ihit],
				      trin->CSI_hits_r_fY[ihit], 
				      trin->CSI_hits_edep[ihit]);
      hisEdepZ->Fill( trin->CSI_hits_r_fX[ihit],
		      trin->CSI_hits_r_fZ[ihit],
		      trin->CSI_hits_edep[ihit]);

      TotalEnergy += trin->CSI_hits_edep[ihit];
      Double_t Radius = 300;//mm
      Double_t Theta  = 1./6*TMath::Pi();//rad
      Double_t CsiPosX,CsiPosY;
      ConvertPosition(Radius, Theta , trin->CSI_hits_r_fX[ihit], trin->CSI_hits_r_fY[ihit] ,CsiPosX, CsiPosY );
      CsIEne->Fill( CsiPosX, CsiPosY, trin->CSI_hits_edep[ihit] );
    }
    for( int ibin = 0; ibin < 2716; ibin++){
      if( CsIEne->GetBinContent( ibin +1 ) > 0.5 ){
	CsIEneT->Fill( ibin, (double)CsIEne->GetBinContent( ibin +1 ));
      }
    }

    can->cd(1);
    gPad->SetLogz();
    hisEdep->Draw("colz");
    can->cd(2);
    gPad->SetLogz();
    CsIEne->Draw("colz L");
    can->cd(3);
    gPad->SetLogz();
    CsIEneT->Draw("colz L");
    can->Update();
    can->Modified();
    getchar();

  }
  app->Run();
  app->Terminate();
  return 0; 
}
