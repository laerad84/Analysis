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
#include "TRandom.h"

#include "CsI_Module.h"
#include "CsIPoly.h"
#include "IDHandler.h"
#include "PulseGenerator.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include "gamma/Gamma.h"
#include "cluster/Cluster.h"

#include "gamma/GammaFinder.h"
#include "cluster/ClusterFinder.h"

#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "User_Function.h"


void RotationTheta( Double_t  Theta, Double_t  x, Double_t  y, Double_t & nx,Double_t & ny){
  nx = x*TMath::Cos( Theta ) - y*TMath::Sin( Theta );
  ny = x*TMath::Sin( Theta ) + y*TMath::Cos( Theta );
} 

void ConvertPosition( Double_t Radius, Double_t Theta, Double_t x, Double_t y, Double_t& nx, Double_t& ny ){
  nx = (x+Radius)*TMath::Cos( Theta ) - y*TMath::Sin( Theta );
  ny = (x+Radius)*TMath::Sin( Theta ) + y*TMath::Cos( Theta );
}

int main( int argc , char** argv ){

  TApplication* app = new TApplication("app",&argc, argv);

  CsIPoly* CsIEne = new CsIPoly("CsIEnergy","CsIEnergy");
  CsIPoly* CsIPos = new CsIPoly("CsIPosition","CsIPosition");
  CsIPoly* CsIEneT = new CsIPoly("CsIEnergyT","CsIEnergyT");
  PulseGenerator* gen = new PulseGenerator();
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
  
  Int_t nDigi;
  Int_t ID[2716];
  Double_t Energy[2716];
  Double_t SignalTime[2716];

  for( int ievent  = 0; ievent < trin->fChain->GetEntries(); ievent++ ){
    trin->GetEntry(ievent);
    hisEdep->Reset();
    hisEdepZ->Reset();
    CsIEne->Reset();
    CsIPos->Reset();
    CsIEneT->Reset();
    double TotalEnergy = 0;
    std::vector<int> CrystalIDVec;
    std::vector<double> CrystalEVec;
    std::vector<int>::iterator    itIDVec;
    std::vector<double>::iterator itEVec;
    std::vector<int>    hitIDVec;
    std::vector<double> hitZVec;
    std::vector<double> hitEVec;
    std::vector<double> hitTVec;

    Int_t  CsIID = -1;
    const double Speed_of_Signal = 100;

    for( int ihit = 0; ihit < trin->CSI_hits_; ihit++){

      std::cout<< trin->CSI_hits_r_fX[ihit] <<  " : " 
	       << trin->CSI_hits_r_fY[ihit] << std::endl;
      if( trin->CSI_hits_edep[ihit] == 0 ){ continue; }
      Int_t RTNIndex = hisEdep->Fill( trin->CSI_hits_r_fX[ihit],
				      trin->CSI_hits_r_fY[ihit], 
				      trin->CSI_hits_edep[ihit]);
      hisEdepZ->Fill( trin->CSI_hits_r_fX[ihit],
		      trin->CSI_hits_r_fZ[ihit],
		      trin->CSI_hits_edep[ihit]);

      TotalEnergy += trin->CSI_hits_edep[ihit];
      Double_t Radius = 300;//mm
      Double_t Theta  = 2./6*TMath::Pi();//rad
      //Double_t Theta  = 0;
      Double_t CsiPosX,CsiPosY;
      ConvertPosition(Radius, Theta , trin->CSI_hits_r_fX[ihit], trin->CSI_hits_r_fY[ihit] ,CsiPosX, CsiPosY );
      CsIID  = CsIEne->Fill( CsiPosX, CsiPosY, trin->CSI_hits_edep[ihit] ) -1;      
      if( CsIID <0 ){ continue; }

      hitIDVec.push_back( CsIID );
      hitZVec.push_back( trin->CSI_hits_r_fZ[ihit] );
      hitEVec.push_back( trin->CSI_hits_edep[ihit] );
      hitTVec.push_back( trin->CSI_hits_time[ihit] + (500 - trin->CSI_hits_r_fZ[ihit])/Speed_of_Signal );

      for( itIDVec = CrystalIDVec.begin(); itIDVec != CrystalIDVec.end(); itIDVec++){
	if( (*itIDVec) == CsIID ){ break; }
      }
      if( itIDVec == CrystalIDVec.end() ){
	CrystalIDVec.push_back( CsIID );
      }
    }    
    /*
    for( int index = 0; index < hitIDVec.size() ; index++){
      std::cout<< hitIDVec[index] << " : " << hitZVec[index] << " : " << hitEVec[index] << std::endl;
    }
    */

    for( itIDVec = CrystalIDVec.begin(); itIDVec != CrystalIDVec.end() ; itIDVec++){
      Double_t Energy = CsIEne->GetBinContent( (*itIDVec) + 1 );
      CrystalEVec.push_back(Energy);
      //std::cout<<  *itIDVec << " : " << Energy << " : " << TotalEnergy  << std::endl;
    }

    itIDVec = CrystalIDVec.begin();
    itEVec  = CrystalEVec.begin();

    nDigi = 0;
    for( int i = 0; i< 2716; i++){
      ID[i] = -1;
      Energy[i] = 0; 
      SignalTime[i] = 0; 
    }

    for(;itIDVec != CrystalIDVec.end() && itEVec != CrystalEVec.end() ; itIDVec++, itEVec++ ){
      // Making Waveform // 

      if( *itEVec < 3 ){ continue; }

      std::vector<double> EVec;
      std::vector<double> TVec;
      for( int ihit = 0; ihit < hitZVec.size(); ihit++){
	if( hitIDVec[ihit] == (*itIDVec)){
	  EVec.push_back( hitEVec[ihit] );
	  TVec.push_back( hitTVec[ihit] );
	}
	/// Get Waveform // 
      }

      TF1* func = gen->GetWaveform( EVec, TVec , 150./13.7);
      //std::cout<< *itEVec << " MeV" << " : " <<  func->GetMaximum( 0, 500)/(150./13.7) <<  std::endl;
      ID[nDigi]     = (*itIDVec);
      Energy[nDigi] = func->GetMaximum( 0, 500 )/ ( 150./13.7);
      SignalTime[nDigi] = func->GetMaximumX( 0, 500 );
      nDigi++;
      delete func;
    }    

    /// prepared to clustering /// 
    




    
    
    for( int ibin = 0; ibin < 2716; ibin++){
      if( CsIEne->GetBinContent( ibin +1 ) > 3 ){
	CsIEneT->Fill( ibin, (double)CsIEne->GetBinContent( ibin +1 ));
      }
    }



    /*
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
    */
  }
  app->Run();
  app->Terminate();
  return 0; 
}
