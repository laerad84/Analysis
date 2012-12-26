#include <cstdio>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
Double_t QuadSum( double a, double b){
  return TMath::Sqrt( a*a + b*b );
}

void DrawVelocity() {
  gStyle->SetOptFit(11111111);

  TFile* tfCosmic = new TFile("CosmicOut_V4.root");
  TTree* trCosmic = (TTree*)tfCosmic->Get("trOut");
  
  Double_t FitP1[2];
  Double_t FitP0[2];
  Double_t FitChisq[2];
  Int_t    CosmicTrigUp;
  Int_t    CosmicTrigDn;
  Double_t Theta;
  Double_t Roh;

  trCosmic->SetBranchAddress("FitP1"       ,FitP1 );
  trCosmic->SetBranchAddress("FitP0"       ,FitP0 );
  trCosmic->SetBranchAddress("FitChisq"    ,FitChisq);
  trCosmic->SetBranchAddress("CosmicTrigUp",&CosmicTrigUp);
  trCosmic->SetBranchAddress("CosmicTrigDn",&CosmicTrigDn);
  trCosmic->SetBranchAddress("Theta"       ,&Theta);
  trCosmic->SetBranchAddress("Roh"         ,&Roh);
  
  std::cout<< "Output" << std::endl;
  //TFile* tfHist = new TFile("CosmicOut_V4_TrigHist.root","recreate");
  TH1D* CosmicP1UpTrig[5];
  TH1D* CosmicP1Trig[5][5];
  TH1D* CosmicP1TrigLarge[5][5];

  for( int i = 0; i< 5; i++){
    CosmicP1UpTrig[i] = new TH1D( Form("CosmicP1UpTrig%d",i),
				  Form("CosmicP1UpTrig%d",i),
				  100,-0.01,0.01);
    for( int j = 0; j< 5; j++){
      CosmicP1Trig[i][j] = new TH1D( Form("CosmicP1_%d_%d",i,j),
				     Form("CosmicP1_%d_%d",i,j),
				     100,-0.01,0.01);
      CosmicP1TrigLarge[i][j] = new TH1D( Form("CosmicP1Large_%d_%d",i,j),
					  Form("CosmicP1Large_%d_%d",i,j),
					  100,-0.01,0.01);
    }
  }

  std::cout<< "LOOP" << std::endl;
  for( int ievent  = 0; ievent < trCosmic->GetEntries(); ievent++){
    trCosmic->GetEntry(ievent);
    if( TMath::Abs(Theta) < 5 ){ continue; }
    if( FitChisq[1] > 10 ){ continue; }
    CosmicP1UpTrig[ CosmicTrigUp ]->Fill( FitP1[1]*TMath::Cos(Theta*TMath::Pi()/180.));
    CosmicP1Trig[CosmicTrigUp][CosmicTrigDn]->Fill( FitP1[1]*TMath::Cos(Theta*TMath::Pi()/180.) );
    if( Roh > 600){
      CosmicP1TrigLarge[CosmicTrigUp][CosmicTrigDn]->Fill( FitP1[1]*TMath::Cos(Theta*TMath::Pi()/180.) );
    }
  }

  std::cout <<  "Draw" << std::endl;
  TGraphErrors* grSlope = new TGraphErrors();
  TGraphErrors* grDelta = new TGraphErrors();
  TGraphErrors* grSlopeLarge = new TGraphErrors();
  grSlope->SetMarkerStyle( 7 );
  grDelta->SetMarkerStyle( 7 );
  grSlopeLarge->SetMarkerStyle( 7 );
  TCanvas* can = new TCanvas("can","",1200,800);
  can->Divide(3,2);

  for( int i = 0; i< 5; i++){
    for( int j = 0; j< 5; j++){
      CosmicP1Trig[i][j]->Fit("gaus","Q","",
			      CosmicP1Trig[i][j]->GetMean()-CosmicP1Trig[i][j]->GetRMS(), 
			      CosmicP1Trig[i][j]->GetMean()+CosmicP1Trig[i][j]->GetRMS());
      TF1* func = CosmicP1Trig[i][j]->GetFunction("gaus");
      grSlope->SetPoint( grSlope->GetN(),i*5+j+1, func->GetParameter(1));
      grSlope->SetPointError( grSlope->GetN()-1, 0, func->GetParameter(2));
      CosmicP1TrigLarge[i][j]->Fit("gaus","Q","",
				   CosmicP1Trig[i][j]->GetMean()-1.5*CosmicP1Trig[i][j]->GetRMS(), 
				   CosmicP1Trig[i][j]->GetMean()+1.5*CosmicP1Trig[i][j]->GetRMS());
      func = CosmicP1TrigLarge[i][j]->GetFunction("gaus");
      grSlopeLarge->SetPoint( grSlopeLarge->GetN(), i*5+j+1, func->GetParameter(1));
      grSlopeLarge->SetPointError( grSlopeLarge->GetN()-1,0,func->GetParameter(2));
    }
  }

  for( int i = 0; i< 5; i++){
    can->cd( i +1 );
    CosmicP1UpTrig[i]->Draw();
    for( int j = 0; j < 5; j++){
      CosmicP1Trig[i][j]->SetLineColor(j+1);
      CosmicP1Trig[i][j]->Draw("same");
      CosmicP1TrigLarge[i][j]->Draw("same");
    }
  }
  double sol = 299.97;//[mm/ns]
  TLine* solLine = new TLine( -1, -1./sol, 25, -1./sol);
  //can->cd( 6 );


  TCanvas* can1 = new TCanvas("can1","",800,400);
  can1->Divide( 2,1);
  can1->cd(1);
  grSlope->GetXaxis()->SetRangeUser( -1, 25);  
  grSlope->Draw("AP");
  solLine->SetLineStyle(2);
  solLine->Draw("same");
  TCanvas* can2 = new TCanvas("can2","",1200,800);
  can2->Divide(3,2);
  for( int i = 0; i< 5; i++){
    can2->cd( i+1 );
    CosmicP1TrigLarge[i][0]->Draw();
    for( int j = 0; j< 5; j++){
      CosmicP1TrigLarge[i][j]->Draw("same");
    }
  }
  can2->cd(6);
  grSlopeLarge->Draw("AP");

  
  Double_t H_trig = 2800;
  Double_t R      = 1200;
  TGraph* grDelta1 = new TGraph();
  TGraph* grDelta2 = new TGraph();

  TGraphErrors * grVCluster = new TGraphErrors();

  for( int i = 0; i< 5; i++){
    for( int j  =0; j< 5; j++){
      Double_t iv_0 = TMath::Abs(grSlope->GetY()[ i*5 +j ]);
      Double_t ive_0 = grSlope->GetEY()[ i*5 + j];
      

      Double_t deltaZ = (i-j)*100*R/H_trig;//mm            

      Double_t EdeltaZ = TMath::Sqrt(2)*50*R/H_trig;

      Double_t l  = TMath::Sqrt( R*R + deltaZ*deltaZ);
      Double_t El = deltaZ*EdeltaZ/TMath::Sqrt( R*R + deltaZ*deltaZ); 
      Double_t delta_t0 = l/sol;
      Double_t delta_t1 = iv_0*R - delta_t0;
      Double_t Edelta_t1 = QuadSum(ive_0*R, El/sol);
      
      Double_t V_CsI = deltaZ / delta_t1;
      Double_t EV_CsI = QuadSum( EdeltaZ / delta_t1, Edelta_t1/delta_t1);

      Double_t delta_t1_adj = deltaZ/70.;
      Double_t delta_t1_adj1 = deltaZ/80.;
      Double_t delta_t1_adj2 = deltaZ/90.;
      
      grDelta->SetPoint(grDelta->GetN(), (delta_t0+delta_t1_adj)/R , iv_0);
      grDelta1->SetPoint(grDelta1->GetN(), (delta_t0+delta_t1_adj1)/R , iv_0);
      grDelta2->SetPoint(grDelta2->GetN(), (delta_t0+delta_t1_adj2)/R , iv_0);
      
      if( TMath::Abs(deltaZ) <= 2*EdeltaZ ){ continue; }
      grVCluster->SetPoint( grVCluster->GetN(), i*5+j+1 , V_CsI);
      grVCluster->SetPointError( grVCluster->GetN()-1 , 0, EV_CsI);
      
    }
  }
  can1->cd(2);
  grVCluster->SetMarkerStyle( 6);
  grVCluster->Draw("AP");
  
  grVCluster->Fit("pol0","","",0,26);
  /*
  grDelta1->SetMarkerColor(2);
  grDelta2->SetMarkerColor(3);
  grDelta1->SetMarkerStyle(6);
  grDelta2->SetMarkerStyle(6);


  can1->cd(2);
  grDelta->Draw("AP");
  grDelta1->Draw("P");
  grDelta2->Draw("P");
  
  TF1* func = new TF1("func","x",0,0.1);
  func->Draw("same");

  */


  /*
  for( int i = 0; i< 5; i++){
    CosmicP1UpTrig[i]->Write();
    for( int j  =0 ;j < 5; j++){
      CosmicP1Trig[i][j]->Write();
    }
  }

  tfHist->Close();
  */
}

  


  

  
