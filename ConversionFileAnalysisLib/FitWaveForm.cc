#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TPDF.h"
#include "TPostScript.h"
#include "TApplication.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "WaveformFitter.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TPDF.h"

int main( int argc, char** argv ){
  std::cout<< "Start" << std::endl;
  std::string InputFilename = argv[1];

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  int Integer=0;
  char** Files=NULL;
  std::cout<< "Start" << std::endl;
  TApplication* app = new TApplication( "app", &Integer, Files );
  TFile * inputFile = new TFile(InputFilename.c_str());
  TTree * inputTree = (TTree*)inputFile->Get("WFTree");
  Int_t Data[48]; 
  Int_t Time[48]; 
  Int_t ID;
  Int_t EventNumber;
  inputTree->SetBranchAddress("Data"       ,Data);
  inputTree->SetBranchAddress("Time"       ,Time);
  inputTree->SetBranchAddress("ID"         ,&ID);
  inputTree->SetBranchAddress("EventNumber",&EventNumber);
  std::cout<< __LINE__ << std::endl;

  TGraph* gr     = new TGraph();
  TGraph* grTemp = new TGraph();

  TFile* outputFile = new TFile("Template_Cosmic.root","RECREATE");
  
  TH2D* hisTemp[2716];
  TH1D* hisAspect[2716]; 
  for( int i = 0; i< 2716; i++){
    hisTemp[i]   = new TH2D(Form("hisTemp%d",i),"", 400,-200,200,450,-0.325,1.8);
    hisAspect[i] = new TH1D(Form("hisAspect%d",i), "",100,0,1); 
  }

  gr->SetMarkerStyle(6);
  TCanvas* can = new TCanvas("can","",0,0,1600,800);
  WaveformFitter* wav = new WaveformFitter( 48, kFALSE );
  TH1D* his = new TH1D("hisChsq","", 200,0,20);
  std::cout<< "Loop" << std::endl;
  for( int ievent  = 0; ievent < inputTree->GetEntries(); ievent++){    
    inputTree->GetEntry(ievent);
    gr->Set(0);
    wav->Clear();
    Double_t Min=20000;
    Double_t Max=0;
    Int_t    MinIndex = 0; 
    Int_t    MaxIndex = 0;
    for( int index = 0; index < 48; index++){
      if( Data[index] > Max ) { Max = Data[index] ; MaxIndex = index;}
      if( Data[index] < Min ) { Min = Data[index] ; MinIndex = index;}
    }
    
    //if( Max - Min < 200 ){continue; }

    for( int ipoint  = 0; ipoint < 48; ipoint++){
      gr->SetPoint( ipoint , ipoint*8 , Data[ipoint] );
    }
    Double_t gndRMS[42];
    Double_t gndMean[42]; 
    Double_t currentMinimumRMS=999999;
    Double_t currentGND=0;
    Int_t    currentPosition=0;
    for( int iCombination  = 0; iCombination < 48 -7 +1 ; iCombination++){ 
      Double_t MaxValue=0;
      for( int ipoint = iCombination; ipoint < iCombination + 7; ipoint++){
	if( MaxValue < Data[ipoint] ){ MaxValue = Data[ipoint];}
      }
      gndRMS[iCombination]  = 0.;
      gndMean[iCombination] = 0.;
      for( int ipoint = iCombination; ipoint < iCombination + 7; ipoint++){
	gndRMS[iCombination]+=pow( Data[ipoint] - MaxValue,2);
	gndMean[iCombination]+=Data[ipoint];
      }
      gndRMS[iCombination]  = sqrt(gndRMS[iCombination])/(7.-1.);
      gndMean[iCombination]/=7.;
      if( gndRMS[iCombination] < currentMinimumRMS ){ 
	currentMinimumRMS = gndRMS[iCombination];
	currentGND        = gndMean[iCombination];
	currentPosition   = 8*iCombination;
      }
    }

    wav->Fit(gr);
    TF1* func = wav->GetFunction();
    //if(func->GetChisquare()/func->GetNDF()<1000 ){ continue; }
    //if(func->GetChisquare()/func->GetNDF()/func->GetParameter(0) >= 10 ){ continue; }
    //if( func->GetParameter(0)> 2000 || func->GetParameter(0) < 1000 ){ continue; }
    //if( func->GetParameter(1)< 80   || func->GetParameter(1) >40*8 ){ continue; }
    Double_t Sum =0;
    if( func->GetParameter(0) != 0 ){
      for( int i = 0; i< 48; i++){
	Sum += (Data[i]-func->GetParameter(4))/(func->GetParameter(0)*func->GetParameter(2));
      }
      hisAspect[ID]->Fill(Sum);
    }
    
    for( int i = 0; i< 48; i++){
      grTemp->SetPoint( i, Time[i]- func->GetParameter(1),
			(Data[i]- func->GetParameter(4))/func->GetParameter(0) );
      hisTemp[ID]->Fill( Time[i]- func->GetParameter(1),
			 (Data[i]- func->GetParameter(4))/func->GetParameter(0) );
    }
    
    his->Fill(func->GetChisquare()/func->GetNDF()/func->GetParameter(0));

    //gr->Fit(func,"","",MaxIndex*8-48, MaxIndex*8+48);
    //can->Update();
    //can->Modified();
    //getchar();

  }
  std::cout<< "End Loop" << std::endl;
  
  //can->Divide( 5,4 );
  TPDF* pdf = new TPDF("dist.pdf",111);
  for( int i = 0; i< 2716; i++){
    //can->cd( i+1 );
    hisTemp[i]->Draw("col");
    //prof[i] = new TProfile(Form("prof%d",i),"",250,-100,150);
    //prof[i] = hisTemp[i]->ProfileX();
    can->Update();
  }
  
  for( int i = 0; i< 2716; i++){
    gPad->SetLogy();
    hisAspect[i]->Draw();    
    //prof[i] = new TProfile(Form("prof%d",i),"",250,-100,150);
    //prof[i] = hisTemp[i]->ProfileX();
    can->Update();
  }

  pdf->Close();

  TProfile* prof[2716];
  for( int i = 0; i< 2716; i++){
    hisTemp[i]->Write();
    hisAspect[i]->Write();
    //prof[i] = new TProfile(Form("prof%d",i),"",250,-100,150);
    //prof[i] = hisTemp[i]->ProfileX();
  }
  outputFile->Close();


  /*
  //his->Draw();
  can->Divide(2,5);
  for( int i = 0; i< 10; i++){
    can->cd(i+1);
    gPad->SetLogz();
    //hisTemp[i]->Draw("col");
    prof[i]->Draw("same");
  }
  can->Update();
  can->Modified();
  */

  //app->Run();
  return 0;
}
    
