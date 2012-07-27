#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>

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
#include "TSpectrum.h"

int main( int argc, char** argv ){

  std::string InputFilename = argv[1];
  TApplication* app = new TApplication( "app", &argc, argv );

  TFile * inputFile = new TFile(InputFilename.c_str());
  TTree * inputTree = (TTree*)inputFile->Get("WFTree");
  Int_t Data[48]; 
  Int_t Time[48]; 
  Int_t ID;
  Int_t EventNumber;
  inputTree->SetBranchAddress("Data",Data);
  inputTree->SetBranchAddress("Time",Time);
  inputTree->SetBranchAddress("ID",&ID);
  inputTree->SetBranchAddress("EventNumber",&EventNumber);
  TCanvas* can = new TCanvas("can","",0,0,600,800);
  can->Divide(1,2);
  TGraph* gr = new TGraph(); 
  TGraph* grRMS = new TGraph();
  gr->SetMarkerStyle(6);
  grRMS->SetMarkerStyle(6);
  TH1D* his = new TH1D("his","",48,0,48);
  //for( int i = 0; i< inputTree->GetEntries(); i++){
  TSpectrum* spec = new TSpectrum(3);
  for( int ievent = 0; ievent< 10000; ievent++){

    inputTree->GetEntry(ievent);
    gr->Set(0);
    Double_t Min = 20000;
    Double_t Max = 0;
    for( int index = 0; index < 48; index++){
      if( Data[index] > Max ) { Max = Data[index] ; }
      if( Data[index] < Min ) { Min = Data[index] ; }
    }
   // 8 point minimum RMS // 
    Double_t Mean8point[45]={0};
    Double_t RMS8point[45] ={0};    
    for( int index  = 0; index < 41; index ++){
      for( int subindex = 0 ; subindex < 8; subindex ++){
	Mean8point[index] += Data[index + subindex];
	RMS8point[index]  += TMath::Power( Data[index + subindex] , 2 );
      }
      RMS8point[index] = RMS8point[index]/8.;
      Mean8point[index] /= 8.;
      RMS8point[index] = TMath::Sqrt( RMS8point[index] -Mean8point[index]*Mean8point[index] );
    }
    Double_t minRMS   = 999999999;
    Int_t    minIndex = 0;
    for( int index = 0; index < 41; index++){
      if( minRMS > RMS8point[index] ){ minRMS = RMS8point[index] ; minIndex = index ;}      
    }
    
    for( int index = 0; index < 48; index++){
      gr->SetPoint( index , index, Data[index] - Mean8point[minIndex]);
    }
    for( int index = 0; index < 41; index++){
      grRMS->SetPoint( index , index, RMS8point[index] );
    }
   
    for( int ipoint  = 0; ipoint < 48; ipoint++){
      his->SetBinContent( ipoint , Data[ipoint]-Mean8point[minIndex] );
    }
    
    if( Max - Min < 20 ){continue;}
    int npeaks = spec->Search(his);
    if( npeaks == 0 ){ continue; }
    //TH1* h = spec->Background(his,18, "same");
    can->cd(1);
    his->Draw();
    //h->Draw("same");
    can->Update();
    can->Modified();
    getchar();
  }

  app->Run();
}
