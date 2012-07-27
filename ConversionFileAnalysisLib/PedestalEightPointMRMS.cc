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
  //for( int i = 0; i< inputTree->GetEntries(); i++){
  for( int ievent = 0; ievent< 10000; ievent++){
    inputTree->GetEntry(ievent);
    gr->Set(0);

    Double_t RMS8[48-16+1]={0};
    Double_t RMSAfter[48-16+1]={0};
    Double_t Min=20000;
    Double_t Max=0;
    for( int index = 0; index < 48; index++){
      if( Data[index] > Max ) { Max = Data[index] ; }
      if( Data[index] < Min ) { Min = Data[index] ; }
    }
    
    for( int index  = 0 ; index < 33; index++ ){
      for( int subindex = 0; subindex < 16; subindex++){
	RMS8[index] += (Data[ index+subindex ])*(Data[ index+subindex ]);
      }
      RMS8[index] = TMath::Sqrt( RMS8[index] )/8;
    }
    Int_t minimumIndex  = 0;
    Double_t MinimumRMS = 999999;
    for( int index = 0; index < 33 ; index++){
      if( RMS8[index] < MinimumRMS ){ MinimumRMS = RMS8[index] ; minimumIndex = index ; }
    }
    Double_t MeanFirstIt = 0;
    for( int index = minimumIndex; index < minimumIndex +16; index++){
      MeanFirstIt += Data[index];
    }
    MeanFirstIt = MeanFirstIt / 16.;

    for( int index  = 0 ; index < 33; index++ ){
      for( int subindex = 0; subindex < 16; subindex++){
	RMSAfter[index] += (Data[index+subindex]-MeanFirstIt)*(Data[index+subindex] -MeanFirstIt );
      }
      RMSAfter[index] = TMath::Sqrt( RMSAfter[index] )/16;
    }

    minimumIndex  = 0;
    MinimumRMS = 999999;
    for( int index = 0; index < 33 ; index++){
      if( RMSAfter[index] < MinimumRMS ){ MinimumRMS = RMSAfter[index] ; minimumIndex = index ; }
    }
    Double_t MeanSecondIt = 0;
    for( int index = minimumIndex; index < minimumIndex +16; index++){
      MeanSecondIt += Data[index];
    }
    MeanSecondIt = MeanSecondIt/16.;

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
   

    if( Max - Min < 100 ){continue;}
    std::cout<< minIndex <<std::endl;
    can->cd(1);
    gr->Draw("AP");
    gPad->SetGridx();
    gPad->SetGridy();
    can->cd(2);
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    grRMS->Draw("AP");
    can->Update();
    can->Modified();
    getchar();
  }

  app->Run();
}
