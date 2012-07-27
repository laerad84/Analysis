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
  
  TH1D* histRMS1 = new TH1D("histRMS1","",100,0,100);
  TH1D* histRMS2 = new TH1D("histRMS2","",100,0,4000);
  TH1D* histRMS3 = new TH1D("histRMS3","",100,0,16000);
  TH2D* hisHeightRMS = new TH2D("hisHeightRMS","",48,0,48*8,200,0,0.5);
  TCanvas* can = new TCanvas("can","",1200,800);
  can->Divide(3,2);
  TGraph* gr  = new TGraph();
  TGraph* grHeightRMS = new TGraph();
  TH2D* hisRMSDist8 = new TH2D("hisRMSDist8","",48,0,48,500,0.1,10);
  TH2D* hisRMSDist = new TH2D("hisRMSDist","",48,0,48,500,0.1,100);

  TGraph* grRMS = new TGraph();
  TGraph* grRMS8f = new TGraph();
  TGraph* grRMS8r = new TGraph();
  TGraph* grRatiofr = new TGraph();
  TGraph* grRatiorf = new TGraph();

  gr->SetMarkerStyle(6);
  grRMS->SetMarkerStyle(6);
  grRMS8f->SetMarkerStyle(6);
  grRMS8r->SetMarkerStyle(6);
  grRatiorf->SetMarkerStyle(6);
  grRatiofr->SetMarkerStyle(6);
  //for( int i = 0; i< inputTree->GetEntries(); i++){
  for( int i = 0; i< 10000; i++){
    gr->Set(0);
    grRMS->Set(0);
    grRMS8f->Set(0);
    grRMS8r->Set(0);
    grRatiofr->Set(0);
    grRatiorf->Set(0);

    inputTree->GetEntry(i);
    Double_t RMS_16point[33]={0};
    Double_t RMS_8point[33][2]={{0}};
    Double_t RMS[48]={0};
    Double_t max = 0;
    Double_t min = 18000;
    Double_t Sigsq = 0;
    Double_t PeakTime = 0;
    Double_t Sum=0;
    for( int index = 0; index < 48; index++ ){
      gr->SetPoint( index , Time[index], Data[index] );
      if( Data[index] < min ){ min = Data[index] ;}
      if( Data[index] > max ){ max = Data[index] ; PeakTime = Time[index]; }      
    }
    for( int index = 0; index < 48; index++ ){      
      Sigsq      += ( Data[index] - min )*( Data[index] - min );
      RMS[index]  = ( Data[index] - min )*( Data[index] - min );
      Sum += (Data[index] -min );
    }  
    Sigsq = TMath::Sqrt( Sigsq ) / 48. ;
    for( int index = 0; index < 48-16+1 ; index++ ){
      for( int subindex = 0; subindex < 16; subindex++ ){
	RMS_16point[index] += RMS[ index +subindex ];
	if( subindex < 8 ){
	  RMS_8point[index][0] += RMS[ index + subindex ];
	}else{ 
	  RMS_8point[index][1] += RMS[ index + subindex ];
	}
      }
      RMS_16point[index]   = TMath::Sqrt( RMS_16point[index]   ) / 16. ;
      RMS_8point[index][0] = TMath::Sqrt( RMS_8point[index][0] ) / 8.  ;
      RMS_8point[index][1] = TMath::Sqrt( RMS_8point[index][1] ) / 8.  ;
      grRMS    ->SetPoint( index , ( index + 7.5 ) * 8 , RMS_16point[index]   );
      grRMS8f  ->SetPoint( index , ( index + 1.5 ) * 8 , RMS_8point[index][0] );
      grRMS8r  ->SetPoint( index , ( index + 5.5 ) * 8 , RMS_8point[index][1] );
      grRatiofr->SetPoint( index , ( index + 7.5 ) * 8 , RMS_8point[index][0] / RMS_8point[index][1] );
      grRatiorf->SetPoint( index , ( index + 7.5 ) * 8 , RMS_8point[index][1] / RMS_8point[index][0] );
    }

    
    if( Sum > 1000 ){
      for( int index  = 0; index < 48-16+1; index++){
	hisRMSDist->Fill( index, RMS_16point[index] );
	hisRMSDist8->Fill( index , RMS_8point[index][0]/RMS_8point[index][1]);
      } 
    }
    grHeightRMS->SetPoint( grHeightRMS->GetN(), PeakTime, TMath::Sqrt( Sigsq )/( max-min ));
    hisHeightRMS->Fill(PeakTime,TMath::Sqrt( Sigsq )/( max-min ));
    can->cd(1);
    gr->Draw("AP");
    can->cd(2);
    grRMS->Draw("AP");
    can->cd(3);
    grRMS8f->Draw("AP");
    can->cd(4);
    grRMS8r->Draw("AP");
    can->cd(5);
    grRatiofr->Draw("AP");
    can->cd(6);
    grRatiorf->Draw("AP");

    std::cout << gr->GetRMS(2) << std::endl;
    //gr->GetYaxis()->SetRangeUser(0.1, 16000);
    histRMS1->Fill(gr->GetRMS(2));
    histRMS2->Fill(gr->GetRMS(2));
    histRMS3->Fill(gr->GetRMS(2));
    can->Update();
    can->Modified();
    getchar();         
  }

  TCanvas* can1 = new TCanvas("can1","",800,800);
  can1->Divide(2,2);
  can1->cd(1);
  histRMS1->Draw();
  can1->cd(2);
  //histRMS2->Draw();
  hisRMSDist8->Draw("col");
  can1->cd(3);
  //histRMS3->Draw();
  hisRMSDist->Draw("col");
  can1->cd(4);
  //grHeightRMS->Draw("AP");
  hisHeightRMS->Draw("col");
  app->Run();
}
