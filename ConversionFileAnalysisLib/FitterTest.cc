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
#include "TGraphErrors.h"
#include "TSpline.h"

#include "E14WaveFitterMinimum.h"

int main( int argc, char** argv ){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  std::string InputFilename = argv[1];
  TApplication* app = new TApplication( "app", &argc, argv );
  TFile* tempFile = new TFile("Template_Cosmic.root");
  TH2D* hisTemp[20];
  TProfile* prof[20];
  TSpline3*  spls[20];
  TGraph* grErr[20];

  for( int i = 0; i< 20; i++){
    hisTemp[i] = (TH2D*)tempFile->Get(Form("hisTemp%d",i));
    prof[i] = hisTemp[i]->ProfileX();
    grErr[i] = new TGraph();
    for( int index = 0; index< prof[i]->GetNbinsX(); index++){
      if(prof[i]->GetBinEntries(index) < 10 ){continue; }
      grErr[i]->SetPoint( grErr[i]->GetN() , prof[i]->GetBinCenter( index ) ,prof[i]->GetBinContent( index ));
      //grErr[i]->SetPointError( index , 0.5, prof[i]->GetBinError( index ));
    }
    spls[i]    = new TSpline3(Form("spls%d",i), grErr[i] );

  }
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

  TGraph* gr     = new TGraph();
  TGraph* grFront = new TGraph();
  TGraph* grBack  = new TGraph();

  TGraph* grTemp      = new TGraph();
  TGraph* grTempFront = new TGraph();
  TGraph* grTempBack  = new TGraph();
    
  gr         ->SetMarkerStyle(6);
  grFront    ->SetMarkerStyle(6);
  grBack     ->SetMarkerStyle(6);
  grTemp     ->SetMarkerStyle(6);
  grTempFront->SetMarkerStyle(6);
  grTempBack ->SetMarkerStyle(6);
  TCanvas* can = new TCanvas("can","",0,0,800,1140);
  can->Divide(2,3);
  WaveformFitter* wav = new WaveformFitter( 48, kFALSE );

  TH1D* his           = new TH1D("hisChsq"      ,"", 200,0,20);
  TH1D* hisFitResult  = new TH1D("hisFitResult" ,"",1000,0,10);
  TH1D* hisNormDist   = new TH1D("hisNormDist"  ,"",1000,0,100);
  TH1D* hisSqrtDist   = new TH1D("hisSqrtDist"  ,"",1000,0,100);
  TH1D* hisAspectDist = new TH1D("hisAspectDist","",1000,0,2);
  TH1D* hisNormDistCUT   = new TH1D("hisNormDistCUT"  ,"",1000,0,100);
  TH1D* hisSqrtDistCUT   = new TH1D("hisSqrtDistCUT"  ,"",1000,0,100);
  TH1D* hisAspectDistCUT = new TH1D("hisAspectDistCUT","",1000,0,2);
  TH1D* hisCorrelation   = new TH1D("hisCorrelation"  ,"",4000,-1,1);


  E14WaveFitterMinimum* fitter = new E14WaveFitterMinimum();
  
  for( int ievent  = 0; ievent < inputTree->GetEntries(); ievent++){    
    inputTree->GetEntry(ievent);
    if( ID ==7 ){continue;}
    gr->Set(0);
    grFront->Set(0);
    grBack->Set(0);    
    grTemp->Set(0);
    grTempFront->Set(0);
    grTempBack->Set(0);

    Double_t Min=20000;
    Double_t Max=0;
    Int_t    MinIndex = 0; 
    Int_t    MaxIndex = 0;
    for( int index = 0; index < 48; index++){
      if( Data[index] > Max ) { Max = Data[index] ; MaxIndex = index;}
      if( Data[index] < Min ) { Min = Data[index] ; MinIndex = index;}
    }
    
    if( Max - Min < 100 ){continue; }
    Double_t gnd = 0;
    for( int ipoint  = 0; ipoint < 48; ipoint++){
      gr->SetPoint( ipoint , ipoint*8 , Data[ipoint] );
      grFront->SetPoint( ipoint , ipoint*8 , Data[ipoint] );
      grBack->SetPoint( ipoint , ipoint*8 , Data[ipoint] );
      if( ipoint < 8){ gnd += Data[ipoint]; }
    }
    gnd = gnd / 8.;
    // Test Ground Search // 
    
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


    //std::cout <<ID << " : " << gnd << " : " << currentGND << " : " << currentPosition << std::endl;
    E14WaveFitterMinimum::m_spl = spls[ID];


    fitter->m_FitFunc->SetParLimits( 2 , currentGND, currentGND); 
    fitter->m_FitFunc->SetParLimits( 1 , MaxIndex*8-8, (MaxIndex)*8 + 8);
    fitter->m_FitFunc->SetParLimits( 0 , (Max-Min)*0.9, (Max-Min)*5);
    fitter->m_FitFunc->SetParameter( 0 , Max-Min);
    fitter->m_FitFunc->SetParameter( 1 , MaxIndex*8 );
    //fitter->m_FitFunc->SetParameter( 2 , (getGnd[0] + getGnd[1])/2);
    fitter->m_FitFunc->SetParameter( 2 , currentGND);
    fitter->Fit( gr, "","",MaxIndex*8-150,MaxIndex*8+75);


    can->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    gr->Draw("AP");        
    can->cd(2);
    gPad->SetGridx();/////////////////////////////////////////////////////////
    gPad->SetGridy();
    grTemp->Draw("AP");
    //spls[ID]->Draw();
    can->cd(3);
    grFront->Draw("AP");
    can->cd(4);
    grTempFront->Draw("AP");
    can->cd(5);
    grBack->Draw("AP");
    can->cd(6);
    grTempBack->Draw("AP");
    can->Update();
    can->Modified();
    getchar();
  }
  app->Run();
}
  
