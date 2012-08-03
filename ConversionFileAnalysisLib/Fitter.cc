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
TProfile* FitProfile;
TSpline*  spl;

double templateFunction( double* x , double* par ){
  double t=x[0];
  double height = par[0];
  double mean   = par[1];
  double ped    = par[2];
  /*
  double b      = par[2];
  double a      = par[3];
  */
  double t_fcn  = 0;
  //std::cout<< t << " : " << spl->Eval(t+mean) << std::endl;
  if( t - mean  < -100 || t- mean >150 ){ return ped;}
  t_fcn         = height*spl->Eval(t - mean) + ped;
  return t_fcn;
}



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

  TF1* func      = new TF1("func"     ,templateFunction,0,500,3);
  TF1* funcFront = new TF1("funcFront",templateFunction, 0, 500, 3);
  TF1* funcBack  = new TF1("funcBack" ,templateFunction, 0, 500, 3);
  func     ->SetLineWidth(1);
  funcFront->SetLineWidth(1);
  funcBack ->SetLineWidth(1);
  funcFront->SetLineColor(2);
  funcBack ->SetLineColor(3);

  TH1D* his           = new TH1D("hisChsq"      ,"", 200,0,20);
  TH1D* hisFitResult  = new TH1D("hisFitResult" ,"",1000,0,10);
  TH1D* hisNormDist   = new TH1D("hisNormDist"  ,"",1000,0,100);
  TH1D* hisSqrtDist   = new TH1D("hisSqrtDist"  ,"",1000,0,100);
  TH1D* hisAspectDist = new TH1D("hisAspectDist","",1000,0,2);
  TH1D* hisNormDistCUT   = new TH1D("hisNormDistCUT"  ,"",1000,0,100);
  TH1D* hisSqrtDistCUT   = new TH1D("hisSqrtDistCUT"  ,"",1000,0,100);
  TH1D* hisAspectDistCUT = new TH1D("hisAspectDistCUT","",1000,0,2);
  TH1D* hisCorrelation   = new TH1D("hisCorrelation"  ,"",4000,-1,1);
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
    spl = spls[ID];

    //funcFront->SetParLimits( 2 , gnd-500, gnd+500);
    funcFront->SetParLimits( 2 , currentGND,currentGND);
    funcFront->SetParLimits( 1 , MaxIndex*8-8, (MaxIndex)*8 + 8);
    funcFront->SetParLimits( 0 , (Max-Min)*0.9, (Max-Min)*5);
    funcFront->SetParameter( 0 , Max-Min);
    funcFront->SetParameter( 1 , MaxIndex*8 );
    funcFront->SetParameter( 2 , currentGND);
    
    funcBack->SetParLimits( 2 , currentGND, currentGND);
    funcBack->SetParLimits( 1 , MaxIndex*8-8, (MaxIndex)*8 + 8);
    funcBack->SetParLimits( 0 , (Max-Min)*0.9, (Max-Min)*5);
    funcBack->SetParameter( 0 , Max-Min);
    funcBack->SetParameter( 1 , MaxIndex*8 );
    funcBack->SetParameter( 2 , gnd);
    
    grFront->Fit(funcFront,"Q","",MaxIndex*8-150,MaxIndex*8-50);
    grBack ->Fit(funcBack ,"Q","",MaxIndex*8+50  ,MaxIndex*8+150);

    Double_t getGnd[2];
    getGnd[0] = funcFront->GetParameter( 2 );
    getGnd[1] = funcBack->GetParameter( 2 );

    
    //    func->SetParLimits( 2 , ((getGnd[0]+getGnd[1])-TMath::Abs(getGnd[0]-getGnd[1]))/2,
    //			((getGnd[0]+getGnd[1])-TMath::Abs(getGnd[0]-getGnd[1]))/2 );

    func->SetParLimits( 2 , currentGND, currentGND); 
    func->SetParLimits( 1 , MaxIndex*8-8, (MaxIndex)*8 + 8);
    func->SetParLimits( 0 , (Max-Min)*0.9, (Max-Min)*5);
    func->SetParameter( 0 , Max-Min);
    func->SetParameter( 1 , MaxIndex*8 );
    //func->SetParameter( 2 , (getGnd[0] + getGnd[1])/2);
    func->SetParameter( 2 , currentGND);
    
    gr     ->Fit(func     ,"Q","",MaxIndex*8-150,MaxIndex*8+75);

    Int_t  nPoint              = 0;
    double NormDist            = 0.;
    double SqSumDist           = 0.;
    double aspectHeight        = 0.;
    double NormalizedSqSumDist = 0.;
    double FitAreaWidth       = 0.;
    double SigAreaWidth       = 0.;
    double Correlation        = 0.;
    TGraph* grCorr = new TGraph(); 
    for( int index = 0; index < 48; index++){
      
      grTemp->SetPoint( index, index*8,func->Eval(index*8)-Data[index]);
      grTempFront->SetPoint( index, index*8,funcFront->Eval(index*8)-Data[index]);
      grTempBack->SetPoint( index, index*8,funcBack->Eval(index*8)-Data[index]);

      if( index*8 > MaxIndex*8-100 && index*8 < MaxIndex*8 +100 ){

	NormDist += abs( func->Eval(index*8) - Data[index] );
	SqSumDist+= pow( func->Eval(index*8) - Data[index] ,2 );
	FitAreaWidth  += func->Eval(index*8) - func->GetParameter(2);
	SigAreaWidth  += Data[index]         - func->GetParameter(2);
	nPoint++;
	grCorr->SetPoint( grCorr->GetN(), 
			  func->Eval(index*8) - func->GetParameter(2),
			  Data[index]         - func->GetParameter(2)); 
      }
    }

    NormDist = NormDist         / nPoint; 
    SqSumDist= sqrt(SqSumDist)  / nPoint;
    FitAreaWidth = FitAreaWidth / nPoint;
    SigAreaWidth = SigAreaWidth / nPoint;
    hisNormDist->Fill(NormDist);
    hisSqrtDist->Fill(SqSumDist);
    hisAspectDist->Fill(FitAreaWidth/SigAreaWidth);
    hisCorrelation->Fill(grCorr->GetCorrelationFactor());
    if( grCorr->GetCorrelationFactor() > 0.998){
      hisNormDistCUT->Fill(NormDist);
      hisSqrtDistCUT->Fill(SqSumDist);
      hisAspectDistCUT->Fill(FitAreaWidth/SigAreaWidth);
    }
    /*
    if( abs(FitAreaWidth/SigAreaWidth -1 ) < 0.025){
      if( SqSumDist < 5 ){
	hisNormDistCUT->Fill(NormDist);
      }

      hisSqrtDistCUT->Fill(SqSumDist);
    }
    
    if( SqSumDist < 5 ){
      hisAspectDistCUT->Fill(FitAreaWidth/SigAreaWidth);
    }
    */
    if( grCorr->GetCorrelationFactor() <0.998 ){ 
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
    
  }

  can->cd(1);
  hisNormDist->Draw();
  hisNormDistCUT->SetLineColor(2);
  hisNormDistCUT->Draw("same");
  can->cd(3);
  hisSqrtDist->Draw();
  hisSqrtDistCUT->SetLineColor(2);
  hisSqrtDistCUT->Draw("same");
  can->cd(5);
  hisAspectDist->Draw();
  hisAspectDistCUT->SetLineColor(2);
  hisAspectDistCUT->Draw("same");
  can->cd(6);
  hisCorrelation->Draw();
  can->Update();
  can->Modified();
  
  app->Run();
}
  
