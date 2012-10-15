#include <cstdlib>
#include <string>
#include "TSpline.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TTree.h"

TSpline3* spl;
Double_t splFunc( double* x, double* par ){
  double t=x[0];
  double height = par[0];
  double mean   = par[1];
  double ped    = par[2];

  double t_fcn  = 0;
  if( t - mean  <= -150 || t- mean >= 300 ){ return ped;}
  //if( t -mean < -80 || t-mean >= 90 ){ return ped ;}
  if( t == 1./0 ){ return ped; }
  t_fcn         = height*spl->Eval(t - mean) + ped;
  return t_fcn;
}

Double_t splDoubleFunc( double* x, double* par ){
  double t =x[0];
  double height1 = par[0];
  double mean1   = par[1];
  double ped     = par[2];
  double height2 = par[3];
  double mean2   = par[4];
  double t_fcn = 0;
  
  if( mean1 < mean2 ){
    if( t - mean1  < -150 ||
	t - mean2 >=  300 ){ 
      return ped;
    }
  }else{
    if( t - mean2 < -150 || t-mean1 >=300){ return ped; }
  }
  if( t == 1./0){ return ped ;}
  t_fcn = height1*spl->Eval(t-mean1)+height2*spl->Eval(t-mean2)+ped;
  return t_fcn;
}



void Test2Peak(){

  TF1* func = new TF1("func",splFunc,0.,500.,3);
  TF1* func1 = new TF1("func1",splFunc,0.,500.,3);
  TF1* func2 = new TF1("func2",splDoubleFunc,0.,500.,5);
  TF1* func3 = new TF1("func3",splFunc,0.,500.,3);
  
  TFile* tfWaveform = new TFile("TEMPLATE_COMPLETE_GRAPH_3pi0RunList_200_400.root");
  
  std::string Filename   = Form("/Volume0/ExpData/2012_Feb_Beam/RootFile_wav/TEMPLATE_FIT_RESULT_DUMP_%d.root",4202);
  
  TFile* tf = new TFile(Filename.c_str());
  TTree* tr = (TTree*)tf->Get("Waveform");
  
  Int_t    EventNumber;
  Int_t    ModuleNumber;
  Double_t Waveform[48];
  Short_t  TimeInfo[48];
  /*
  Double_t ChisqNDF;
  Double_t PeakTime;
  Double_t HHTime;
  Double_t Height;
  Double_t Pedestal;
  */
  {
    tr->SetBranchAddress("EventNumber",&EventNumber);
    tr->SetBranchAddress("ModuleNumber",&ModuleNumber);
    tr->SetBranchAddress("Waveform",Waveform);
    tr->SetBranchAddress("TimeInfo",TimeInfo);
    /*
    tr->SetBranchAddress("PeakTime",&PeakTime);
    tr->SetBranchAddress("HHTime",&HHTime);
    tr->SetBranchAddress("Height",&Height);
    tr->SetBranchAddress("Pedestal",&Pedestal);
    */
  }

  //for( int ievent = 0; ievent < tr->GetEntries(); ievent++){
  //tr->GetEntry(ievent);
  TGraph* grTemplate = (TGraph*) tfWaveform->Get("Template_1090");
  spl = new TSpline3("spl",grTemplate);
  tr->GetEntry ( 2060 );
  
  TGraph* grWave = new TGraph();
  TGraph* grWaveDelta = new TGraph();
  TGraph* grWaveDeltaDelta = new TGraph();
  grWave->SetMarkerStyle(6);
  grWaveDelta->SetMarkerStyle(6);
  grWaveDeltaDelta->SetMarkerStyle(6);
  for( int ipoint = 0; ipoint < 48; ipoint++ ){
    grWave->SetPoint( ipoint , ipoint*8, Waveform[ipoint]);
  }

  Double_t Ped     = 510;
  Double_t Mean1   = 150;
  Double_t Height1 = -1*Ped+3600;
  Double_t Mean2   = 200;
  Double_t Height2 = 120;

  Double_t FitRegion[4][2]={{0}};

  FitRegion[0][0] = Mean1-100;
  FitRegion[0][1] = Mean1+8;
  FitRegion[1][0] = Mean1-8;
  FitRegion[1][1] = Mean1+150;

  func->SetParameter(2, Ped);
  func->SetParameter(1, Mean1);
  func->SetParameter(0, Height1);
  func->SetParLimits(0, Height1*0.9,Height1*1.1);
  func->SetParLimits(1, Mean1-16, Mean1+16);
  //  func->SetParLimits(2, Ped-4, Ped+4);
  func->SetParLimits( 2, Ped, Ped);
  func1->SetParameter(2, Ped);
  func1->SetParameter(1, Mean1);
  func1->SetParameter(0, Height1);
  func1->SetParLimits(0, Height1*0.9,Height1*1.1);
  func1->SetParLimits(1, Mean1-16, Mean1+16);
  //func1->SetParLimits(2, Ped-4, Ped+4);
  func1->SetParLimits(2, Ped,Ped);


  grWave->Fit(func ,""  ,"",FitRegion[0][0],FitRegion[0][1]);
  grWave->Fit(func1,"R+","",FitRegion[1][0],FitRegion[1][1]);

  Int_t frIndex = 0;
  if( func->GetChisquare()/func->GetNDF() < func1->GetChisquare()/func1->GetNDF()){
    frIndex = 1; 
    std::cout<< "/////////////////////////////////////\n";
    std::cout<< "Fit With Front" << std::endl;
    std::cout<< "/////////////////////////////////////\n";
    for( int i = 0; i<grWave->GetN(); i++){
      Double_t Delta  = grWave->GetY()[i]- func->Eval( grWave->GetX()[i] );
      grWaveDelta->SetPoint(i, grWave->GetX()[i], Delta );
    }
  }else{
    frIndex = 2;
    std::cout<< "/////////////////////////////////////\n";
    std::cout<< "Fit With Rear " << std::endl;
    std::cout<< "/////////////////////////////////////\n";
    for( int i = 0; i<grWave->GetN(); i++){
      Double_t Delta  = grWave->GetY()[i]- func1->Eval( grWave->GetX()[i] );
      grWaveDelta->SetPoint(i, grWave->GetX()[i], Delta );
    }
  }

  func3->SetParameter(2, 0);
  func3->SetParameter(1, Mean2);
  func3->SetParameter(0, Height2);
  func3->SetParLimits(0, Height2*0.9,Height2*1.1);
  func3->SetParLimits(1, Mean2 - 16, Mean2+ 16 );
  func3->SetParLimits(2, -0.1,0.1);
  FitRegion[2][0] = Mean2 - 50;
  FitRegion[2][1] = Mean2 + 50;
  grWaveDelta->Fit( func3, "","",FitRegion[2][0],FitRegion[2][1]);
  if( frIndex  == 1){
    std::cout<< "FRONT" << std::endl;
    func2->SetParameter(0, func->GetParameter(0));
    func2->SetParameter(1, func->GetParameter(1));
    func2->SetParLimits(0, func->GetParameter(0)*0.8,func->GetParameter(0)*1.2);
    func2->SetParLimits(1, func->GetParameter(1)-16,func->GetParameter(1)+16);
  }else{
    std::cout<< "REAR" << std::endl;
    func2->SetParameter(0, func1->GetParameter(0));
    func2->SetParameter(1, func1->GetParameter(1));
    func2->SetParLimits(0, func1->GetParameter(0)*0.8,func->GetParameter(0)*1.2);
    func2->SetParLimits(1, func1->GetParameter(1)-16,func->GetParameter(1)+16);
  }
  func2->SetParameter(3, func3->GetParameter(0));
  func2->SetParameter(4, func3->GetParameter(1));
  func2->SetParLimits(3, func3->GetParameter(0)*0.8, func3->GetParameter(0)*1.2);
  func2->SetParLimits(4, func3->GetParameter(1)-16,func3->GetParameter(1)+16);
  //func2->SetParameter(2, func->GetParameter(2));
  //func2->SetParLimits(2, func->GetParameter(2)-4,func->GetParameter(2)+4);
  func2->SetParameter(2, Ped);
  func2->SetParLimits(2, Ped,Ped);
  if( frIndex  ==1 ){
    if( func->GetParameter(1) < func3->GetParameter(1) ){
      FitRegion[3][0] = func->GetParameter(1)  - 100;
      FitRegion[3][1] = func3->GetParameter(1) + 80;
    }else{
      FitRegion[3][0] = func3->GetParameter(1) - 100;
      FitRegion[3][1] = func->GetParameter(1)  + 80;
    }
  }else{
    if( func1->GetParameter(1) < func3->GetParameter(1) ){
      FitRegion[3][0] = func1->GetParameter(1) - 100;
      FitRegion[3][1] = func3->GetParameter(1) + 80;
    }else{
      FitRegion[3][0] = func3->GetParameter(1) - 100;
      FitRegion[3][1] = func1->GetParameter(1) + 80;
    }
  }
  grWave->Fit(func2, "R+","",FitRegion[3][0],FitRegion[3][1]);
  //grWave->Fit(func2, "R+","",0,250);

  for( int i = 0; i<grWave->GetN(); i++){
    Double_t DeltaDelta  = grWave->GetY()[i]- func2->Eval( grWaveDelta->GetX()[i] );
    grWaveDeltaDelta->SetPoint(i, grWaveDelta->GetX()[i], DeltaDelta );
  }

  
  TCanvas* can = new TCanvas("can1","",800,800);
  can->Divide(2,2);
  can->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  grWave->Draw("AP");  
  can->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  grWaveDelta->Draw("AP");
  can->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  grWave->Draw("AP");
  func2->Draw("same");
  can->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  grWaveDeltaDelta->Draw("AP");
    //}
}
  
  
