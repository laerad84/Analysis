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
  if( t - mean1 < -200 || t-mean1 >=300){ return ped; }
  //if( t == 1./0){ return ped ;}
  t_fcn = height1*spl->Eval(t-mean1)+height2*spl->Eval(t-mean2)+ped;
  return t_fcn;
}



void Test2Peak(){

  TF1* func = new TF1("func",splFunc,0.,500.,3);
  TF1* func1 = new TF1("func1",splFunc,0.,500.,3);
  TF1* func2 = new TF1("func2",splDoubleFunc,0.,500.,5);
  TF1* func3 = new TF1("func3",splFunc,0.,500.,3);
  
  TFile* tfWaveform = new TFile("TEMPLATE_COMPLETE_GRAPH_3pi0RunList_200_400.root");
  TGraph* grTemplate = (TGraph*) tfWaveform->Get("Template_1220");
  spl = new TSpline3("spl",grTemplate);
  
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
  tr->GetEntry ( 2097965 );
  
  TGraph* grWave = new TGraph();
  TGraph* grWaveDelta = new TGraph();
  TGraph* grWaveDeltaDelta = new TGraph();
  grWave->SetMarkerStyle(6);
  grWaveDelta->SetMarkerStyle(6);
  grWaveDeltaDelta->SetMarkerStyle(6);
  for( int ipoint = 0; ipoint < 48; ipoint++ ){
    grWave->SetPoint( ipoint , ipoint*8, Waveform[ipoint]);
  }
  func->SetParameter(2, 340);
  func->SetParameter(1, 95);
  func->SetParameter(0, 800);
  func->SetParLimits(0, 750,850);
  func->SetParLimits(1, 90,100);
  func->SetParLimits(2, 335,345);

  func1->SetParameter(2, 340);
  func1->SetParameter(1, 95);
  func1->SetParameter(0, 800);
  func1->SetParLimits(0, 750,850);
  func1->SetParLimits(1, 90,100);
  func1->SetParLimits(2, 335,345);
  
  TCanvas* can = new TCanvas("can1","",800,800);
  can->Divide(2,2);
  can->cd(1);
  grWave->Draw("AP");
  grWave->Fit(func,"","",0,100);
  grWave->Fit(func1,"R+","",95,300);
  for( int i = 0; i<grWave->GetN(); i++){
    Double_t Delta  = grWave->GetY()[i]- func->Eval( grWave->GetX()[i] );
    grWaveDelta->SetPoint(i, grWave->GetX()[i], Delta );
  }
  can->cd(2);
  func3->SetParameter(2, 0);
  func3->SetParameter(1, 150);
  func3->SetParameter(0, 350);
  func3->SetParLimits(0, 300,375);
  func3->SetParLimits(1, 130,170);
  func3->SetParLimits(2, -4,4);
  grWaveDelta->Fit( func3, "","",0,300);
  grWaveDelta->Draw("AP");

  can->cd(3);

  func2->SetParameter(0, func->GetParameter(0));
  func2->SetParameter(1, func->GetParameter(1));
  func2->SetParameter(2, func->GetParameter(2));
  func2->SetParameter(3, func3->GetParameter(0));
  func2->SetParameter(4, func3->GetParameter(1));
  func2->SetParLimits(0, func->GetParameter(0)*0.95,func->GetParameter(0)*1.05);
  func2->SetParLimits(1, func->GetParameter(1)-4,func->GetParameter(1)+4);
  func2->SetParLimits(2, func->GetParameter(2)-4,func->GetParameter(2)+4);
  func2->SetParLimits(3, func3->GetParameter(0)*0.9, func3->GetParameter(0)*1.1);
  func2->SetParLimits(4, func3->GetParameter(1)-4,func3->GetParameter(1)+4);

  grWave->Fit(func2, "R+","",0,250);
  grWave->Draw("AP");
  func2->Draw("same");
  for( int i = 0; i<grWave->GetN(); i++){
    Double_t DeltaDelta  = grWave->GetY()[i]- func2->Eval( grWaveDelta->GetX()[i] );
    grWaveDeltaDelta->SetPoint(i, grWaveDelta->GetX()[i], DeltaDelta );
  }
  can->cd(4);
  grWaveDeltaDelta->Draw("AP");
    //}
}
  
  
