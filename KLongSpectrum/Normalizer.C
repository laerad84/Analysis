#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>
#include "TMath.h"
#include "TVector2.h"

const double KLMass = 497.648;//MeV                                                                                

double KLSpectrum(double* x,double*p){
  //double klp   = TMath::Sqrt( x[0]*x[0]-KLMass*KLMass);                                                          
  double klp = x[0];
  double sigma = p[2]*(1-(p[3]+p[4]*klp)*(klp-p[1]));
  double value = p[0]*TMath::Exp(-0.5*(klp-p[1])*(klp-p[1])/(sigma*sigma));
  return value;
}

void Normalizer(){

  TF1* soltFunc  = new TF1("soltfunc",KLSpectrum,0,12,5);
  TF1* sugarFunc = new TF1("sugarfunc",KLSpectrum,0,12,5);
  const double soltPar[5] = { 1, 1.37,7.48e-01,-2.9e-01,1.68e-04};
  const double sugarPar[5] = {1,1.41991,0.810237,-0.301413,0.017092};
  soltFunc->SetParameters(soltPar);
  sugarFunc->SetParameters(sugarPar);

  TFile* tf = new TFile("DistributionTest.root");
  TH1D* his[2];
  his[0] = (TH1D*)tf->Get("hisKLP_0");
  his[1] = (TH1D*)tf->Get("hisKLP_1");

  TH1D* hisNorm = new TH1D("hisKLEN","hisKLEN",100,0,10000);
  TH1D* hisAspect = new TH1D("hisFactor","hisFactor",100,0,10000);
  TH1D* hisRatio  = new TH1D("hisRatio","hisRatio",100,0,10000);


  for( int ibin = 1; ibin <= his[0]->GetNbinsX(); ibin++){
    double KLP    = his[0]->GetBinCenter(ibin)/1000;
    double weight = sugarFunc->Eval(KLP)/soltFunc->Eval(KLP);
    hisNorm->SetBinContent( ibin, his[0]->GetBinContent(ibin)*weight);
    hisAspect->SetBinContent( ibin, weight);
  }
  double ScaleFactor=  his[1]->Integral()/his[0]->Integral();
  double ScaleFactor_N = his[1]->Integral()/hisNorm->Integral();
  hisNorm->Scale(his[1]->Integral()/hisNorm->Integral());

  for( int ibin  =1; ibin < his[0]->GetNbinsX();ibin++){
    double data = his[1]->GetBinContent(ibin);
    double sim  = hisNorm->GetBinContent(ibin);
    double edata = his[1]->GetBinError( ibin );
    double esim  = hisNorm->GetBinError( ibin );
    double error = TMath::Sqrt( esim*esim +esim*esim );
    if( sim != 0 ){
      hisRatio->SetBinContent(ibin,data/sim);
      hisRatio->SetBinError(ibin,error/sim);
    }
  }


  TCanvas * can = new TCanvas("can","can",800,800);
  TPad* pad0 = new TPad("pad0","pad0",0,0.5,1,1);
  TPad* pad1 = new TPad("pad1","pad1",0,0.25,1,0.5);
  TPad* pad2 = new TPad("pad2","pad2",0,0,1,0.25);  
  can->Draw();
  pad0->Draw();
  pad1->Draw();
  pad2->Draw();
  pad0->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  his[0]->SetLineWidth(2);
  his[1]->SetLineWidth(2);
  hisNorm->SetLineWidth(2);
  hisNorm->SetLineColor(4);
  his[1]->SetMarkerStyle(21);
  his[1]->Draw("EP");
  his[0]->Scale(ScaleFactor);
  his[0]->Draw("same");
  hisNorm->Draw("same");
  pad1->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  hisAspect->SetMarkerStyle(21);
  hisAspect->Draw("ep");
  pad2->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  hisRatio->SetMarkerStyle(21);
  hisRatio->Draw("ep");

  /*
  hisNorm->Scale( his[1]->Integral()/hisNorm->Integral());
  hisNorm->SetLineColor(4);
  hisNorm->Draw("same");
  */
}
