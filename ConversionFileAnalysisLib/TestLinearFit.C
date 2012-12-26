#include <TSpline.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TF1.h>

TSpline3* spl; 
double templeteFunction( double* x, double* par){
  double t = x[0]; 
  double height = par[0]; 
  double mean   = par[1];
  double ped    = par[2]; 
  double t_fcn =0;
  if( t- mean < -100 || t -mean > 200 ) { return ped ; }
  t_fcn = height*spl->Eval(t-mean) +ped;
  return t_fcn;
}
  

void TestLinearFit(){
  
  TFile* tf = new TFile("TEMPLETE_OUT_HEIGHT_0.root");
  TGraphErrors* gre = (TGraphErrors*)tf->Get("Waveform_Height_1200_0");
  gre->Draw();
  spl          = new TSpline3("spl",(TGraph*)gre);
  TF1* waveform = new TF1("waveform",templeteFunction,-150.,200.,3);    
  waveform->SetParameter(0, 1);
  waveform->SetParameter(1, 0);
  waveform->SetParameter(2, 0);
  double hht = waveform->GetX( 0.5, -150, 0); 
  TGraphErrors* grFit = new TGraphErrors();
  for( int i = -120; i< 0 ; i++){
    gre->Fit("pol1","","",i,i+24 );
    TF1* func = gre->GetFunction("pol1");
    double chi = func->GetChisquare();
    double ndf = func->GetNDF();
    grFit->SetPoint(grFit->GetN(), i, chi/ndf);
  }

  gre->Fit("pol1","","",hht-12,hht+12);
  TF1* tfunc = gre->GetFunction("pol1");
  TGraph* grDelta = new TGraph();
  for( int i = 0; i< gre->GetN(); i++){
    if( gre->GetX()[i] < hht-16; || gre->GetX()[i] > hht+16){continue;}
    grDelta->SetPoint( grDelta->GetN(), gre->GetX()[i], gre->GetY()[i]- tfunc->Eval(gre->GetX()[i]));
    
  }
  TCanvas* can = new TCanvas("can","",0,0,800,800);
  can->Divide(2,2);
  can->cd(1);
  grFit->Draw("AP");
  can->cd(2);
  std::cout<< tfunc->GetChisquare()/ tfunc->GetNDF() << std::endl ;
  gre->Draw("AP");
  std::cout<< hht << std::endl;
  can->cd(3);
  grDelta->Draw("AP");
}

