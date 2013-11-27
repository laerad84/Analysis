#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
//#include "myplot/Plotter.h"
#include "math.h"
#include <iostream>
#include <string>
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
/*
double asymGausExpo(double *x,double *par){
  double X = -(x[0]-1);
  double height = par[0];
  double peak = -(par[1]-1);
  double sigma = fabs(par[2]);
  double sigma2=fabs(par[3]);
  double alpha = fabs(par[4]);
  X-=peak;
  //  double offset = fabs(par[6]);
  double f;
  if(X<0){
    X/=sigma;
    f = height*exp(-X*X/2);
  }else{
    X/=sigma2;
    if(X<alpha){
      f = height*exp(-X*X/2);
    }else{
      f = exp(-alpha*(X-alpha))*height*exp(-alpha*alpha/2);
    }
  }
  // return f+offset;
  return f;
}
*/

int main(){
  //TCanvas *canvas = Plotter::getPlotter()->canvas();
  TCanvas* canvas = new TCanvas("can","can",800,600);
  std::string psfile = "plot.ps";
  canvas->Print((psfile+"[").c_str());
  TChain chain("tro");
  chain.Add("data/run1*");
  chain.Add("data/run2*");
  std::cout<<"# of entry:"<<chain.GetEntries()<<std::endl;
  TFile* tfout = new TFile("tfile.root","recreate");
  
  int const ncal = 10;
  double cal[ncal] = {0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05};
  double x[12];
  for(int i=0;i<12;i++) x[i] = 100*(i+1);
  //  TF1 func("myfunc",asymGausExpo,0.5,1.3,5);
  double rms[ncal][12],sig[ncal][12];
  TGraph *gra[ncal]={0};
  for(int ical=0;ical<ncal;ical++){
    std::cout<<"ical:"<<ical<<std::endl;
    chain.Draw(Form("eclus[%d]/etot:genE>>h(12,50,1250,100,0.5,1.3)",ical),"genTheta<25","colz");
    canvas->Print(psfile.c_str());
    TH2F* h2d = (TH2F*)gDirectory->Get("h");
    for(int ibin=0;ibin<12;ibin++){
      TH1D* h1d = h2d->ProjectionY("hprj",ibin+1,ibin+1);
      rms[ical][ibin] = h1d->GetRMS();
      //    h1d->Fit(&func);
      //      sig[ical][ibin] = sqrt(func.Variables(0.5,1.3));
    }
    gra[ical] = new TGraph(12,x,rms[ical]);
    //    graS[ical] = new TGraph(12,x,sig[ical]);
    gra[ical]->SetMarkerColor(ical%5+1);
    //    graS[ical]->SetMarkerColor(ical/6+1);
    gra[ical]->SetLineColor(ical%5+1);
    //    graS[ical]->SetLineColor(ical/6+1);
    gra[ical]->SetNameTitle(Form("gr_%d",ical),Form("gr_%d",ical));
  }
  
  {
    TH2D frame("frame","",100,0,1300,100,0,0.1);
    frame.SetStats(0);
    frame.Draw();
    for(int i=0;i<5;i++) gra[i]->Draw("*l");
    canvas->Print(psfile.c_str());
    frame.Draw();
    for(int i=5;i<ncal;i++) gra[i]->Draw("*l");
    canvas->Print(psfile.c_str());
  }
  double param[4][ncal];
  double dparam[4][ncal];
  TF1 resoFunc("resof","sqrt(pow([0],2)+pow([1]/sqrt(x/1000.),2)+pow([2]/(x/1000.),2))",100,1200);
  resoFunc.SetParameters(0.01,0.02,0.005);
  for(int i=0;i<ncal;i++){
    gra[i]->Fit(&resoFunc);
    for(int j=0;j<3;j++){
      param[j][i] = fabs(resoFunc.GetParameter(j));
      if(i>0){
	dparam[j][i] = (param[j][i]>param[j][0])
	  ? sqrt(pow(param[j][i],2)-pow(param[j][0],2))
	  : -sqrt(pow(param[j][0],2)-pow(param[j][i],2));
      }else{
	dparam[j][i] = 0;
      }
    }
    param[3][i] = sqrt(pow(param[0][i],2)+pow(param[1][i],2)+pow(param[2][i],2));
    dparam[3][i] = (param[3][i]>param[3][0])
      ? sqrt(pow(param[3][i],2)-pow(param[3][0],2))
      : -sqrt(pow(param[3][0],2)-pow(param[3][i],2));

  }

  {
  TH2D frame("frame","",100,0,0.06,100,0,0.1);
  frame.SetStats(0);
  frame.Draw();
  TGraph *resoGra[4];
  for(int j=0;j<4;j++){
    resoGra[j] = new TGraph(ncal,cal,param[j]);
    resoGra[j]->SetLineColor(j+1);
    resoGra[j]->SetMarkerColor(j+1);
    resoGra[j]->Draw("*l");
    resoGra[j]->SetNameTitle(Form("gr_resolution_%d",j),Form("gr_resolution_%d"));
    resoGra[j]->Write();
  }

  canvas->Print(psfile.c_str());
  }

  {
    TH2D frame("frame","",100,0,0.06,100,-0.035,0.065);
    frame.SetStats(0);
    frame.Draw();
    TGraph *resoGra[4];
    for(int j=0;j<4;j++){
      //      delete resoGra[j];
      resoGra[j] = new TGraph(ncal,cal,dparam[j]);
      resoGra[j]->SetLineColor(j+1);
      resoGra[j]->SetMarkerColor(j+1);
      resoGra[j]->Draw("*l");
    }
    canvas->Print(psfile.c_str());
  }
  

  canvas->Print((psfile+"]").c_str());
  for( int i = 0; i< ncal; i++){
    gra[i]->Write();
  }

  tfout->Close();

}
