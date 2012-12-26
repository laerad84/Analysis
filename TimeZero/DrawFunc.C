#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPad.h"

double fcn( double *x ,double *par ){
  double x0     = x[0];
  double L      = 2622;
  double l      = 100;
  double deltaT = (l*TMath::Sqrt(L*L+x0*x0) - 3*l*L)/300./TMath::Sqrt(L*L+x0*x0);
  //double deltaT = (l*TMath::Sqrt(L*L+x0*x0) - 3*l*L)/300./TMath::Sqrt(L*L+x0*x0);
  return deltaT;
}

double func0( double *x ,double *par ){

  double x0     = x[0];
  double L      = 2622;
  double l      = 100;
  double deltaT = (L*L +x0*x0 + l*TMath::Sqrt(L*L+x0*x0) - 3*l*L)/300./TMath::Sqrt(L*L+x0*x0);
  return deltaT;
}

void DrawFunc(){

  TF1* fcnt   = new TF1("fcn",fcn, 0,1000,0);
  TF1* func  = new TF1("func",func0, 0, 1000,0);
  
  std::cout <<func->Eval( 0) << std::endl;
  TCanvas* can  =new TCanvas("can","",800,800);
  can->Divide(2,2);
  can->cd(1);
  fcnt->Draw();  
  can->cd(2);
  func->Draw()
;
}
