#include "GeneralFunctions.h"
#include <iostream>

double* GenLogArray( int n, double min, double max ){
  if( min < 0 || max < 0 ){
    std::cerr << __LINE__ << " " << __func__ << " : Value is negative" << std::endl;
    return NULL;
  }
  double *x  = new double[n+1];  
  double minlog = TMath::Log(min);
  double maxlog = TMath::Log(max); 
  for( int i = 0; i<= n; i++){
    x[i] = TMath::Exp( minlog+(maxlog-minlog)*i/n );
  }
  return x; 
}
double* GenLinArray( int n, double min, double max ){
  if( min > max ){
    std::cerr << __LINE__ << " " << __func__ << " : Min value is bigger than max value" << std::endl;
  }
  double *x = new double[n+1];
  for( int i = 0; i<= n; i++){
    x[i] = min +(max - min)*i/n;
  }
  return x; 
}
void    DrawRatioPad(TPad* pad){
  pad->cd();
  TPad* pad1 = new TPad(Form("%s_1",pad->GetName()),Form("%s_1",pad->GetName()),
			0.,0.33,1,1);
  TPad* pad2 = new TPad(Form("%s_2",pad->GetName()),Form("%s_2",pad->GetName()),
			0.,0.,1,0.33);
  pad1->SetNumber(1);
  pad2->SetNumber(2);
  pad1->SetMargin(0.15,0.1,0,0.1);
  pad2->SetMargin(0.15,0.1,0.2,0);
  pad1->SetGridx();
  pad1->SetGridy();
  pad2->SetGridx();
  pad2->SetGridy();
  pad1->Draw();
  pad2->Draw();
}
TH1D*   GenRatioHist( TH1D* h1, TH1D* h2 ){
  TH1D* hisRatio = new TH1D(Form("%s_%s_Ratio",h1->GetName(),h2->GetName()),
			    Form("%s_%s_Ratio;%s;Ratio",h1->GetName(),h2->GetName(),h1->GetXaxis()->GetTitle()),
			    h1->GetNbinsX(),h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
  for( int ibin = 1; ibin <= h1->GetNbinsX(); ibin++){
    double e[2];
    double n[2]={0,0};
    n[0]= h1->GetBinContent(ibin);
    n[1]= h2->GetBinContent(ibin);
    e[0]= h1->GetBinError(ibin);
    e[1]= h2->GetBinError(ibin);
    if( n[0] ==0 || n[1] == 0 ){
      continue; 
    }
    double Error=TMath::Sqrt( TMath::Power(e[0]*n[1]/n[0],2)+ TMath::Power(e[1],2))/n[1];
    hisRatio->SetBinContent(ibin,n[0]/n[1]);
    hisRatio->SetBinError(ibin,Error);
  }
  return hisRatio;
}
