#include "TMath.h"
#include "TPad.h"
#include "TH1.h"


double  ConvertRadToDeg( double rad){
  return rad/TMath::Pi()*180;
}
double ConvertDegToRad( double deg ){
  return deg*TMath::Pi()/180;
}

double* GenLogArray(int n,double xmin, double xmax);
double* GenLinArray(int n,double xmin, double xmax);
void    DrawRatioPad(TPad* pad);
TH1D*   GenRatioHist( TH1D* h1, TH1D* h2);
