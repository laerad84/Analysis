#ifndef  GENERALFUNCTIONS_H_
#define  GENERALFUNCTIONS_H_

#include "TMath.h"
#include "TPad.h"
#include "TH1.h"
#include "klong/Klong.h"
#include "gamma/Gamma.h"

double  ConvertRadToDeg( double rad);
double  ConvertDegToRad( double deg );

double* GenLogArray(int n,double xmin, double xmax);
double* GenLinArray(int n,double xmin, double xmax);
void    DrawRatioPad(TPad* pad);
TH1D*   GenRatioHist( TH1D* h1, TH1D* h2);
double  THCorrectionFunction( double* x, double* p);

double  THCorrectionFunction( double* x, double *p){
  Double_t value(0);
  value = p[0]+p1[0]*TMath::Log(1+p[2]*TMath::Exp(x[0]/2000));
  return value;
}



#endif  //GENERALFUNCTIONS_H_
