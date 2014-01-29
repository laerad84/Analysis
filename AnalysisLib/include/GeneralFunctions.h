#ifndef  GENERALFUNCTIONS_H_
#define  GENERALFUNCTIONS_H_

#include "TMath.h"
#include "TPad.h"
#include "TH1.h"


double  ConvertRadToDeg( double rad);
double  ConvertDegToRad( double deg );

double* GenLogArray(int n,double xmin, double xmax);
double* GenLinArray(int n,double xmin, double xmax);
void    DrawRatioPad(TPad* pad);
TH1D*   GenRatioHist( TH1D* h1, TH1D* h2);

#endif  //GENERALFUNCTIONS_H_
