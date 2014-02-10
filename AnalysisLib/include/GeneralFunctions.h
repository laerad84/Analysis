#ifndef  GENERALFUNCTIONS_H_
#define  GENERALFUNCTIONS_H_

#include "TMath.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "klong/Klong.h"
#include "gamma/Gamma.h"

double  ConvertRadToDeg( double rad);
double  ConvertDegToRad( double deg );

double* GenLogArray(int n,double xmin, double xmax);
double* GenLinArray(int n,double xmin, double xmax);
void    DrawRatioPad(TPad* pad);
TH1D*   GenRatioHist( TH1D* h1, TH1D* h2);
double  THCorrectionFunction( double* x, double* p);

//////////////////////////////////////////////////////
// Utility Functions 
//////////////////////////////////////////////////////

void DrawSlices(TH2D* his, TH1D** hisArr, int cord=0);

#endif  //GENERALFUNCTIONS_H_
