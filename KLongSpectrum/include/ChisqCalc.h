#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

double GetTimeResSq( double Ene );
void RegionTimeChisq( TH2D* TimeShapeHist, TH2D* RefHist, int n, Double_t *R, Double_t *D, Double_t *E, Double_t *T, Double_t Offset, Double_t& chisq , Double_t& NDF);
void TimeChisq( TH2D* TimeShapeHist, TH2D* RefHist, int n, Double_t *R, Double_t *D, Double_t *E, Double_t *T, Double_t Offset, Double_t& chisq , Double_t& NDF);
