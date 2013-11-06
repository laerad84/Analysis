#include "ChisqCalc.h"


double GetTimeResSq( double Ene ){
  double value = 0;
  value = 0.3*0.3+3.592*3.592/Ene + 12.24*12.24/Ene/Ene;
  //value /=value/2;
  //value = 1/Ene;
  return value;
}

void RegionTimeChisq( TH2D* TimeShapeHist, TH2D* RefHist, int n, Double_t *R, Double_t *D, Double_t *E, Double_t *T, Double_t Offset, Double_t& chisq , Double_t& NDF){
  int mndf = 0; 
  double mchisq=0;
  for( int i = 0; i< n ;i++){
    if( TMath::Abs( D[i] ) > 12.5 * 5 ){ continue; }
    if( R[i] < -12.5*5 || R[i] > 12.5 * 7 ){ continue; }
    int ibinCount = RefHist->Fill( R[i],D[i] );
    double binerror = TimeShapeHist->GetBinError(ibinCount);
    if( binerror == 0 ){ continue; }
    double sigma = GetTimeResSq(E[i]);
    double BinCenter = TimeShapeHist->GetBinContent( ibinCount );
    double fracChisq = TMath::Power( BinCenter - (T[i]- Offset),2 )/sigma;
    mchisq += fracChisq;
    mndf++;
  }
  if( mndf > 1 ){
    chisq = mchisq/(mndf -1 );
    NDF = mndf;
  }else{
    chisq = 0xFFFF;
    NDF   = 0;
  } 
  return;
}
void TimeChisq( TH2D* TimeShapeHist, TH2D* RefHist, int n, Double_t *R, Double_t *D, Double_t *E, Double_t *T, Double_t Offset, Double_t& chisq , Double_t& NDF){
  int mndf = 0; 
  double mchisq=0;
  for( int i = 0; i< n ;i++){
    if( TMath::Abs( D[i] ) > 12.5 * 5 ){ continue; }
    if( R[i] < -12.5*5 || R[i] > 12.5 * 7 ){ continue; }
    int ibinCount = RefHist->Fill( R[i],D[i] );
    double binerror = TimeShapeHist->GetBinError(ibinCount);
    if( binerror == 0 ){ continue; }
    double sigma = GetTimeResSq(E[i]);
    double BinCenter = TimeShapeHist->GetBinContent( ibinCount );
    double fracChisq = TMath::Power( BinCenter - (T[i]- Offset),2 )/sigma;
    mchisq += fracChisq;
    mndf++;
  }
  if( mndf > 1 ){
    chisq = mchisq/(mndf -1 );
    NDF = mndf;
  }else{
    chisq = 0xFFFF;
    NDF  = 0;
  } 
  return;
}


