#ifndef E14COSMIC_CHISQ_H_
#define E14COSMIC_CHISQ_H_

#include <cmath>

#include "TMath.h"
#include "TObject.h"
#include "TF2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"

class E14CosmicChisq{
  
 public:

  Bool_t   m_Cosmic;
  Double_t m_roh;
  Double_t m_theta;
  Double_t m_chisq;
  Double_t m_init_roh;
  Double_t m_init_theta;

  TF2*     m_func;
  
  E14CosmicChisq();
  virtual ~E14CosmicChisq();
  
  void Reset();
  
  Bool_t   SetFunction(TGraph* gr);
  Bool_t   SetRange( double roh, double theta);
  
  Double_t CalChisq();

  void GetRohTheta(double& roh, double& theta) const ;
  Double_t GetRoh()   const ;
  Double_t GetTheta() const ;
  Double_t GetChisq() const ;
  Double_t CalDistance(double x, double y) const ;
  Double_t CalWeightedDistance( double x, double y ) const; 
  TF2*     GetFunction() const;
  TLine*        GetLine();
  Double_t GetCalibrationFactor();
};
#endif //E14COSMIC_CHISQ_H_
