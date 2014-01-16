#ifndef E14COSMIC_HOUGH_H_
#define E14COSMIC_HOUGH_H_
#include <cmath>

#include <TROOT.h>
#include <TObject.h>
#include <TMath.h>

#include <TF1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TLine.h>

class E14CosmicHough{
public:
	static const Int_t    nBinx;
	static const Int_t    nBiny;
	static const Double_t thetaMax;
	static const Double_t rohMax;
	static const Int_t    nDiv;
	static const Double_t UnitDeg;
	static const Double_t DegreeToPi;
  //  static       Int_t    nHoughCsI;
public:
  Bool_t   m_Cosmic;
  Double_t m_roh;
  Double_t m_theta;
  TH2D*    m_hisHough;
  
  E14CosmicHough();
  virtual ~E14CosmicHough();

  void          Reset();
  void          GetRohTheta(double& roh, double& theta);
  Bool_t        CosmicJudgment(TGraphErrors* grEvent);
  Bool_t        CosmicJudgment(TGraph* grEvent);
  Double_t      GetRoh()   const ;
  Double_t      GetTheta() const ;
  Double_t      GetTangentXbase() const ;
  Double_t      GetTangentYbase() const ;
  Double_t      GetOffsetXbase()  const ;
  Double_t      GetOffsetYbase()  const ;
  TH2D*         GetHisHough();
  Bool_t        IsCosmic();
  TF1*          GetFunction();
  Double_t      CalDistance(Double_t x , Double_t y) const ;
  Double_t      CalWeightedDistance( double x, double y ) const ;
  TLine*        GetLine();
  Double_t      GetCalibrationFactor();

  //TGraphErrors* GetRemains();
private:
  void          Init();

  ClassDef(E14CosmicHough,1)
};

#endif
