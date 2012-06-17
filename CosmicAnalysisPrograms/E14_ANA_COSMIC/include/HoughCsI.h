#ifndef HOUGHCSI_H_
#define HOUGHCSI_H_
#include <cmath>

#include <TROOT.h>
#include <TObject.h>
#include <TMath.h>

#include <TF1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TLine.h>

class HoughCsI{
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
  
  HoughCsI();
  virtual ~HoughCsI();

  void          Reset();
  void          GetRohTheta(double& roh, double& theta);
  Bool_t        CosmicJudgment(TGraphErrors* grEvent);
  Bool_t        CosmicJudgment(TGraph* grEvent);
  Double_t      GetRoh();
  Double_t      GetTheta();
  Double_t      GetTangentXbase();
  Double_t      GetTangentYbase();
  Double_t      GetOffsetXbase();
  Double_t      GetOffsetYbase();
  TF1*          GetFunction();
  TH2D*         GetHisHough();
  Bool_t        IsCosmic();
  Double_t      GetCalibrationFactor();
  TLine*        GetLine();
  Double_t      CalDistance(Double_t x , Double_t y);

  //TGraphErrors* GetRemains();
private:
  void          Init();

  ClassDef(HoughCsI,1)
};

#endif
