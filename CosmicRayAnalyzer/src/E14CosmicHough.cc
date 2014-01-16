#include <TMath.h>
#include <iostream>
#include "E14CosmicHough.h"
#if !defined(__CINT__)
ClassImp(E14CosmicHough)
#endif

const Double_t E14CosmicHough::rohMax   = 1000.;
const Double_t E14CosmicHough::thetaMax = TMath::Pi();
const Int_t    E14CosmicHough::nDiv     = 10;
const Int_t    E14CosmicHough::nBinx    = 180*E14CosmicHough::nDiv;
const Int_t    E14CosmicHough::nBiny    = 120;
const Double_t E14CosmicHough::UnitDeg  = 1.0/E14CosmicHough::nDiv;
const Double_t E14CosmicHough::DegreeToPi = TMath::Pi()/180.;
E14CosmicHough::E14CosmicHough(){
	Init();
}

E14CosmicHough::~E14CosmicHough(){
  delete m_hisHough;
}

void E14CosmicHough::Init(){
  std::cout << "Init" << std::endl;
  m_hisHough = new TH2D("hisHough","",
			nBinx+1,
			-(90+E14CosmicHough::UnitDeg/2)*E14CosmicHough::DegreeToPi,
			( 90+E14CosmicHough::UnitDeg/2)*E14CosmicHough::DegreeToPi,
			E14CosmicHough::nBiny+1,
			-(1+0.5/E14CosmicHough::nBiny)*E14CosmicHough::rohMax,
			(1+0.5/E14CosmicHough::nBiny)*E14CosmicHough::rohMax);
  m_roh = 0;
  m_theta = 0;
  m_Cosmic = kFALSE;
}

void E14CosmicHough::Reset(){
  m_hisHough->Reset();
  m_roh = 0;
  m_theta = 0;
  m_Cosmic = kFALSE;
}

Bool_t E14CosmicHough::CosmicJudgment(TGraphErrors* grEvent){
  Double_t* Xarr = grEvent->GetX();
  Double_t* Yarr = grEvent->GetY();
  Double_t* Earr = grEvent->GetEY();

  for( Int_t indexTheta = 0; indexTheta<nBinx+1; indexTheta++){
    for( Int_t index  = 0; index< grEvent->GetN(); index++){
      Double_t theta = ((double)indexTheta/nBinx)*TMath::Pi();
      Double_t roh = Xarr[index]*TMath::Cos(theta) + Yarr[index]*TMath::Sin(theta);
      m_hisHough->Fill( theta, roh, Earr[index]);
    }
  }
  
  Int_t x,y,z;
  m_hisHough->GetMaximumBin(x,y,z);
  //std::cout << x  << " ; " << y << " ; " << z << std::endl;
  Double_t OveralWeight=0;
  Double_t WeightedX=0;
  Double_t WeightedY=0;
  Double_t OveralRMS[3]={};
  
  m_theta = x - 90;
  m_roh   = -1*rohMax+2*rohMax/E14CosmicHough::nBiny*((double)y-0.5);
  m_Cosmic = kTRUE;
  return m_Cosmic;
}


Bool_t E14CosmicHough::CosmicJudgment(TGraph* grEvent){
  
  Double_t* Xarr = grEvent->GetX();
  Double_t* Yarr = grEvent->GetY();
  
  if( grEvent->GetN() >200  || grEvent->GetN() < 10){
    return kFALSE;
  }

  for( Int_t indexTheta = 0; indexTheta<nBinx+1; indexTheta++){
    for( Int_t index  = 0; index< grEvent->GetN(); index++){
      Double_t theta = ((double)indexTheta-90.)*E14CosmicHough::DegreeToPi;
      Double_t roh = Xarr[index]*TMath::Cos(theta) + Yarr[index]*TMath::Sin(theta);
      m_hisHough->Fill( theta, roh);
    }
  }

  Int_t x,y,z;
  m_hisHough->GetMaximumBin(x,y,z);
  //std::cout << x  << " ; " << y << " ; " << z << std::endl;
  Double_t OveralWeight=0;
  Double_t WeightedX=0;
  Double_t WeightedY=0;
  Double_t OveralRMS[3]={};
  
  m_theta = x*E14CosmicHough::UnitDeg - 90;
  m_roh   = -1*rohMax+2*rohMax/E14CosmicHough::nBiny*((double)y-0.5);
  m_Cosmic = kTRUE;

  //  std::cout << "theta:" << m_theta << " roh:" << m_roh << std::endl;
  return m_Cosmic;
}

void     E14CosmicHough::GetRohTheta(double &roh,double &theta){
  roh   = m_roh;
  theta = m_theta;
}

Double_t E14CosmicHough::GetRoh() const {
  return m_roh;
}

Double_t E14CosmicHough::GetTheta() const {
  return m_theta;
}

Double_t E14CosmicHough::GetTangentXbase() const {
  if(TMath::Sin(m_theta) == 0){
    return 999999;
  }else{
    return -1*TMath::Cos(m_theta)/TMath::Sin(m_theta);
  }
}

Double_t E14CosmicHough::GetOffsetXbase() const {
  if(TMath::Sin(m_theta) == 0){
    return 999999;
  }else{
    return m_roh/TMath::Sin(m_theta);
	}
}

Double_t E14CosmicHough::GetTangentYbase() const {
  if(TMath::Cos(m_theta) == 0){
    return 999999;
  }else{
    return -1*TMath::Sin(m_theta)/TMath::Cos(m_theta);
  }
}
Double_t E14CosmicHough::GetOffsetYbase() const {
  if(TMath::Cos(m_theta) == 0){
    return 999999;
  }else{
    return m_roh/TMath::Cos(m_theta);
  }
}

TH2D*    E14CosmicHough::GetHisHough(){
  return m_hisHough;
}


TF1*     E14CosmicHough::GetFunction(){
  if( TMath::Sin(m_theta) == 0 ){
    m_theta = 0.000000001;
  }
  TF1* func = new TF1("func",Form("%lf*x+%lf",
				  -1*TMath::Cos(m_theta*E14CosmicHough::DegreeToPi)/TMath::Sin(m_theta*E14CosmicHough::DegreeToPi),
				  m_roh/TMath::Sin(m_theta*E14CosmicHough::DegreeToPi)),
		      -1000.,1000.);
  //  std::cout<< m_theta << " : " 
  //	   << m_roh   << " : " 
  //	   << -1*TMath::Cos(m_theta*E14CosmicHough::DegreeToPi)/TMath::Sin(m_theta*E14CosmicHough::DegreeToPi) << " : " 
  //	   << func->GetParameter(0) << " : "
  //	   << func->GetParameter(1) <<  std::endl;
  return func;  
}

Bool_t   E14CosmicHough::IsCosmic(){
  return m_Cosmic;
}

Double_t E14CosmicHough::GetCalibrationFactor(){
  Double_t factor = TMath::Cos(m_theta*E14CosmicHough::DegreeToPi);
  return factor; 
}

TLine* E14CosmicHough::GetLine(){
  Double_t cos  = TMath::Cos(m_theta*E14CosmicHough::DegreeToPi);
  Double_t sin  = TMath::Sin(m_theta*E14CosmicHough::DegreeToPi);
  Double_t LineX[2];
  Double_t LineY[2];

  Double_t x[2], y[2];
  if( TMath::Sin(m_theta) == 0){
    m_theta = 0.000000001;
  }

  double b = TMath::Sqrt( 950*950 - m_roh*m_roh );
  
  LineX[0] =(m_roh*cos +b*sin);
  LineX[1] =(m_roh*cos -b*sin);
  LineY[0] =(m_roh*sin -b*cos);
  LineY[1] =(m_roh*sin +b*cos);
  /*
  LineX[0] =0;
  LineX[1] = m_roh*cos;
  LineY[0] = 0;
  LineY[1] = m_roh*sin;
  */
  TLine* line  = new TLine(LineX[0],LineY[0],LineX[1],LineY[1]);
  return line;
}

Double_t E14CosmicHough::CalDistance(Double_t x , Double_t y) const {
  Double_t cos = TMath::Cos(m_theta*E14CosmicHough::DegreeToPi);
  Double_t sin = TMath::Sin(m_theta*E14CosmicHough::DegreeToPi);
  Double_t Distance = TMath::Abs((x-(1./cos*(m_roh - y*sin)))*cos);
  //  std::cout << cos << " : " << sin << " : " << Distance << std::endl;

  return Distance;
}

Double_t E14CosmicHough::CalWeightedDistance(Double_t x , Double_t y) const {
  Double_t cos = TMath::Cos(m_theta*E14CosmicHough::DegreeToPi);
  Double_t sin = TMath::Sin(m_theta*E14CosmicHough::DegreeToPi);
  Double_t Distance = TMath::Abs((x-(1./cos*(m_roh - y*sin)))*cos);
  Double_t Weight = 0;
  if( TMath::Abs( x ) > 600 || TMath::Abs( y ) > 600 ){
    Weight = 50;
  }else{
    Weight = 25;
  }
  //  std::cout << cos << " : " << sin << " : " << Distance << std::endl;
  return Distance/Weight;
}
  
