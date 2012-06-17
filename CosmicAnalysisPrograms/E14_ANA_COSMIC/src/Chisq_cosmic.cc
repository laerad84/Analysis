#include "Chisq_cosmic.h"

Chisq_cosmic::Chisq_cosmic(){
  
  m_Cosmic = false;
  m_roh    = 0;
  m_theta  = -90;
  m_func   = NULL;
  
}

Chisq_cosmic::~Chisq_cosmic(){
  Reset();
}

void Chisq_cosmic::Reset(){
  if( m_func != NULL ){
    m_func->Clear();
    m_func = NULL;
  }
  m_Cosmic = false;
  m_roh    = 0;
  m_theta  = -90;
  m_func   = NULL;
  m_init_roh = 0;
  m_init_theta = -90;
}

void Chisq_cosmic::GetRohTheta(double& roh, double &theta){
  roh = GetRoh();
  theta = GetTheta();
}

Double_t Chisq_cosmic::GetRoh() const {
  return m_roh;
}
Double_t Chisq_cosmic::GetTheta() const {
  return m_theta;
}

Bool_t Chisq_cosmic::SetFunction(TGraph* gr){
  Double_t x;
  Double_t y;
  Double_t xi_sq=0;
  Double_t yi_sq=0;
  Double_t xi_yi=0;
  Double_t xi=0;
  Double_t yi=0;
  Int_t nPoint  = gr->GetN();
  Double_t crystalWidth = 100;
  Double_t corr = 0;
  for( int ipoint =0 ;ipoint < nPoint; ipoint++){
    x = gr->GetX()[ipoint];
    y = gr->GetY()[ipoint];
    /*
    if( abs(x) >600 || abs(y) >600){//small or large
      crystalWidth = 50.;
    }else{
      crystalWidth = 25.;
    }
    */
    crystalWidth = 25;
    

    /*
    corr  = 1;//(crystalWidth*crystalWidth);    
    xi_sq += x*x/(crystalWidth*crystalWidth);    
    yi_sq += y*y/(crystalWidth*crystalWidth);
    xi_yi += x*y/(crystalWidth*crystalWidth);
    xi    += x/(crystalWidth*crystalWidth);
    yi    += y/(crystalWidth*crystalWidth);
    */
    corr  += 1;
    xi_sq += x*x;
    yi_sq += y*y;
    xi_yi += x*y;
    xi    += x;
    yi    += y;
    
    
  }
  
  m_func = new TF2("Func",Form("%lf*cos(x)*cos(x)+%lf*sin(x)*sin(x)+2*%lf*cos(x)*sin(x)-2*y*(%lf*cos(x)+%lf*sin(x))+%lf*y*y",
			       xi_sq,yi_sq,xi_yi,xi,yi,corr),
		   -TMath::Pi()/2,TMath::Pi()/2,
		   -1000,1000);
  //-900,900);
		   /*
		   (m_init_theta-10)*TMath::Pi()/180,
		   (m_init_theta+10)*TMath::Pi()/180,
		   m_init_roh-100,m_init_roh+100);
		   */
  //std::cout << xi_sq << "\t" << yi_sq << std::endl;
  if( m_func == NULL ){
    return false;
  }else{
    return true;
  }
}

Bool_t Chisq_cosmic::SetRange( double roh, double theta){
  m_init_roh = roh;
  m_init_theta = theta;
  return true;
}

Double_t Chisq_cosmic::CalChisq(){
  if( m_func == NULL ){
    return -1;
  } 
  //std::cout<< m_init_roh << "\t" << m_init_theta << std::endl;
  Double_t m_theta_rad;
  m_func->GetMinimumXY(m_theta_rad,m_roh);
  m_theta = m_theta_rad/TMath::Pi()*180;
  //std::cout<< "theta " << m_theta << " Roh:" << m_roh << std::endl;
  m_chisq = m_func->Eval(m_theta_rad,m_roh);
  return  m_chisq;
}

Double_t Chisq_cosmic::GetChisq() const {
  return m_chisq;
}

Double_t Chisq_cosmic::GetDistance(double x, double y) const {
  Double_t cos = TMath::Cos(m_theta*TMath::Pi()/180);
  Double_t sin = TMath::Sin(m_theta*TMath::Pi()/180);
  Double_t Distance = TMath::Abs((x-(1./cos*(m_roh - y*sin)))*cos);
  return Distance;
}

TF2*     Chisq_cosmic::GetFunction() const {
  return m_func;
}
  

TLine* Chisq_cosmic::GetLine(){
  Double_t cos  = TMath::Cos(m_theta*TMath::Pi()/180.);
  Double_t sin  = TMath::Sin(m_theta*TMath::Pi()/180.);
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

Double_t Chisq_cosmic::GetCalibrationFactor(){
  Double_t factor = TMath::Cos(m_theta*TMath::Pi()/180.);
  return factor; 
}
