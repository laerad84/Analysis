#include "E14CosmicAnalyzer.h"

E14CosmicAnalyzer::E14CosmicAnalyzer(): mb_chisq(false), mb_hough(false) {
  Init();
}

E14CosmicAnalyzer::~E14CosmicAnalyzer(){

}

bool E14CosmicAnalyzer::Init(){
  mc_hough      = new HoughCsI();
  mc_chi2Cosmic = new Chisq_cosmic();
  m_gr          = new TGraph();
  mb_chisq      = true;
  mb_hough      = true; 
  mc_hough->Reset();
  mc_chi2Cosmic->Reset();
  return true;
  
}

bool E14CosmicAnalyzer::Reset(){
  md_houghTheta = 90;
  md_houghRoh   = 0;
  md_chi2Theta  = 90;
  md_chi2Roh    = 0; 
  
  if( mb_chisq && mb_hough ) {
    mc_hough->Reset();
    mc_chi2Cosmic->Reset();
    return true;
  }else{
    return false;
  }
}

bool E14CosmicAnalyzer::GetResult(TGraph* gr, double& roh, double& theta){

  if( gr == NULL ){ return false;}
  if( !(mc_hough->CosmicJudgment(gr)) ){ return false;}
  m_gr->Set(0);
  md_houghRoh   =  mc_hough->GetRoh();
  md_houghTheta =  mc_hough->GetTheta(); 
  double *x = gr->GetX();
  double *y = gr->GetY();
  int index  =0; 
  for( int igr  =0; igr <  gr->GetN(); igr++){
    if( mc_hough->CalDistance( x[igr], y[igr] ) <= 50) {
      m_gr->SetPoint( index, x[igr], y[igr] );
      index++;
    }
  }
  mc_chi2Cosmic->SetFunction( m_gr );
  mc_chi2Cosmic->SetRange( md_houghRoh, md_houghTheta );
  mc_chi2Cosmic->CalChisq();
  md_chi2Roh   = mc_chi2Cosmic->GetRoh();
  md_chi2Theta = mc_chi2Cosmic->GetTheta();  
  roh   = md_chi2Roh;
  theta = md_chi2Theta;

  return true; 
}


bool E14CosmicAnalyzer::GetResult(TGraph* gr, double& roh, double& theta, double& roh1, double& theta1){

  if( gr == NULL ){ return false;}
  if( !(mc_hough->CosmicJudgment(gr)) ){ return false;}
  m_gr->Set(0);
  md_houghRoh   =  mc_hough->GetRoh();
  md_houghTheta =  mc_hough->GetTheta(); 
  double *x = gr->GetX();
  double *y = gr->GetY();
  int index  =0; 
  for( int igr  =0; igr <  gr->GetN(); igr++){
    if( mc_hough->CalDistance( x[igr], y[igr] ) <= 50) {
      m_gr->SetPoint( index, x[igr], y[igr] );
      index++;
    }
  }
  mc_chi2Cosmic->SetFunction( m_gr );
  mc_chi2Cosmic->SetRange( md_houghRoh, md_houghTheta );
  mc_chi2Cosmic->CalChisq();
  md_chi2Roh   = mc_chi2Cosmic->GetRoh();
  md_chi2Theta = mc_chi2Cosmic->GetTheta();  
  roh   = md_chi2Roh;
  theta = md_chi2Theta;
  roh1  = md_houghRoh;
  theta1= md_houghTheta;

  return true; 
}



