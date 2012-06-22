/* *********************************************************

class WaveformFitter : FADC waveform fitter class
  author    : Eito IWAI ( iwai # champ.hep.sci.osaka-u.ac.jp )

history :
   v0.0.2    add Approx method to evaluate the channel to do fitting
   v0.0.1b   use TGraph::GetMean instead of TGraph::Fit("pol0") to calculate mean values
   v0.0.1a   development version
   
********************************************************** */
#include "WaveformFitter.h"
#if !defined(__CINT__)
ClassImp(WaveformFitter)
#endif


WaveformFitter::WaveformFitter(int Nsamples, bool fixed, int PedSmpl):m_Nsamples(Nsamples),m_fixed(fixed),m_pedsmpl(PedSmpl),m_width(30.),m_asymm(0.1){
  m_fitfunc=new TF1("m_fitfunc",AsymmetricGaussian,16.,500.,5);
}

WaveformFitter::~WaveformFitter(void){
  delete m_fitfunc;
}

bool WaveformFitter::Approx(TGraph* gr){
  
  double gnd=0., height=0., mean=0.;
  double ttx[3],tty[4];
  
  tty[3]=0.;
  for ( int i=1;i<m_Nsamples;i++ ){
    for ( int j=0;j<2;j++ ){
      gr->GetPoint(i+j,ttx[j],tty[j]);
    }
    tty[2]=tty[0]+tty[1];
    if ( tty[3]<=tty[2] ){
      tty[3]=tty[2];
      height=(tty[0]+tty[1])/2.;
      mean=(ttx[0]+ttx[1])/2.;
    }
  }
  const int m = m_pedsmpl-1;
  double tx[m], ty[m];
  //double tx[7], ty[7];
  if ( mean>8*m_pedsmpl+50. ){
    for ( int i=1;i<m_pedsmpl;i++ ){ gr->GetPoint(i,tx[i-1],ty[i-1]); }
  }else{
    for ( int i=1;i<m_pedsmpl;i++ ){ gr->GetPoint(m_Nsamples-i,tx[i-1],ty[i-1]); }
  }
  
  TGraph tgr(m_pedsmpl-1,tx,ty);
  /*
    tgr.Fit("pol0","Q");
    gnd=tgr.GetFunction("pol0")->GetParameter(0);
  */
  gnd=tgr.GetMean(2);
  height-=gnd;
  
  //gr->Clear();

  m_fitfunc->SetRange(mean-m_width*5.,mean+m_width*f_sigma);
  m_fitfunc->SetParLimits(0,0.9*height,2.0*height);
  m_fitfunc->SetParLimits(1,mean-16.,mean+16.);
  m_fitfunc->SetParLimits(4,gnd,gnd);
  if ( m_fixed ){
    m_fitfunc->SetParLimits(2,m_width,m_width);
    m_fitfunc->SetParLimits(3,m_asymm,m_asymm);
  }else{
    m_fitfunc->SetParLimits(2,10.,50.);
    m_fitfunc->SetParLimits(3,0.,1.);
  }
  m_fitfunc->SetParameters(height,mean,m_width,m_asymm,gnd);
  
  if ( height<10 ){
    return kFALSE;
  }else{
    return kTRUE;
  }
}

bool WaveformFitter::Fit(TGraph* gr){
  if ( !Approx(gr) ){ return kFALSE; }
  gr->Fit(m_fitfunc,"QR");
  return kTRUE;
}

bool WaveformFitter::Clear(){
  m_fitfunc->SetParameters(0,0,0,0,0);
  return kTRUE;
}

double WaveformFitter::GetParameter(int parNum){
  return m_fitfunc->GetParameter(parNum);
}
  
double AsymmetricGaussian(double* x, double* par){
  double t=x[0];
  double height=par[0];
  double mean=par[1];
  double b=par[2];
  double a=par[3];
  double ped=par[4];
  double t_fcn=0.;
  
  // asynmetric gaussian
  double sigma=a*(t-mean)+b;
  if ( sigma<0 ){ return ped; }
  t_fcn=height*TMath::Gaus(t,mean,sigma)+ped;
  
  return t_fcn;
}
double LinearGaussian(double* x, double* par){
	double t=x[0];
	double height = par[0];
	double mean   = par[1];
	double b      = par[2];
	double a      = par[3];
	double ped    = par[4];
	double t_fcn  = 0;
	double sigma  = a*(t-mean);
	if( sigma < 0 ){return ped;}
	t_fcn         = height*(1+t)/TMath::Sqrt(1+t*t)*TMath::Gaus(t,mean,sigma)+ped;
}
