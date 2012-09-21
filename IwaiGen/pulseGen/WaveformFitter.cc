/* *********************************************************

class WaveformFitter : FADC waveform fitter class
  author    : Eito IWAI ( iwai # champ.hep.sci.osaka-u.ac.jp )

history :
   v0.0.1a   development version
   
********************************************************** */
/***********************************************************
 JWLee : Divide to Waveform.h and Waveform.cc
 Arrange Parameter
 [0] : pedestal
 [1] : height
 [2] : mean
 ***********************************************************/

#include "WaveformFitter.h"

WaveformFitter::WaveformFitter(int Nsamples, bool fixed, int PedSmpl)
  :m_Nsamples(Nsamples),m_fixed(fixed),m_pedsmpl(PedSmpl),m_fitType(0){
  m_fitfunc=new TF1("m_fitfunc",AsymmetricGaussian,16.,500.,8);
  //m_fitfunc=new TF1("m_fitfunc",ScintiFunction,16.,500.,8);
  //m_fitfunc=new TF1("m_fitfunc",LnsFunction,16.,500.,6);
}

WaveformFitter::WaveformFitter(int Nsamples, bool fixed, int PedSmpl, int fitType)
 :m_Nsamples(Nsamples),m_fixed(fixed),m_pedsmpl(PedSmpl){
  m_fitType=fitType;
  switch(m_fitType){
  case 0:{
    m_fitfunc=new TF1("m_fitfunc1",AsymmetricGaussian,16.,500.,8);
    break;
  }
  case 1:{
    m_fitfunc=new TF1("m_fitfunc2",ScintiFunction,16.,500.,5);
    break;
  }
  case 2:{
    m_fitfunc=new TF1(Form("m_fitfunc3_%d",m_fixed),TypicalFunction,16.,500.,6);
    
    m_f=new TFile("combined.root");
    int t_nSeg=(m_fixed)? 4 : 0;
    for ( int i=t_nSeg;i<5;i++ ){
      for ( int y=0;y<12;y++ ){
	for ( int x=0;x<12;x++ ){
	  m_typShape[i][y][x]=(TH1D*)m_f->Get(Form("typShape%d_%d_%d_pfx_px",x,y,i));
	}
      }
    }
    
    m_tf=new TF1("m_tf","pol2",0.,500.);
  }
  }
}

WaveformFitter::~WaveformFitter(void){
  delete m_fitfunc;
}

bool WaveformFitter::Approx(TGraph* gr){
  double gnd=0., height=0., mean=0.;
  double ttx[2],tty[4];
  
  // fast search for peak height/position
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
  
  // pedestal calculation
  double tx[m_pedsmpl-1], ty[m_pedsmpl-1];
  if ( mean>8*m_pedsmpl+50. ){
    for ( int i=1;i<m_pedsmpl;i++ ){ gr->GetPoint(i,tx[i-1],ty[i-1]); }
  }else{
    for ( int i=1;i<m_pedsmpl;i++ ){ gr->GetPoint(m_Nsamples-i,tx[i-1],ty[i-1]); }
  }
  TGraph tgr(m_pedsmpl-1,tx,ty);
  gnd=tgr.GetMean(2);
  height-=gnd;
  
  if ( m_fitType==0 && m_fixed ){
    if ( m_fitPar[0]<35 ){
      m_fitfunc->SetRange(mean-30*5.,mean+30.*f_sigma); // 09/09/2010
    }else{
      m_fitfunc->SetRange(mean-40*5.,mean+50.*f_sigma); // 09/09/2010
    }
  }else{
    m_fitfunc->SetRange(mean-30*5.,mean+30.*f_sigma); // 09/09/2010
  }
  
  m_fitfunc->FixParameter(0,gnd);
  m_fitfunc->SetParLimits(1,0.5*height,1.5*height);
  m_fitfunc->SetParLimits(2,mean-32.,mean+32.);
  
  m_fitfunc->SetParameters(gnd,height,mean);
  if ( height>=0 ){
    switch(m_fitType){
    case 0:{ // asymmetric gaussian
      if ( m_fixed ){
	for ( int i=0;i<5;i++ ){ m_fitfunc->FixParameter(3+i,m_fitPar[i]); }
      }else{
	m_fitfunc->FixParameter(5,0.);
	m_fitfunc->FixParameter(6,0.);
	m_fitfunc->FixParameter(7,0.);
	m_fitfunc->SetParLimits(3,10.,70.);
	m_fitfunc->SetParLimits(4,-0.5,1.);
      }
      break;
    }
    case 1:{ // scinti funciton
      if ( m_fixed ){
	for ( int i=0;i<2;i++ ){ m_fitfunc->FixParameter(3+i,m_fitPar[i]); }
      }else{
	m_fitfunc->SetParameter(3,20.);
	m_fitfunc->SetParameter(4,30.);
      }
      break;
    }
    case 2:{ // typical function
      for ( int i=0;i<2;i++ ){ m_fitfunc->FixParameter(3+i,m_fitPar[i]); }
      if ( m_fixed ){ // only use the typical shape of larger pulses
	m_fitfunc->FixParameter(5,4);
	m_fitPar[2]=4;
      }else{
	int t_idx=5*TMath::Log10(height/30.)/TMath::Log10(5000./30);
	if ( t_idx>=5 ){ t_idx=4; }
	if ( t_idx<0 ){ t_idx=0; }
	m_fitfunc->FixParameter(5,t_idx);
	m_fitPar[2]=t_idx;
      }
      break;
    }
    default:
      ;
    }
  }
  
  if ( height<0 ){
    return kFALSE;
  }else{
    return kTRUE;
  }
}

bool WaveformFitter::Fit(TGraph* gr){
  if ( !Approx(gr) ){ return kFALSE; }

  if ( m_fitType==2 ){ // typical function
    double height=0., mean=0., sigma=0.;
    double gnd=m_fitfunc->GetParameter(0);
    double t_h=m_fitfunc->GetParameter(1);
    double t_t=m_fitfunc->GetParameter(2);
    
    height=t_h;
    mean=t_t;
    
    // better search for the peak
    double ttx[4], tty[5];
    if ( t_h>20 ){
      const int step=4;
      tty[4]=0.;
      for ( int i=step;i<m_Nsamples-step;i++ ){
	tty[3]=0.;
	for ( int j=0;j<3;j++ ){
	  gr->GetPoint(i+j-1,ttx[j],tty[j]);
	  tty[3]+=tty[j];
	}
	if ( tty[3]>=tty[4] ){
	  tty[4]=tty[3];
	  ttx[3]=i;
	}
      }
      double x[3], y[3];
      for ( int j=0;j<3;j++ ){
	gr->GetPoint(ttx[3]+step*(j-1),x[j],y[j]);
	y[j]-=gnd;
      }
      
      // gaussian correction
      double pars[2][3];
      for ( int j=0;j<3;j++ ){
	for ( int i=0;i<2;i++ ){
	  pars[i][0]=x[(i+j)%3]-x[(i+1+j)%3];
	  pars[i][1]=-x[(i+j)%3]*x[(i+j)%3]+x[(i+1+j)%3]*x[(i+1+j)%3];
	  pars[i][2]=TMath::Log(y[(i+1+j)%3]/y[(i+j)%3]);
	}
	if ( pars[0][0]*pars[1][2]-pars[1][0]*pars[0][2]==0. ){ continue; }
	mean=(pars[1][1]*pars[0][2]-pars[0][1]*pars[1][2])/(2*(pars[0][0]*pars[1][2]-pars[1][0]*pars[0][2]));
	for ( int i=0;i<2;i++ ){
	  pars[i][0]=(x[(i+j+1)%3]-mean)*(x[(i+j+1)%3]-mean);
	  pars[i][1]=TMath::Log(y[(i+j+1)%3]);
	}
	if ( pars[0][0]-pars[1][0]==0. ){ continue; }
	height=(pars[1][1]*pars[0][0]-pars[0][1]*pars[1][0])/(pars[0][0]-pars[1][0]);
	if ( height-pars[0][1]==0. ){ continue; }
	if ( 0.5*pars[0][0]/(height-pars[0][1])<0 ){ continue; }
	sigma=TMath::Sqrt(0.5*pars[0][0]/(height-pars[0][1]));
	height=TMath::Exp(height);
	
	// correction, when step=4
	mean-=0.6;
      }
      
      /*
      // pol2 correction
      double a=0., b=0., c=0.;
      for ( int j=0;j<3;j++ ){
      a=y[(0+j)%3]*(x[(2+j)%3]-x[(0+j)%3])-y[(1+j)%3]*(x[(2+j)%3]-x[(1+j)%3]);
      if ( x[(1+j)%3]-x[(0+j)%3] ){ continue; }
      a/=(x[(1+j)%3]-x[(0+j)%3]);
      a+=(y[(2+j)%3]-y[(1+j)%3]-y[(0+j)%3]);
      if ( x[(2+j)%3]*(x[(2+j)%3]-x[(1+j)%3]-x[(0+j)%3])+x[(0+j)%3]*x[(1+j)%3] ){ continue; }
      a/=(x[(2+j)%3]*(x[(2+j)%3]-x[(1+j)%3]-x[(0+j)%3])+x[(0+j)%3]*x[(1+j)%3]);
      b=(y[(0+j)%3]-y[(1+j)%3])/(x[(0+j)%3]-x[(1+j)%3])-a*(x[(0+j)%3]+x[(1+j)%3]);
      c=(y[(0+j)%3]*x[(1+j)%3]-y[(1+j)%3]*x[(0+j)%3])/(x[(1+j)%3]-x[(0+j)%3])+a*x[(0+j)%3]*x[(1+j)%3];
      
      mean=-0.5*b/a;
      height=c-0.25*b*b/a;
      
      if ( a>=0. || TMath::Abs(mean-t_t)>12. ||
      TMath::Abs(height/t_h-1)>0.2 || height<10 || mean<0 || mean>1000 ){
      std::cerr << t_h << "\t" << t_t << std::endl;
      std::cerr << height << "\t" << mean << " : " << a << std::endl << std::endl;
      height=-100.;
      }
      }
      */
      
      // spike candidate
      if ( height>80 && (sigma<20.|| TMath::Abs(height/t_h-1)>0.2) ){
	tty[4]=0.;
	for ( int i=step;i<m_Nsamples-step;i++ ){
	  if ( 8*(i+step+2)>mean ){ break; } // spike position
	  tty[3]=0.;
	  for ( int j=0;j<3;j++ ){
	    gr->GetPoint(i+j-1,ttx[j],tty[j]);
	    tty[3]+=tty[j];
	  }
	  if ( tty[3]>=tty[4] ){
	    tty[4]=tty[3];
	    ttx[3]=i;
	  }
	}
	double x[3], y[3];
	for ( int j=0;j<3;j++ ){
	  gr->GetPoint(ttx[3]+step*(j-1),x[j],y[j]);
	  y[j]-=gnd;
	}
	
	// gaussian correction
	double pars[2][3];
	for ( int j=0;j<3;j++ ){
	  for ( int i=0;i<2;i++ ){
	    pars[i][0]=x[(i+j)%3]-x[(i+1+j)%3];
	    pars[i][1]=-x[(i+j)%3]*x[(i+j)%3]+x[(i+1+j)%3]*x[(i+1+j)%3];
	    pars[i][2]=TMath::Log(y[(i+1+j)%3]/y[(i+j)%3]);
	  }
	  if ( pars[0][0]*pars[1][2]-pars[1][0]*pars[0][2]==0. ){ continue; }
	  mean=(pars[1][1]*pars[0][2]-pars[0][1]*pars[1][2])/(2*(pars[0][0]*pars[1][2]-pars[1][0]*pars[0][2]));
	  for ( int i=0;i<2;i++ ){
	    pars[i][0]=(x[(i+j+1)%3]-mean)*(x[(i+j+1)%3]-mean);
	    pars[i][1]=TMath::Log(y[(i+j+1)%3]);
	  }
	  if ( pars[0][0]-pars[1][0]==0. ){ continue; }
	  height=(pars[1][1]*pars[0][0]-pars[0][1]*pars[1][0])/(pars[0][0]-pars[1][0]);
	  if ( height-pars[0][1]==0. ){ continue; }
	  if ( 0.5*pars[0][0]/(height-pars[0][1])<0 ){ continue; }
	  sigma=TMath::Sqrt(0.5*pars[0][0]/(height-pars[0][1]));
	  height=TMath::Exp(height);
	  
	  // correction, when step=4
	  mean-=0.6;
	}
	//std::cerr << "[Warning] WaveformFitter::Fit >> Detect spike, search for another peak." << std::endl;
	//std::cerr << height << ", " << mean << ", " << sigma << std::endl;
	
	if ( !m_fixed ){
	  int t_idx=5*TMath::Log10(height/30.)/TMath::Log10(5000./30);
	  if ( t_idx>=5 ){ t_idx=4; }
	  if ( t_idx<0 ){ t_idx=0; }
	  m_fitfunc->FixParameter(5,t_idx);
	  m_fitPar[2]=t_idx;
	}
      }
      
      // invalid smaller pulse : use default value
      if ( sigma==0. || height<10 || mean<0 || mean>1000 ){
	height=t_h;
	mean=t_t;
      }
    }
    
    m_fitfunc->SetParameter(1,height);
    //m_fitfunc->SetParameter(1,mean);
    m_fitfunc->SetParameter(2,mean-20.5);
    if ( height>100 ){ m_fitfunc->SetParLimits(0,0.5*height,1.5*height); }
    else{ m_fitfunc->SetParLimits(0,2.,5.*height); }
    m_fitfunc->SetParLimits(2,mean-32.,mean+32.);
    m_fitfunc->SetRange(mean-30.*5,mean+30.*f_sigma);
  }
  
  if ( m_fitType==2 && !m_fixed && std::atoi(m_typShape[(int)m_fitPar[2]][(int)m_fitPar[1]][(int)m_fitPar[0]]->GetTitle())<(int)m_fitPar[2] ){
    std::cerr << "[Warning] WaveformFitter::Fit >> Fitter is about to use the typical shape of the smaller pulses" << std::endl;
    std::cerr << "(x,y)=(" << (int)m_fitPar[0] << "," << (int)m_fitPar[1] << "), idx=" << std::atoi(m_typShape[(int)m_fitPar[2]][(int)m_fitPar[1]][(int)m_fitPar[0]]->GetTitle()) << " < " << (int)m_fitPar[2] << std::endl;
  }
  gr->Fit(m_fitfunc,"QR+");
  return kTRUE;
}

double AsymmetricGaussian(double* x, double* par){

  double t      = x[0];
  double ped    = par[0];
  double height = par[1];
  double mean   = par[2];
  
  double t_par[5];
  for ( int i=0;i<5;i++ ){
    t_par[i]=par[3+i];
  }
  
  // asynmetric gaussian
  double sigma=0.;
  for ( int i=0;i<5;i++ ){
    sigma=sigma*(t-mean)+t_par[4-i];
    //sigma=sigma*(t-mean)+par[7-i];
  }
  if ( sigma<0 ){ return ped; }
  return height*TMath::Gaus(t,mean,sigma)+ped;
}

double ScintiFunction(double* x, double* par){
  double t      = x[0];
  double ped    = par[0];
  double height = par[1];
  double mean   = par[2];
  double tauR   = par[3];
  double tau    = par[4];
  
  return 2*height*TMath::Freq((t-mean)/tauR)*TMath::Exp(-(t-mean)/tau)+ped;
}

double TypicalFunction(double* x, double* par){
  double t     = x[0];
  double ped   = par[0];
  double height= par[1];
  double mean  = par[2];
  int tx       = (int)par[3];
  int ty       = (int)par[4];
  int idx      = (int)par[5];
  
  double p[3], q[3];
  double tt=t+(140.-mean);
  int t_mid=(int)(tt+1);
  if ( t_mid<1 || 48*8<t_mid ){ return ped; }
  if ( m_typShape[idx][ty][tx]->GetBinContent(t_mid)==0. ){
    int i=1;
    while ( 1 ){
      if ( m_typShape[idx][ty][tx]->GetBinContent(t_mid-TMath::Power(-1,i)*i)!=0. ){
	t_mid-=(TMath::Power(-1,i)*i);
	break;
      }
      ++i;
    }
  }
  p[1]=m_typShape[idx][ty][tx]->GetBinCenter(t_mid);
  q[1]=m_typShape[idx][ty][tx]->GetBinContent(t_mid);
  
  int i=1;
  while ( 1 ){
    if ( m_typShape[idx][ty][tx]->GetBinContent(t_mid+i)!=0. ){ break; }
    ++i;
  }
  p[2]=m_typShape[idx][ty][tx]->GetBinCenter(t_mid+i);
  q[2]=m_typShape[idx][ty][tx]->GetBinContent(t_mid+i);
  
  i=1;
  while ( 1 ){
    if ( m_typShape[idx][ty][tx]->GetBinContent(t_mid-i)!=0. ){ break; }
    ++i;
  }
  p[0]=m_typShape[idx][ty][tx]->GetBinCenter(t_mid-i);
  q[0]=m_typShape[idx][ty][tx]->GetBinContent(t_mid-i);
  
  double tmp=q[0]*(p[2]-p[0])-q[1]*(p[2]-p[1]);
  tmp/=(p[1]-p[0]);
  tmp+=(q[2]-q[1]-q[0]);
  tmp/=(p[2]*(p[2]-p[1]-p[0])+p[0]*p[1]);
  m_tf->SetParameter(2,tmp);
  
  tmp=(q[0]-q[1])/(p[0]-p[1])-m_tf->GetParameter(2)*(p[0]+p[1]);
  m_tf->SetParameter(1,tmp);
  
  tmp=(q[0]*p[1]-q[1]*p[0])/(p[1]-p[0])+m_tf->GetParameter(2)*p[0]*p[1];
  m_tf->SetParameter(0,tmp);
  
  //std::cerr << t << ", " << par[0] << ", " << par[1] << ", " << height*m_tf->Eval(tt)+ped << std::endl;
  
  return height*m_tf->Eval(tt)+ped;
}

double LnsFunction(double* x, double* par){
  double t     = x[0];
  double ped   = par[0];
  double height= par[1];
  double mean  = par[2];
  double TauR  = par[3];
  double Tau   = par[4];  
  double tmp   = (1-TMath::Exp(-(t-mean)/TauR))*TMath::Exp(-(t-mean)/Tau);
  if ( tmp<0 ){
    return ped;
  }else{
    return height*tmp+ped;
  }
}


#endif // Waveform_h
