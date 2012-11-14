#include "PulseGenerator.h"

const double PulseGenerator::pdfPar[5]={3.28113, 0.541526, -0.00582465, 4.53179e-05, -1.16776e-07};  
const double PulseGenerator::corPar[5]={-0.197912, 0.0695365, 0.161864, 0.193783, 0.232639};

PulseGenerator::PulseGenerator(){
  MakeFunction();
}
PulseGenerator::~PulseGenerator(){
  if( pePDF    != NULL ) { delete pePDF; pePDF = NULL; }
  if( pfPDF    != NULL ) { delete pfPDF; pfPDF = NULL; }
  if( Waveform != NULL ) { delete Waveform; Waveform = NULL; }
}
void PulseGenerator::MakeFunction(){
  
  pfPDF = new TF1("pfPDF", agaus, 0., 8.*nSample, 4 );
  pfPDF->SetParameter( 0, 1.652);//height
  // 1.652 = 240[cnt] / 13.7[MeV] / 12.7[P.E./MeV] * 16000[P.E] / 13360[cnt];
  // E_to_PE = 12.7;//for iwai crystal;
  //pfPDF->SetParameter( 0, normF);     //height
  pfPDF->SetParameter( 1, 4.*nSample);//mean
  pfPDF->SetParameter( 2, 0.44321);   //a
  pfPDF->SetParameter( 3, 23.71 );    //b

  pePDF = new TF1("pePDF", Asymmetric, 0., 140+110., 8 );
  pePDF->SetParameter( 0, 1.);     //height
  pePDF->SetParameter( 1, peMean );//mean
  pePDF->SetParameter( 2, 0.);     //ped
  for( int i = 0; i< 5; i++){
    pePDF->SetParameter( 3+i , pdfPar[i]* (1+corPar[i] ));// Sigma of Asymmetric gaussian
  }
  
  for( int pointIndex = 0; pointIndex < fSeg*8*nSample; pointIndex++){
    pFuncArray[pointIndex] = pfPDF->Eval(1./fSeg*pointIndex );
  }
}
void PulseGenerator::SetLightYield( double E_to_PE ){
  double NormalizeFactor = E_to_PE/12.7*1.652;
  pfPDF->SetParameter( 0, NormalizeFactor );
  for( int pointIndex = 0; pointIndex < fSeg*8*nSample; pointIndex++){
    pFuncArray[pointIndex] = pfPDF->Eval(1./fSeg*pointIndex );
  }
}
TF1* PulseGenerator::GetWaveform( double Energy, double SignalTime, double E_to_PE){
  SetLightYield( E_to_PE);
  if( Waveform != NULL ){ delete Waveform ; Waveform = NULL;}
  //Double_t t_npe = 0;
  Int_t nPE = 0 ;
  while( nPE == 0){
    nPE = gRandom->Poisson( Energy*E_to_PE );
  }
  for( int iPE = 0; iPE < nPE; iPE++){
    Waveform->SetParameter( 2*iPE+2, gRandom->PoissonD( lambdaPMT )/ lambdaPMT );
    Waveform->SetParameter( 2*iPE+1, pePDF->GetRandom() + SignalTime );
  }

  return Waveform;
}
TF1* PulseGenerator::GetWaveform( std::vector<double> EnergyArr, std::vector<double> SignalTimeArr , Double_t E_to_PE){
  SetLightYield( E_to_PE );
  if( Waveform != NULL ){ delete Waveform; Waveform = NULL; }  
  Double_t t_npe  = 0;
  std::vector<double> nPEArrInt;
  for( std::vector<double>::iterator it = EnergyArr.begin();it != EnergyArr.end(); it++){
    Int_t nPE = 0; 
    while ( nPE == 0 ){
      nPE = gRandom->Poisson( (*it)*E_to_PE );
    }
    t_npe += nPE;
    nPEArrInt.push_back( nPE );
  }
  
  Waveform = new TF1("Waveform",pulseFunc,0.,8.*nSample, 2*t_npe+1);
  Waveform->SetParameter(0,t_npe);
  std::vector<double>::iterator itPE = nPEArrInt.begin();
  std::vector<double>::iterator itTime   = SignalTimeArr.begin();
  for(; itPE != nPEArrInt.end() || itTime != SignalTimeArr.end() ; itPE++, itTime++){
    for( int iPE = 0; iPE < *itPE; iPE++){
      Waveform->SetParameter( 2*iPE+2, gRandom->PoissonD( lambdaPMT )/lambdaPMT );
      Waveform->SetParameter( 2*iPE+1, pePDF->GetRandom() + *itTime );
    }
  }
  return Waveform;
}
double PulseGenerator::agaus(double* x, double* par){
    double t=x[0];
    double height=par[0];
    double mean=par[1];
    double a=par[2];
    double b=par[3];
    double t_fcn=0.;    
    // asynmetric gaussian
    double sigma=a*(t-mean)+b;
    if ( sigma<0 ){ return 0.; }
    t_fcn=height*TMath::Gaus(t,mean,sigma);    
    return t_fcn;
}
double PulseGenerator::pulseFunc(double* x, double* par){
  double t=x[0];
  double t_fcn=0.;
  
  int t_npe=(int)par[0];
  for ( int i=0;i<t_npe;i++ ){
    /*
      pfPDF->SetParameter(1,(double)par[2*i+1]);
      t_fcn+=par[2*i+2]*pfPDF->Eval(t);
    */
    // to boost generation speed
    //int t_int=(int)(t+4*nSample-par[2*i+1])*fSeg;
    int t_int=(int)(t-par[2*i+1])*fSeg;
    if ( t_int>=0 && t_int<8*nSample*fSeg ){
      t_fcn+=par[2*i+2]*pFuncArray[t_int];
    }
  }    
  return t_fcn;
}
double PulseGenerator::Asymmetric( double *x, double *par ){
  
  double t      = x[0];
  double height = par[0];
  double mean   = par[1];
  double ped    = par[2];

  double t_par[5];
  for ( int i=0;i<5;i++ ){
    t_par[i]=par[3+i];
  }
  
  // asynmetric gaussian
  double sigma=0.;
  // Sigma is pol4 of (t-mean) 
  for ( int i=0;i<5;i++ ){
    sigma=sigma*(t-mean)+t_par[4-i];
    //sigma=sigma*(t-mean)+par[7-i];
  }
  if ( sigma<0 ){ return ped; }
  return height*TMath::Gaus(t,mean,sigma)+ped;
}
