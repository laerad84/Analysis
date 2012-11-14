#include <iostream>
#include <string>
#include "TGraph.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStrig.h"
#include "TH2.h"

class PulseGenerator {
 public:

  static const int nSample = 64;;
  static const int peMean  = 160;
  static const double absLY= 12.7;
  static const double gndSgm = 2.4; 
  static const double normF = 1.652;
  static const double lambdaPMT = 14.; 
  static const int    fSeg = 10; 
  // p3-p7, valid in range [p1-8,p1+160] (p1:peak time)
  static const double pdfPar[5]={3.28113, 0.541526, -0.00582465, 4.53179e-05, -1.16776e-07};  
  // correction factor
  static const double corPar[5]={-0.197912, 0.0695365, 0.161864, 0.193783, 0.232639};

  double pFuncArray[8*nSample*fSeg];  
  double agaus( double *x, double* par );
  double pulseFunc( double *x, double *par );
  double DecayFunc( double *x, double *par );

  TF1* pePDF;
  TF1* pfPDF;
  TF1* Waveform;

  void MakeFunction();
  TF1* GetWaveForm( std::vector<double>  Energy, std::vector<double> time,
		    double NormalizeFactor = 1.652 , double LightYield = 12.7);  
  PulseGenerator();
  virtual ~PulseGenerator();
};

void PulseGenerator::MakeFunction(){
  
  pfPDF = new TF1("pfPDF", agaus, 0., 8.*nSameple, 4 );
  pfPDF->SetParameter( 0, normF);     //height
  pfPDf->SetParameter( 1, 4.*nSample);//mean
  pfPDF->SetParameter( 2, 0.44321);   //a
  pfPDF->SetParameter( 3, 23.71 );    //b

  pePDF = new TF1("pePDF", this->AsymmetricGaussian, 0., 140+110., 8 );
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

void SetNormalizeFactor( double NormalizeFactor ){
  pfPDF->SetParameter( 0, NormalizeFactor );
  for( int pointIndex = 0; pointIndex < fSeg*8*nSample; pointIndex++){
    pFuncArray[pointIndex] = pfPDF->Eval(1./fSeg*pointIndex );
  }
}

TF1* PulseGenerator::GetWaveForm( std::vector<double> Energy, std::vector<double> time, double NormalizeFactor, double LightYield ){
  SetGainFactor( NormalizeFactor );
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


double PulseGenerator::AsymmetricGaussian( double *x, double *par ){
  
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
