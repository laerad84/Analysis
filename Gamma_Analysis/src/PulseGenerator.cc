#include "PulseGenerator.h"
#if !defined(__CINT__)
ClassImp(PulseGenerator)
#endif

//const double PulseGenerator::pdfPar[5]={3.28113, 0.541526, -0.00582465, 4.53179e-05, -1.16776e-07};  
//const double PulseGenerator::corPar[5]={-0.197912, 0.0695365, 0.161864, 0.193783, 0.232639};
const int PulseGenerator::nSample = 64;
const int PulseGenerator::fSeg    = 10;
double* PulseGenerator::pFuncArray = new double[PulseGenerator::nSample*PulseGenerator::fSeg*8];
PulseGenerator::PulseGenerator():
  peMean(160) , absLY(12.7) , absGain(17.52), conversionConstant(1.2),
  gndSgm( 2.4), normF(1.652), lambdaPMT(14.)
{  
  pdfPar = new double[5];
  corPar = new double[5];
  pdfPar[0] = 3.28113;
  pdfPar[1] = 0.541526;
  pdfPar[2] = -0.00582465;
  pdfPar[3] = 4.53179e-05;
  pdfPar[4] = -1.16776e-07;
  corPar[0] = -0.197912;
  corPar[1] = 0.0695265;
  corPar[2] = 0.161864;
  corPar[3] = 0.193783;
  corPar[4] = 0.232639;

  MakeFunction();
}
PulseGenerator::~PulseGenerator(){
  delete pdfPar;
  delete corPar;
  /*
  if( pePDF    != NULL ) { delete pePDF; pePDF = NULL; }
  if( pfPDF    != NULL ) { delete pfPDF; pfPDF = NULL; }
  if( Waveform != NULL ) { delete Waveform; Waveform = NULL; }
  */
}

void PulseGenerator::MakeFunction(){
  normF = 17.52/181.8;//one photon height.when E_to_CNT= 1//
  // Shape of one photon signal //
  pfPDF = new TF1("pfPDF", agaus, 0., 8.*nSample, 4 );
  pfPDF->SetParameter( 0, normF);//height
  // 1.652 = 240[cnt] / 13.7[MeV] / 12.7[P.E./MeV] * 16000[P.E] / 13360[cnt];
  // E_to_PE = 12.7;//for iwai crystal;  
  pfPDF->SetParameter( 1, 4.*nSample);//mean
  pfPDF->SetParameter( 2, 0.044321);   //a
  pfPDF->SetParameter( 3, 23.71 );    //b

  // Decay func // 
  pePDF = new TF1("pePDF", Asymmetric, 0., 300., 8 );
  pePDF->SetParameter( 0, 1.);     //height
  pePDF->SetParameter( 1, peMean );//mean
  pePDF->SetParameter( 2, 0.);     //ped
  for( int i = 0; i< 5; i++){
    pePDF->SetParameter( 3+i , pdfPar[i]*(1+corPar[i] ));// Sigma of Asymmetric gaussian
  }  
  for( int pointIndex = 0; pointIndex < fSeg*8*nSample; pointIndex++){
    pFuncArray[pointIndex] = pfPDF->Eval(1./PulseGenerator::fSeg*pointIndex );
  }
}

void PulseGenerator::SetLightYield(  double E_to_CNT , double E_to_PE ){
  //std::cout<< __PRETTY_FUNCTION__ << std::endl;
}

TF1* PulseGenerator::GetWaveform( double Energy, double SignalTime, double E_to_CNT, double E_to_PE){
  //std::cout<< __PRETTY_FUNCTION__ << std::endl;
  //SetLightYield( E_to_CNT, E_to_PE);
  Double_t t_npe = 0;
  Int_t nPE = 0 ;
  //while( nPE != 0){
  //nPE = gRandom->Poisson( Energy*E_to_PE );
  //}
  nPE = gRandom->Poisson(Energy*E_to_PE);
  t_npe += nPE;
  Waveform = new TF1("Waveform",pulseFunc,0.,8.*nSample, 2*t_npe+2);
  Waveform->SetParameter(0,t_npe);
  Waveform->SetParameter(1,E_to_CNT);
  for( int iPE = 1; iPE < nPE; iPE++){
    Waveform->SetParameter( 2*iPE, gRandom->PoissonD( lambdaPMT )/ lambdaPMT );//one photion output deviation // 
    Waveform->SetParameter( 2*iPE+1, pePDF->GetRandom() + SignalTime );//Add Signal Time // 
  }  
  TF1* newFunc = (TF1*)Waveform->IsA()->New();
  Waveform->Copy(*newFunc);
  Reset();
 return newFunc;
}
void PulseGenerator::Reset(){
  if( Waveform != NULL ){ delete Waveform ; Waveform = NULL;}
}

TGraph* PulseGenerator::GenSignal( double Energy, double SignalTime, double E_to_CNT, double E_to_PE, double gnd){
  //std::cout<< __PRETTY_FUNCTION__ << std::endl;
  //SetLightYield( E_to_CNT, E_to_PE);
  Double_t t_npe = 0;
  Int_t nPE = 0 ;
  //while( nPE != 0){
  //nPE = gRandom->Poisson( Energy*E_to_PE );
  //}
  nPE = gRandom->Poisson(Energy*E_to_PE);
  t_npe += nPE;
  Waveform = new TF1("Waveform",pulseFunc,0.,8.*nSample, 2*t_npe+2);
  Waveform->SetParameter(0,t_npe);
  Waveform->SetParameter(1,E_to_CNT);
  for( int iPE = 1; iPE < nPE; iPE++){
    Waveform->SetParameter( 2*iPE, gRandom->PoissonD( lambdaPMT )/ lambdaPMT );//one photion output deviation // 
    Waveform->SetParameter( 2*iPE+1, pePDF->GetRandom() + SignalTime );//Add Signal Time // 
  }  

  TGraph* gr = new TGraph();
  for( int ipoint = 0 ; ipoint < PulseGenerator::nSample; ipoint++){
    double offset = 0;//8*(gRandom->Rndm()*-0.5);
    double tsum   = 0.;
    int d_tmp  = (int)(Waveform->Eval(ipoint*8.-offset)+gRandom->Gaus(0.,gndSgm)+gnd);
    gr->SetPoint(ipoint,ipoint*8,d_tmp);
  }      
  Reset();
  return gr;
}

TF1* PulseGenerator::GetWaveform( std::vector<double> EnergyArr, std::vector<double> SignalTimeArr , Double_t E_to_CNT, Double_t E_to_PE){
  //std::cout << __PRETTY_FUNCTION__ << std::endl;
  //SetLightYield( E_to_CNT, E_to_PE );
  Double_t t_npe  = 0;
  //std::cout<< t_npe << std::endl;
  std::vector<double> nPEArrInt;
  for( std::vector<double>::iterator it = EnergyArr.begin();it != EnergyArr.end(); it++){
    Int_t nPE = 0; 
    
    //while ( nPE == 0 ){
    //  if( *it==0 ){ break; }
    //  nPE = gRandom->Poisson( (*it)*E_to_PE );
    //}
    nPE = gRandom->Poisson((*it)*E_to_PE );
    t_npe += nPE;
    nPEArrInt.push_back( nPE );
  }
  Waveform = new TF1("Waveform",pulseFunc,0.,8.*nSample, 2*t_npe+2);
  Waveform->SetParameter(0,t_npe);
  std::vector<double>::iterator itPE = nPEArrInt.begin();
  std::vector<double>::iterator itTime   = SignalTimeArr.begin();
  int ipe = 0;
  for(; itPE != nPEArrInt.end() && itTime != SignalTimeArr.end() ; itPE++, itTime++){
    //std::cout<< t_npe << " : " << *itPE << " : " << *itTime << std::endl;
    Waveform->SetParameter(1,E_to_CNT);
    for( int iPE = 1; iPE < (*itPE); iPE++){
      Waveform->SetParameter( 2*iPE, gRandom->PoissonD( lambdaPMT )/lambdaPMT );//one photon output deviation//
      Waveform->SetParameter( 2*iPE+1, (pePDF->GetRandom()) + (*itTime) );
      ipe++;
    }
  }
  return Waveform;
  
}
TGraph* PulseGenerator::GenSignal( std::vector<double> EnergyArr, std::vector<double> SignalTimeArr , Double_t E_to_CNT, Double_t E_to_PE, double gnd){
  //std::cout << __PRETTY_FUNCTION__ << std::endl;
  //SetLightYield( E_to_CNT, E_to_PE );
  Double_t t_npe  = 0;
  //std::cout<< t_npe << std::endl;
  std::vector<double> nPEArrInt;
  for( std::vector<double>::iterator it = EnergyArr.begin();it != EnergyArr.end(); it++){
    Int_t nPE = 0; 
    
    //while ( nPE == 0 ){
    //  if( *it==0 ){ break; }
    //  nPE = gRandom->Poisson( (*it)*E_to_PE );
    //}
    nPE = gRandom->Poisson((*it)*E_to_PE );
    t_npe += nPE;
    nPEArrInt.push_back( nPE );
  }
  Waveform = new TF1("Waveform",pulseFunc,0.,8.*nSample, 2*t_npe+2);
  Waveform->SetParameter(0,t_npe);
  std::vector<double>::iterator itPE = nPEArrInt.begin();
  std::vector<double>::iterator itTime   = SignalTimeArr.begin();
  int ipe = 0;
  for(; itPE != nPEArrInt.end() && itTime != SignalTimeArr.end() ; itPE++, itTime++){
    //std::cout<< t_npe << " : " << *itPE << " : " << *itTime << std::endl;
    Waveform->SetParameter(1,E_to_CNT);
    for( int iPE = 1; iPE < (*itPE); iPE++){
      Waveform->SetParameter( 2*iPE, gRandom->PoissonD( lambdaPMT )/lambdaPMT );//one photon output deviation//
      Waveform->SetParameter( 2*iPE+1, (pePDF->GetRandom()) + (*itTime) );
      ipe++;
    }
  }
  //return Waveform;
  TGraph* gr = new TGraph();
  for( int ipoint = 0 ; ipoint < PulseGenerator::nSample; ipoint++){
    double offset = 0;//8*(gRandom->Rndm()*-0.5);
    double tsum   = 0.;
    int    d_tmp  = (int)(Waveform->Eval(ipoint*8.-offset)+gRandom->Gaus(0.,gndSgm)+gnd);
    gr->SetPoint(ipoint,ipoint*8,d_tmp);
  }      
  return gr;
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
  for ( int i=1;i<t_npe;i++ ){

  //pfPDF->SetParameter(1,(double)par[2*i+1]);
  //t_fcn+=par[2*i+2]*pfPDF->Eval(t);
    // to boost generation speed
    //int t_int=(int)(t+4*nSample-par[2*i+1])*fSeg;
    int t_int=(int)(t-par[2*i+1]+4*PulseGenerator::nSample)*fSeg;
    if ( t_int>=0 && t_int<8*nSample*fSeg ){
      t_fcn+=par[2*i]*pFuncArray[t_int];
    }
  }
  t_fcn = par[1]*t_fcn;  
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

