#include <iostream>
#include <string>
#include "TGraph.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TH2F.h"

#include "TUtil.C"
#include "WaveformFitter.cc"

const int nEvent=12;
const int nSample=64;
//const double peMean=160.; // ns, peak of the pePDF
const double peMean = 250;
const int nnpe=8;
const int npe[nnpe]={40,80,160,320,640,1200,2000,4000};

const double absLY=12.7; // pe/MeV
const double gndSgm=2.4; // [cnt] raw value

// MIP 13.7MeV-240cnt, 16000pe-height13360
// 240/13.7/12.7*16000/13360 = 1.38*1.2 =1.652
const double normF=1.652;
											 


//const double lambdaPMT=3.952; // RV580UV@1800V
//const double lambdaPMT=8.283; // H7195@2100V
const double lambdaPMT=14.; // R4125@1400V
//const double lambdaPMT=17.; // R4125@1900V

const int fSeg=10; // ex. fSeg=10 means 1/10=0.1ns
double pFuncArray[8*nSample*fSeg];

// p3-p7, valid in range [p1-8,p1+160] (p1:peak time)
const double pdfPar[5]={3.28113, 0.541526, -0.00582465, 4.53179e-05, -1.16776e-07};

// correction factor
const double corPar[5]={-0.197912, 0.0695365, 0.161864, 0.193783, 0.232639};


double agaus(double* x, double* par){
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

TF1* pfPDF;
double pulseFunc(double* x, double* par){
    double t=x[0];
    double t_fcn=0.;
  
    int t_npe=(int)par[0];
    for ( int i=0;i<t_npe;i++ ){
      /*
	pfPDF->SetParameter(1,(double)par[2*i+1]);
	t_fcn+=par[2*i+2]*pfPDF->Eval(t);
      */
      // to boost generation speed
      int t_int=(int)(t+4*nSample-par[2*i+1])*fSeg;
      if ( t_int>=0 && t_int<8*nSample*fSeg ){
	t_fcn+=par[2*i+2]*pFuncArray[t_int];
      }
    }
    
    return t_fcn;
}

void PulseGenerator(void){
	
  // single photoelectron's filter output shape
  TCanvas* can = new TCanvas("can","",800,800);
  can->Divide(2,2);
  pfPDF=new TF1("pfPDF",agaus,0.,8.*nSample,4);
  pfPDF->SetParameters(normF,4.*nSample+100,0.04432,23.71);
  can->cd(1);
  pfPDF->Draw("");

  // to boost generation speed
  for ( int i=0;i<fSeg*8*nSample;i++ ){
    pFuncArray[i]=pfPDF->Eval(1./fSeg*i);
  }
	
	
  TF1* pePDF=new TF1("pePDF",AsymmetricGaussian,0.,140+110.,8);
  pePDF->SetParameter(0,1.);
  pePDF->SetParameter(1,peMean);
  pePDF->SetParameter(2,0.);  
  for ( int i=0;i<5;i++ ){ pePDF->SetParameter(3+i,pdfPar[i]*(1+corPar[i])); }
  can->cd(2);
  pePDF->Draw();

  TGraph* gr=new TGraph(nSample);
  for ( int i=0;i<nnpe;i++ ){
    for ( int j=0;j<nEvent;j++ ){
      gr->Set(0);
      double t_npe=0.;
      while ( t_npe==0. ){ t_npe=gRandom->Poisson(npe[i]); }
      
      TF1* waveform=new TF1("waveform",pulseFunc,-120,8.*nSample,2*t_npe+1);
      waveform->SetParameter(0,t_npe);
      
      for ( int k=0;k<t_npe;k++ ){
	waveform->SetParameter(2*k+2,gRandom->PoissonD(lambdaPMT)/lambdaPMT);
	waveform->SetParameter(2*k+1,pePDF->GetRandom());
      }
      double offset=8*(gRandom->Rndm()-0.5);
      for ( int k=0;k<nSample;k++ ){
	double d_tmp=waveform->Eval(k*8.-offset)+gRandom->Gaus(0.,gndSgm);
	gr->SetPoint(k,k*8,d_tmp);
      }
      can->cd(3);
      waveform->Draw();
      if ( j<12 ){ gr->Draw("PL"); }      
      can->Update();
      can->Modified();
      getchar();
      delete waveform;
    }
  }  
}
