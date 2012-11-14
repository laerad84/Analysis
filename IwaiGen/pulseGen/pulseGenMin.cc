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


const double flt_a=0.07214;
const double flt_b=27.8;
const double absLY=12.7; // pe/MeV
const double gndSgm=2.4; // [cnt] raw value
//const bool fit_fixed=kTRUE;
const bool fit_fixed=kFALSE;

// p3-p7, valid in range [p1-8,p1+160] (p1:peak time)
const double pdfPar[5]={3.28113, 0.541526, -0.00582465, 4.53179e-05, -1.16776e-07};

// correction factor
const double corPar[5]={-0.197912, 0.0695365, 0.161864, 0.193783, 0.232639};

const int fSeg=10; // ex. fSeg=10 means 1/10=0.1ns
double pFuncArray[8*nSample*fSeg];

//const double lambdaPMT=3.952; // RV580UV@1800V
//const double lambdaPMT=8.283; // H7195@2100V
const double lambdaPMT=14.; // R4125@1400V
//const double lambdaPMT=17.; // R4125@1900V

const int nnpe=8;
const int npe[nnpe]={40,80,160,320,640,1200,2000,4000};

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

void pulseGenMin(void){
  TCanvas* can = new TCanvas("can","",800,800);
  can->Divide(2,2);
	
  // single photoelectron's filter output shape
  pfPDF=new TF1("pfPDF",agaus,0.,8.*nSample,4);
	// MIP 13.7MeV-240cnt, 16000pe-height13360
	// 240/13.7/12.7*16000/13360 = 1.38*1.2 =1.652
  pfPDF->SetParameters(1.652,4.*nSample,0.04432,23.71);
  can->cd(1);
  pfPDF->Draw();
  // to boost generation speed
  for ( int i=0;i<fSeg*8*nSample;i++ ){
    pFuncArray[i]=pfPDF->Eval(1./fSeg*i);
  }
  
  //TUtil* u=new TUtil;

  WaveformFitter* wf=new WaveformFitter(nSample-2,fit_fixed,8);
  TF1* f1=new TF1("f1","pol1",0.,600.);
  f1->SetLineColor(4);
  
  TH1D *time[nnpe], *time_cfm[nnpe], *energy[nnpe], *sum[nnpe];
  for ( int i=0;i<nnpe;i++ ){
    time[i]=new TH1D(Form("time_%d",npe[i]),"",100,172-120/TMath::Sqrt(npe[i]),172+120/TMath::Sqrt(npe[i]));
    time_cfm[i]=new TH1D(Form("time_cfm_%d",npe[i]),"",100,172-120/TMath::Sqrt(npe[i]),172+120/TMath::Sqrt(npe[i]));
    energy[i]=new TH1D(Form("energy_%d",i),"",100,(1.38-1*8/TMath::Sqrt(npe[i]))*npe[i],(1.38+1*8/TMath::Sqrt(npe[i]))*npe[i]);
    sum[i]=new TH1D(Form("sum_%d",i),"",100,9.*(1.38-1*8/TMath::Sqrt(npe[i]))*npe[i],9.*(1.38+1*8/TMath::Sqrt(npe[i]))*npe[i]);
  }
  
  TH1F* ahist=new TH1F("ahist","",100,0.06,0.09);
  TH1F* bhist=new TH1F("bhist","",100,27.,29.);
  
  TF1* pePDF=new TF1("pePDF",AsymmetricGaussian,0.,140+110.,8);
  pePDF->SetParameter(0,1.);
  pePDF->SetParameter(1,160.);
  pePDF->SetParameter(2,0.);
  can->cd(2);
  pePDF->Draw();

  for ( int i=0;i<5;i++ ){ pePDF->SetParameter(3+i,pdfPar[i]*(1+corPar[i])); }
  TGraph* gr=new TGraph(nSample);
  for ( int i=0;i<nnpe;i++ ){
    for ( int j=0;j<nEvent;j++ ){
      gr->Set(0);
      double t_npe=0.;
      while ( t_npe==0. ){ t_npe=gRandom->Poisson(npe[i]); }
      
      TF1* waveform=new TF1("waveform",pulseFunc,0.,8.*nSample,2*t_npe+1);
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
      
      //if ( j<12 ){ u->Draw(gr,"AP"); }
      if( j < 12 ){ 
	can->cd( 3);
	waveform->Draw();
	gr->SetMarkerStyle(6);
	gr->Draw("P");
      }
      can->Update();
      can->Modified();
      getchar();
      delete waveform;
    }
  }
  
  //delete u;
}
