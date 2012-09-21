/* *********************************************************

class WaveformFitter : FADC waveform fitter class
  author    : Eito IWAI ( iwai # champ.hep.sci.osaka-u.ac.jp )

history :
   v0.0.1a   development version
   
********************************************************** */
/***********************************************************
 JWLee : Divide to Waveform.h and Waveform.cc
 ***********************************************************/

#ifndef Waveform_h
#define Waveform_h

#define WAVEFORM_VERSION "WAVEFORM v0.0.1a_20081012"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>

#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>

const double f_sigma=1.2;
//const double f_sigma=2.;

double AsymmetricGaussian(double* x, double* par);
double ScintiFunction(double* x, double* par);
double TypicalFunction(double* x, double* par);
double LnsFunction(double* x, double* par);

static TH1D* m_typShape[5][12][12];
static TF1* m_tf;

class WaveformFitter{
public:
	WaveformFitter(int Nsamples=64, bool fixed=false, int PedSmpl=8);
	WaveformFitter(int Nsamples, bool fixed, int PedSmpl, int fitType);
	virtual ~WaveformFitter(void);
  
	TF1* GetFunction(void){ return m_fitfunc; }
	void SetParameters(double p0, double p1, double p2=0., double p3=0., double p4=0.){ m_fitPar[0]=p0; m_fitPar[1]=p1; m_fitPar[2]=p2; m_fitPar[3]=p3; m_fitPar[4]=p4; }
	bool Approx(TGraph* gr);
	bool Fit(TGraph* gr);
  
protected:

	int    m_Nsamples;
	bool   m_fixed;
	TF1*   m_fitfunc;
	int    m_pedsmpl;
	
	double m_fitPar[5];
	int    m_fitType;
	TFile* m_f;
	
	//ClassDef(WaveformFitter,1)
};

#endif // Waveform_h
