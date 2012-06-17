#ifndef SemiOnlinePlot_h
#define SemiOnlinePlot_h

#include <iostream>
#include <fstream>
#include <math.h>
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphErrors.h"

#include "EventTree.h"
#include "CsIImage.h"
#include "E14Mapper.h"

using namespace std;

class SemiOnlinePlot
{
 public:
  SemiOnlinePlot();
  virtual ~SemiOnlinePlot();

 public:
  virtual void Init( int runID );
  virtual void InitValue();

  virtual void Fill( int event );

  virtual void RunDisplay( );
  virtual void EventDisplay( int event );
  virtual void DrawAll();

  virtual void SetIntegratedADC( int idet, int ich, Double_t Ene );
  virtual void SetPeakHeight( int idet, int ich, Short_t Peak );
  virtual void SetDetEne( int idet, int ich, Double_t Ene );
  virtual void SetEtSum( int icrate, int ifadc, Int_t *ADC );

 private:
  int NDet; // Should be same for DetName, map, calib
  std::string DetName[7];
  E14Mapper     map[7];

  int DetNumber[7];
  int NOfIntegratedADC[7];
  int NOfDetEne[7];
  Double_t DetIntegratedADC[7][2716];
  Short_t  DetPeakHeight[7][2716];
  Double_t DetEne[7][2716];
  Int_t    EtSum[10][16][48];  // 10 crates, 16 FADC, 48 samples.

  CsIImage *adcDistribution;
  CsIImage *adcDistributionLaser;
  CsIImage *energyDistribution;
  CsIImage *hitPosition;
  CsIImage *eventDisplay;

  TFile *hfile;
  TH1F *hIntegratedADCSum;
  TH1F *hEneSum;
  TH1F *hPinDiode[4];
  TH1F *hLaserRef[2716];       // Laser / PINdiode
  TH1F *hLaserRefNorm[2716];   // Laser / PINdiode / normalized
  TH1F *hPedNoise[2716];       // Pedestal noise
  TH2D *hEtSum;
  TH2D *hEtSumLaser;
  TH1F *hCV[10];

  TGraphErrors *gLaserStabilityS[35];
  TGraphErrors *gLaserStabilityL[4];

 private:
  virtual void ReadLaserRefValue();

  int RunNum;
  int threshold;
  int EneThreshold;
  int TargetCH;

  char psFileName1[256];  
  char psFileName2[256];  
  IDHandler* handler;
  TCanvas *c1;
  TCanvas *c2;
  TLatex *ttex;
  TLine *tline;

  float LaserRefMean[2716];   // Laser / PINdiode
  float LaserRefSigma[2716];
 
};

#endif
