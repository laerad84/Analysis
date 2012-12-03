#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <string>

#include "EventTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TApplication.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TProfile.h"
#include "TRandom.h"

#include "CsIPoly.h"
#include "IDHandler.h"
#include "PulseGenerator.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include "gamma/Gamma.h"
#include "cluster/Cluster.h"

#include "gamma/GammaFinder.h"
#include "cluster/ClusterFinder.h"

#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "User_Function.h"

class EDepositAnalysis {
 public:
  Int_t nDigi;
  Int_t ID[2716];//nDigi
  Double_t Energy[2716];//nDigi
  Double_t SignalTime[2716];//nDigi

  CsIPoly*        CsIEne;
  PulseGenerator* gen;
  TChain*         chain; 
  EventTree*      trin;
  TFile*          OutputFile;
  TTree*          OutputTree;
  
  std::string m_InputFilename;
  std::string m_OutputFilename;


  EDepositAnalysis(const  char* InputFilename,const  char* OutputFilename, Int_t InjectionDirection = 0);
  virtual ~EDepositAnalysis();
  void Init();
  int  EventProcess(int ievent );
  int  Loop( int Entries = 0);
  
  static const double Speed_of_Signal;
  void Export();
  void Close();

  E14GNAnaDataContainer* data;
  ClusterFinder*         clusterFinder;
  void DrawEvent();
 private:
  int  m_ZDirection;
  void PrepareIO();
  int  PrepareDump();
  int  ResetDump();
  void ConvertPosition(const  Double_t Radius, const Double_t Theta, const Double_t x, const Double_t y, Double_t& nx , Double_t& ny);

};
