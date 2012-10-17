#ifndef E14EVENTBUILDER__H__
#define E14EVENTBUILDER__H__
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>

#include "E14ConvReader.h"
#include "E14ConvWriter.h"
#include "E14IDHandler.h"
#include "E14WaveFitter.h"
#include "IDHandler.h"
#include "E14WaveformAnalyzer.h"
#include "WaveformFitter.h"
#include "EnergyConverter.h"

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "Structs.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TText.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "CsI_Module.h"

static const int nCrateFeb     = 11; 
static const int nCVModule     = 10;
static const int nCosmicModule = 20; 
static const int CosmicArr[20] = {4 ,5 ,2 ,3 ,6 ,7 ,0 ,1 ,12,13,10,11,14,15,8 ,9 ,16,17,18,19};
  

class E14EventBuilder_V0 {
 public:
  E14EventBuilder_V0(TTree* trout, Int_t RunNumber);
  ~E14EventBuilder_V0();
  bool Init();
  bool InitEnvironment(); //Written
  bool InitIOFile();      //Written
  bool InitTemplate();    //Written
  bool InitTimeOffset();  //Written
  bool InitTrigger();     //Written

  int  TriggerDicision();
  int  AnalyzeCsIData();
  void Clear();
  void DrawEvent(TCanvas* canvas);
  void DrawResult(TCanvas* canvas);
  int  EventProcess(int ievent);
  int  LoopAll();
  
  // Utility Classes // 
  WaveformFitter*      wavFitter; //for Trigger waveform Fitting
  E14WaveFitter*       Fitter;    //for CsI  waveform Fitting 
  E14WaveformAnalyzer* wavAnalyzer;
  EnergyConverter*     Converter;
  IDHandler*           idHandler;
  // IO //
  E14ConvReader* conv[nCrateFeb];
  TFile*         tf[nCrateFeb];
  E14ConvWriter* wConv;
  // Template // 
  TGraph*       tempGr[2716];
  TSpline3*     tempSpl[2716]; 
  
  TTree* m_trOut;
  long   m_Entries;
  TH1D*  m_TimeClusterHist;
  TH2D*  m_EnergyTimeDistrib;

  int m_RunNumber;
  int m_EventNumber;

  TGraph*     grTrigger;
  TGraph*     grCsI;
  std::string ANALIBDIR;
  std::string CONVFILEDIR;
  std::string WAVEFILEDIR;
  std::string SUMFILEDIR;

  double TimeOffset[2716];

  /// Trigger;

  int CsiModuleID;
  int CC03ModuleID;
  int OEVModuleID;
  int CVModuleID;
  int CosmicModuleID;
  int LaserModuleID;
  
  int    TotalTriggerFlag;
  double COSMIC_THRESHOLD[20];
  double CosmicSignal[20];
  double CosmicTime[20];
  double CVSignal[10];
  double CVTime[10];
  int    LaserCFC[3];
  int    CVCFC[10][3];  
  int    CosmicCFC[20][3];

};

#endif 
