#ifndef E14VMESum_h
#define E14VMESum_h

#include "EventTree.h"
#include "E14Mapper.h"
#include "E14Calibrator.h"
#include "SemiOnlinePlot.h"

#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>

#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"
#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <sstream>
#include <iomanip>

class E14VMESum
{
 public:

  E14VMESum();
  virtual ~E14VMESum();

  // Fill raw data from *_conv.root.
  virtual void  Fill();

  virtual void  SetDebugMode() { DebugMode = true; }
  virtual int  SetRunNumber(  int num );

 private:
  
  virtual void  InitTree();
  virtual void  InitStorage();
  virtual void  WriteTree();
  virtual void  PrintParameters();  
  virtual bool  CheckTimeStampSync();

 private:

  SemiOnlinePlot *sol;

  int NDet; // Should be same for DetName, map, calib
  std::string DetName[7];
  E14Mapper     map[7];
  E14Calibrator calib[7];
  EventTree *eventTree[12];

  TFile *hfile;
  TTree *tree;
  TTree *calibTree;
  TTree *mapTree;

  // Data container for 2012 run
  Int_t EventNo;
  Short_t SpillNo;
  Int_t TimeStamp;
  Short_t Error;
  Int_t DetNumber[7];
  Int_t DetModID[7][2716];
  Double_t DetEne[7][2716];
  Double_t DetIntegratedADC[7][2716];
  Double_t DetTime[7][2716];  

  int nCrate;
  int ConvCrateID[12];
  char InFileName[256];
  bool INFILEJUDGE[12];
  char ConvCompName[12][256];

  int RunNumber;
  float EnergyThreshold;


  int readend;
  bool DebugMode;

  ClassDef(E14VMESum,1)
    };

#endif // E14VMESum_h
