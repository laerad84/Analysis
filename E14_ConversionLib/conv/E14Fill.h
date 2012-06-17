#ifndef E14Fill_h
#define E14Fill_h

#include "E14RawData.h"
#include "E14DataReadVME.h"
#include "SemiOfflineHist.h"

#include "E14Mapper.h"
#include "E14Calibrator.h"

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
#include "TObject.h"
#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <sstream>
#include <iomanip>

class E14Fill : public TObject
{
 public:

  E14Fill();
  virtual ~E14Fill();

  // Fill raw data from *.dat.
  virtual void  Fill();
  // Fill raw data from *_conv.dat for 2010 semi online production.

  virtual void  SetDebugMode() { DebugMode = true; }
  virtual void  SetInputFile(  std::string s) { i_rootFile = s; }
  virtual void  SetOutputFile(  std::string s) { o_rootFile = s; }
  virtual void  SetInputDirectory(  std::string s) { id_rootFile = s; }
  virtual void  SetInputCompName(  std::string s) { ic_rootFile = s; }
  virtual void  SetCrateID(  int num) { CrateID = num; }

  virtual void  SetRunNumber(  int num ) { RunNumber = num; }
  virtual void  SetnCrate(  int num ) { nCrate = num; }

 private:
  
  virtual void  PrintParameters();

 private:

  TFile *hfile;
  //  TTree *tree[64];  
  TTree *tree;  
  E14RawData *rawdata[64];
  E14DataReadVME *e14data[64];
  SemiOfflineHist *sol;

  int NDet; // Should be same for DetName, map, calib
  std::string DetName[4];
  E14Mapper     map[4];
  E14Calibrator calib[4];

  TTree *calibTree;
  TTree *mapTree;

  char InFileName[256];
  char OutFileName[256];
  std::string i_rootFile;
  std::string id_rootFile;
  std::string ic_rootFile;
  std::string o_rootFile;
  int CrateID;

  int RunNumber;
  int ADCThreshold;

  int nCrate;
  bool INFILEJUDGE[64];

  int readend;
  bool DebugMode;

  ClassDef(E14Fill,1)
    };

#endif // E14Fill_h
