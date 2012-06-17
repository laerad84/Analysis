#ifndef E14Fill_h
#define E14Fill_h

#include <iostream>
#include <sstream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>

#include "E14RawData.h"
#include "E14DataRead2010.h"
#include "E14Mapper.h"
#include "E14Calibrator.h"

#include <GsimData/GsimEventData.h>
#include <GsimData/GsimGenParticleData.h>
#include <GsimData/GsimTrackData.h>
#include <GsimData/GsimDetectorEventData.h>
#include <GsimData/GsimDetectorHitData.h>

#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"


/**
 *  @class GsimDigiData
 *  @brief Digitization data.
 *  This class provides raw information.
 */
class E14Fill
{
 public:

  E14Fill();
  virtual ~E14Fill();

  // Fill raw data from *.dat.
  virtual void  Fill();
  // Fill raw data from *_conv.dat for 2010 semi online production.
  virtual void  FillFromConv();

  virtual void  SetDebugMode() { DebugMode = true; }
  virtual void  SetInputFile(  std::string s) { i_rootFile = s; }
  virtual void  SetInputDirectory(  std::string s) { id_rootFile = s; }
  virtual void  SetOutputFile(  std::string s) { o_rootFile = s; }
  virtual void  SetRunNumber(  int num ) { RunNumber = num; }

  virtual void  SetnCrate(  int num ) { nCrate = num; }

 private:
  
  virtual void  InitTree();
  virtual void  WriteTree();
  virtual void  PrintParameters();

 private:

  TFile*          hfile;
  TTree*          tree[64];  
  E14RawData      rawdata[64];
  E14DataRead2010 e14data[64];

  // Should be same for DetName, map, calib, gsim
  int             NDet;
  std::string     DetName[4];
  E14Mapper       map[4];
  E14Calibrator   calib[4];

  TFile *gsimFile;
  TTree *gsimTree;
  TTree *calibTree;
  TTree *mapTree;
  GsimDetectorEventData gsim[4];

  char InFileName[256];
  std::string i_rootFile;
  std::string id_rootFile;
  std::string o_rootFile;

  int  RunNumber;
  int  ADCThreshold;

  int  nCrate;
  bool INFILEJUDGE[64];

  int  readend;
  bool DebugMode;

  ClassDef(E14Fill,1)
};

#endif // E14Fill_h
