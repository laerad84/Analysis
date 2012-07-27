#ifndef E14CONVWRITER__H__
#define E14CONVWRITER__H__

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <cstring>

#include "TTree.h"
#include "TGraph.h"
#include "GeneralTypes.h"
#include "E14MapReader.h"
#include "Structs.h"
#include "E14ConvWriterModule.h"


class E14ConvWriter {
 private:
  static const int nMaxModule = 32; 
  TTree* m_tr;
  int    m_nModule;
  bool   bInitialize;    
  std::string m_mapFilename;
  std::list<std::string> m_modList;
  int    m_RunNo;
  int    m_EventNo;
  int    m_TrigFlag;
  int    m_CosmicTrig;
  int    m_LaserTrig;
  int    m_CVTrig;
  int    m_CosmicTrigFlagUp;
  int    m_CosmicTrigFlagDn;
  int    m_CVTrigFlag;
  int    m_LaserTrigFlag;
  int    m_ModID_Cocmic;
  int    m_ModID_CV;
  int    m_ModID_Laser;

  TGraph* m_gr;
 public:
  E14MapReader*        map;
  struct MapStruct     ModMap[32];   
  E14ConvWriterModule* mod[32];
  
  E14ConvWriter( char* , TTree* );
  ~E14ConvWriter();

  bool  AddModule(char*);
  bool  Set();
  bool  SetMap();
  bool  Branch();
  bool  SetBranchAddress();  
  bool  InitTriggerFlag();
  bool  TrigJudgement();
  bool  ScanMod(char*);
  int   GetNmodule();
  int   GetModuleID( char* );
  int   Fit( int );
  int   FitAll( );
};

#endif
