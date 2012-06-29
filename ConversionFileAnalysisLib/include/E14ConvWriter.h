#ifndef E14CONVWRITER__H__
#define E14CONVWRITER__H__

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>

#include "TTree.h"
#include "TGraph.h"
#include "GeneralTypes.h"
#include "E14MapReader.h"
#include "Structs.h"
#include "E14ConvWriterModule.h"


class E14ConvWriter {
 private:
  static const int nMaxModule = 32; 
  bool bInitialize;    
  TTree* m_tr;
  int m_nModule;
  std::string m_mapFilename;
  
 public:
  E14MapReader*        map;
  struct MapStruct     ModMap[32];   
  E14ConvWriterModule* mod[32];
  
  E14ConvWriter( char* , TTree* );
  ~E14ConvWriter();

  bool  AddModule(char*);
  bool  SetMap();
  bool  Set();
  int   GetNmodule();
  bool  Branch();
  bool  SetBranchAddress();  
  bool  InitData();
};

#endif
