#ifndef E14CONVWRITER__H__
#define E14CONVWRITER__H__

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>

#include "TGraph.h"
#include "GeneralTypes.h"
#include "E14MapReader.h"
#include "Structs.h"

class E14ConvWriter {
 private:
  static const nMaxModule = 32; 
  bool bInitialize;    
  TTree* m_tr;
  TFile* m_tf;
  int m_nModule;
  std::string m_mapFilename;
  std::string m_outFilename;

 public:
  E14MapReader* map;
  struct MapStruct ModMap[32];   

  E14ConvWriter(char*,TFile*);
  ~E14ConvWriter();

  bool  InitMap();
  bool  Branch(char*);
  bool  MakeBranchAll();
  bool  OpenMapFile();
  bool  ReadMap();
  bool  Write();
  

};

#endif
