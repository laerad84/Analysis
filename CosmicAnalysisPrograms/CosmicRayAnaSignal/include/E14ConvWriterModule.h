#ifndef CONVWRITEMODULE__H__
#define CONVWRITEMODULE__H__

#include "TTree.h"
#include "Structs.h"
#include "GeneralTypes.h"

class E14ConvWriterModule {
 public:
  TTree* m_Tree;
  char   m_DetectorName[128];
  /*
  int    m_nDigi;
  int    m_Fit[4096];
  int    m_ID[4096];
  */
  int    m_nDigi;
  int    m_Fit[4096];
  short  m_ID[4096];

  double m_Pedestal[4096];
  double m_Signal[4096];
  double m_Timing[4096];
  double m_FitTiming[4096];
  double m_HHTiming[4096];
  double m_SplTiming[4096];
  double m_ParA[4096];
  double m_ParB[4096];
  double m_Ene[4096];
 public:
  E14ConvWriterModule(TTree*, char*);
  ~E14ConvWriterModule();
  bool InitData();
  bool SetBranchAddress();
  bool Branch();
  bool Close();
};

#endif
