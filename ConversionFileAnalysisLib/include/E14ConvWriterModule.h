#ifndef CONVWRITEMODULE__H__
#define CONVWRITEMODULE__H__

#include "TTree.h"
#include "Structs.h"
#include "GeneralTypes.h"

class E14ConvWriterModule {
 public:
  TTree* m_Tree;
  char   m_DetectorName[128];
  int    m_DetectorID;

  int    m_nDigi;
  short  m_ID[4096];
  short  m_FitHeight[4096];
  short  m_FitTime[4096];
  double m_FitShape[4096];
  double m_DeltaDiff[4096];

  double m_Pedestal[4096];
  double m_Signal[4096];
  double m_Time[4096];
  double m_HHTime[4096];

  double m_ParA[4096];
  double m_ParB[4096];

  double m_Chisq[4096];
  short  m_NDF[4096];

 public:
  E14ConvWriterModule(TTree*, char*);
  ~E14ConvWriterModule();
  bool InitData();
  bool SetBranchAddress();
  bool Branch();
  bool Close();
};

#endif
