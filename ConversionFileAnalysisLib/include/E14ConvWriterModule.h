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

  short  m_nDigi;
  short  m_ID[4096];//m_nDigi
  short  m_FitHeight[4096];//m_nDigi
  short  m_FitTime[4096];//m_nDigi
  double m_FitShape[4096];//m_nDigi
  double m_DeltaDiff[4096];//m_nDigi

  double m_Pedestal[4096];//m_nDigi
  double m_Signal[4096];//m_nDigi
  double m_Time[4096];//m_nDigi
  double m_HHTime[4096];//m_nDigi

  double m_ParA[4096];//m_nDigi
  double m_ParB[4096];//m_nDigi

  double m_Chisq[4096];//m_nDigi
  short  m_NDF[4096];//m_nDigi

 public:
  E14ConvWriterModule(TTree*, char*);
  ~E14ConvWriterModule();
  bool InitData();
  bool SetBranchAddress();
  bool Branch();
  bool Close();
};

#endif
