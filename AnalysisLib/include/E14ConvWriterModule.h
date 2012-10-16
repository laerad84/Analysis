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
  // Overall Information

  int    m_nOverFlow;
  int    m_nUnderFlow;
  int    m_IDOverFlow[2716];
  int    m_IDUnderFlow[2716];

  //////////////////////////////////////////////////////
  // Channel Information 
  /////////////////////////////////////////////////////
  double m_TotalEnergy;
  int    m_nTimeCluster;
  int    m_TotalEnergyInTimeCluster[32];
  double m_TimeClusterHead[32];
  double m_TimeClusterTail[32];

  int    m_nDigi;  
  short  m_ID[2716];  
  int    m_TimeClusterID[2716];
  short  m_FitHeight[2716];
  short  m_FitTime[2716];
  double m_FitShape[2716];
  double m_DeltaDiff[2716];
  double m_ADC[2716];
  double m_FitADC[2716];

  double m_Pedestal[2716];
  double m_Signal[2716];
  double m_Time[2716];
  double m_HHTime[2716];
  double m_Energy[2716];

  double m_ParA[2716];
  double m_ParB[2716];

  double m_Chisq[2716];
  short  m_NDF[2716];

  // For WavAnalysis Parameter // 
  double m_wav_SlopeDelta[2716];
  double m_wav_Height[2716];
  double m_wav_PeakTime[2716];
  double m_wav_Pedestal[2716];
  double m_wav_Width[2716];
  double m_wav_FrontHalfTime[2716];
  double m_wav_RearHalfTime[2716];

  // For Result of Fitting //
  double m_Fit_Pedestal[2716];
  double m_Fit_Time[2716];
  double m_Fit_Height[2716];
  double m_Fit_HHTime[2716];
  double m_Fit_ChisqNDF[2716];
  double m_Fit_ChisqPed[2716];
  double m_Fit_ChisqFront[2716];
  double m_Fit_ChisqRear[2716];
  double m_Fit_ChisqTail[2716];
  // m_Fit_ChisqNDF   : ChisqNDF of Fit Region( peak -120 ~ peak + 50 )
  // m_Fit_ChisqPed   : ChisqNDF of peak - 144 ~ peak - 96
  // m_Fit_ChisqFront : ChisqNDF of peak - 96  ~ peak 
  // m_Fit_ChisqRear  : ChisqNDF of peak       ~ peak + 96
  // m_Fit_ChisqTail  : ChisqNDF of peak + 96  ~ peak + 240 

  // For Converted Data
  double m_Conv_Energy[2716];
  double m_Conv_Time[2716];
  // m_Conv_Energy : Height Adjusted Energy ( same with normal Energy )
  // m_Conv_Time   : Height Adjusted Time 

 public:
  E14ConvWriterModule(TTree*, char*);
  ~E14ConvWriterModule();
  bool InitData();
  bool SetBranchAddress();
  bool Branch();
  bool Close();

};

#endif
