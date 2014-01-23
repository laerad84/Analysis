#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "TH1.h"
#include "CsIPoly.h"
#include "TPad.h"

const int nMaxCrate = 20;
class L1TrigCounter {  
 private:
  
  int    m_nTrigCrate;
  double m_Threshold;
  double m_CrateCount[nMaxCrate];
  int    m_L1Map[2716];

  double m_Et;
  int    m_EtFlag;
  double m_EtThreshold;

 public:
  
  L1TrigCounter();
  ~L1TrigCounter();
  TH1D* hisTrigCounter; 
  CsIPoly* CsiL1Map;

  void   Reset(){ m_Et = 0; m_EtFlag = 0; hisTrigCounter->Reset();}
  int    Fill( int id, double value );
  void   SetThreshold(double threshold){ m_Threshold = threshold; }
  double GetThreshold(){ return m_Threshold;}
  bool   ReadMapFile( char* filename = "ch_map_CsI_L1.txt"); 
  int    GetCrateCnt();
  std::vector<int> GetTriggedCrate();
  std::vector<double> GetCount();
  double GetEtCount();
  void DrawTrig();
  void DrawTriggerMap();
};
