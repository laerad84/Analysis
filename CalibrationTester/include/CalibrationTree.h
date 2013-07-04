#ifndef CALIBRATIONTREE__H__
#define CALIBRATIONTREE__H__
#include "Structs.h"
#include "TChain.h"

class CalibrationTree {
 public:
  TTree* m_Tree;
  int FlagKL_prefit;
  int FlagKL[6];
  int FlagCluster[6];
  int FlagCalibrated[6];
  int CorrID[6];
  double Corr[6];
  double CorrE[6];
  double Ratio[6];
  double SecondRatio[6];
  double GammaEnergy[6];
  double GammaSigma[6];
  double chisq[6];
  double LeadingEnergy[6];
  double LeadingHeight[6];
  int    LeadingChID[6];
  int    nCalibrated;
  
 public:  
  struct CalibrationData m_calibData;
  CalibrationTree();
  ~CalibrationTree();
  void Branch(TTree* tr);
  void SetBranchAddress(TTree* tr);
  int GetEntries();
  int GetEntry(int ientry);
  int InitValue();
};
#endif
