#ifndef CALIBRATION__H__
#define CALIBRATION__H__
#include <cstdio>
#include <cstdlib>
#include <list>
#include <vector>
#include <iostream>
#include <iomanip>

#include "TH1.h"
#include "TROOT.h"
#include "TTree.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "gamma/Gamma.h"
#include "pi0/Pi0.h"
#include "klong/Klong.h"
#include "klong/RecKlong.h"

#include "Structs.h"

#include "CalibrationTree.h"

#define MAX_ITERACTIONS 30
#define MASS_KL         497.648
#define MASS_PI0        134.9766

//const int N_CSI = 2716;
//const int N_EDGE_CSI = 148;


class Calibration {
 private:
  int    m_FlagKL_prefit; // Flag of KL_prefit
  int    m_FlagKL[6];     // Flag of KL(Number of Gamma)
  int    m_FlagCluster[6];// Flag of Cluster(Number of Gamma)
  int    m_FlagCalibrated[6];

  int    m_CorrID[6];     // Calibrated Channel ID of EachGamma(Number of Gamma, leading Channel)
  double m_CorrE[6];      // Calibrated Channel Energy of Each Gamma( Number of Gamma )
  double m_Corr[6];       // Calibration Factor of Each Gamma(Number of Gamma), Valid only if m_FlagXX == 0 
  double m_Ratio[6];      // Energy Raito leading channel / Energy of Gamma 
  double m_SecondRatio[6];// Energy Ratio second leading Channel / Energy of Gamma
  double m_GammaEnergy[6];     // Energy of Gamma
  double m_GammaSigma[6]; // Sigma of Gamma Energy
  double m_chisq[6];      // Chisquare of Gamma
  int    m_nCalibrated;   // Number of Calibrated KL in 6 Gamma 
  Klong  KL;
  Klong  KL_prefit;  
  std::vector<Klong>  KL_Vec;

 public:
  const static double klmassdiff = 10.;//5. initial
  const static double pimassdiff = 6.;//3. initial
  const static double threshold  = 0.2;
  
  Calibration();
  ~Calibration();

  int InitValue();
  int InitValue(Klong &KL_init);
  int InitValue(std::vector<Klong> &KL_VecInit);
  int CalEnergy_idv(std::vector<Klong> &KL_VecInit);
  int Calibration_Gamma();
  int SortGamma(int idx , int GammaIndex);// idx 0~2, GammaIndex 0~1
  
  // return Error by int number; 0 is Good.
  int SelectKL_prefit();
  int SelectKL(int idx, int GammaIndex);
  int SelectCluster( int idx, int GammaIndex);// idx 0~2, GammaIndex 0-1
  double Calibrate();//return chisq
  void GetResult(CalibrationTree &calib);
  void Print();
};

#endif
