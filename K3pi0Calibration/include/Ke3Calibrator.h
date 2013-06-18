#ifndef Ke3Calibrator_h__
#define Ke3Calibrator_h__ 

#include "gamma/Gamma.h"

class Ke3Calibrator{
 public:
  Ke3Calibrator(int numRequest);
  ~Ke3Calibrator();

  void fillGamma(Gamma const& gam,double referenceEnergy,double sigmaUser=-1);
  bool gaussianElimination();

  /// get functions 
  double getCalibrationFactor(int id){return m_calibFactor[id];}
  void getCalibrationFactors(double *calib){
    for(int i=0;i<s_maxArrSize;i++){
      calib[i] = m_calibFactor[i];
    }
  }

  static int const s_maxArrSize = 3000;

 private:
  int const m_numRequest;
  int m_counter[s_maxArrSize];
  double m_calibFactor[s_maxArrSize];
  double *m_matrix[s_maxArrSize];
};

#endif  //Ke3Calibrator_h__
