#ifndef CALIBRATION_STRUCTS__H__
#define CALIBRATION_STRUCTS__H__
#include <string>

struct CalibrationData {
  int FlagKL_prefit;
  int FlagKL[6];
  int FlagCluster[6];
  int FlagCalibrated[6];
  
  int CorrID[6];
  double Corr[6];
  double CorrE[6];
  double SecondRatio[6];
  double Ratio[6];
  double GammaEnergy[6];
  double GammaSigma[6];
  double chisq[6];
  double nCalibrated;

};

struct MapStruct {
  int modID;
  int nMod;
  int Map[4096][3];
  std::string name;
};


#endif //CALIBRATION_STRUCTS__H__
