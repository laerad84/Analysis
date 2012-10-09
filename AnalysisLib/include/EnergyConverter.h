#ifndef ENERGYCONVERTER__H__
#define ENERGYCONVERTER__H__
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include "PeakCompensater.h"
#include "TFile.h"
#include "TTree.h"

class EnergyConverter {
 private:
  double m_CalibrationConstant[4096];
  int    m_nChannel;
  bool   m_fCalibrationConstant[4096];
  bool   m_fInit;
  bool   m_fRead;
  std::string m_DetectorName;
  std::string m_CalFilename; 

 public:
  PeakCompensater* m_Compensater;  
  EnergyConverter();
  EnergyConverter( const char* DetectorName, int nChannel = 2716 ); 
  ~EnergyConverter();  
  virtual bool   ReadCalibrationTextFile( const char* );
  virtual bool   ReadCalibrationRootFile( const char* );
  virtual bool   IsGoodChannel         ( int )         const ;
  virtual double GetCalibrationConstant( int )         const ;
  virtual double ConvertToEnergy       ( int, double ) const ; 
  
 private:
  virtual void Init();
  virtual void Reset();
  
};

#endif

