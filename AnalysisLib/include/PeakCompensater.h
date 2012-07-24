#ifndef PEAKCOMPENSATER__H__
#define PEAKCOMPENSATER__H__
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include "TSpline.h"
#include "TGraph.h"

class PeakCompensater
{
 private:
  TSpline3* m_spl[3];
  TGraph*   m_gr[3];
 public:
  PeakCompensater();
  ~PeakCompensater();
  double Compensate( int, double );
 private:
    bool Init();
};
#endif //PEAKCOMPENSATER__H__
