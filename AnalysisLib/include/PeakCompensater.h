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
  ;
 public:
  TSpline3* m_spl[3];
  TGraph*   m_gr[3];

  PeakCompensater();
  ~PeakCompensater();
  double Compensate( int, double );
  void Draw(int ,char*);
 private:
    bool Init();
};
#endif //PEAKCOMPENSATER__H__
