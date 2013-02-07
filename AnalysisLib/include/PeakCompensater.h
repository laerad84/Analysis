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
  TSpline3* m_splInv[3];
  TGraph*   m_grInv[3];
  PeakCompensater();
  ~PeakCompensater();
  virtual double Compensate( int, double );
  virtual double InvCompensate( int,double );
  virtual void   Draw(int ,char*);
 private:
  virtual bool Init();
};
#endif //PEAKCOMPENSATER__H__
