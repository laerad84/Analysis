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
#include "TFile.h"
#include <fstream>
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
  PeakCompensater( int );
  ~PeakCompensater();
  int     m_version;
  double  m_map[2716][4];
  virtual void   SetMap();
  virtual double Compensate( int, double );
  virtual double InvCompensate( int,double );
  virtual void   Draw(int ,char*);

 private:
  virtual bool Init( int version =0);
};
#endif //PEAKCOMPENSATER__H__
