#ifndef CSIPOLY__H__
#define CSIPOLY__H__

#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TObject.h"
#include "TPolyLine.h"
#include "TBox.h"
#include "TH2Poly.h"

class CsIPoly : public TH2Poly {
 public:
  static const int numberOfCsI = 2716;
  static const int numberOfSmallCsI = 2240;
  static const int numberOfLargeCsI = 476;
  
  CsIPoly();
  CsIPoly( const char* name, const char* title );
  virtual ~CsIPoly();
  virtual void  Reset(){
    TH2Poly::Reset("");
  }
  void Init();
  int Fill( int ibin );
  int Fill( int ibin, double w );
  int Fill( double x, double y );
  int Fill( double x, double y, double w );
  
 private:
  Double_t Dep[numberOfCsI];
  Double_t CsI_xx[numberOfCsI];
  Double_t CsI_yy[numberOfCsI];
  
  ClassDef(CsIPoly,0)
};
#endif
