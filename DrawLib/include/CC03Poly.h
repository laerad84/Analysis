#ifndef CC03POLY__H__
#define CC03POLY__H__
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TObject.h"
#include "TH2Poly.h"
#include "TPolyLine.h"
#include "TBox.h"

class CC03Poly : public TH2Poly {
 public:
  static const int numberOfCC03 = 32;

  CC03Poly();
  CC03Poly( const char* name, const char* title );
  virtual ~CC03Poly();
  virtual void Reset( void ){
    TH2Poly::Reset("");
  }
  void Init();

  int Fill(int ibin, Double_t w =1);
  int Fill(double x, double y, double w = 1);


 private:

  Double_t Dep[numberOfCC03];
  Double_t CC03_xx[numberOfCC03];
  Double_t CC03_yy[numberOfCC03];

  ClassDef(CC03Poly,0);
};
#endif
