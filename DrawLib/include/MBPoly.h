#ifndef MBPoly__H__
#define MBPoly__H__
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TObject.h>
#include <TH2Poly.h>
#include <TPolyLine.h>
class MBPoly : public TH2Poly {

 public:
  static const int numberOfMB = 32;

  MBPoly();
  MBPoly(const char* name, const char* title);
  virtual ~MBPoly();
  void Reset( void ){
    TH2Poly::Reset("");    
  }
  void Init( void );
  int Fill(int ibin, Double_t w = 1);
  int Fill(double x, double y, double w =1);

 private:

  Double_t MB_xx[numberOfMB];
  Double_t MB_yy[numberOfMB];

  ClassDef(MBPoly,0);
};
#endif
