#ifndef FBPoly__H__
#define FBPoly__H__
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TObject.h>
#include <TH2Poly.h>
#include <TPolyLine.h>

class FBPoly : public TH2Poly {

 public:
  static const int numberOfFB = 16;

  FBPoly( );
  FBPoly( const char* name , const  char* title );
  virtual ~FBPoly();
  void Init( void );
  void Reset( void ){
    TH2Poly::Reset("");
  }
  int Fill(int ibin, Double_t w = 1);
  int Fill(double x, double y, double w = 1 );

 private:

  Double_t FB_xx[numberOfFB];
  Double_t FB_yy[numberOfFB];

  ClassDef(FBPoly,0);
};
#endif
