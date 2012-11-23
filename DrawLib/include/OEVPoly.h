#ifndef OEVPoly__H__
#define OEVPoly__H__
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TObject.h>
#include <TH2Poly.h>
#include <TPolyLine.h>

class OEVPoly : public TH2Poly {

 public:
  static const int numberOfOEV = 44;

  OEVPoly();
  OEVPoly(const char* name, const char* title);
  virtual ~OEVPoly();
  void Reset( void ){ TH2Poly::Reset("");}
  void Init(  );
  int Fill(int ibin, Double_t w =1 );
  int Fill(double x, double y, double w =1 );

 private:

  Double_t OEV_xx[numberOfOEV];
  Double_t OEV_yy[numberOfOEV];  

  ClassDef(OEVPoly,0);

};
#endif
