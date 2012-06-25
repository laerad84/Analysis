#ifndef CC03_MODULE__H__
#define CC03_MODULE__H__
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TObject.h>
#include <TH2Poly.h>
#include <TPolyLine.h>
#include <TBox.h>

class CC03_Module : public TObject {
 public:
  static const int numberOfCC03 = 32;

  CC03_Module(const char*);
  virtual ~CC03_Module();
  void Init( const char* );
  void Reset( void );
  int Fill(int ibin);
  int Fill(int ibin, Double_t w);
  int Fill(double x, double y, double w );
  void SetTitle(const char*);  
  void Draw(const char*);
  void DrawWithRange( double, double, const char* );

 private:

  TH2Poly* CC03;
  TBox* CC03_Box[numberOfCC03];
  Double_t Dep[numberOfCC03];
  Double_t CC03_xx[numberOfCC03];
  Double_t CC03_yy[numberOfCC03];


  ClassDef(CC03_Module,0);
};
#endif
