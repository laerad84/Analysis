#ifndef CsI_MODULE__H__
#define CsI_MODULE__H__
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TObject.h>
#include <TH2Poly.h>
#include <TPolyLine.h>
#include <TBox.h>

class CsI_Module : public TObject {
 public:
  static const int numberOfCsI = 2716;
  static const int numberOfSmallCsI = 2240;
  static const int numberOfLargeCsI = 576;

  CsI_Module(const char*);
  virtual ~CsI_Module();
  void Init( const char* );
  void Reset( void );
  int Fill(int ibin);
  int Fill(int ibin, Double_t w);
  int Fill(double x, double y, double w);
  void SetTitle(const char*);  
  void Draw(const char*);
  void DrawWithRange( double , double , const char*);

  TH2Poly* CsI;
 private:
  TBox* CsI_Box[numberOfCsI];
  Double_t Dep[numberOfCsI];
  Double_t CsI_xx[numberOfCsI];
  Double_t CsI_yy[numberOfCsI];


  ClassDef(CsI_Module,0);
};
#endif
