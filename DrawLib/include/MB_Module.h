#ifndef MB_MODULE__H__
#define MB_MODULE__H__
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TObject.h>
#include <TH2poly.h>
#include <TPolyLine.h>
class MB_Module : public TObject {

 public:
  static const int numberOfMB = 32;

  MB_Module(const char*);
  virtual ~MB_Module();
  void Init( const char* );
  void Reset( void );
  int Fill(int ibin);
  int Fill(int ibin, Double_t w);
  int Fill(double x, double y, double w);
  void ResetContent( void );
  void UpdateValue(void);
  void SetTitle(const char*);  
  void Draw(const char*);
  void DrawWithRange( double, double, const char*);

 private:

  TH2Poly* MB;
  TPolyLine* MB_Poly[numberOfMB];
  Double_t Dep[numberOfMB];
  Double_t MB_xx[numberOfMB];
  Double_t MB_yy[numberOfMB];


  ClassDef(MB_Module,0);
};
#endif
