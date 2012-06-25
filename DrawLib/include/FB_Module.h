#ifndef FB_MODULE__H__
#define FB_MODULE__H__
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TObject.h>
#include <TH2poly.h>
#include <TPolyLine.h>

class FB_Module : public TObject {

 public:
  static const int numberOfFB = 16;

  FB_Module(const char*);
  virtual ~FB_Module();
  void Init( const char* );
  void Reset( void );
  int Fill(int ibin);
  int Fill(int ibin, Double_t w);
  int Fill(double x, double y, double w);
  void ResetContent( void );
  void UpdateValue(void);
  void SetTitle(const char*);  
  void Draw(const char*);
  void DrawWithRange( double , double , const char* );

 private:

  TH2Poly* FB;
  TPolyLine* FB_Poly[numberOfFB];
  Double_t Dep[numberOfFB];
  Double_t FB_xx[numberOfFB];
  Double_t FB_yy[numberOfFB];

  ClassDef(FB_Module,0);
};
#endif
