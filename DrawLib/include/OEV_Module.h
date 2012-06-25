#ifndef OEV_MODULE__H__
#define OEV_MODULE__H__
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TObject.h>
#include <TH2Poly.h>
#include <TPolyLine.h>

class OEV_Module : public TObject {

 public:
  static const int numberOfOEV = 44;

  OEV_Module(const char*);
  virtual ~OEV_Module();
  void Init( const char* );
  void Reset( void );
  int Fill(int ibin);
  int Fill(int ibin, Double_t w);
  int Fill(double x, double y, double w);
  void ResetContent( void );
  void UpdateValue(void);
  void SetTitle(const char*);  
  void Draw(const char*);
  void DrawWithRange( double, double , const char* );
 private:

  Double_t OEV_xx[numberOfOEV];
  Double_t OEV_yy[numberOfOEV];  

  TH2Poly* OEV;
  TPolyLine* OEV_Poly[numberOfOEV];
  Double_t Dep[numberOfOEV];

  ClassDef(OEV_Module,0);

};
#endif
