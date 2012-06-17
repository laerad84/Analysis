#ifndef CSIIMAGE__H__
#define CSIIMAGE__H__
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TObject.h>
#include <TH2D.h>
#include <TBox.h>
#include "IDHandler.h"

class CsIImage : public TObject {
public:
  CsIImage( IDHandler* );
  virtual ~CsIImage();
  
  void FillRandom( const char*, Int_t );
  const double  GetCrystalContent( int );
  void Draw(char* = "COLZ");
  void DrawWithRange(char* = "COLZ",Double_t = 0, Double_t = 1);

  bool GetBin( int, int&, int& );
  bool Fill( int, double = 1 );
  bool SetFillColor( int, Color_t );
  void SetTitle( const char* );
  void Reset( void );

private:
  TH2D* smallHist;
  TH2D* largeHist;
  TH2D* CC03Hist[4];


  static const int numberOfSmall   = 2240;
  static const int numberOfLarge   = 476;
  static const int numberOfCrystals = 2716;
  static const int numberOfCC03    = 32;

  IDHandler* handler;
  
  TBox* box[numberOfCrystals+numberOfCC03];

  ClassDef(CsIImage,0)
};
#endif
