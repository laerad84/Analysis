#ifndef CsIImage_h
#define CsIImage_h

//============================================================================
// Name        : CsIImage.cpp
// Author      : Takahiko Masuda
// Version     : 1.0
// Copyright   : taka
// Date        : 2010/10/25
//============================================================================

#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TH2D.h"
#include "TBox.h"
#include "IDHandler.h"


class CsIImage {
public:
  CsIImage( IDHandler* );
  virtual ~CsIImage();
  
  void FillRandom( const char*, Int_t );
  const double  GetCrystalContent( int );
  void Draw();

  bool GetBin( int, int&, int& );
  bool Fill( int, double = 1 );
  bool SetFillColor( int, Color_t );
  void SetTitle( const char* );
  void Reset( void );

	
private:
  TH2D* smallHist;
  TH2D* largeHist;
  TH2D* CC03Hist[4];

  int numberOfSmall;
  int numberOfLarge;
  int numberOfCrystals;
  int numberOfCC03;

  IDHandler* handler;
  
  TBox* box[2716*3];
  
};


#endif
