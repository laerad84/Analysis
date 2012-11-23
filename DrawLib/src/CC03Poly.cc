#ifndef CC03POLY__H__
#include "CC03Poly.h"
#endif

#include "TVector2.h"
#include "TMath.h"
#if !defined(__CINT__)
ClassImp(CC03Poly)
#endif

CC03Poly::CC03Poly(){
  Initialize( -100,100,-100,100, 25,25 );
  SetName("NoName");
  SetName("NoTitle");
  SetFloat();
  Init();
}
CC03Poly::CC03Poly( const char* name, const char* title){
  Initialize( -100,100,-100,100, 25,25 );
  SetName( name );
  SetName( title );
  SetFloat(kFALSE);
  Init();
}
CC03Poly::~CC03Poly(){
  ;
}
void    CC03Poly::Init( ){
  Double_t CC03_initx[numberOfCC03] = {  90.0  ,  90.0,  90.0,  90.0,  90.0,   90.0,   90.0,  90.0,
					 88.75 , 66.25, 43.75, 21.25, -1.25, -23.75, -46.25, -68.75,
					-90.0  , -90.0, -90.0, -90.0, -90.0,  -90.0,  -90.0,  -90.0,
					-88.75 ,-66.25,-43.75,-21.25,  1.25,  23.75,  46.25,  68.75};
  
  Double_t CC03_inity[numberOfCC03] = {  88.75 ,66.25 ,43.75  ,21.25 ,-1.25 ,-23.75,-46.25,-68.75,
					-90.0 ,-90.0 ,-90.0  ,-90.0 ,-90.0 ,-90.0 ,-90.0 ,-90.0 ,
					-88.75,-66.25,-43.75 ,-21.25,1.25  ,23.75 ,46.25 ,68.75 ,
					90.0  ,90.0  ,90.0   ,90.0  ,90.0  ,90.0  ,90.0  ,90.0  };
  for( int i = 0; i < numberOfCC03; i++){
    CC03_xx[i] = CC03_initx[(i+numberOfCC03-4)%numberOfCC03];
    CC03_yy[i] = CC03_inity[(i+numberOfCC03-4)%numberOfCC03];
  }

  for( int i = 0; i< numberOfCC03; i++ ){
  //for( int k = 0; k < numberOfCC03; k++){
    int k = (i+numberOfCC03-4)%numberOfCC03;
    Double_t x[5]={0};
    Double_t y[5]={0};
    switch( k/8 ){
    case 0:
    case 2:
      x[0] = CC03_xx[i]-10;
      x[1] = CC03_xx[i]+10;
      x[2] = CC03_xx[i]+10;
      x[3] = CC03_xx[i]-10;
      x[4] = CC03_xx[i]-10;
      y[0] = CC03_yy[i]-11.25;
      y[1] = CC03_yy[i]-11.25;
      y[2] = CC03_yy[i]+11.25;
      y[3] = CC03_yy[i]+11.25;
      y[4] = CC03_yy[i]-11.25;      
      break;
    case 1:
    case 3:
      x[0] = CC03_xx[i]-11.25;
      x[1] = CC03_xx[i]+11.25;
      x[2] = CC03_xx[i]+11.25;
      x[3] = CC03_xx[i]-11.25;
      x[4] = CC03_xx[i]-11.25;
      y[0] = CC03_yy[i]-10;
      y[1] = CC03_yy[i]-10;
      y[2] = CC03_yy[i]+10;
      y[3] = CC03_yy[i]+10;
      y[4] = CC03_yy[i]-10;
      break;
    default:
      std::cerr << "Wrong index Please confirm:" << __FUNCTION__ << " : " << __LINE__ << std::endl;
    }
    this->AddBin(5,x,y);
  }
}
int     CC03Poly::Fill(int ibin, double w){
  return TH2Poly::Fill(CC03_xx[ibin],CC03_yy[ibin],w);
}
int     CC03Poly::Fill(double x, double y, double w){
  return TH2Poly::Fill(x,y,w);
}

