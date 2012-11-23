#ifndef FBPoly__H__
#include "FBPoly.h"
#endif

#include "TVector2.h"
#include "TMath.h"
#if !defined(__CINT__)
ClassImp(FBPoly)
#endif

FBPoly::FBPoly(){
  Initialize( -1500,1500,-1500,1500,25,25);
  SetName("NoName");
  SetTitle("NoTitle");
  SetFloat();
  Init();
}
FBPoly::FBPoly( const char* name, const char* title ){
  Initialize( -1500,1500,-1500,1500,25,25);
  SetName( name );
  SetTitle( title );
  SetFloat( kFALSE );
  Init();
}

FBPoly::~FBPoly(){
}
void    FBPoly::Init(){

  Double_t FB_x[5] = { -310,   -825, -785.5, -310, -310};
  Double_t FB_y[5] = {-60.8,  -60.8,  259.5, 60.8,-60.8};

  TVector2 vecCenter[numberOfFB];
  vecCenter[0].Set((-310-825-785.5-310)/4,(-60.8-60.8+259.5+60.8)/4);
  for( int i = 1; i< numberOfFB; i++){
    vecCenter[i] = vecCenter[0].Rotate(-2*TMath::Pi()/numberOfFB*i);
  }
  for( int i = 0; i< numberOfFB; i++){
    FB_xx[i] = vecCenter[i].X();
    FB_yy[i] = vecCenter[i].Y();
  }

  TVector2 vec[5][numberOfFB];
  for( int i = 0; i< 5; i++){
    vec[i][0].Set(FB_x[i],FB_y[i]);
  }
  for( int i = 1; i< numberOfFB; i++){
    for( int j = 0; j< 5; j++){
      vec[j][i] = vec[j][0].Rotate(-2*TMath::Pi()/numberOfFB*i);
    }
  }

  for( int i = 0; i< numberOfFB; i++){
    Double_t FB_Posx[5]={0};
    Double_t FB_Posy[5]={0};
    for( int j = 0; j < 5; j++){
      FB_Posx[j] = vec[j][i].X();
      FB_Posy[j] = vec[j][i].Y();
    }        
    this->AddBin(5,FB_Posx,FB_Posy);
  }     
}
int     FBPoly::Fill(int ibin, double w){
  return TH2Poly::Fill(FB_xx[ibin],FB_yy[ibin],w);
}
int     FBPoly::Fill(double x, double y, double w ){
  return TH2Poly::Fill(x,y,w);
}

