#ifndef MBPOLY__H__
#include "MBPoly.h"
#endif

#include "TVector2.h"
#include "TMath.h"
#if !defined(__CINT__)
ClassImp(MBPoly)
#endif

MBPoly::MBPoly( ){
  Initialize( -1500,1500,-1500,1500 ,25,25);
  SetName("NoName");
  SetTitle("NoTitle");
  SetFloat();
  Init();
}
MBPoly::MBPoly(const char* name , const char* title ){
  Initialize( -1500,1500,-1500,1500 ,25,25);
  SetName( name );
  SetTitle( title );
  SetFloat();
  Init();
}
MBPoly::~MBPoly(){
  ;
}
void    MBPoly::Init( ){

  Double_t MB_x[5] = {-1003,-1349.5,-1349.5,-1003,-1003};
  Double_t MB_y[5] = { -100,   -100,  168.5,  100, -100};
  TVector2 vecCenter[numberOfMB];
  vecCenter[0].Set((-1003-1349.5-1349.5-1003)/4 , (-100-100,168.5,100)/4);
  for( int i = 1 ; i <numberOfMB; i++){
    vecCenter[i] = vecCenter[0].Rotate(-2*TMath::Pi()/numberOfMB*i);
  } 
  for( int i = 0; i< numberOfMB; i++){
    MB_xx[i] = vecCenter[i].X();
    MB_yy[i] = vecCenter[i].Y();
  }

  TVector2 vec[5][numberOfMB];
  for( int i = 0; i< 5; i++){
    vec[i][0].Set(MB_x[i],MB_y[i]);
  }
  for( int i = 1; i< numberOfMB; i++){
    for( int j = 0; j< 5; j++){
      vec[j][i] = vec[j][0].Rotate(-2*TMath::Pi()/numberOfMB*i);
    }
  }
  for( int i = 0; i< numberOfMB; i++){
    Double_t MB_Posx[5]={0};
    Double_t MB_Posy[5]={0};
    for( int j = 0; j < 5; j++){
      MB_Posx[j] = vec[j][i].X();
      MB_Posy[j] = vec[j][i].Y();
    }    
    this->AddBin(5,MB_Posx,MB_Posy);
  }
}
int     MBPoly::Fill(int ibin, Double_t w ){
  return TH2Poly::Fill(MB_xx[ibin],MB_yy[ibin],w);
}
int     MBPoly::Fill(double x, double y, double w ){
  return TH2Poly::Fill(x,y,w);
}


