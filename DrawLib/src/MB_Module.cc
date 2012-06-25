#ifndef MB_MODULE__H__
#include "MB_Module.h"
#endif

#include "TVector2.h"
#include "TMath.h"
#if !defined(__CINT__)
ClassImp(MB_Module)
#endif


MB_Module::MB_Module(const char* name){
  Init(name);
  Reset();
}
MB_Module::~MB_Module(){
  delete MB;
}

void    MB_Module::Init(const char* name ){
  MB = new TH2Poly(name, "" ,-1500,1500,-1500,1500);
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
    MB_Poly[i] = new TPolyLine(5, MB_Posx, MB_Posy);
    MB->AddBin(5,MB_Posx,MB_Posy);
  }
  for( int i = 0; i< numberOfMB; i++){
    Dep[numberOfMB]=0;
  }
}

int     MB_Module::Fill(int ibin){
  Dep[ibin]++;  
  return MB->Fill(MB_xx[ibin],MB_yy[ibin]);
}

int     MB_Module::Fill(int ibin, Double_t w){
  Dep[ibin] += w;
  return MB->Fill(MB_xx[ibin],MB_yy[ibin],w);
}
int     MB_Module::Fill(double x, double y, double w = 1){
  return MB->Fill(x,y,w);
}

void    MB_Module::SetTitle(const char* title){
  MB->SetNameTitle("MB",title);
}
void    MB_Module::Reset( void ){
  ResetContent();
  UpdateValue();
  MB->Reset();
}

void    MB_Module::ResetContent( void ){
  for ( int i = 0; i< numberOfMB; i++){
    Dep[i] = 0; 
  }
}

void    MB_Module::UpdateValue( void ){
  for( int i = 0; i<numberOfMB; i++){
    MB->SetBinContent( i, Dep[i] );
  }
}

void    MB_Module::Draw(const char* option = "colz"){
  MB->Draw(option);
  for( int i = 0; i< numberOfMB ; i++){
    MB_Poly[i]->Draw();
  }
}

void MB_Module::DrawWithRange( double minuser, double maxuser, const char* option = "colz"){
  MB->SetMaximum( maxuser );
  MB->SetMimimum( minuser );
  MB->Draw(option);
  for( int i  =0; i< numberOfMB; i++){
    MB_Poly[i]->Draw();
  }
}
