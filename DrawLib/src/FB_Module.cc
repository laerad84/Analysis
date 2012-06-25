#ifndef FB_MODULE__H__
#include "FB_Module.h"
#endif

#include "TVector2.h"
#include "TMath.h"
#if !defined(__CINT__)
ClassImp(FB_Module)
#endif


FB_Module::FB_Module(const char* name){
  Init(name);
  Reset();
}
FB_Module::~FB_Module(){
  delete FB;
}

void    FB_Module::Init(const char* name ){
  FB = new TH2Poly(name, "" ,-1500,1500,-1500,1500);
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
    FB_Poly[i] = new TPolyLine(5,FB_Posx,FB_Posy);
    FB->AddBin(5,FB_Posx,FB_Posy);
  }  
   
   for( int i = 0; i< numberOfFB; i++){
     Dep[numberOfFB]=0;
   }
}

int     FB_Module::Fill(int ibin){
  return FB->Fill(FB_xx[ibin],FB_yy[ibin]);
}

int     FB_Module::Fill(int ibin, Double_t w){
  return FB->Fill(FB_xx[ibin],FB_yy[ibin],w);
}

int     FB_Module::Fill(double x, double y, double w = 1 ){
  return FB->Fill(x,y,w);
}
void    FB_Module::SetTitle(const char* title){
  FB->SetNameTitle(FB->GetName(),title);
}
void    FB_Module::Reset( void ){
  ResetContent();
  UpdateValue();
  FB->Reset();
}

void    FB_Module::ResetContent( void ){
  for ( int i = 0; i< numberOfFB; i++){
    Dep[i] = 0; 
  }
}

void    FB_Module::UpdateValue( void ){
  for( int i = 0; i<numberOfFB; i++){
    FB->SetBinContent( i, 0 );
  }
}

void    FB_Module::Draw(const char* option = "colz"){
  FB->Draw(option);
  for( int i = 0; i< numberOfFB; i++){
    FB_Poly[i]->Draw();
  }
}

void FB_Module::DrawWithRange( double minuser , double maxuser, const char* option = "colz"){
  FB->SetMaximum( maxuser );
  FB->SetMinimum( minuser );
  FB->Draw(option);
  for( int i = 0; i< numberOfFB; i++){
    FB_Poly[i]->Draw();
  }
}
