#ifndef CC03_MODULE__H__
#include "CC03_Module.h"
#endif

#include "TVector2.h"
#include "TMath.h"
#if !defined(__CINT__)
ClassImp(CC03_Module)
#endif


CC03_Module::CC03_Module(const char* name){
  Init(name);
  Reset();
}
CC03_Module::~CC03_Module(){
  delete CC03;
}

void    CC03_Module::Init(const char* name ){
  CC03 = new TH2Poly(name, "" ,-100,100,-100,100);
  Double_t CC03_initx[numberOfCC03] = {  90.0  ,  90.0,  90.0,  90.0,  90.0,   90.0,   90.0,  90.0,
					88.75 , 66.25, 43.75, 21.25, -1.25, -23.75, -46.25, -68.75,
					-90.0 , -90.0, -90.0, -90.0, -90.0,  -90.0,  -90.0,  -90.0,
					-88.75,-66.25,-43.75,-21.25,  1.25,  23.75,  46.25,  68.75};
  
  Double_t CC03_inity[numberOfCC03] = {  88.75 ,66.25 ,43.75  ,21.25 ,-1.25 ,-23.75,-46.25,-68.75,
					-90.0 ,-90.0 ,-90.0  ,-90.0 ,-90.0 ,-90.0 ,-90.0 ,-90.0 ,
					-88.75,-66.25,-43.75 ,-21.25,1.25  ,23.75 ,46.25 ,68.75 ,
					90.0  ,90.0  ,90.0   ,90.0  ,90.0  ,90.0  ,90.0  ,90.0  };
  for( int i = 0; i < numberOfCC03; i++){
    CC03_xx[i] = CC03_initx[i];
    CC03_yy[i] = CC03_inity[i];
  }

  for( int i = 0; i< numberOfCC03; i++ ){
    Double_t x[5]={0};
    Double_t y[5]={0};
    switch( i/8 ){
    case 0:
    case 2:
      CC03_Box[i] = new TBox( CC03_xx[i] -10, CC03_yy[i]-11.25, CC03_xx[i] +10, CC03_yy[i]+11.25);
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
      CC03_Box[i] = new TBox( CC03_xx[i] -11.25, CC03_yy[i]-10, CC03_xx[i] +11.25, CC03_yy[i]+10);
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
    CC03_Box[i]->SetFillColor(0);
    CC03_Box[i]->SetFillStyle(0);
    CC03->AddBin(5,x,y);
  }
}

int     CC03_Module::Fill(int ibin){
  return CC03->Fill(CC03_xx[ibin],CC03_yy[ibin]);
}

int     CC03_Module::Fill(int ibin, Double_t w){
  return CC03->Fill(CC03_xx[ibin],CC03_yy[ibin],w);
}
int     CC03_Module::Fill(double x, double y, double w = 1){
  return CC03->Fill(x,y,w);
}

void    CC03_Module::SetTitle(const char* title){
  CC03->SetNameTitle("CC03",title);
}

void    CC03_Module::Reset( void ){
  CC03->Reset();
}

void    CC03_Module::Draw(const char* option = "colz"){
  CC03->Draw(option);
  for( int i = 0; i< numberOfCC03; i++){
    CC03_Box[i]->Draw();
  }
}

void   CC03_Module::DrawWithRange( double minuser, double maxuser, const char* option = "colz"){
  CC03->SetMaximum(maxuser);
  CC03->SetMinimum(minuser);
  CC03->Draw(option);
  for( int i = 0; i< numberOfCC03; i++){
    CC03_Box[i]->Draw();
  }
}
