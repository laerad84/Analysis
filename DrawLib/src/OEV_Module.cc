#ifndef OEV_MODULE__H__
#include "OEV_Module.h"
#endif

#include "TVector2.h"
#include "TMath.h"

#if !defined(__CINT__)
ClassImp(OEV_Module)
#endif


OEV_Module::OEV_Module(const char* name){
  Init(name);
  Reset();
}
OEV_Module::~OEV_Module(){
  //delete OEV;
}

void    OEV_Module::Init(const char* name ){
  OEV = new TH2Poly(name, "" ,-1000,1000,-1000,1000);

  std::vector<Double_t> OEVx[6];
  std::vector<Double_t> OEVy[6];

  Double_t OEV_x[44][8];
  Double_t OEV_y[44][8];
  for( int i = 0; i< 44; i++){
    for( int j = 0; j< 8; j++){
      OEV_x[i][j] = 0;
      OEV_y[i][j] = 0;
    }
  }

  Double_t OEV_arr2[6]={ 925, 925, 875, 825, 751, 651};
  Double_t OEV_arr1[6]={  50, 150, 350, 450, 551, 651};

  for( int i = 0; i< numberOfOEV; i++){
    OEV_xx[i] = 0;
    OEV_yy[i] = 0;
  }

  const Double_t R=950.;
  Double_t xzero=0;
  Double_t yzero=0;
  Double_t yfirst=0;
  Double_t ysecond=0;
  Double_t x_7[8] = {600,650,650,700,700,TMath::Sqrt(950*950-700*700),600,600};
  Double_t y_7[8] = {650,650,600,600,TMath::Sqrt(950*950-700*700),700,700,650};
  for( int index = 0; index < 6 ; ++index){
    switch(index){
    case 0:
      (OEVx[index]).push_back(100*index);
      (OEVx[index]).push_back(100*(index+1));
      (OEVx[index]).push_back(100*(index+1));
      (OEVx[index]).push_back(100*index);
      (OEVx[index]).push_back(100*index);
      
      (OEVy[index]).push_back(900);
      (OEVy[index]).push_back(900);
      (OEVy[index]).push_back(TMath::Sqrt(R*R-100*(index+1)*100*(index+1)));
      (OEVy[index]).push_back(TMath::Sqrt(R*R-100*(index)*100*(index)));
      (OEVy[index]).push_back(900);
      break;
    case 1:
      (OEVx[index]).push_back(100);
      (OEVx[index]).push_back(TMath::Sqrt(R*R-900*900));
      (OEVx[index]).push_back(200);
      (OEVx[index]).push_back(100);
      (OEVx[index]).push_back(100);
      
      (OEVy[index]).push_back(900);
      (OEVy[index]).push_back(900);
      (OEVy[index]).push_back(TMath::Sqrt(R*R-200*200));
      (OEVy[index]).push_back(TMath::Sqrt(R*R-100*100));
      (OEVy[index]).push_back(900);
      break;
    case 2:
      xzero = 100*3;
      yzero = 850;
      yfirst= yzero+50;
      (OEVx[index]).push_back(xzero);
      (OEVx[index]).push_back(TMath::Sqrt(R*R-yzero*yzero));
      (OEVx[index]).push_back(TMath::Sqrt(R*R-yfirst*yfirst));
      (OEVx[index]).push_back(xzero);
      (OEVx[index]).push_back(xzero);
      
      (OEVy[index]).push_back(yzero);
      (OEVy[index]).push_back(yzero);
      (OEVy[index]).push_back(yfirst);
      (OEVy[index]).push_back(yfirst);
      (OEVy[index]).push_back(yzero);
      break;
    case 3:
      xzero = 400;
      yzero = 800;
      yfirst= yzero+50;

      (OEVx[index]).push_back(xzero);
      (OEVx[index]).push_back(TMath::Sqrt(R*R-yzero*yzero));
      (OEVx[index]).push_back(TMath::Sqrt(R*R-yfirst*yfirst));
      (OEVx[index]).push_back(xzero);
      (OEVx[index]).push_back(xzero);
      
      (OEVy[index]).push_back(yzero);
      (OEVy[index]).push_back(yzero);
      (OEVy[index]).push_back(yfirst);
      (OEVy[index]).push_back(yfirst);
      (OEVy[index]).push_back(yzero);
      
      break;
    case 4:
      xzero  = 550;
      yzero  = 700;
      yfirst = 750;
      ysecond= 800;
      
      (OEVx[index]).push_back(xzero);
      (OEVx[index]).push_back(TMath::Sqrt(R*R-yzero*yzero));
      (OEVx[index]).push_back(TMath::Sqrt(R*R-yfirst*yfirst));
      (OEVx[index]).push_back(TMath::Sqrt(R*R-ysecond*ysecond));
      (OEVx[index]).push_back(xzero-50);
      (OEVx[index]).push_back(xzero-50);
      (OEVx[index]).push_back(xzero);
      (OEVx[index]).push_back(xzero);
      
      (OEVy[index]).push_back(yzero); 
      (OEVy[index]).push_back(yzero); 
      (OEVy[index]).push_back(yfirst); 
      (OEVy[index]).push_back(ysecond); 
      (OEVy[index]).push_back(ysecond); 
      (OEVy[index]).push_back(yfirst); 
      (OEVy[index]).push_back(yfirst); 
      (OEVy[index]).push_back(yzero); 
      break;
    case 5:
      for( int ipoint = 0; ipoint< 8; ++ipoint){
	(OEVx[index]).push_back(x_7[ipoint]);
	(OEVy[index]).push_back(y_7[ipoint]);
      }
      break;
    default: 
      std::cerr << "Out Of Range" << std::endl;
    }
  }

  //std::cout << "SetNumber" << std::endl;
  int nForthOEV=11;
  int nMiddleOEV=5;
  for( int i = 0; i< numberOfOEV; ++i){
    int VectorIndex = nMiddleOEV-TMath::Abs(i%nForthOEV-nMiddleOEV);
    
    int VectorxySign = 0;
    int VectorpmxSign = 1;
    int VectorpmySign = 1;
    int position = i/nForthOEV;
    
    if( position%2 == 0){	
      if( (i%nForthOEV) < nMiddleOEV ){
	VectorxySign = 1;
      }else{
	VectorxySign = 0;
      }
    }else{
      if( (i%nForthOEV) < nMiddleOEV ){
	VectorxySign = 0;
      }else{
	VectorxySign = 1;
      }
    }

    switch( position ){
    case 0:	
      VectorpmxSign = -1;
      break;
      case 1:
	break;
    case 2:
      VectorpmySign =-1;
      break;
    case 3:
      VectorpmxSign =-1;
      VectorpmySign =-1;
      break;
    default:
      std::cerr <<" Out Of Range" << std::endl;
    }

    for( unsigned int  vecIndex = 0; vecIndex < OEVx[VectorIndex].size(); ++vecIndex){

      if(VectorxySign == 0){
	OEV_x[i][vecIndex] = OEVx[VectorIndex][vecIndex]*VectorpmxSign;
	OEV_y[i][vecIndex] = OEVy[VectorIndex][vecIndex]*VectorpmySign;
      }else{
	OEV_x[i][vecIndex] = OEVy[VectorIndex][vecIndex]*VectorpmxSign;
	OEV_y[i][vecIndex] = OEVx[VectorIndex][vecIndex]*VectorpmySign;
      }
    }
    if(VectorxySign ==0){
      OEV_xx[i] = OEV_arr1[VectorIndex]*VectorpmxSign;
      OEV_yy[i] = OEV_arr2[VectorIndex]*VectorpmySign;
    }else{    
      OEV_xx[i] = OEV_arr2[VectorIndex]*VectorpmxSign;
      OEV_yy[i] = OEV_arr1[VectorIndex]*VectorpmySign;
    }

    OEV->AddBin(OEVx[VectorIndex].size(),OEV_x[i],OEV_y[i]);
    OEV_Poly[i] = new TPolyLine(OEVx[VectorIndex].size(),OEV_x[i],OEV_y[i]);
  }
  OEV->SetLineColor(1);
  for( int i = 0; i< numberOfOEV; i++){
    Dep[numberOfOEV]=0;
  }
}

int     OEV_Module::Fill(int ibin){
  return OEV->Fill(OEV_xx[ibin],OEV_yy[ibin]);
}

int     OEV_Module::Fill(int ibin, Double_t w){
  return OEV->Fill(OEV_xx[ibin],OEV_yy[ibin],w);
}
int     OEV_Module::Fill(double x, double y, double w = 1){
  return OEV->Fill(x,y,w);
}
void    OEV_Module::SetTitle(const char* title){
  OEV->SetNameTitle(OEV->GetName(),title);
}
void    OEV_Module::Reset( void ){
  ResetContent();
  UpdateValue();
  OEV->Reset("");
}

void    OEV_Module::ResetContent( void ){
  for ( int i = 0; i< numberOfOEV; i++){
    Dep[i] = 0; 
  }
}

void    OEV_Module::UpdateValue( void ){
  for( int i = 0; i<numberOfOEV; i++){
    OEV->SetBinContent( i, Dep[i] );
  }
}

void    OEV_Module::Draw(const char* option = "COLZ"){
  OEV->Draw(option);
  for( int i = 0; i< numberOfOEV; i++){
    OEV_Poly[i]->Draw();
  }
}


void   OEV_Module::DrawWithRange( double minuser, double maxuser, const char* option = "colz"){
  double maxValue = 0, minValue = 1000000000,content = 0;
  OEV->SetMaximum( maxuser );
  OEV->SetMinimum( minuser );
  OEV->Draw(option);
  for( int i = 0; i< numberOfOEV; i++){
    OEV_Poly[i]->Draw();
  }
}