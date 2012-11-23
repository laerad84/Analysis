#ifndef OEVPoly__H__
#include "OEVPoly.h"
#endif

#include "TVector2.h"
#include "TMath.h"
#if !defined(__CINT__)
ClassImp(OEVPoly)
#endif


OEVPoly::OEVPoly( ){
  Initialize( -1000,1000,-1000,1000,25,25);
  SetName("NoName");
  SetTitle("NoTitle");
  SetFloat();
  Init();
}
OEVPoly::OEVPoly( const char* name , const char* title ){
  Initialize( -1000,1000,-1000,1000,25,25);
  SetName( name );
  SetTitle( title );
  SetFloat( kFALSE );
  Init();
}
OEVPoly::~OEVPoly(){  
  ;
}
void    OEVPoly::Init( ){

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

    this->AddBin(OEVx[VectorIndex].size(),OEV_x[i],OEV_y[i]);
  }
  this->SetLineColor(1);
}
int     OEVPoly::Fill(int ibin, Double_t w){
  return TH2Poly::Fill(OEV_xx[ibin],OEV_yy[ibin],w);
}
int     OEVPoly::Fill(double x, double y, double w){
  return TH2Poly::Fill(x,y,w);
}


