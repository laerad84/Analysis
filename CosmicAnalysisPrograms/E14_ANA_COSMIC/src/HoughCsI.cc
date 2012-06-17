#include <TMath.h>
#include <iostream>
#include "HoughCsI.h"
#if !defined(__CINT__)
ClassImp(HoughCsI)
#endif

const Double_t HoughCsI::rohMax   = 1000.;
const Double_t HoughCsI::thetaMax = TMath::Pi();
const Int_t    HoughCsI::nDiv     = 10;
const Int_t    HoughCsI::nBinx    = 180*HoughCsI::nDiv;
const Int_t    HoughCsI::nBiny    = 120;
const Double_t HoughCsI::UnitDeg  = 1.0/HoughCsI::nDiv;
const Double_t HoughCsI::DegreeToPi = TMath::Pi()/180.;
HoughCsI::HoughCsI(){
	Init();
}

HoughCsI::~HoughCsI(){
  delete m_hisHough;
}

void HoughCsI::Init(){
  std::cout << "Init" << std::endl;
  m_hisHough = new TH2D("hisHough","",
			nBinx+1,
			-(90+HoughCsI::UnitDeg/2)*HoughCsI::DegreeToPi,
			( 90+HoughCsI::UnitDeg/2)*HoughCsI::DegreeToPi,
			HoughCsI::nBiny+1,
			-(1+0.5/HoughCsI::nBiny)*HoughCsI::rohMax,
			(1+0.5/HoughCsI::nBiny)*HoughCsI::rohMax);
  m_roh = 0;
  m_theta = 0;
  m_Cosmic = kFALSE;
}

void HoughCsI::Reset(){
  m_hisHough->Reset();
  m_roh = 0;
  m_theta = 0;
  m_Cosmic = kFALSE;
}

Bool_t HoughCsI::CosmicJudgment(TGraphErrors* grEvent){
  Double_t* Xarr = grEvent->GetX();
  Double_t* Yarr = grEvent->GetY();
  Double_t* Earr = grEvent->GetEY();

  for( Int_t indexTheta = 0; indexTheta<nBinx+1; indexTheta++){
    for( Int_t index  = 0; index< grEvent->GetN(); index++){
      Double_t theta = ((double)indexTheta/nBinx)*TMath::Pi();
      Double_t roh = Xarr[index]*TMath::Cos(theta) + Yarr[index]*TMath::Sin(theta);
      m_hisHough->Fill( theta, roh, Earr[index]);
    }
  }
  
  Int_t x,y,z;
  m_hisHough->GetMaximumBin(x,y,z);
  //std::cout << x  << " ; " << y << " ; " << z << std::endl;
  Double_t OveralWeight=0;
  Double_t WeightedX=0;
  Double_t WeightedY=0;
  Double_t OveralRMS[3]={};
  
  /*
    if(x!= 1 && x!=HoughCsI::nBinx+1){
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y  ),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y+1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y  ),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[0] = OveralRMS[0]/TMath::Power(m_hisHough->GetBinContent(x-1,y),2)/8 - 1;
    
    
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[1]=OveralRMS[1]/TMath::Power(m_hisHough->GetBinContent(x,y),2)/8 - 1;
    
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y  ),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y  ),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y+1),2);
    OveralRMS[2]=OveralRMS[2]/TMath::Power(m_hisHough->GetBinContent(x+1,y),2)/8 - 1;
    
    if( OveralRMS[0] < OveralRMS[1] ){
    if(OveralRMS[0] < OveralRMS[2]){
    x = x-1;
    }else{
    x = x+1;
    }
    }else if (OveralRMS[1] < OveralRMS[2]){
    x = x;
    }else{
    x = x+1;
    }
    OveralWeight +=m_hisHough->GetBinContent(x-1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x-1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x-1,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y-1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y  );
    WeightedX    -=m_hisHough->GetBinContent(x-1,y+1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y-1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y  );
    WeightedX    +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    /=OveralWeight;
    WeightedX    +=x;
    }else if( x ==1){
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[1]=OveralRMS[1]/TMath::Power(m_hisHough->GetBinContent(x,y),2)/8 - 1;
    
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y  ),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y  ),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y+1),2);
    OveralRMS[2]=OveralRMS[2]/TMath::Power(m_hisHough->GetBinContent(x+1,y),2)/8 - 1;
    
    if( OveralRMS[1] < OveralRMS[2]){
    x = x;
    }else{
    x = x+1;
    }
    OveralWeight +=m_hisHough->GetBinContent(x-1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x-1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x-1,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y-1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y  );
    WeightedX    -=m_hisHough->GetBinContent(x-1,y+1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y-1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y  );
    WeightedX    +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    /=OveralWeight;
    WeightedX    +=x;
    }else if( x == HoughCsI::nBinx+1){
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y  ),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y+1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y  ),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[0]=OveralRMS[0]/TMath::Power(m_hisHough->GetBinContent(x-1,y),2)/8 - 1;
    
    
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[1]=OveralRMS[1]/TMath::Power(m_hisHough->GetBinContent(x,y),2)/8 - 1;
    
    if( OveralRMS[0] < OveralRMS[1]){
    x = x-1;
    }else{
    x = x;
    }
    OveralWeight +=m_hisHough->GetBinContent(x-1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x-1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x-1,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y-1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y  );
    WeightedX    -=m_hisHough->GetBinContent(x-1,y+1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y-1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y  );
    WeightedX    +=m_hisHough->GetBinContent(x+1,y+1);
		WeightedX    /=OveralWeight;
		WeightedX    +=x;
		}
  */
  
  m_theta = x - 90;
  m_roh   = -1*rohMax+2*rohMax/HoughCsI::nBiny*((double)y-0.5);
  m_Cosmic = kTRUE;
  return m_Cosmic;
}


Bool_t HoughCsI::CosmicJudgment(TGraph* grEvent){
  
  Double_t* Xarr = grEvent->GetX();
  Double_t* Yarr = grEvent->GetY();
  
  if( grEvent->GetN() >200  || grEvent->GetN() < 10){
    return kFALSE;
  }

  for( Int_t indexTheta = 0; indexTheta<nBinx+1; indexTheta++){
    for( Int_t index  = 0; index< grEvent->GetN(); index++){
      Double_t theta = ((double)indexTheta-90.)*HoughCsI::DegreeToPi;
      Double_t roh = Xarr[index]*TMath::Cos(theta) + Yarr[index]*TMath::Sin(theta);
      m_hisHough->Fill( theta, roh);
    }
  }

  Int_t x,y,z;
  m_hisHough->GetMaximumBin(x,y,z);
  //std::cout << x  << " ; " << y << " ; " << z << std::endl;
  Double_t OveralWeight=0;
  Double_t WeightedX=0;
  Double_t WeightedY=0;
  Double_t OveralRMS[3]={};
  
  /*
  if(x!= 1 && x!=HoughCsI::nBinx+1){
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y  ),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y+1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y  ),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[0] = OveralRMS[0]/TMath::Power(m_hisHough->GetBinContent(x-1,y),2)/8 - 1;
    
    
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[1]=OveralRMS[1]/TMath::Power(m_hisHough->GetBinContent(x,y),2)/8 - 1;
    
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y  ),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y  ),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y+1),2);
    OveralRMS[2]=OveralRMS[2]/TMath::Power(m_hisHough->GetBinContent(x+1,y),2)/8 - 1;
    
    if( OveralRMS[0] < OveralRMS[1] ){
      if(OveralRMS[0] < OveralRMS[2]){
	x = x-1;
      }else{
	x = x+1;
      }
    }else if (OveralRMS[1] < OveralRMS[2]){
      x = x;
    }else{
      x = x+1;
    }
    OveralWeight +=m_hisHough->GetBinContent(x-1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x-1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x-1,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y-1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y  );
    WeightedX    -=m_hisHough->GetBinContent(x-1,y+1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y-1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y  );
    WeightedX    +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    /=OveralWeight;
    WeightedX    +=x;
  }else if( x ==1){
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[1]=OveralRMS[1]/TMath::Power(m_hisHough->GetBinContent(x,y),2)/8 - 1;
    
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y  ),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y-1),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y  ),2);
    OveralRMS[2]+=TMath::Power(m_hisHough->GetBinContent(x+2,y+1),2);
    OveralRMS[2]=OveralRMS[2]/TMath::Power(m_hisHough->GetBinContent(x+1,y),2)/8 - 1;
    
    if( OveralRMS[1] < OveralRMS[2]){
			x = x;
    }else{
      x = x+1;
    }
    OveralWeight +=m_hisHough->GetBinContent(x-1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x-1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x-1,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y-1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y  );
    WeightedX    -=m_hisHough->GetBinContent(x-1,y+1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y-1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y  );
    WeightedX    +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    /=OveralWeight;
    WeightedX    +=x;
  }else if( x == HoughCsI::nBinx+1){
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y  ),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-2,y+1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y  ),2);
    OveralRMS[0]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[0]=OveralRMS[0]/TMath::Power(m_hisHough->GetBinContent(x-1,y),2)/8 - 1;
    
    
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x-1,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x  ,y+1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y-1),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y  ),2);
    OveralRMS[1]+=TMath::Power(m_hisHough->GetBinContent(x+1,y+1),2);
    OveralRMS[1]=OveralRMS[1]/TMath::Power(m_hisHough->GetBinContent(x,y),2)/8 - 1;
    
    if( OveralRMS[0] < OveralRMS[1]){
      x = x-1;
    }else{
      x = x;
    }
    OveralWeight +=m_hisHough->GetBinContent(x-1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x-1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x-1,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x  ,y+1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y-1);
    OveralWeight +=m_hisHough->GetBinContent(x+1,y  );
    OveralWeight +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y-1);
    WeightedX    -=m_hisHough->GetBinContent(x-1,y  );
    WeightedX    -=m_hisHough->GetBinContent(x-1,y+1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y-1);
    WeightedX    +=m_hisHough->GetBinContent(x+1,y  );
    WeightedX    +=m_hisHough->GetBinContent(x+1,y+1);
    WeightedX    /=OveralWeight;
    WeightedX    +=x;
  }
   */
  m_theta = x*HoughCsI::UnitDeg - 90;
  m_roh   = -1*rohMax+2*rohMax/HoughCsI::nBiny*((double)y-0.5);
  m_Cosmic = kTRUE;

  //  std::cout << "theta:" << m_theta << " roh:" << m_roh << std::endl;
  return m_Cosmic;
}

void     HoughCsI::GetRohTheta(double &roh,double &theta){
  roh   = m_roh;
  theta = m_theta;
}

Double_t HoughCsI::GetRoh(){
  return m_roh;
}

Double_t HoughCsI::GetTheta(){
  return m_theta;
}

Double_t HoughCsI::GetTangentXbase(){
  if(TMath::Sin(m_theta) == 0){
    return 999999;
  }else{
    return -1*TMath::Cos(m_theta)/TMath::Sin(m_theta);
  }
}

Double_t HoughCsI::GetOffsetXbase(){
  if(TMath::Sin(m_theta) == 0){
    return 999999;
  }else{
    return m_roh/TMath::Sin(m_theta);
	}
}

Double_t HoughCsI::GetTangentYbase(){
  if(TMath::Cos(m_theta) == 0){
    return 999999;
  }else{
    return -1*TMath::Sin(m_theta)/TMath::Cos(m_theta);
  }
}
Double_t HoughCsI::GetOffsetYbase(){
  if(TMath::Cos(m_theta) == 0){
    return 999999;
  }else{
    return m_roh/TMath::Cos(m_theta);
  }
}

TH2D*    HoughCsI::GetHisHough(){
  return m_hisHough;
}


TF1*     HoughCsI::GetFunction(){
  if( TMath::Sin(m_theta) == 0 ){
    m_theta = 0.000000001;
  }
  TF1* func = new TF1("func",Form("%lf*x+%lf",
				  -1*TMath::Cos(m_theta*HoughCsI::DegreeToPi)/TMath::Sin(m_theta*HoughCsI::DegreeToPi),
				  m_roh/TMath::Sin(m_theta*HoughCsI::DegreeToPi)),
		      -1000.,1000.);
  //  std::cout<< m_theta << " : " 
  //	   << m_roh   << " : " 
  //	   << -1*TMath::Cos(m_theta*HoughCsI::DegreeToPi)/TMath::Sin(m_theta*HoughCsI::DegreeToPi) << " : " 
  //	   << func->GetParameter(0) << " : "
  //	   << func->GetParameter(1) <<  std::endl;
  return func;  
}

Bool_t   HoughCsI::IsCosmic(){
  return m_Cosmic;
}

Double_t HoughCsI::GetCalibrationFactor(){
  Double_t factor = TMath::Cos(m_theta*HoughCsI::DegreeToPi);
  return factor; 
}

TLine* HoughCsI::GetLine(){
  Double_t cos  = TMath::Cos(m_theta*HoughCsI::DegreeToPi);
  Double_t sin  = TMath::Sin(m_theta*HoughCsI::DegreeToPi);
  Double_t LineX[2];
  Double_t LineY[2];

  Double_t x[2], y[2];
  if( TMath::Sin(m_theta) == 0){
    m_theta = 0.000000001;
  }

  double b = TMath::Sqrt( 950*950 - m_roh*m_roh );
  
  LineX[0] =(m_roh*cos +b*sin);
  LineX[1] =(m_roh*cos -b*sin);
  LineY[0] =(m_roh*sin -b*cos);
  LineY[1] =(m_roh*sin +b*cos);
  /*
  LineX[0] =0;
  LineX[1] = m_roh*cos;
  LineY[0] = 0;
  LineY[1] = m_roh*sin;
  */
  TLine* line  = new TLine(LineX[0],LineY[0],LineX[1],LineY[1]);
  return line;
}

Double_t HoughCsI::CalDistance(Double_t x , Double_t y){
  Double_t cos = TMath::Cos(m_theta*HoughCsI::DegreeToPi);
  Double_t sin = TMath::Sin(m_theta*HoughCsI::DegreeToPi);
  Double_t Distance = TMath::Abs((x-(1./cos*(m_roh - y*sin)))*cos);
  //  std::cout << cos << " : " << sin << " : " << Distance << std::endl;

  return Distance;
}
  
