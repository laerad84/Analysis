#include "KlongParameter.h"


KlongParameter::KlongParameter(){
  ;
}
KlongParameter::~KlongParameter(){
  ;
}
void KlongParameter::Branch( TTree* trout){
  trout->Branch("GMaxR",&GMaxR,"GMaxR/D");
  trout->Branch("GMinX",&GMinX,"GMinX/D");
  trout->Branch("GMinY",&GMinY,"GMinY/D");
  trout->Branch("GMinE",&GMinE,"GMinE/D");
  trout->Branch("GMaxTDelta",&GMaxTDelta,"GMaxTDelta/D");

  trout->Branch("klChisqZ",&klChisqZ,"klChisqZ/D");
  trout->Branch("klpt" ,&klpt,"klpt/D");
  trout->Branch("klpos",klpos,"klpos[3]/D");
  trout->Branch("klE"  ,&klE ,"klE/D");

  trout->Branch("pi0ptMax",&pi0ptMax,"pi0ptMax/D");
  trout->Branch("pi0pt",pi0pt,"pi0pt[3]/D");
  
}
void KlongParameter::InitPar(){
  GMaxR=0;
  GMinX=0;
  GMinY=0;
  GMinE=100000;
  klChisqZ=0;
  klpt=0;
  klE=0;
  pi0ptMax=0;
  for( int i = 0; i< 3; i++){
    klpos[i]=0;
    pi0pt[i]=0;
  }
}
void KlongParameter::Judge(Klong kl){

  double GMaxR;
  double GMinX;
  double GMinY;
  double GMinE;
  double GMaxTDelta;

  double klChisqZ;
  double klpt;
  double klpos[3];
  double klE;
  double pi0ptMax;
  double pi0pt[3];
  InitPar();
  
  for( int i = 0; i< kl.pi0().size(); i++){
    double x;
    double y;
    double e;
    double r;
    double ppt;

    ppt = kl.pi0()[i].p3().perp();
    pi0pt[i] = ppt;
    if( ppt > pi0ptMax ){
      pi0ptMax = ppt;
    }

    x = kl.pi0()[i].g1().x();
    y = kl.pi0()[i].g1().y();
    e = kl.pi0()[i].g1().e();
    r = TMath::Sqrt( x*x + y*y );
    if( r > GMaxR ){
      GMaxR = r;
    }
    if( abs(x) < GMinX ){
      GMinX = x;
    }
    if( abs(y) < GMinY ){
      GMinY = y;
    }
    if( e < GMinE ){
      GMinE = e;
    } 

    x = kl.pi0()[i].g2().x();
    y = kl.pi0()[i].g2().y();
    e = kl.pi0()[i].g2().e();
    r = TMath::Sqrt( x*x + y*y );
    if( r > GMaxR ){
      GMaxR = r;
    }
    if( abs(x) < GMinX ){
      GMinX = x;
    }
    if( abs(y) < GMinY ){
      GMinY = y;
    }
    if( e < GMinE ){
      GMinE = e;
    } 


  }
}
