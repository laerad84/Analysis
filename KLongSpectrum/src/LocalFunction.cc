#ifndef  _LOCAL_FUNCTION_H_
#include "LocalFunction.h"
#endif   _LOCAL_FUNCTION_H_

void SetGammaTime(Gamma &g){
  Double_t t = GetTiming(g);
  g.setTime(t);
}

void SetGammaTime(std::list<Gamma> glist){
  std::list<Gamma>::iterator git = glist.begin();
  for( ; git != glist.end(); git++){
    Double_t t = GetTiming( (*git) );
    (*git).setTime(t);
  }
}
double GetWeight( Gamma g){
  Double_t w(0);
  for( int i  =0; i< g.clusterEVec().size(); i++){
    w += TMath::Sqrt( g.clusterEVec()[i]);
  }
  return w;
}
double GetTiming( Gamma g){
  Double_t t(0);
  Double_t w(0);
  for( int i = 0; i< g.clusterTimeVec().size(); i++){
    t += g.clusterTimeVec()[i]*TMath::Sqrt( g.clusterEVec()[i]);
    w += TMath::Sqrt(g.clusterEVec()[i]);
  }
  t = t/w;
  return t;
}
double GetClusterTSigma( Gamma g){
  double tMean = GetTiming( g );
  double tWeight = GetWeight( g );
  Double_t Sig(0);
  Double_t wg(0);
  for( int i = 0; i< g.clusterTimeVec().size(); i++){
    Sig += TMath::Power( g.clusterTimeVec()[i] - tMean , 2)*TMath::Sqrt(g.clusterEVec()[i]);
    wg  += TMath::Sqrt(g.clusterEVec()[i]);
  }
  Double_t SigtCluster = TMath::Sqrt( Sig/wg );
  return SigtCluster;
}
double CalGammaTOF( Klong kl, Gamma g ){
  double length = 0;
  length = sqrt( pow(g.x()-kl.vx(),2)+ pow(g.y()-kl.vy(),2)+pow(g.z()-kl.vz(),2));
  return length/SpeedOfLight;  
}
