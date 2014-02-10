#include "User_Functions.h"

bool User_RecG2(std::list<Gamma> const &glist, std::list<Pi0>& piList){
  static Rec2g rec2g;      

  // reconstruction 
  piList = rec2g.recPi0withConstM(glist);
  if(piList.size()!=1)   return false;
  
  // position correction for angle dependency
  E14GNAnaFunction::getFunction()->correctPosition(piList.front().g1());  
  E14GNAnaFunction::getFunction()->correctPosition(piList.front().g2());  
  
  // re-reconstruction with corrected gamma
  std::list<Gamma> glist2;
  glist2.push_back(piList.front().g1());
  glist2.push_back(piList.front().g2());
  piList = rec2g.recPi0withConstM(glist2);
  if(piList.size()!=1)   return false;

  // shape chi2 evaluation 
  E14GNAnaFunction::getFunction()->shapeChi2(piList.front().g1());  
  E14GNAnaFunction::getFunction()->shapeChi2(piList.front().g2());  

  return true;
}
bool User_RecG4(std::list<Gamma> const &glist, std::vector<Klong> &klVec){
  static RecKlong recKl;      
  //// reconstruction
  //  klVec =  recKl.recK2pi0(glist,VERTEX_FIX_XYZERO) ;
  klVec =  recKl.recK2pi0(glist) ;
  
  if(klVec.size()==0) {
    return false;
  }

  //// gamma position correction for angle dependency
  std::list<Gamma> glist2;
  for( std::vector<Pi0>::iterator it=klVec[0].pi0().begin();
       it!=klVec[0].pi0().end(); it++){
    E14GNAnaFunction::getFunction()->correctPosition(it->g1());  
    glist2.push_back(it->g1());
    E14GNAnaFunction::getFunction()->correctPosition(it->g2());  
    glist2.push_back(it->g2());
  }
  
  //// re-reconstruction with corrected gamma
  //  klVec =  recKl.recK2pi0(glist2,VERTEX_FIX_XYZERO) ;
  klVec =  recKl.recK2pi0(glist2) ;
  if(klVec.size()==0)    return false;

  //// shape chi2 evaluation
  for( std::vector<Pi0>::iterator it=klVec[0].pi0().begin();
       it!=klVec[0].pi0().end(); it++){
    E14GNAnaFunction::getFunction()->shapeChi2(it->g1());  
    E14GNAnaFunction::getFunction()->shapeChi2(it->g2());  
  }

  return true;
}
bool User_RecG6(std::list<Gamma> const &glist, std::vector<Klong> &klVec){
  static RecKlong recKl;      
  //// reconstruction
  //  klVec =  recKl.recK3pi0(glist,VERTEX_FIX_XYZERO) ;
  klVec =  recKl.recK3pi0(glist) ;
  
  if(klVec.size()==0) {
    return false;
  }

  //// gamma position correction for angle dependency
  std::list<Gamma> glist2;
  for( std::vector<Pi0>::iterator it=klVec[0].pi0().begin();
       it!=klVec[0].pi0().end(); it++){
    E14GNAnaFunction::getFunction()->correctPosition(it->g1());  
    glist2.push_back(it->g1());
    E14GNAnaFunction::getFunction()->correctPosition(it->g2());  
    glist2.push_back(it->g2());
  }
  
  //// re-reconstruction with corrected gamma
  //  klVec =  recKl.recK3pi0(glist2,VERTEX_FIX_XYZERO) ;
  klVec =  recKl.recK3pi0(glist2) ;
  if(klVec.size()==0)    return false;

  //// shape chi2 evaluation
  for( std::vector<Pi0>::iterator it=klVec[0].pi0().begin();
       it!=klVec[0].pi0().end(); it++){
    E14GNAnaFunction::getFunction()->shapeChi2(it->g1());  
    E14GNAnaFunction::getFunction()->shapeChi2(it->g2());  
  }

  return true;
}

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


GammaCut::GammaCut(){
  ;
}
GammaCut::~GammaCut(){
  ;
}
void GammaCut::Reset(){
  GMinEnergy = 10000;
  GMinDist   = 10000;
  GMinY      = 10000;
  GMinX      = 10000;
  GMaxR      = 0;
  GMaxY      = 0;
  GMaxChi    = 0;
  GMaxDeltaT = 0;
  GMaxDeltaTID = -1;
}
void GammaCut::Decision(Klong kl){
  Reset();
  int gammaIndex(0);
  double tmpE[6];
  double tmpX[6];
  double tmpY[6];
  double tmpR[6];
  double tmpT[6];
  double tmpChi[6];
  std::list<Gamma> glist;
  for( int i = 0; i< kl.pi0().size();i++){
    glist.push_back(kl.pi0()[i].g1());
    glist.push_back(kl.pi0()[i].g2());
  }
  /*
  for( int i = 0; i< kl.pi0().size(); i++){
    tmpE[gammaIndex] = kl.pi0()[i].g1().e();
    tmpX[gammaIndex] = kl.pi0()[i].g1().x();
    tmpY[gammaIndex] = kl.pi0()[i].g1().y();
    tmpR[gammaIndex] = TMath::Sqrt( tmpX[gammaIndex]*tmpX[gammaIndex] + tmpY[gammaIndex]*tmpY[gammaIndex]);
    tmpT[gammaIndex] = kl.pi0()[i].g1().t();
    tmpChi[gammaIndex] = kl.pi0()[i].g1().chisq();
    gammaIndex++;
    tmpE[gammaIndex] = kl.pi0()[i].g2().e();
    tmpX[gammaIndex] = kl.pi0()[i].g2().x();
    tmpY[gammaIndex] = kl.pi0()[i].g2().y();
    tmpR[gammaIndex] = TMath::Sqrt( tmpX[gammaIndex]*tmpX[gammaIndex] + tmpY[gammaIndex]*tmpY[gammaIndex]);
    tmpT[gammaIndex] = kl.pi0()[i].g2().t();
    tmpChi[gammaIndex] = kl.pi0()[i].g2().chisq();
    gammaIndex++;
  }
  for( int i = 0; i< gammaIndex; i++){
    if( tmpE[i] < GMinEnergy ){
      GMinEnergy = tmpE[i];
    }
    if( TMath::Abs(tmpY[i]) > GMaxY ){
      GMaxY = TMath::Abs(tmpY[i]);
    }
    if( TMath::Abs(tmpY[i]) < GMinY ){
      GMinY = TMath::Abs(tmpY[i]);
    }
    if( TMath::Abs(tmpX[i]) < GMinX ){
      GMinX = TMath::Abs(tmpX[i]);
    }
    if( tmpChi[i] > GMaxChi ){
      GMaxChi = tmpChi[i];
    }
    if( i < gammaIndex -1 ){
      for( int j = i+1; j< gammaIndex; j++){
	double R = TMath::Sqrt( (tmpX[i]-tmpX[j])*(tmpX[i]-tmpX[j])+(tmpY[i]-tmpY[j])*(tmpY[i]-tmpY[j]));
	if( R < GMinDist ){
	  GMinDist = R;
	}
      }
    }
  }
  */
  Decision( glist );
}
void GammaCut::Decision( std::vector<Klong> kl ){
  Decision( kl[0] );
}
void GammaCut::Decision( std::list<Gamma> g ){
  Reset();
  std::list<Gamma>::iterator git = g.begin();
  int    gammaIndex =0;
  double tmpE[16];
  double tmpX[16];
  double tmpY[16];
  double tmpR[16];
  double tmpT[16];
  double tmpChi[16];
  double tmpMeanDelta[16];

  for( ; git != g.end(); git++){    
    if( gammaIndex == 16 ){ break;}
    tmpE[gammaIndex] = (*git).e();
    tmpX[gammaIndex] = (*git).x();
    tmpY[gammaIndex] = (*git).y();
    tmpR[gammaIndex] = TMath::Sqrt( tmpX[gammaIndex]*tmpX[gammaIndex] + tmpY[gammaIndex]*tmpY[gammaIndex]);
    tmpT[gammaIndex] = (*git).t();
    tmpChi[gammaIndex] = (*git).chisq();
    gammaIndex++;
  }

  for( int i = 0; i < gammaIndex; i++){
    tmpMeanDelta[i] = 0;
    for( int j = 0; j< gammaIndex; j++){
      if( i==j ){ continue; }
      tmpMeanDelta[i] = tmpT[j] - tmpT[i];
    }
    tmpMeanDelta[i] = tmpMeanDelta[i]/( gammaIndex - 1 );
  }

  for( int i = 0; i< gammaIndex; i++){
    if( tmpE[i] < GMinEnergy ){
      GMinEnergy = tmpE[i];
    }
    if( TMath::Abs(tmpY[i]) > GMaxY ){
      GMaxY = TMath::Abs(tmpY[i]);
    }
    if( TMath::Abs(tmpY[i]) < GMinY ){
      GMinY = TMath::Abs(tmpY[i]);
    }
    if( TMath::Abs(tmpX[i]) < GMinX ){
      GMinX = TMath::Abs(tmpX[i]);
    }
    if( tmpChi[i] > GMaxChi ){
      GMaxChi = tmpChi[i];
    }
    if( TMath::Abs(tmpMeanDelta[i]) > TMath::Abs(GMaxDeltaT) ){
      GMaxDeltaTID = i;
      GMaxDeltaT = tmpMeanDelta[i];
    }
    if( i < gammaIndex -1 ){
      for( int j = i+1; j< gammaIndex; j++){
	double R = TMath::Sqrt( (tmpX[i]-tmpX[j])*(tmpX[i]-tmpX[j])+(tmpY[i]-tmpY[j])*(tmpY[i]-tmpY[j]));
	if( R < GMinDist ){
	  GMinDist = R;
	}
      }
    }
  }

}
void GammaCut::SetBranchAddress( TTree* tr ){
  tr->SetBranchAddress("GMinEnergy",&GMinEnergy);
  tr->SetBranchAddress("GMinDist",&GMinDist);
  tr->SetBranchAddress("GMinX",&GMinX);
  tr->SetBranchAddress("GMinY",&GMinY);
  tr->SetBranchAddress("GMaxR",&GMaxR);
  tr->SetBranchAddress("GMaxY",&GMaxR);
  tr->SetBranchAddress("GMaxDeltaT",&GMaxDeltaT);
  tr->SetBranchAddress("GMaxDeltaTID",&GMaxDeltaTID);
}
void GammaCut::Branch( TTree* tr ){
  tr->Branch("GMinEnergy"  ,&GMinEnergy  ,"GMinEnergy/D");
  tr->Branch("GMinDist"    ,&GMinDist    ,"GMinDist/D");
  tr->Branch("GMinX"       ,&GMinX       ,"GMinX/D");
  tr->Branch("GMinY"       ,&GMinY       ,"GMinY/D");
  tr->Branch("GMaxR"       ,&GMaxR       ,"GMaxR/D");
  tr->Branch("GMaxY"       ,&GMaxY       ,"GMaxY/D");
  tr->Branch("GMaxDeltaT"  ,&GMaxDeltaT  ,"GMaxDeltaT/D");
  tr->Branch("GMaxDeltaTID",&GMaxDeltaTID,"GMaxDeltaTID/I");
}

CsiCut::CsiCut(){
  l1         = new L1TrigCounter();
  CIDHandler = new CrateIDHandler();
  map        = CsiMap::getCsiMap();
  hisCsiTime = new TH1D("hisCsiTime","hisCsiTime",500,0,500); 
  mWTimeWindow = 10;//ns, initial value.
  /*
  CsiL1TrigCountThreshold={1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
			   1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};
  */
  l1->ReadMapFile();
  l1->SetThreshold(1800);
  l1->Reset();
}
CsiCut::~CsiCut(){
  delete hisCsiTime;
}
void CsiCut::Reset(){
  hisCsiTime->Reset();
  l1->Reset();
  CsiEventTime  = 0;
  CsiEventTimeSigma = 0;
  CsiNumber = 0;
  for( int i = 0; i< 2716; i++){
    CsiID[i] = -1;
    CsiEne[i] = 0;
    CsiTime[i] = 0;
    CsiSignal[i] = 0;
    CsiChisq[i] = 0;
    CsiNDF[i] = 0;
    CsiCrate[i] = -1;
    CsiL1[i] = -1;
    CsiGB[i] = 0;
    CsiPosID[i] = 0;
  }
  CsiL1nTrig = 0;
  for( int i = 0; i< 20; i++){
    CsiL1TrigCount[i] = 0;
  }

}

void CsiCut::Decision( int ICsiNumber, int* ICsiID, double* ICsiEne, double* ICsiTime , double* ICsiSignal, double* ICsiChisq, short* ICsiNDF){
  Reset();
  for( int i = 0; i< ICsiNumber; i++){
    hisCsiTime->Fill( ICsiTime[i],ICsiEne[i] );    
  }
  CsiEventTime = hisCsiTime->GetBinCenter(hisCsiTime->GetMaximumBin());
  CsiEventTimeSigma= hisCsiTime->GetRMS();
  for( int i = 0; i< ICsiNumber; i++){
    if( ICsiSignal[i] < 5 || ICsiEne[i] < 0.5){ continue; }
    if( TMath::Abs( ICsiTime[i] - CsiEventTime ) < mWTimeWindow ){
      CsiID[CsiNumber]= ICsiID[i];
      CsiEne[CsiNumber]= ICsiEne[i];
      CsiTime[CsiNumber]=ICsiTime[i];
      CsiSignal[CsiNumber]=ICsiSignal[i];
      CsiChisq[CsiNumber] =ICsiChisq[i];
      CsiNDF[CsiNumber]   =ICsiNDF[i];
      //std::cout<< CsiNumber << "\t" << CsiID[CsiNumber] << std::endl;
      double x(0),y(0);
      x = map->getX( CsiID[CsiNumber] );
      y = map->getY( CsiID[CsiNumber] );
      if( x > 0 ){
	if( y > 0 ){ CsiPosID[CsiNumber]  = 0; }
	else{ CsiPosID[CsiNumber] = 1; }
      }else{
	if( y > 0 ){ CsiPosID[CsiNumber] = 2; }
	else{ CsiPosID[CsiNumber]  =3; }
      }
      if( y > 0 ){
	if( x < -100 ){
	  CsiGB[CsiNumber] = -1;
	}else{
	  CsiGB[CsiNumber] = 1;
	}
      }else{
	if( x < 75 ){
	  CsiGB[CsiNumber] = -1;
	}else{
	  CsiGB[CsiNumber] = 1;
	}
      }
      
      CsiCrate[CsiNumber] = CIDHandler->GetCrate(CsiID[CsiNumber]);
      CsiL1[CsiNumber]    = CIDHandler->GetL1(CsiID[CsiNumber]);
      CsiNumber++;
    }
  }
  ///////////////////////////////////////////////////////////////////////
  //L1 Trigger Counter
  for( int i = 0; i< CsiNumber; i++){
    l1->Fill( CsiID[i], CsiSignal[i]);
  }
  
  std::vector<double> vecCount = l1->GetCount();
  if( vecCount.size() != 20 ){ std::cout<< "VectorSize Error" << std::endl;}
  for( int i = 0; i< 20; i++){
    CsiL1TrigCount[i] = vecCount.at(i);
    if( CsiL1TrigCount[i] > CsiL1TrigCountThreshold[i] ){
      CsiL1nTrig++;
    }
  }
}


void CsiCut::SetCutValue( double wTimeWindow ){
  mWTimeWindow = wTimeWindow;
}
void CsiCut::Branch( TTree* trout){
    trout->Branch("CsiNumber"  ,&CsiNumber   ,"CsiNumber/I");
    trout->Branch("CsiModID"   ,CsiID        ,"CsiModID[CsiNumber]/I");//CsiNumber
    trout->Branch("CsiEne"     ,CsiEne       ,"CsiEne[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiTime"    ,CsiTime      ,"CsiTime[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiHHTime"  ,CsiHHTime    ,"CsiHHTime[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiSignal"  ,CsiSignal    ,"CsiSignal[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiChisq"   ,CsiChisq     ,"CsiChisq[CsiNumber]/D");//CsiNumber
    trout->Branch("CsiNDF"     ,CsiNDF       ,"CsiNDF[CsiNumber]/S");//CsiNumber
    trout->Branch("CsiPosID"   ,CsiPosID     ,"CsiPosID[CsiNumber]/S");//CsiNumber
    trout->Branch("CsiGB"      ,CsiGB        ,"CsiGB[CsiNumber]/S");//CsiNumber
    trout->Branch("CsiCrate"   ,CsiCrate     ,"CsiCrate[CsiNumber]/S");//CsiNumber
    
    trout->Branch("CsiL1nTrig" ,&CsiL1nTrig  ,"CsiL1nTrig/I");
    trout->Branch("CsiL1TrigCount",CsiL1TrigCount,"CsiL1TrigCount[20]/D");
    
    trout->Branch("CsiEventTime",&CsiEventTime,"CsiEventTime/D");
    trout->Branch("CsiEventTimeSigma",&CsiEventTimeSigma,"CsiEventTimeSigma/D");
}

void CsiCut::SetBranchAddress( TTree* tr ){

}
