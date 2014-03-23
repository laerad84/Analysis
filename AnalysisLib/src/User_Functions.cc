#include "User_Functions.h"
double Rec_Mass2g( std::list<Gamma> const &glist, double recPosition ){
  double mass=0;
  CLHEP::Hep3Vector Vtx(0,0,recPosition);
  int id = 0;
  if( glist.size() != 2 ){ return 0; }
  Gamma g1 = glist.front();
  Gamma g2 = glist.back();
  CLHEP::Hep3Vector GPos[2];
  GPos[0] = g1.pos() -Vtx;
  GPos[1] = g2.pos() -Vtx;
  GPos[0].setMag(g1.e());
  GPos[1].setMag(g2.e());
  double MassSq = TMath::Power((g1.e() + g2.e()),2) - (GPos[0]-GPos[1]).mag2();// (E1+E2)^2 - G1*G2 
  mass = TMath::Sqrt( MassSq );
  return mass;
}
void RecVtx_ConstM( const Gamma& g1, const Gamma& g2, double Mass, double* recZ, double* recZsig2 ){
  //std::cout << "rec2g::recVtxWithConstM() : Mass=" << Mass << std::endl;

  double r1 = g1.pos().perp2();
  double r2 = g2.pos().perp2();
  double ip = g1.x()*g2.x() + g1.y()*g2.y();
  double ct = 1. - Mass*Mass/(2.*g1.e()*g2.e());

  double A = 1. - ct*ct;
  double B = 2.*ip - (r1 + r2)*ct*ct;
  double C = ip*ip - r1*r2*ct*ct;

  double D = B*B - 4.*A*C;



  if( D > 0 ) {
    double dz1 = (-B + sqrt(D))/(2.*A);
    double dz2 = (-B - sqrt(D))/(2.*A);

    //    std::cout<<"DEBUG : "<<"E : "<<g1.e()<<"\t sigmaE : "<<g1.sigmaE()<<std::endl; 
    //  getchar();
    //    std::cout<<"DEBUG : "<<"sigamE : "<<g2.sigmaE()<<"\t sigmaX : "<<g2.sigmaX()<<std::endl; 
    //
    double r12 = (g1.pos() - g2.pos()).perp();

    double rerr_e1   = g1.sigmaE()/g1.e();
    double rerr_e2   = g2.sigmaE()/g2.e();

    //
    double sig2r1 = 
      4.*( g1.x()*g1.x()*g1.sigmaX()*g1.sigmaX() +
	   g1.y()*g1.y()*g1.sigmaY()*g1.sigmaY() );
    double sig2r2 = 
      4.*( g2.x()*g2.x()*g2.sigmaX()*g2.sigmaX() +
	   g2.y()*g2.y()*g2.sigmaY()*g2.sigmaY() );
    double sig2ip = 
      g2.x()*g2.x()*g1.sigmaX()*g1.sigmaX() +
      g1.x()*g1.x()*g2.sigmaX()*g2.sigmaX() +
      g2.y()*g2.y()*g1.sigmaY()*g1.sigmaY() +
      g1.y()*g1.y()*g2.sigmaY()*g2.sigmaY();
    double sig2ct = 
      ( (Mass*Mass)/(2*g1.e()*g2.e()) )*( (Mass*Mass)/(2*g1.e()*g2.e()) )*
      ( rerr_e1*rerr_e1 + rerr_e2*rerr_e2 );


    double DA_Dct = -2.*ct;
    double sig2_A = DA_Dct*DA_Dct*sig2ct;
    
    double DB_Dip = 2;
    double DB_Dr1 = -1.*ct*ct;
    double DB_Dr2 = -1.*ct*ct;
    double DB_Dct = -2.*(r1 + r2)*ct;
    double sig2_B = 
      DB_Dip*DB_Dip*sig2ip + DB_Dr1*DB_Dr1*sig2r1 + DB_Dr2*DB_Dr2*sig2r2 + DB_Dct*DB_Dct*sig2ct;

    double DC_Dip = 2*ip;
    double DC_Dr1 = -1.*r2*ct*ct;
    double DC_Dr2 = -1.*r1*ct*ct;
    double DC_Dct = -2.*r1*r2*ct;
    double sig2_C = 
      DC_Dip*DC_Dip*sig2ip + DC_Dr1*DC_Dr1*sig2r1 + DC_Dr2*DC_Dr2*sig2r2 + DC_Dct*DC_Dct*sig2ct;
      
    
    double DD_DB  =  2.*B;
    double DD_DA  = -4.*C;
    double DD_DC  = -4.*A;
    double sig2_D = DD_DB*DD_DB*sig2_B + DD_DA*DD_DA*sig2_A + DD_DC*DD_DC*sig2_C;


    if( dz1 > 0 ) {
      double Ddz1_DB  = -1./(2.*A);
      double Ddz1_DA  = -1.*(-1.*B + sqrt(D))/(2.*A*A);
      double Ddz1_DD  = +1./(4.*A*sqrt(D));
      double sig2_dz1 = Ddz1_DB*Ddz1_DB*sig2_B + Ddz1_DA*Ddz1_DA*sig2_A + Ddz1_DD*Ddz1_DD*sig2_D;

      recZ[0]     = g1.z() - sqrt(dz1);
      recZsig2[0] = (-1./(2.*sqrt(dz1)))*(-1./(2.*sqrt(dz1)))*sig2_dz1;

    }
    if( dz2 > 0 ) {
      double Ddz2_DB  = -1./(2.*A);
      double Ddz2_DA  = -1.*(-1.*B - sqrt(D))/(2.*A*A);
      double Ddz2_DD  = -1./(4.*A*sqrt(D));
      double sig2_dz2 = Ddz2_DB*Ddz2_DB*sig2_B + Ddz2_DA*Ddz2_DA*sig2_A + Ddz2_DD*Ddz2_DD*sig2_D;

      recZ[1]     = g1.z() - sqrt(dz2);
      recZsig2[1] = (-1./(2.*sqrt(dz2)))*(-1./(2.*sqrt(dz2)))*sig2_dz2;
    }
  }
}

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

void GammaTimeDeltaCut( std::list<Gamma> glist, std::list<Gamma>& glistOut , double TimeThreshold ){
  int nGamma=0;
  double GTime[20];
  double GTimeDeltaMean[20];
  for( int i = 0; i<20 ; i++){
    GTime[i] = 0;
    GTimeDeltaMean[i] = 0;
  }
  std::list<Gamma>::iterator git = glist.begin();
  for( ; git != glist.end(); git++){
    GTime[nGamma] = (*git).t();
    nGamma++;
  }
  //std::cout<< nGamma << std::endl;

  for( int i = 0; i< nGamma; i++){    
    for( int j = 0; j< nGamma; j++){
      if( i == j){ continue; }
      GTimeDeltaMean[i] += GTime[j]-GTime[i];
    }
    GTimeDeltaMean[i] = GTimeDeltaMean[i]/(nGamma-1);
  }
  git = glist.begin();
  int gIndex = 0;
  for( ; git != glist.end(); git++){
    if( abs( GTimeDeltaMean[gIndex] ) < TimeThreshold ){ 
      glistOut.push_back((*git));
    }
    gIndex++;
  }  
}

void GammaTimeDeltaCutEventTime( std::list<Gamma> glist, std::list<Gamma>& glistOut , double EventTime, double TimeThreshold ){
  int nGamma=0;
  double GTime[20];
  double GTimeDeltaMean[20];
  double MinDeltaTime=500;
  double SecMinDeltaTime=500;
  for( int i = 0; i<20 ; i++){
    GTime[i] = 0;
    GTimeDeltaMean[i] = 0;
  }

  std::list<Gamma>::iterator git = glist.begin();
  for( ; git != glist.end(); git++){
    GTime[nGamma] = (*git).t();
    nGamma++;
  }
  //std::cout<< nGamma << std::endl;

  for( int i = 0; i< nGamma; i++){    
    for( int j = 0; j< nGamma; j++){
      if( i == j){ continue; }
      GTimeDeltaMean[i] += GTime[j]-GTime[i];
    }
    GTimeDeltaMean[i] = GTimeDeltaMean[i]/(nGamma-1);
  }
  git = glist.begin();
  int gIndex = 0;
  for( ; git != glist.end(); git++){
    if( abs( GTime[gIndex]-EventTime ) < TimeThreshold ){ 
      glistOut.push_back((*git));
    }
    gIndex++;
  }  
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
  GMinXY      = 10000;
  GMaxR      = 0;
  GMaxY      = 0;
  GMaxChi    = 0;
  GMaxDeltaT = 0;
  GMaxDeltaTID = -1;
}
void GammaCut::Decision(Klong kl){
  Reset();
  std::list<Gamma> glist;
  for( int i = 0; i< kl.pi0().size();i++){
    glist.push_back(kl.pi0()[i].g1());
    glist.push_back(kl.pi0()[i].g2());
  }
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
      tmpMeanDelta[i] += tmpT[j] - tmpT[i];
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

    /*
    if( TMath::Abs(tmpY[i]) < GMinY ){
      GMinY = TMath::Abs(tmpY[i]);
    }
    if( TMath::Abs(tmpX[i]) < GMinX ){
      GMinX = TMath::Abs(tmpX[i]);
    }
    */

    if( TMath::Abs( tmpY[i] ) > TMath::Abs( tmpX[i] ) ){
      if( TMath::Abs( tmpY[i] )< GMinXY ){
	GMinXY  =TMath::Abs( tmpY[i] );
      }
    }else{
      if( TMath::Abs( tmpX[i] ) < GMinXY ){
	GMinXY  =TMath::Abs( tmpX[i] );
      }
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
  tr->SetBranchAddress("GMinXY",&GMinXY);
  //tr->SetBranchAddress("GMinX",&GMinX);
  //tr->SetBranchAddress("GMinY",&GMinY);
  tr->SetBranchAddress("GMaxR",&GMaxR);
  tr->SetBranchAddress("GMaxY",&GMaxR);
  tr->SetBranchAddress("GMaxDeltaT",&GMaxDeltaT);
  tr->SetBranchAddress("GMaxDeltaTID",&GMaxDeltaTID);
}
void GammaCut::Branch( TTree* tr ){
  tr->Branch("GMinEnergy"  ,&GMinEnergy  ,"GMinEnergy/D");
  tr->Branch("GMinDist"    ,&GMinDist    ,"GMinDist/D");
  tr->Branch("GMinXY"      ,&GMinXY      ,"GMinXY/D");
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
void CsiCut::DecisionForPi0Run( int ICsiNumber, int* ICsiID, double* ICsiEne, double* ICsiTime , double* ICsiSignal, double* ICsiChisq, short* ICsiNDF){
  Reset();
  for( int i = 0; i< ICsiNumber; i++){
    hisCsiTime->Fill( ICsiTime[i],ICsiEne[i] );    
  }
  CsiEventTime = hisCsiTime->GetBinCenter(hisCsiTime->GetMaximumBin());
  CsiEventTimeSigma= hisCsiTime->GetRMS();
  for( int i = 0; i< ICsiNumber; i++){
    if( ICsiSignal[i] < 5 || ICsiEne[i] < 0.5){ continue; }
    if( TMath::Abs( ICsiTime[i] - CsiEventTime ) < mWTimeWindow ){
      CsiID[CsiNumber]    = ICsiID[i];
      CsiEne[CsiNumber]   = ICsiEne[i];
      CsiTime[CsiNumber]  = ICsiTime[i];
      CsiSignal[CsiNumber]= ICsiSignal[i];
      CsiChisq[CsiNumber] = ICsiChisq[i];
      CsiNDF[CsiNumber]   = ICsiNDF[i];
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
    if( CsiL1TrigCount[i] > CsiL1TrigCountThresholdPi0[i] ){
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
  tr->SetBranchAddress("CsiNumber"  ,&CsiNumber);
  tr->SetBranchAddress("CsiModID"   ,CsiID        );//CsiNumber
  tr->SetBranchAddress("CsiEne"     ,CsiEne       );//CsiNumber
  tr->SetBranchAddress("CsiTime"    ,CsiTime      );//CsiNumber
  tr->SetBranchAddress("CsiHHTime"  ,CsiHHTime    );//CsiNumber
  tr->SetBranchAddress("CsiSignal"  ,CsiSignal    );//CsiNumber
  tr->SetBranchAddress("CsiChisq"   ,CsiChisq     );//CsiNumber
  tr->SetBranchAddress("CsiNDF"     ,CsiNDF       );//CsiNumber
  tr->SetBranchAddress("CsiPosID"   ,CsiPosID     );//CsiNumber
  tr->SetBranchAddress("CsiGB"      ,CsiGB        );//CsiNumber
  tr->SetBranchAddress("CsiCrate"   ,CsiCrate     );//CsiNumber
  
  tr->SetBranchAddress("CsiL1nTrig" ,&CsiL1nTrig  );
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  
  tr->SetBranchAddress("CsiEventTime",&CsiEventTime);
  tr->SetBranchAddress("CsiEventTimeSigma",&CsiEventTimeSigma);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// KLCut
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

KLCut::KLCut(){
  ;
}
KLCut::~KLCut(){
  ;
}
void KLCut::Branch(TTree* tr){
  tr->Branch("Pi0PtMax",&Pi0PtMax,"Pi0PtMax/D");  
}
void KLCut::SetBranchAddress(TTree* tr ){
  tr->SetBranchAddress("Pi0PtMax",&Pi0PtMax);
}
void KLCut::Reset(){
  Pi0PtMax = 0;
}
void KLCut::Decision(Klong kl){
  Reset();
  for( int i = 0; i< kl.pi0().size(); i++){
    double pi0pt = TMath::Sqrt(kl.pi0()[i].p3()[0]*kl.pi0()[i].p3()[0]+kl.pi0()[i].p3()[1]*kl.pi0()[i].p3()[1]);
    if( pi0pt > Pi0PtMax ){
      Pi0PtMax = pi0pt;
    }
  }
}
void KLCut::Decision( std::vector<Klong> kl){
  Decision(kl[0]);
}
