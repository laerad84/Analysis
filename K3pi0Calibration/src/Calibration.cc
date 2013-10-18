#include "Calibration.h"
////

Calibration::Calibration()
{
  InitValue();
}

Calibration::~Calibration(){
  std::cout << "Calibration is end" << std::endl; 
}

int Calibration::InitValue(){
  KL_Vec = std::vector<Klong> ();
  m_FlagKL_prefit = -1; 
  for( int i = 0; i< 6; ++i){
    m_FlagKL[i]      = -1; 
    m_FlagCluster[i] = -1;
    m_FlagCalibrated[i] = -1; 
    m_CorrID[i]      = 0;
    m_CorrE[i]       = 0;
    m_Corr[i]        = 0;
    m_Ratio[i]       = 0; 
    m_SecondRatio[i] = 0; 
    m_GammaEnergy[i] = 0; 
    m_GammaSigma[i]  = 0;
    m_LeadingEnergy[i] = 0;
    m_LeadingHeight[i] = 0;
    m_LeadingChID[i] = -1;
    m_chisq[i]       = 0;
 }
  m_nCalibrated = 0;
  return  0;
}

int Calibration::InitValue(Klong &KL_init){
  KL_Vec = std::vector<Klong> ();
  m_FlagKL_prefit = -1; 
  for( int i = 0; i< 6; ++i){
    m_FlagKL[i]      = -1; 
    m_FlagCluster[i] = -1;
    m_FlagCalibrated[i] = -1; 
    m_CorrID[i]      = 0;
    m_CorrE[i]       = 0;
    m_Corr[i]        = 0;
    m_Ratio[i]       = 0; 
    m_SecondRatio[i] = 0; 
    m_GammaEnergy[i] = 0; 
    m_GammaSigma[i]  = 0;
    m_LeadingEnergy[i] = 0;
    m_LeadingHeight[i] = 0;
    m_LeadingChID[i] = -1;
    m_chisq[i]       = 0;
 }
  m_nCalibrated = 0;
  KL_prefit = KL_init;
  return  0;
}

int Calibration::InitValue(std::vector<Klong> &KL_VecInit){
  KL_Vec    = KL_VecInit;
  KL_prefit = KL_Vec[0];
  m_FlagKL_prefit = -1; 
  for( int i = 0; i< 6; ++i){
    m_FlagKL[i]      = -1; 
    m_FlagCluster[i] = -1;
    m_FlagCalibrated[i] = -1; 
    m_CorrID[i]      = 0;
    m_CorrE[i]       = 0;
    m_Corr[i]        = 0;
    m_Ratio[i]       = 0; 
    m_SecondRatio[i] = 0; 
    m_GammaEnergy[i] = 0; 
    m_GammaSigma[i]  = 0;
    m_LeadingEnergy[i] = 0;
    m_LeadingChID[i] = -1;
    m_chisq[i]       = 0;
  }
  m_nCalibrated = 0;
  return  0;
}

void Calibration::GetResult(CalibrationTree &calib){
  calib.FlagKL_prefit = m_FlagKL_prefit;
  for( int i = 0; i< 6; ++i ){
    calib.FlagKL[i]         = m_FlagKL[i];
    calib.FlagCluster[i]    = m_FlagCluster[i];
    calib.FlagCalibrated[i] = m_FlagCalibrated[i];
    calib.CorrID[i]         = m_CorrID[i];
    calib.CorrE[i]          = m_CorrE[i];
    calib.Corr[i]           = m_Corr[i];
    calib.Ratio[i]          = m_Ratio[i];
    calib.SecondRatio[i]    = m_SecondRatio[i];
    calib.GammaEnergy[i]    = m_GammaEnergy[i];
    calib.GammaSigma[i]     = m_GammaSigma[i];
    calib.LeadingChID[i]    = m_LeadingChID[i];
    calib.LeadingEnergy[i]  = m_LeadingEnergy[i];
    calib.chisq[i]          = m_chisq[i];
  }
  calib.nCalibrated         = m_nCalibrated;
}

void Calibration::Print(){
  std::cout<< m_nCalibrated << "xGamma used for Calibration" << "\n";
  std::cout<< KL_prefit.pi0().size() <<std::endl;
  std::cout<< KL.pi0().size()        << std::endl;
  std::cout<< "FlagKL_prefit : " <<  m_FlagKL_prefit  << "\n";
  std::cout<< "FlagKL:FlagCluster:FlagClaibrated:CorrID:Corr:Chisq\n"; 
  for( int i = 0 ; i < 6; ++i){
    std::cout<< m_FlagKL[i]         << " : " 
	     << m_FlagCluster[i]    << " : " 
	     << m_FlagCalibrated[i] << " : "
	     << m_CorrID[i]         << " : "
	     << m_Corr[i]           << " : "
	     << m_GammaEnergy[i]    << " : "
	     << m_GammaSigma[i]     << " : " 
	     << m_chisq[i]          << std::endl;
  }
}


//int Calibration::CalEnergy_idv(Klong &KL_init){
int Calibration::CalEnergy_idv(std::vector<Klong> &KL_VecInit){
  int rst = -1;
  rst = InitValue(KL_VecInit);
  if( rst != 0 ){ return m_nCalibrated;}
  rst = SelectKL_prefit();
  if( rst != 0 ){ return m_nCalibrated;}

  //std::cout <<"RST:"<< rst << std::endl;
  //if( rst != 0 ){ return rst; }
  //////////////////////////////////////////////
  // loop for all pi0s
  for( unsigned int idx = 0; idx < KL_prefit.pi0().size() ; idx++ ){
    for( int GammaIndex = 0; GammaIndex <2; ++GammaIndex){
      //Calibration//
      // Sort Gamma -> Calibrate -> Select Vaild Event // 
      rst = SortGamma(idx, GammaIndex);
      if( rst != 0){ continue; }
      m_chisq[idx*2+GammaIndex] = Calibrate();      
      rst = SelectKL(idx,GammaIndex);
      if( rst != 0){ continue; }
      rst =  SelectCluster(idx,GammaIndex);
      if( rst != 0){ continue; }
      
      ++m_nCalibrated;
      m_FlagCalibrated[idx*2+GammaIndex] = 0;      
    }
  }
  return m_nCalibrated;
}


int Calibration::SortGamma( int idx ,int GammaIndex){
  KL = KL_prefit;
  // std::cout<< "SortGamma : " << KL.pi0().size() << std::endl;
  KL.pi0().clear();
  for( unsigned int i = 0; i< KL_prefit.pi0().size() ; ++i ){
    if( i!= idx ){
      KL.pi0().push_back( KL_prefit.pi0()[i] );
    }
  }
  KL.pi0().push_back( KL_prefit.pi0()[idx] );
  if( GammaIndex == 0){//first Gamma
    KL.pi0()[KL.pi0().size() -1].setGamma( KL_prefit.pi0()[idx].g2(),KL_prefit.pi0()[idx].g1());
  }else{
    //if GammaIndex == 1
    //no needs for exchange gamma
    KL.pi0()[KL.pi0().size() -1].setGamma( KL_prefit.pi0()[idx].g1(),KL_prefit.pi0()[idx].g2());
  }      
  if( GammaIndex != 0 && GammaIndex != 1){ 
    return -1; 
  }
  return 0;
}

int Calibration::SelectKL_prefit(){  

  m_FlagKL_prefit = 0;
  CLHEP::Hep3Vector pos[6];
  //Test Klong
  if( KL_prefit.chisqZ() > chisqCut ){ m_FlagKL_prefit |= 1;}//Cut Value initial = 10;
  if( KL_Vec.size() > 1 ){
    if( KL_Vec[1].chisqZ() - KL_Vec[0].chisqZ() < deltachisqCut){
      m_FlagKL_prefit |= 2; 
    }
  }
  
  if( fabs( KL_prefit.m() - MASS_KL ) > KLMassCut ){
    m_FlagKL_prefit |=4;
  }//Cut Value initial=10;
  if( KL_prefit.vz() > 5000 || KL_prefit.vz() <3000){
    m_FlagKL_prefit |= 8;
  }
  if( KL_prefit.pi0().size() != 3){ m_FlagKL_prefit |=16;}
  if( KL_prefit.pi0()[0].status() != 1 ||
      KL_prefit.pi0()[1].status() != 1 ||
      KL_prefit.pi0()[2].status() != 1 ){
    m_FlagKL_prefit |= 32 ;
  }

  for( unsigned int idx  = 0; idx < KL_prefit.pi0().size(); ++idx){
    if( fabs( KL_prefit.pi0()[idx].m() -MASS_PI0 ) > pimassCut){ m_FlagKL_prefit |= 64;}
    pos[idx*2]   = KL_prefit.pi0()[idx].g1().pos();
    pos[idx*2+1] = KL_prefit.pi0()[idx].g2().pos();    
  }

  for( int i = 0; i < 6; ++i ){
    for( int j = i+1; j < 6; ++j){
      if(( pos[i] - pos[j] ).perp() < gammadistCut){ m_FlagKL_prefit |= 128;}
    }
  }  
  
  for( unsigned int idx = 0; idx < KL_prefit.pi0().size(); ++idx){
    if( KL_prefit.pi0()[idx].g1().e() < gammaEneCut ){
      m_FlagKL_prefit |= 256;
    }
    if( KL_prefit.pi0()[idx].g2().e() < gammaEneCut ){
      m_FlagKL_prefit |= 256;
    }
  }

  return m_FlagKL_prefit;  
}

int Calibration::SelectKL(int idx, int GammaIndex){
  int Flag = 0;
  if( m_chisq[idx*2 + GammaIndex]/2. > gammaChisqCut ){ Flag |= 1; }
  if( fabs( KL.m() - MASS_KL) >  KLMassCut ){ Flag |= 2; }
  if( fabs( KL.pi0()[KL.pi0().size() -1 ].m() - MASS_PI0 ) > pimassCut ){ Flag |= 4; }
  if( fabs( KL.vz() - KL_prefit.vz() ) > deltaCalibrationZCut ){ Flag |= 8;}
  m_FlagKL[idx*2 + GammaIndex] = Flag;
  return   m_FlagKL[idx*2 + GammaIndex];
  
} 

int Calibration::SelectCluster( int idx , int GammaIndex){
  std::vector<int>::const_iterator i_id;
  std::vector<int>::const_iterator j_id;
  std::vector<int>::const_iterator ii_id;
  std::vector<int>::const_iterator jj_id;
  std::vector<double>::const_iterator i_e;
  std::vector<double>::const_iterator j_e;
  std::vector<double>::const_iterator ii_e;
  std::vector<double>::const_iterator jj_e;
  
  if( GammaIndex != 0 && GammaIndex != 1){
    m_FlagCluster[2*idx+GammaIndex] = 1;
    return -1;
  }
  
  m_FlagCluster[2*idx+GammaIndex] = 0; 
  
  double clusterNoCorrE = 0;
  double CorrFactor=0;
  double Corr = 0.;
  double err  = 0.;
  int    Flag = 0;
  
  // If there is the channel which is more than 0.2*Cluster Energy even if not leading channel.
  // Ignore the event 
  // leading channel is more than 0.2*Cluster Energy 
  // second leading channel is less than 0.2*Cluster Energy
  //  std::cout<< "Select Cluster" << std::endl;
  
  if( GammaIndex == 0){//Gamma1
    clusterNoCorrE = 0;
    for( jj_e = KL_prefit.pi0()[idx].g1().clusterEVec().begin();
	 jj_e != KL_prefit.pi0()[idx].g1().clusterEVec().end();
	 ++jj_e){
      clusterNoCorrE += *jj_e;
    }
    CorrFactor = KL_prefit.pi0()[idx].g1().e() - clusterNoCorrE;
    Corr = (KL.pi0()[KL.pi0().size() -1 ].g2().e() - KL_prefit.pi0()[idx].g1().e())/clusterNoCorrE + 1;
    //Corr = (KL.pi0()[KL.pi0().size() -1 ].g2().e() -CorrFactor)/(KL_prefit.pi0()[idx].g1().e() - CorrFactor);
    err = KL.pi0()[ KL.pi0().size() -1 ].g2().sigmaE() / KL_prefit.pi0()[idx].g1().e();

    /*
      std::cout<< "Check"<< std::endl;
      std::cout<< clusterNoCorrE                << std::endl;
      std::cout<< KL_prefit.pi0()[idx].g1().e() << std::endl;
      std::cout<< KL_prefit.pi0()[KL.pi0().size() -1 ].g2().e() << std::endl;
    */

    for( i_id = KL_prefit.pi0()[idx].g1().clusterIdVec().begin(),
	   i_e = KL_prefit.pi0()[idx].g1().clusterEVec().begin();
	 i_id != KL_prefit.pi0()[idx].g1().clusterIdVec().end();
	 ++i_id,++i_e){
      
      if( *i_id != *KL_prefit.pi0()[idx].g1().clusterIdVec().begin() &&
	  *i_e  >= KL_prefit.pi0()[idx].g1().e()*threshold){
	Flag |= 2;
	
      }

      // Check leading crystal Energy domination
      
      if( *i_id == *KL_prefit.pi0()[idx].g1().clusterIdVec().begin()){
	if( *i_e  >= KL_prefit.pi0()[idx].g1().e()*threshold ){
	  m_Corr[idx*2   + GammaIndex] = Corr;
	  m_CorrID[idx*2 + GammaIndex] = *i_id;
	  //m_CorrE[idx*2  + GammaIndex] = *i_e;
	  m_CorrE[idx*2 + GammaIndex] = KL.pi0()[KL.pi0().size() -1].g2().e();
	  m_Ratio[idx*2  + GammaIndex] = *i_e / KL_prefit.pi0()[idx].g1().e();
	  m_GammaEnergy[idx*2 + GammaIndex] = KL_prefit.pi0()[idx].g1().e();
	  m_GammaSigma[idx*2 +GammaIndex]   = KL_prefit.pi0()[idx].g1().sigmaE();
	  m_LeadingChID[idx*2+GammaIndex]   = (*i_id);
	  m_LeadingEnergy[idx*2+GammaIndex] = (*i_e);

	}else{
	  Flag |= 4;
	}
      }

      //std::cout<< *i_id << " : " << *i_e << std::endl; 
    }
    if( Flag == 0 ){ 
      i_id = KL_prefit.pi0()[idx].g1().clusterIdVec().begin();
      i_e  = KL_prefit.pi0()[idx].g1().clusterEVec().begin();
      ++i_id;
      ++i_e;
      m_SecondRatio[idx*2+GammaIndex] =  *i_e / KL_prefit.pi0()[idx].g1().e();  
    }
    
  }else{//Gamma2
    clusterNoCorrE = 0;
    for( jj_e  = KL_prefit.pi0()[idx].g2().clusterEVec().begin();
	 jj_e != KL_prefit.pi0()[idx].g2().clusterEVec().end();
	 ++jj_e){
      clusterNoCorrE += *jj_e;
    }

    CorrFactor = KL_prefit.pi0()[idx].g2().e() - clusterNoCorrE;
    Corr = (KL.pi0()[KL.pi0().size() -1 ].g2().e() - KL_prefit.pi0()[idx].g2().e())/clusterNoCorrE + 1;
    //Corr = (KL.pi0()[KL.pi0().size() -1 ].g2().e() -CorrFactor)/(KL_prefit.pi0()[idx].g2().e() - CorrFactor); 
    err = KL.pi0()[ KL.pi0().size() -1 ].g2().sigmaE() / KL_prefit.pi0()[idx].g2().e();
    
    /*
    std::cout<< "Check"<< std::endl;
    std::cout<< clusterNoCorrE                << std::endl;
    std::cout<< KL_prefit.pi0()[idx].g1().e() << std::endl;
    std::cout<< KL_prefit.pi0()[KL.pi0().size() -1 ].g2().e() << std::endl;
    */

    for( i_id = KL_prefit.pi0()[idx].g2().clusterIdVec().begin(),
	   i_e = KL_prefit.pi0()[idx].g2().clusterEVec().begin();
	 i_id != KL_prefit.pi0()[idx].g2().clusterIdVec().end();
	 ++i_id,++i_e){
      
      if( *i_id != *KL_prefit.pi0()[idx].g2().clusterIdVec().begin() &&
	  *i_e  >= KL_prefit.pi0()[idx].g2().e()*threshold){
	Flag |= 2;
      }
      //Check leading crystal Energy domination
      if( *i_id == *KL_prefit.pi0()[idx].g2().clusterIdVec().begin()){
	if( *i_e  >= KL_prefit.pi0()[idx].g2().e()*threshold ){
	  
	  m_Corr[idx*2   + GammaIndex] = Corr;
	  m_CorrID[idx*2 + GammaIndex] = *i_id;
	  //m_CorrE[idx*2  + GammaIndex] = *i_e;
	  m_CorrE[idx*2 + GammaIndex] = KL.pi0()[KL.pi0().size() -1].g2().e();
	  m_Ratio[idx*2  + GammaIndex] = *i_e / KL_prefit.pi0()[idx].g2().e();
	  //m_GammaEnergy[idx*2 + GammaIndex] = KL_prefit.pi0()[idx].g2().e();
	  m_GammaEnergy[idx*2 + GammaIndex] = KL_prefit.pi0()[idx].g2().e();
	  m_GammaSigma[idx*2 +GammaIndex]   = KL_prefit.pi0()[idx].g2().sigmaE();
	  m_LeadingChID[idx*2+GammaIndex] = (*i_id);
	  m_LeadingEnergy[idx*2+GammaIndex] = (*i_e);
	  //std::cout<< *i_id << " : " << *i_e << std::endl; 

	}else{
	  Flag |= 4;
	}
      }      
    }
    if( Flag == 0 ){ 
      i_id = KL_prefit.pi0()[idx].g2().clusterIdVec().begin();
      i_e  = KL_prefit.pi0()[idx].g2().clusterEVec().begin();
      ++i_id;
      ++i_e;
      m_SecondRatio[idx*2+GammaIndex] =  *i_e / KL_prefit.pi0()[idx].g2().e();  
    }
  }
  m_FlagCluster[idx*2 +GammaIndex] = Flag; 
  return m_FlagCluster[idx*2 + GammaIndex];
}

double Calibration::Calibrate(){
  
  double chisq_keep   = HUGE;
  double chisq_keep_v = HUGE;
  double chisq        = -1.0;
  
  int ifail         = 0;
  int error_flag    = 0;
  int fit_updated   = 0;
  const int maxPi0  = 3;
  const int npi0s   = KL.pi0().size();

  double px1[maxPi0];
  double py1[maxPi0];
  double pz1[maxPi0];
  double E1[maxPi0];
  double x1[maxPi0];
  double y1[maxPi0];
  double z1[maxPi0];
  double R1[maxPi0];

  double px2[maxPi0];
  double py2[maxPi0];
  double pz2[maxPi0];
  double E2[maxPi0];
  double x2[maxPi0];
  double y2[maxPi0];
  double z2[maxPi0];
  double R2[maxPi0];

  double klpx, klpy, klpz, klE, klE0 = 0.;
  double vx, vy, vz;

  // Prepare initial values
  int idx;
  for( idx=0; idx<npi0s; idx++) {
    x1[idx] = KL.pi0()[idx].g1().x();
    y1[idx] = KL.pi0()[idx].g1().y();
    z1[idx] = KL.pi0()[idx].g1().z();
    E1[idx] = KL.pi0()[idx].g1().e();
    x2[idx] = KL.pi0()[idx].g2().x();
    y2[idx] = KL.pi0()[idx].g2().y();
    z2[idx] = KL.pi0()[idx].g2().z();
    E2[idx] = KL.pi0()[idx].g2().e();
    
    klE0 += E1[idx] + E2[idx];
  }
  
  vx = KL.vx();
  vy = KL.vy();
  vz = KL.vz();
  
  // a() & v()
  CLHEP::HepVector a0(npi0s*6-1),a(npi0s*6-1);
  CLHEP::HepVector v0(4),v(4);
  
  for( idx=0; idx<npi0s; idx++) {
    a0(idx*6 + 1)  = x1[idx];
    a0(idx*6 + 2)  = y1[idx];
    a0(idx*6 + 3)  = E1[idx];
    a0(idx*6 + 4)  = x2[idx];
    a0(idx*6 + 5)  = y2[idx];
    if (idx!=npi0s-1){
      a0(idx*6 + 6)  = E2[idx];
    }
  }

  v0(1) = vx;
  v0(2) = vy;
  v0(3) = vz;
  v0(4) = E2[npi0s-1];

  CLHEP::HepSymMatrix Va0(npi0s*6-1,0);  // (npi0s*6-1) x (npi0s*6-1)
  
  for( idx=0; idx<npi0s; idx++) {
    Va0(idx*6+1 ,idx*6+1) = pow(KL.pi0()[idx].g1().sigmaX(),2);
    Va0(idx*6+2 ,idx*6+2) = pow(KL.pi0()[idx].g1().sigmaY(),2);
    Va0(idx*6+3 ,idx*6+3) = pow(KL.pi0()[idx].g1().sigmaE(),2);
    Va0(idx*6+4 ,idx*6+4) = pow(KL.pi0()[idx].g2().sigmaX(),2);
    Va0(idx*6+5 ,idx*6+5) = pow(KL.pi0()[idx].g2().sigmaY(),2);
    if (idx!=npi0s-1)
      Va0(idx*6+6 ,idx*6+6) = pow(KL.pi0()[idx].g2().sigmaE(),2);
  }
  //   for( idx=0; idx<npi0s; idx++) {
  //     Va0(idx*6+1 ,idx*6+1) = 0.6*0.6;
  //     Va0(idx*6+2 ,idx*6+2) = 0.6*0.6;
  //     Va0(idx*6+3 ,idx*6+3) = 0.02*0.02/KL.pi0()[idx].g1().e()+0.02*0.02;
  //     Va0(idx*6+4 ,idx*6+4) = 0.6*0.6;
  //     Va0(idx*6+5 ,idx*6+5) = 0.6*0.6;
  //     if (idx!=npi0s-1)
  //       Va0(idx*6+6 ,idx*6+6) = 0.02*0.02/KL.pi0()[idx].g2().e()+0.02*0.02;
  //   }
  
  a = a0;
  v = v0;
  
  CLHEP::HepMatrix D(npi0s+3,npi0s*6-1,0); // = del H/del alpha(alpha_A)
  CLHEP::HepMatrix E(npi0s+3,4,0);      // = del H/del alpha(v)

  CLHEP::HepVector d(npi0s+3); // = H(alpha_A)
  
  CLHEP::HepVector da0(npi0s*6-1);
  CLHEP::HepVector dv0(4);

  CLHEP::HepMatrix VD(npi0s+3,npi0s+3);
  CLHEP::HepMatrix VE(4,4);

  CLHEP::HepMatrix Lambda0(npi0s+3,1);
  CLHEP::HepMatrix Lambda(npi0s+3,1);

  CLHEP::HepVector v_keep    		= v;
  CLHEP::HepMatrix VE_keep    		= VE;
  CLHEP::HepMatrix VD_keep    		= VD;
  CLHEP::HepMatrix Lambda0_keep	        = Lambda0;
  CLHEP::HepMatrix E_keep      	        = E;
  CLHEP::HepMatrix D_keep      	        = D;
  CLHEP::HepVector a_keep 		= a;

  CLHEP::HepVector v_keep_v    	        = v;
  CLHEP::HepMatrix VE_keep_v    	= VE;
  CLHEP::HepMatrix VD_keep_v    	= VD;
  CLHEP::HepMatrix Lambda0_keep_v	= Lambda0;
  CLHEP::HepMatrix E_keep_v      	= E;
  CLHEP::HepMatrix D_keep_v      	= D;

  //=============================================================
  // start minizing loop

  int it,it_v;
  for(it=0;it<MAX_ITERACTIONS;it++) {		// Global minizing loop

    chisq_keep_v = HUGE;				// reset chisq_keep_v
    for(it_v=0;it_v<MAX_ITERACTIONS;it_v++) {	// KL vertex (v) minizing loop

      v0 = v; 	// set initial vertex to be expanding point
      error_flag = 0;	// reset error flag
      
      double dpx1_dx1[maxPi0], dpx2_dx2[maxPi0];
      double dpx1_dy1[maxPi0], dpx2_dy2[maxPi0];
      double dpx1_dE1[maxPi0], dpx2_dE2[maxPi0];
      
      double dpy1_dx1[maxPi0], dpy2_dx2[maxPi0];
      double dpy1_dy1[maxPi0], dpy2_dy2[maxPi0];
      double dpy1_dE1[maxPi0], dpy2_dE2[maxPi0];
      
      double dpz1_dx1[maxPi0], dpz2_dx2[maxPi0];
      double dpz1_dy1[maxPi0], dpz2_dy2[maxPi0];
      double dpz1_dE1[maxPi0], dpz2_dE2[maxPi0];
      
      double dR1_dx1[maxPi0],  dR2_dx2[maxPi0];
      double dR1_dy1[maxPi0],  dR2_dy2[maxPi0];
      double dR1_dE1[maxPi0],  dR2_dE2[maxPi0];
      
      double dpx1_dvx[maxPi0], dpx2_dvx[maxPi0];
      double dpx1_dvy[maxPi0], dpx2_dvy[maxPi0];
      double dpx1_dvz[maxPi0], dpx2_dvz[maxPi0];
      
      double dpy1_dvx[maxPi0], dpy2_dvx[maxPi0];
      double dpy1_dvy[maxPi0], dpy2_dvy[maxPi0];
      double dpy1_dvz[maxPi0], dpy2_dvz[maxPi0];
      
      double dpz1_dvx[maxPi0], dpz2_dvx[maxPi0];
      double dpz1_dvy[maxPi0], dpz2_dvy[maxPi0];
      double dpz1_dvz[maxPi0], dpz2_dvz[maxPi0];
      
      double dR1_dvx[maxPi0],  dR2_dvx[maxPi0];
      double dR1_dvy[maxPi0],  dR2_dvy[maxPi0];
      double dR1_dvz[maxPi0],  dR2_dvz[maxPi0];
      
      double dklpx_dvx, dklpy_dvx, dklpz_dvx;
      double dklpx_dvy, dklpy_dvy, dklpz_dvy;
      double dklpx_dvz, dklpy_dvz, dklpz_dvz;

      klpx = klpy = klpz = klE = 0.;
      dklpx_dvx = dklpy_dvx = dklpz_dvx = 0.;
      dklpx_dvy = dklpy_dvy = dklpz_dvy = 0.;
      dklpx_dvz = dklpy_dvz = dklpz_dvz = 0.;

      for( idx=0; idx<npi0s; idx++) {

    	  R1[idx] = sqrt(pow(x1[idx]-vx,2) + pow(y1[idx]-vy,2) + pow(z1[idx]-vz,2));
    	  R2[idx] = sqrt(pow(x2[idx]-vx,2) + pow(y2[idx]-vy,2) + pow(z2[idx]-vz,2));

    	  px1[idx] = (x1[idx]-vx)/R1[idx] * E1[idx];
    	  py1[idx] = (y1[idx]-vy)/R1[idx] * E1[idx];
    	  pz1[idx] = (z1[idx]-vz)/R1[idx] * E1[idx];

    	  px2[idx] = (x2[idx]-vx)/R2[idx] * E2[idx];
    	  py2[idx] = (y2[idx]-vy)/R2[idx] * E2[idx];
    	  pz2[idx] = (z2[idx]-vz)/R2[idx] * E2[idx];

    	  klpx += px1[idx] + px2[idx];
    	  klpy += py1[idx] + py2[idx];
    	  klpz += pz1[idx] + pz2[idx];
    	  klE  += E1[idx]  + E2[idx];

    	  dR1_dx1[idx] = (x1[idx] - vx) / R1[idx];
    	  dR1_dy1[idx] = (y1[idx] - vy) / R1[idx];
    	  dR1_dE1[idx] = 0;

    	  dpx1_dx1[idx] = (E1[idx] - px1[idx] * dR1_dx1[idx])/R1[idx];
    	  dpx1_dy1[idx] = (        - px1[idx] * dR1_dy1[idx])/R1[idx];
    	  dpx1_dE1[idx] = (x1[idx] - vx) / R1[idx];

    	  dpy1_dx1[idx] = (        - py1[idx] * dR1_dx1[idx])/R1[idx];
    	  dpy1_dy1[idx] = (E1[idx] - py1[idx] * dR1_dy1[idx])/R1[idx];
    	  dpy1_dE1[idx] = (y1[idx] - vy) / R1[idx];

    	  dpz1_dx1[idx] = (        - pz1[idx] * dR1_dx1[idx])/R1[idx];
    	  dpz1_dy1[idx] = (	   - pz1[idx] * dR1_dy1[idx])/R1[idx];
    	  dpz1_dE1[idx] = (z1[idx] - vz) / R1[idx];

    	  dR2_dx2[idx] = (x2[idx] - vx) / R2[idx];
    	  dR2_dy2[idx] = (y2[idx] - vy) / R2[idx];
    	  dR2_dE2[idx] = 0;

    	  dpx2_dx2[idx] = (E2[idx] - px2[idx] * dR2_dx2[idx])/R2[idx];
    	  dpx2_dy2[idx] = (        - px2[idx] * dR2_dy2[idx])/R2[idx];
    	  dpx2_dE2[idx] = (x2[idx] - vx) / R2[idx];

    	  dpy2_dx2[idx] = (        - py2[idx] * dR2_dx2[idx])/R2[idx];
    	  dpy2_dy2[idx] = (E2[idx] - py2[idx] * dR2_dy2[idx])/R2[idx];
    	  dpy2_dE2[idx] = (y2[idx] - vy) / R2[idx];

    	  dpz2_dx2[idx] = (        - pz2[idx] * dR2_dx2[idx])/R2[idx];
    	  dpz2_dy2[idx] = (	   - pz2[idx] * dR2_dy2[idx])/R2[idx];
    	  dpz2_dE2[idx] = (z2[idx] - vz) / R2[idx];

    	  dR1_dvx[idx] = -(x1[idx] - vx) / R1[idx];
    	  dR1_dvy[idx] = -(y1[idx] - vy) / R1[idx];
    	  dR1_dvz[idx] = -(z1[idx] - vz) / R1[idx];

    	  dpx1_dvx[idx] = (-E1[idx] - px1[idx] * dR1_dvx[idx])/R1[idx];
    	  dpx1_dvy[idx] = (         - px1[idx] * dR1_dvy[idx])/R1[idx];
    	  dpx1_dvz[idx] = (         - px1[idx] * dR1_dvz[idx])/R1[idx];

    	  dpy1_dvx[idx] = (         - py1[idx] * dR1_dvx[idx])/R1[idx];
    	  dpy1_dvy[idx] = (-E1[idx] - py1[idx] * dR1_dvy[idx])/R1[idx];
    	  dpy1_dvz[idx] = (         - py1[idx] * dR1_dvz[idx])/R1[idx];
	  
    	  dpz1_dvx[idx] = (         - pz1[idx] * dR1_dvx[idx])/R1[idx];
    	  dpz1_dvy[idx] = (         - pz1[idx] * dR1_dvy[idx])/R1[idx];
    	  dpz1_dvz[idx] = (-E1[idx] - pz1[idx] * dR1_dvz[idx])/R1[idx];
	  
    	  dR2_dvx[idx] = -(x2[idx] - vx) / R2[idx];
    	  dR2_dvy[idx] = -(y2[idx] - vy) / R2[idx];
    	  dR2_dvz[idx] = -(z2[idx] - vz) / R2[idx];
	  
    	  dpx2_dvx[idx] = (-E2[idx] - px2[idx] * dR2_dvx[idx])/R2[idx];
    	  dpx2_dvy[idx] = (         - px2[idx] * dR2_dvy[idx])/R2[idx];
    	  dpx2_dvz[idx] = (         - px2[idx] * dR2_dvz[idx])/R2[idx];
	  
    	  dpy2_dvx[idx] = (         - py2[idx] * dR2_dvx[idx])/R2[idx];
    	  dpy2_dvy[idx] = (-E2[idx] - py2[idx] * dR2_dvy[idx])/R2[idx];
    	  dpy2_dvz[idx] = (         - py2[idx] * dR2_dvz[idx])/R2[idx];
	  
    	  dpz2_dvx[idx] = (         - pz2[idx] * dR2_dvx[idx])/R2[idx];
    	  dpz2_dvy[idx] = (         - pz2[idx] * dR2_dvy[idx])/R2[idx];
    	  dpz2_dvz[idx] = (-E2[idx] - pz2[idx] * dR2_dvz[idx])/R2[idx];
	  
    	  dklpx_dvx += (dpx1_dvx[idx] + dpx2_dvx[idx]);
    	  dklpx_dvy += (dpx1_dvy[idx] + dpx2_dvy[idx]);
    	  dklpx_dvz += (dpx1_dvz[idx] + dpx2_dvz[idx]);
	  
    	  dklpy_dvx += (dpy1_dvx[idx] + dpy2_dvx[idx]);
    	  dklpy_dvy += (dpy1_dvy[idx] + dpy2_dvy[idx]);
    	  dklpy_dvz += (dpy1_dvz[idx] + dpy2_dvz[idx]);
	  
    	  dklpz_dvx += (dpz1_dvx[idx] + dpz2_dvx[idx]);
    	  dklpz_dvy += (dpz1_dvy[idx] + dpz2_dvy[idx]);
    	  dklpz_dvz += (dpz1_dvz[idx] + dpz2_dvz[idx]);
      }
      
      for( idx=0; idx<npi0s; idx++){
	d(idx+1) =   pow(E1[idx]  + E2[idx] ,2)
	  - pow(px1[idx] + px2[idx],2)
	  - pow(py1[idx] + py2[idx],2)
	  - pow(pz1[idx] + pz2[idx],2)
	  - MASS_PI0*MASS_PI0;
      }
      d(npi0s+1) = klE*klE - klpx*klpx - klpy*klpy - klpz*klpz - MASS_KL*MASS_KL;
      
      d(npi0s+2) = - vx*klE;
      d(npi0s+3) = - vy*klE;
      for( idx=0; idx<npi0s; idx++) {
	d(npi0s+2) += E1[idx]*x1[idx] + E2[idx]*x2[idx];
	d(npi0s+3) += E1[idx]*y1[idx] + E2[idx]*y2[idx];
      }
      
      // Reset core matrix D & E
      for(int i=1;i<=npi0s+3;i++) {
	for(int j=1;j<=npi0s*6-1;j++) 	D(i,j) = 0.;
	for(int j=1;j<=4;j++)  		E(i,j) = 0.;
      }
      
      // Calculate core matrix D & E
      for( idx=0; idx<npi0s; idx++) {
	
	//d(idx+1) =   pow(E1[idx]  + E2[idx] ,2)
	//	     - pow(px1[idx] + px2[idx],2)
	//	     - pow(py1[idx] + py2[idx],2)
	//	     - pow(pz1[idx] + pz2[idx],2)
	//	     - MASS_PI0*MASS_PI0;
	
	D(idx+1,idx*6+1) = -2.*(px1[idx] + px2[idx])*dpx1_dx1[idx]
	  -2.*(py1[idx] + py2[idx])*dpy1_dx1[idx]
	  -2.*(pz1[idx] + pz2[idx])*dpz1_dx1[idx]; // del x1
	D(idx+1,idx*6+2) = -2.*(px1[idx] + px2[idx])*dpx1_dy1[idx]
	  -2.*(py1[idx] + py2[idx])*dpy1_dy1[idx]
	  -2.*(pz1[idx] + pz2[idx])*dpz1_dy1[idx]; // del y1
	D(idx+1,idx*6+3) =  2.*(E1[idx]  + E2[idx])
	  -2.*(px1[idx] + px2[idx])*dpx1_dE1[idx]
	  -2.*(py1[idx] + py2[idx])*dpy1_dE1[idx]
	  -2.*(pz1[idx] + pz2[idx])*dpz1_dE1[idx]; // del E1
	
	D(idx+1,idx*6+4) = -2.*(px1[idx] + px2[idx])*dpx2_dx2[idx]
	  -2.*(py1[idx] + py2[idx])*dpy2_dx2[idx]
	  -2.*(pz1[idx] + pz2[idx])*dpz2_dx2[idx]; // del x2
	D(idx+1,idx*6+5) = -2.*(px1[idx] + px2[idx])*dpx2_dy2[idx]
	  -2.*(py1[idx] + py2[idx])*dpy2_dy2[idx]
	  -2.*(pz1[idx] + pz2[idx])*dpz2_dy2[idx]; // del y2
	if (idx!=npi0s-1)
	  D(idx+1,idx*6+6) =  2.*(E1[idx]  + E2[idx])
	    -2.*(px1[idx] + px2[idx])*dpx2_dE2[idx]
	    -2.*(py1[idx] + py2[idx])*dpy2_dE2[idx]
	    -2.*(pz1[idx] + pz2[idx])*dpz2_dE2[idx]; // del E2
	
	E(idx+1,1)  = -2.*(px1[idx] + px2[idx])*(dpx1_dvx[idx] + dpx2_dvx[idx])
	  -2.*(py1[idx] + py2[idx])*(dpy1_dvx[idx] + dpy2_dvx[idx])
	  -2.*(pz1[idx] + pz2[idx])*(dpz1_dvx[idx] + dpz2_dvx[idx]); // del vx
	E(idx+1,2)  = -2.*(px1[idx] + px2[idx])*(dpx1_dvy[idx] + dpx2_dvy[idx])
	  -2.*(py1[idx] + py2[idx])*(dpy1_dvy[idx] + dpy2_dvy[idx])
	  -2.*(pz1[idx] + pz2[idx])*(dpz1_dvy[idx] + dpz2_dvy[idx]); // del vy
	E(idx+1,3)  = -2.*(px1[idx] + px2[idx])*(dpx1_dvz[idx] + dpx2_dvz[idx])
	  -2.*(py1[idx] + py2[idx])*(dpy1_dvz[idx] + dpy2_dvz[idx])
	  -2.*(pz1[idx] + pz2[idx])*(dpz1_dvz[idx] + dpz2_dvz[idx]); // del vz
	
	if (idx==npi0s-1)
	  E(idx+1,4)  =  2.*(E1[idx]  + E2[idx])
	    -2.*(px1[idx] + px2[idx])*dpx2_dE2[idx]
	    -2.*(py1[idx] + py2[idx])*dpy2_dE2[idx]
	    -2.*(pz1[idx] + pz2[idx])*dpz2_dE2[idx]; // del E2[npi0s-1]
	
	// d(npi0s+1) = klE*klE - klpx*klpx - klpy*klpy - klpz*klpz - MASS_KL*MASS_KL;
	
	D(npi0s+1,idx*6+1) = -2.* klpx * dpx1_dx1[idx]
	  -2.* klpy * dpy1_dx1[idx]
	  -2.* klpz * dpz1_dx1[idx]; // del x1
	D(npi0s+1,idx*6+2) = -2.* klpx * dpx1_dy1[idx]
	  -2.* klpy * dpy1_dy1[idx]
	  -2.* klpz * dpz1_dy1[idx]; // del y1
	D(npi0s+1,idx*6+3) =  2.* klE
	  -2.* klpx * dpx1_dE1[idx]
	  -2.* klpy * dpy1_dE1[idx]
	  -2.* klpz * dpz1_dE1[idx]; // del E1
	
	D(npi0s+1,idx*6+4) = -2.* klpx * dpx2_dx2[idx]
	  -2.* klpy * dpy2_dx2[idx]
	  -2.* klpz * dpz2_dx2[idx]; // del x2
	D(npi0s+1,idx*6+5) = -2.* klpx * dpx2_dy2[idx]
	  -2.* klpy * dpy2_dy2[idx]
	  -2.* klpz * dpz2_dy2[idx]; // del y2
	if (idx!=npi0s-1)
	  D(npi0s+1,idx*6+6) =  2.* klE
	    -2.* klpx * dpx2_dE2[idx]
	    -2.* klpy * dpy2_dE2[idx]
	    -2.* klpz * dpz2_dE2[idx]; // del E2
	
	E(npi0s+1,1)  = -2.* klpx * dklpx_dvx
	  -2.* klpy * dklpy_dvx
	  -2.* klpz * dklpz_dvx; // del vx
	E(npi0s+1,2)  = -2.* klpx * dklpx_dvy
	  -2.* klpy * dklpy_dvy
	  -2.* klpz * dklpz_dvy; // del vy
	E(npi0s+1,3)  = -2.* klpx * dklpx_dvz
	  -2.* klpy * dklpy_dvz
	  -2.* klpz * dklpz_dvz; // del vz
	if (idx==npi0s-1)
	  E(npi0s+1,4)  =  2.* klE
	    -2.* klpx * dpx2_dE2[idx]
	    -2.* klpy * dpy2_dE2[idx]
	    -2.* klpz * dpz2_dE2[idx]; // del E2
	
	//	d(npi0s+2) += sum(E1[idx]*x1[idx] + E2[idx]*x2[idx]) - vx*klE;
	//	d(npi0s+3) += sum(E1[idx]*y1[idx] + E2[idx]*y2[idx]) - vy*klE;
	
	D(npi0s+2,idx*6+1) = E1[idx]; 		// del x1
	D(npi0s+2,idx*6+2) = 0.; 		// del y1
	D(npi0s+2,idx*6+3) = x1[idx] - vx; 	// del E1
	
	D(npi0s+2,idx*6+4) = E2[idx]; 		// del x2
	D(npi0s+2,idx*6+5) = 0.; 		// del y2

	if (idx!=npi0s-1)
	  D(npi0s+2,idx*6+6) = x2[idx] - vx; 	// del E2
	
	E(npi0s+2,1)  = -klE; 			// del vx
	E(npi0s+2,2)  = 0.; 			// del vy
	E(npi0s+2,3)  = 0.; 			// del vz
	
	if (idx==npi0s-1)
	  E(npi0s+2,4)  = x2[idx] - vx; 	// del E2
	
	D(npi0s+3,idx*6+1) = 0.; 		// del x1
	D(npi0s+3,idx*6+2) = E1[idx]; 		// del y1
	D(npi0s+3,idx*6+3) = y1[idx] - vy; 	// del E1
	
	D(npi0s+3,idx*6+4) = 0.; 		// del x2
	D(npi0s+3,idx*6+5) = E2[idx]; 		// del y2
	
	if (idx!=npi0s-1)
	  D(npi0s+3,idx*6+6) = y2[idx] - vy; 	// del E2
	
	E(npi0s+3,1)  = 0.; 			// del vx
	E(npi0s+3,2)  = -klE; 			// del vy
	E(npi0s+3,3)  = 0.; 			// del vz
	
	if (idx==npi0s-1)
	  E(npi0s+3,4)  = y2[idx] - vy; 	// del E2
	
      }
      
      da0 = a0 - a;
      dv0 = v0 - v;
      
      VD = (D*Va0*D.T()).inverse(ifail);
      if (ifail!=0) {
	fprintf(stderr,"[ERROR] : Cannot get inverse matrix (VD) : error code = %d\n",ifail);
	error_flag = -1;
      }else {
	VE = (E.T()*VD*E).inverse(ifail);
	if (ifail!=0) {
	  fprintf(stderr,"[ERROR] : Cannot get inverse matrix (VE) : error code = %d\n",ifail);
	  error_flag = -2;
	}else {
	  
	  Lambda0 = VD*(D*da0 + d);
	  chisq = (Lambda0.T() * (D*da0 + E*dv0 + d))(1,1);
	  v = v - VE * E.T() * Lambda0;
	}
      }
      
      if (chisq < chisq_keep_v && error_flag==0 && chisq>=0.) {
	chisq_keep_v 	= chisq;
	v_keep_v    	= v;
	VE_keep_v    	= VE;
	VD_keep_v    	= VD;
	Lambda0_keep_v	= Lambda0;
	E_keep_v      	= E;
	D_keep_v      	= D;
	fit_updated	= 1;
      }else {
	chisq 	= chisq_keep_v;
	v    	= v_keep_v;
	VE    	= VE_keep_v;
	VD    	= VD_keep_v;
	Lambda0	= Lambda0_keep_v;
	E      	= E_keep_v;
	D      	= D_keep_v;
	break;
      }
      
      // fill back KL vertex first
      vx          = v(1);
      vy          = v(2);
      vz          = v(3);
      E2[npi0s-1] = v(4);
      
    } // end of KL vertex minizing loop
    
    if (error_flag==0) {
      Lambda  = Lambda0 - VD * E * VE * E.T() * Lambda0;
      a = a0 - Va0 * D.T() * Lambda;
    }
    
    if ((it == 0 || chisq < chisq_keep) && error_flag==0 && chisq>=0.) {
      chisq_keep 	= chisq;
      v_keep    	= v;
      VE_keep    	= VE;
      VD_keep    	= VD;
      Lambda0_keep	= Lambda0;
      E_keep      	= E;
      D_keep      	= D;
      a_keep 		= a;
      fit_updated	= 1;
    }else {
      chisq 	= chisq_keep;
      v    	= v_keep;
      VE    	= VE_keep;
      VD    	= VD_keep;
      Lambda0	= Lambda0_keep;
      E      	= E_keep;
      D      	= D_keep;
      a	        = a_keep;
      break;
    }
    
    // fill back other parameters
    for( idx=0; idx<npi0s; idx++) {
      x1[idx] = a(idx*6 + 1);
      y1[idx] = a(idx*6 + 2);
      E1[idx] = a(idx*6 + 3);
      x2[idx] = a(idx*6 + 4);
      y2[idx] = a(idx*6 + 5);
      if (idx!=npi0s-1)
	E2[idx] = a(idx*6 + 6);
    }
    
  } // end of global minizing loop
  
  if (fit_updated == 0) { // fit failed.
    fprintf(stderr,"[ERROR] : Fit failed (%d), return error code = %d.\n",ifail,error_flag);
    chisq = -1.;
    if (error_flag!=0) return (double)error_flag;
    else return -3.;
  }
  
  klpx = klpy = klpz = klE = 0.;
  for( idx=0; idx<npi0s; idx++) {
    
    R1[idx] = sqrt(pow(x1[idx]-vx,2) + pow(y1[idx]-vy,2) + pow(z1[idx]-vz,2));
    R2[idx] = sqrt(pow(x2[idx]-vx,2) + pow(y2[idx]-vy,2) + pow(z2[idx]-vz,2));
    
    px1[idx] = (x1[idx]-vx)/R1[idx] * E1[idx];
    py1[idx] = (y1[idx]-vy)/R1[idx] * E1[idx];
    pz1[idx] = (z1[idx]-vz)/R1[idx] * E1[idx];
    
    px2[idx] = (x2[idx]-vx)/R2[idx] * E2[idx];
    py2[idx] = (y2[idx]-vy)/R2[idx] * E2[idx];
    pz2[idx] = (z2[idx]-vz)/R2[idx] * E2[idx];
    
    klpx += px1[idx] + px2[idx];
    klpy += py1[idx] + py2[idx];
    klpz += pz1[idx] + pz2[idx];
    klE  += E1[idx]  + E2[idx];
  }
  
  // Final error matrix
  CLHEP::HepMatrix VDtu = VD - VD * E * VE * E.T() * VD;
  CLHEP::HepMatrix Va = Va0 - Va0 * D.T() * VDtu * D * Va0;
  
  KL.setMass(sqrt(klE*klE-klpx*klpx-klpy*klpy-klpz*klpz));
  KL.setEnergy(klE);
  KL.setP3(klpx,klpy,klpz);
  KL.setVtx(vx,vy,vz);
  
  for( idx=0; idx<npi0s; idx++) {
    KL.pi0()[idx].setEnergy(E1[idx]+E2[idx]);
    KL.pi0()[idx].setMass(sqrt(pow(E1[idx]+E2[idx],2) -
			       pow(px1[idx]+px2[idx],2) -
			       pow(py1[idx]+py2[idx],2) -
			       pow(pz1[idx]+pz2[idx],2)));

    KL.pi0()[idx].setVtx(vx,vy,vz);
    KL.pi0()[idx].setP3(px1[idx]+px2[idx],py1[idx]+py2[idx],pz1[idx]+pz2[idx]);
    KL.pi0()[idx].setRecZ(vz);
    KL.pi0()[idx].setRecZsig2(VE(3,3));
    
    KL.pi0()[idx].g1().setEnergy(E1[idx]);
    KL.pi0()[idx].g1().setPos(x1[idx],y1[idx],z1[idx]);
    KL.pi0()[idx].g1().setP3(px1[idx],py1[idx],pz1[idx]);
    
    KL.pi0()[idx].g2().setEnergy(E2[idx]);
    KL.pi0()[idx].g2().setPos(x2[idx],y2[idx],z2[idx]);
    KL.pi0()[idx].g2().setP3(px2[idx],py2[idx],pz2[idx]);
  }
  
  return chisq;
}


