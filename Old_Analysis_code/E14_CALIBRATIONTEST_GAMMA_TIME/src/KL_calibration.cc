#include "E14_CALIBRATIONTEST_GAMMA/KL_calibration.h"
////
void
user_ana_recg6_init(TH1D** his1,TH1D** his2)
{
  static bool init=false;
  //
  if(!init){
    his1[0] = new TH1D("his_KL_mass_Before" ,"kl_mass_Before" ,800,400,600);
    his1[1] = new TH1D("his_KL_mass_After"  ,"kl_mass_After"  ,800,400,600);
    his1[2] = new TH1D("his_pi0_mass_Before","pi0_mass_Before",800,100,200);
    his1[3] = new TH1D("his_pi0_mass_After" ,"pi0_mass_After" ,800,100,200);
    his1[4] = new TH1D("his_chisq"          ,"chisq"          ,100,-1,10);
    his1[5] = new TH1D("his_chisq_dof"      ,"chisq_dof"      ,100,-1,10);
    his1[6] = new TH1D("his_prob"           ,"prob"           ,100,0,1);
    for(int ich=0;ich<N_CSI;ich++){
      his2[ich]=new TH1D(Form("his_Calibration_%04d",ich),Form("his_Calibration_%04d",ich),500,0,2);
    }
    init=true;
  }
}

/*
////
int
user_ana_recg6( const std::vector<Gamma>& glist , int do_corr )
{
  //
  static RecKlong recklong;
  std::vector<Klong> kl = recklong.recK3pi0( glist );

  if( kl.size() > 0 ) { // KL reconstruction success
    recdata.klstat[0]   = kl[0].status();
    recdata.klmass[0]   = kl[0].m();
    recdata.kle[0]      = kl[0].e();
    recdata.klpp[0]     = kl[0].p3().mag();
    recdata.klpt[0]     = kl[0].p3().perp();
    recdata.klpx[0]     = kl[0].p3().x();
    recdata.klpy[0]     = kl[0].p3().y();
    recdata.klpz[0]     = kl[0].p3().z();
    recdata.klx[0]      = kl[0].vx();
    recdata.kly[0]      = kl[0].vy();
    recdata.klz[0]      = kl[0].vz();
    recdata.klid[0]     = kl[0].id();
    recdata.nkl         = 1;

    recdata2.klchi2z1 = kl[0].chisqZ();
    recdata2.klchi2z2 = (kl.size()>1) ? kl[1].chisqZ() : -1.0;
    recdata2.klmass2  = (kl.size()>1) ? kl[1].m()      :  0.0;
    recdata2.klnrec   = kl.size();

    recdata.klcast[0] = 0;

    his_Fit[0]->Fill(kl[0].m());
    for( int ipi0 = 0; ipi0 < kl[0].pi0.size() ; ipi0++){
      his_Fit[2]->Fill(kl[0].pi0()[ipi0].m());
    }

    if (do_corr) {

      if (recdata2.klnrec == 1 && recdata2.klchi2z1<5.){
	recdata.klcast[0] = CalEnergy_idv(kl[0]);
      }else if (recdata2.klnrec > 1  &&
		recdata2.klchi2z1<5. &&
		( recdata2.klchi2z2-recdata2.klchi2z1>8. || recdata2.klmass2<0.47 || recdata2.klmass2>0.54 )  ){
	recdata.klcast[0] = CalEnergy_idv(kl[0]);
      }

    }

    // Pi0 associated with KL
    for( std::vector<Pi0>::const_iterator p=kl[0].pi0().begin();
	 p!=kl[0].pi0().end(); p++ ) {

      recdata.pi0stat[ recdata.npi0 ]  = p->status();
      recdata.pi0mass[ recdata.npi0 ]  = p->m();
      recdata.pi0e[ recdata.npi0 ]     = p->e();
      recdata.pi0x[ recdata.npi0 ]     = p->vx();
      recdata.pi0y[ recdata.npi0 ]     = p->vy();
      recdata.pi0z[ recdata.npi0 ]     = p->vz();
      recdata.pi0recz[ recdata.npi0 ]  = p->recZ();
      recdata.pi0sig2[ recdata.npi0 ]  = p->recZsig2();
      recdata.pi0pp[ recdata.npi0 ]    = p->p3().mag();
      recdata.pi0pt[ recdata.npi0 ]    = p->p3().perp();
      recdata.pi0px[ recdata.npi0 ]    = p->p3().x();
      recdata.pi0py[ recdata.npi0 ]    = p->p3().y();
      recdata.pi0pz[ recdata.npi0 ]    = p->p3().z();
      recdata.pi0g1id[ recdata.npi0 ]  = p->g1().id();
      recdata.pi0g2id[ recdata.npi0 ]  = p->g2().id();
      recdata.pi0id[ recdata.npi0 ]    = p->id();
      recdata.npi0++;
    }
  }

  return( 0 );
}
*/

int CalEnergy_idv(Klong &KL_prefit,TH1D** his_Fit, TH1D** his_CSI){

  CLHEP::Hep3Vector pos[6];
  int n=0;
  double klmassdiff = 10.;//5 initial
  double pimassdiff = 6.;//3 initial
  double threshold  = 0.2;

  if( fabs(KL_prefit.m() - MASS_KL) > klmassdiff ) return -1;
  if( KL_prefit.pi0().size() != 3 ) return -2; 
  if( KL_prefit.pi0()[0].status() != 1 ||
      KL_prefit.pi0()[1].status() != 1 ||
      KL_prefit.pi0()[2].status() != 1 ) return -3;
  for( unsigned int idx = 0 ; idx < KL_prefit.pi0().size(); idx++){
    if( fabs( KL_prefit.pi0()[idx].m() - MASS_PI0 ) > pimassdiff ){return -4;}
    pos[n] = KL_prefit.pi0()[idx].g1().pos();
    n++;
    pos[n] = KL_prefit.pi0()[idx].g2().pos();
    n++;
  }
  
  for( int i = 0; i< n ; i++){
    for( int j = i+1; j< n; j++){
      if(( pos[i] - pos[j] ).perp() < 140.){ return -6;}
    }
  }
  
  //////////////////////////////////////////////
  // loop for all pi0s
  Klong kl;
  for( unsigned int idx = 0; idx < KL_prefit.pi0().size() ; idx++ ){
    double corr, err, chisq;
    int flag;
    std::vector<int>::const_iterator i_id;
    std::vector<double>::const_iterator i_e;        
    std::vector<int>::const_iterator jj_id; //cluster ID of g1
    std::vector<double>::const_iterator jj_e;// cluster deposit energy of g1
    
    
    ///////**********************************
    // fitting 1st gamma
    // push_back pi0 to kldata,
    // and exchange gamma.

    kl=KL_prefit;
    kl.pi0().clear();
    
    //push pi0s with pi0[idx] - in last
    for( unsigned int i=0;i<KL_prefit.pi0().size();i++) {
      if (i!=idx) 
	kl.pi0().push_back( KL_prefit.pi0()[i] );
    }
    
    kl.pi0().push_back( KL_prefit.pi0()[idx]);
    
    
    //for last pi0- exchange gammas
    kl.pi0()[kl.pi0().size()-1].setGamma(KL_prefit.pi0()[idx].g2(),KL_prefit.pi0()[idx].g1());
    

    //fitting
    his_Fit[0]->Fill(kl.m());
    his_Fit[2]->Fill(kl.pi0()[0].m());
    his_Fit[2]->Fill(kl.pi0()[1].m());
    his_Fit[2]->Fill(kl.pi0()[2].m());
    
    chisq = Refit_KL_freelast(kl);
    
    ////////////    ////////////    ////////////    ////////////    ////////////    ////////////
    //changed Parameter e:: KL_perfit.pi0[idx].g1().e(); -> kl.pi0()[kl.pi0.size() -1].g2();
    ////////////    ////////////    ////////////    ////////////    ////////////    ////////////

    his_Fit[1]->Fill(kl.m());
    his_Fit[3]->Fill(kl.pi0()[0].m());
    his_Fit[3]->Fill(kl.pi0()[1].m());
    his_Fit[3]->Fill(kl.pi0()[2].m());    
    his_Fit[4]->Fill(chisq);
    his_Fit[5]->Fill(chisq/2.);

    double clustnocorr = 0; // energy of Gamma
    for( std::vector<double>::const_iterator jj = KL_prefit.pi0()[idx].g1().clusterEVec().begin();
	 jj != KL_prefit.pi0()[idx].g1().clusterEVec().end();jj++){
      clustnocorr += *jj;
    }

    // Edited by Lee  Change Some Code.
    double corrfactor = KL_prefit.pi0()[idx].g1().e() - clustnocorr;
    //corr = (kl.pi0()[kl.pi0().size() - 1].g2().e() -corrfactor ) / (KL_prefit.pi0()[idx].g1().e() - corrfactor );
    //err  = kl.pi0()[kl.pi0().size() - 1].g2().sigmaE() / KL_prefit.pi0()[idx].g1().e();
    corrfactor = KL_prefit.pi0()[idx].g1().e() - clustnocorr;    
    corr = (kl.pi0()[kl.pi0().size() - 1].g2().e() - KL_prefit.pi0()[idx].g1().e())/clustnocorr + 1;
    err  = kl.pi0()[kl.pi0().size() - 1].g2().sigmaE() / KL_prefit.pi0()[idx].g1().e();
    
    flag = 0;
    if( chisq/2. > 5. ){ flag  =1;}
    if( fabs(kl.m() - MASS_KL ) > klmassdiff ){flag = 1;}
    if( fabs(kl.pi0()[kl.pi0().size() - 1 ].m() - MASS_PI0 ) > pimassdiff ){
      flag = 1;
    }
    if( fabs(kl.vz() -KL_prefit.vz()) > 200. ){
      flag = 1; 
    }

    for( i_id = KL_prefit.pi0()[idx].g1().clusterIdVec().begin(),
	   i_e = KL_prefit.pi0()[idx].g1().clusterEVec().begin();
	 i_id != KL_prefit.pi0()[idx].g1().clusterIdVec().end(); i_id++,i_e++ ) {
      //  Only maximum deposited crystal in cluster  &&  more than 0.2 of cluster energy //
      if( *i_id != *KL_prefit.pi0()[idx].g1().clusterIdVec().begin() && *i_e >=KL_prefit.pi0()[idx].g1().e()*threshold ){
	flag = 1;
      }
    }
    
    // flag = 1 -> failed gamma
    //if( !flag ){ // success gamma
    if (!flag ){

      double size = KL_prefit.pi0()[idx].g1().clusterIdVec().size();
      for( i_id = KL_prefit.pi0()[idx].g1().clusterIdVec().begin(),
	     i_e = KL_prefit.pi0()[idx].g1().clusterEVec().begin();
	   i_id != KL_prefit.pi0()[idx].g1().clusterIdVec().end(); i_id++,i_e++){
	double clustE = 0;// calcluation for weight 
	
	for( jj_id = KL_prefit.pi0()[idx].g1().clusterIdVec().begin(),
	       jj_e = KL_prefit.pi0()[idx].g1().clusterEVec().begin();
	     jj_id != KL_prefit.pi0()[idx].g1().clusterIdVec().end(); 
	     jj_id++,jj_e++ ){
	  if((!( *i_id != *KL_prefit.pi0()[idx].g1().clusterIdVec().begin() && *i_e >= KL_prefit.pi0()[idx].g1().e()*threshold))){
	    clustE += *jj_e;
	  }
	}
	
	int ff = 0; 
	for( int p = 0; p< N_EDGE_CSI; p++){
	  if( *i_id == edge_csi[p] )ff = 1;
	  break;
	}
	
	double weight = csi_energy[*i_id]/clustE;
	csi_cal_factor_sum[*i_id]   += corr*weight;
	csi_cal_factor_count[*i_id] += 1; 
	csi_cal_factor_weight[*i_id]+= weight;
	//if(  *i_id == *KL_prefit.pi0()[idx].g1().clusterIdVec().begin() && *i_e >= KL_prefit.pi0()[idx].g1().e()*0.4){
	if( *i_e >= KL_prefit.pi0()[idx].g1().e()*0.2){
	  his_CSI[*i_id]->Fill(corr);
	}
      } 
    }
    
    ///////////////////////////////////////
    // Fitting 2nd gamma
    kl=KL_prefit;
    kl.pi0().clear();
    
    //set pi0 - pi0[idx] last
    for( unsigned int index=0;index < KL_prefit.pi0().size();index++) {
      if (index != idx ){
	kl.pi0().push_back( KL_prefit.pi0()[index] );
      }
    }
    kl.pi0().push_back( KL_prefit.pi0()[idx]);
    
    // fitting
    his_Fit[0]->Fill(kl.m());
    his_Fit[2]->Fill(kl.pi0()[0].m());
    his_Fit[2]->Fill(kl.pi0()[1].m());
    his_Fit[2]->Fill(kl.pi0()[2].m());
    
    chisq = Refit_KL_freelast(kl);
    
    his_Fit[1]->Fill(kl.m());
    his_Fit[3]->Fill(kl.pi0()[0].m());
    his_Fit[3]->Fill(kl.pi0()[1].m());
    his_Fit[3]->Fill(kl.pi0()[2].m());
    
    his_Fit[4]->Fill(chisq);
    his_Fit[5]->Fill(chisq/2.);
    
    clustnocorr = 0;
    
    for( std::vector<double>::const_iterator jj = KL_prefit.pi0()[idx].g2().clusterEVec().begin();
	 jj!=KL_prefit.pi0()[idx].g2().clusterEVec().end(); jj++ )
      clustnocorr += *jj;
    
    // corr : 
    corrfactor = KL_prefit.pi0()[idx].g2().e() - clustnocorr;    
    corr = (kl.pi0()[kl.pi0().size() - 1].g2().e() - KL_prefit.pi0()[idx].g2().e())/clustnocorr + 1;
    err  = kl.pi0()[kl.pi0().size() - 1].g2().sigmaE() / KL_prefit.pi0()[idx].g2().e();
	    
    //selection
    flag = 0;
    if (chisq/2. > 5.) flag = 1;
    if (fabs(kl.m()-MASS_KL)> klmassdiff) flag = 1;
    if (fabs(kl.pi0()[kl.pi0().size()-1].m()-MASS_PI0)> pimassdiff ) flag = 1;
    if (fabs(kl.vz()-KL_prefit.vz())> 200.) flag = 1;
    

    //clusterIdVec <- arranged by energy?
    
    for( i_id = KL_prefit.pi0()[idx].g2().clusterIdVec().begin(),i_e = KL_prefit.pi0()[idx].g2().clusterEVec().begin();
	 i_id != KL_prefit.pi0()[idx].g2().clusterIdVec().end();
	 i_id++,i_e++ ) {

      /// only maximum deposited crystal && deposited energy > clusterEnergy*0.2 ///
      if( *i_id != *KL_prefit.pi0()[idx].g2().clusterIdVec().begin() && *i_e >=KL_prefit.pi0()[idx].g2().e()*threshold){
	flag = 1;
      }
    }
 
    if (!flag){
      float size=KL_prefit.pi0()[idx].g2().clusterIdVec().size();
      for( i_id = KL_prefit.pi0()[idx].g2().clusterIdVec().begin(),i_e = KL_prefit.pi0()[idx].g2().clusterEVec().begin();
	   i_id!=KL_prefit.pi0()[idx].g2().clusterIdVec().end();
	   i_id++,i_e++ ) {
	
	double clustE = 0;
	for( jj_id=KL_prefit.pi0()[idx].g2().clusterIdVec().begin(),jj_e=KL_prefit.pi0()[idx].g2().clusterEVec().begin();
	     jj_id!=KL_prefit.pi0()[idx].g2().clusterIdVec().end(); 
	     jj_id++,jj_e++ ){
	  
          if (!(*i_id != *KL_prefit.pi0()[idx].g2().clusterIdVec().begin() && *i_e>=KL_prefit.pi0()[idx].g2().e()*threshold)){
	    clustE += *jj_e;
	  }
	}
	
	int ff=0;
	for(int p=0;p<N_EDGE_CSI;p++)
	  if(*i_id==edge_csi[p]) ff=1;
	//	if(i==KL_prefit.pi0()[idx].g2().clusterIdVec().begin() || ff==1)
	{
	  double weight = csi_energy[*i_id]/clustE;
	  csi_cal_factor_sum[*i_id]    += corr*weight;
	  csi_cal_factor_count[*i_id]  += 1;
	  csi_cal_factor_weight[*i_id] += weight;
	  //if(*i_id == *KL_prefit.pi0()[idx].g2().clusterIdVec().begin() && *i_e >= KL_prefit.pi0()[idx].g2().e()*0.4){
	  if( *i_e >= KL_prefit.pi0()[idx].g2().e()*0.2){
	    his_CSI[*i_id]->Fill(corr);
	  }	  
	}
      }
    }
  }
  return 1;
}

double Refit_KL_freelast(Klong &KL){
  
  double chisq_keep   = HUGE;
  double chisq_keep_v = HUGE;
  double chisq        = -1.0;
  
  int ifail         = 0;
  int error_flag    = 0;
  int fit_updated   = 0;
  const int maxPi0 = 3;
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


