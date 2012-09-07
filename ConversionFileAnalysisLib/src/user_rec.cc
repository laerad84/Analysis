#include <iostream>
#include <list>

#include "TMath.h"

#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "rec2g/Rec2g.h"
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include "csimap/CsiMap.h"

#include "user_func.h"

double rec_mass2g(std::list<Gamma> const &glist, double recPosition){
  double mass = 0;
  CLHEP::Hep3Vector Vtx(0,0,recPosition);
  int id = 0;
  if(glist.size() !=2 ){return 0;}
  Gamma g1 = glist.front();
  Gamma g2 = glist.back();
  // index = 0 : positive solution (status = 1)
  // index = 1 : negative solution (status = 2)
  CLHEP::Hep3Vector GPos[2];
  GPos[0] = (g1).pos() -Vtx;
  GPos[1] = (g2).pos() -Vtx;
  GPos[0].setMag((g1).e());
  GPos[1].setMag((g2).e());      
  double MassSq = ((g1).e()+ (g2).e())*((g1).e()+(g2).e()) - (GPos[0]+GPos[1]).mag2();
  mass   = TMath::Sqrt(MassSq);            
  double recZ     = -999.;
  double recZsig2 = -1.;
  return mass;
}

void   recVtxWithConstM( const Gamma& g1, const Gamma& g2, double Mass,
		       double* recZ, double* recZsig2 )
{
  std::cout << "rec2g::recVtxWithConstM() : Mass=" << Mass << std::endl;

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

std::list<Pi0>  recPi0withConstM( std::list<Gamma> glist, double mass )
{
  //
  std::list<Pi0> pi0list;
  int id = 0;
  for( std::list<Gamma>::iterator g1=glist.begin();
       g1!=glist.end(); g1++ ) {
    for( std::list<Gamma>::iterator g2=glist.begin();
	 g2!=glist.end(); g2++ ) {
      //
      if( g1->id() >= g2->id() )
	continue;  // skip this combination

      // index = 0 : positive solution (status = 1)
      // index = 1 : negative solution (status = 2)

      double recZ[2]     = { -999.,-999. };
      double recZsig2[2] = { -1.,-1. };
      recVtxWithConstM( *g1, *g2,  mass, recZ, recZsig2 );
      for( int i=0; i<2; i++ ) {
	if( recZsig2[i] >= 0 ) {
	  Pi0 pi0;
	  
	  pi0.setId( id );
	  pi0.setRecZ( recZ[i] );
	  pi0.setRecZsig2( recZsig2[i] );
	  pi0.setVtx( 0, 0, recZ[i] );
	  pi0.setGamma( *g1, *g2 );
	  pi0.setStatus( i+1 );

	  pi0.updateVars();
	  
	  pi0list.push_back( pi0 );
	  id++;
	}
      }
    }
  }
  return( pi0list );
}

bool   user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList,double &mass , double recPosition){
  //static Rec2g rec2g;        
  // reconstruction 
  
  mass = rec_mass2g( glist,recPosition);  
  std::cout<< "user_rec : " << mass << std::endl;
  piList = recPi0withConstM(glist,mass);
  if(piList.size()!=1)   return false;
  
  // position correction for angle dependency
  E14GNAnaFunction::getFunction()->correctPosition(piList.front().g1());  
  E14GNAnaFunction::getFunction()->correctPosition(piList.front().g2());  
  
  // re-reconstruction with corrected gamma
  std::list<Gamma> glist2;
  glist2.push_back(piList.front().g1());
  glist2.push_back(piList.front().g2());

  mass = rec_mass2g( glist2,recPosition);  
  piList = recPi0withConstM(glist2,mass);
  if(piList.size()!=1)   return false;

  // shape chi2 evaluation 
  E14GNAnaFunction::getFunction()->shapeChi2(piList.front().g1());  
  E14GNAnaFunction::getFunction()->shapeChi2(piList.front().g2());  

  return true;
}

bool   user_rec(std::list<Gamma> const &glist, std::vector<Klong> &klVec){
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

int    kinematicCut( std::vector<Klong> const &klVec){
  int iCut  = 0;

  double Egamma_Min = 0.10 * 1000; //GeV to MeV
  double Egamma_Max = 2.0 * 1000; //GeV to MeV
  
  double XY_Min =150;
  double R_Max =850;

  double gDist_Min =150.; //mm

  double Zvert_Min =3000;
  double Zvert_Max =5000;

  double Pt_Max = sqrt(1.25e-4) * 1000; //GeV to MeV
  
  double BestChi2_Max = 3.0;

  double Chi2_2ndChi2_Min = 4.0;

  double pi0Mass_diff_Max = 5.125; //MeV

  double pi0VertZ_diff_Max = 100; // mm

  int const nPi0 = klVec[0].pi0().size();
  double clusterPos[10][2]={{0}};
  for(int ipi=0;ipi<nPi0;ipi++){ 
    Pi0 const & pi0 = klVec[0].pi0().at(ipi);
    //pi0 mass cut
    double const pi0Mass = 134.976; //MeV
    if(fabs(pi0.m()-pi0Mass)>pi0Mass_diff_Max)
      iCut = (iCut | 1<<PI0_MASS_CUT );

    for(int igam = 0;igam<2;igam++){
      Gamma const &g = (igam==0) ? pi0.g1() : pi0.g2();      
      //Energy Cut
      if( g.e()  < Egamma_Min || g.e() > Egamma_Max )
	iCut = (iCut | 1<<ENERGY_CUT );

      //Fiducial in CsI      
      if( g.pos().perp() > R_Max ||  
	  (fabs(g.x())<XY_Min && fabs(g.y())<XY_Min) ||
	  fabs(g.y()) > 550){
	iCut = ( iCut | 1<<CSI_FIDUCIAL_CUT );
      }
      clusterPos[ipi*2+igam][0] = g.x();
      clusterPos[ipi*2+igam][1] = g.y();
    }
  }

  //Two gamma distance
  for(int i=0; i<nPi0*2; i++){
    for(int j=i+1; j<nPi0*2; j++){
      double dist = sqrt( pow(clusterPos[i][0]-clusterPos[i][1],2)
			  + pow(clusterPos[i][0]-clusterPos[i][1],2) );
      if ( dist < gDist_Min ) iCut = ( iCut | 1<<DISTANCE_CUT );      
    }
  }

  // KL Vertex Cut
  double RecVertZ = klVec[0].vz();
  if ( RecVertZ < Zvert_Min || RecVertZ > Zvert_Max )
    iCut = ( iCut | 1<<KL_VERTEX_CUT );
  
  // KL PT Cut
  double RecPt = klVec[0].p3().perp();
  if ( RecPt > Pt_Max )   iCut = ( iCut | 1<<KL_PT_CUT ); 

  // KL Best Chi2 cut
  double Chi2 = klVec[0].chisqZ();
  if ( Chi2 > BestChi2_Max ) iCut = ( iCut | 1<<KL_CHI2_CUT );
    
  // KL 2nd Best - Best Chi2
  if(klVec.size()>1){
    double Chi2_2nd = klVec[1].chisqZ();
    if( Chi2_2nd - Chi2 < Chi2_2ndChi2_Min ) 
      iCut = ( iCut | 1<<KL_2NDBEST_CHI2_CUT );
  }
  
  if(klVec[0].deltaZ()>pi0VertZ_diff_Max)
    iCut = ( iCut | 1<<PI0_Z_CUT );
  
  return iCut;  
}

int    shapeCut(std::vector<Klong> const &klVec){
  double cutVal = 2.5;
  int iCut = 0;
  for(std::vector<Pi0>::const_iterator it = klVec[0].pi0().begin();
      it!=klVec[0].pi0().end(); it++){
    if(it->g1().chisq() > cutVal || it->g2().chisq() > cutVal )
      iCut = ( iCut | 1<<SHAPE_CHI2_CUT );
  }
  
  return iCut;
}

int    vetoCut(E14GNAnaDataContainer const &data){
  double threshold[100];
  for(int i=0;i<100;i++) threshold[i] = -1;
  threshold[DigiReader::CC01]=3;
  threshold[DigiReader::CC02]=3;
  threshold[DigiReader::CC03]=3;
  threshold[DigiReader::CC04]=3;
  threshold[DigiReader::CC05]=3;
  threshold[DigiReader::CC06]=3;
  threshold[DigiReader::CBAR]=3;
  threshold[DigiReader::FBAR]=3;
  threshold[DigiReader::BHPV]=3;
  threshold[DigiReader::CV]=3;
  threshold[DigiReader::BCV]=3;
  threshold[DigiReader::BHCV]=3;
  
  int iVeto = 0;
  for(int i=0;i<data.VetoNumber;i++){
    int detid = data.VetoDetID[i];
    double thre = threshold[detid];
    if(thre<0){
      static bool isFirstOne = true;
      if(isFirstOne){
	GsimMessage::getInstance()->
	  report("warning",Form("there is the unknown detector (detID==%d)",detid));
	isFirstOne = false;
      }
    }else if(thre<data.VetoEne[i]){
      iVeto = ( iVeto | (1<<detid) );
    }
  }
  return iVeto;
}

int    csiVetoCut(E14GNAnaDataContainer const &data,std::vector<Klong> const &klVec){
  bool isVeto = false;
  int nClusterCsi=0;
  int clusterCsiId[3000];
  int const nPi0 = klVec[0].pi0().size();
  double cx[10]={0},cy[10]={0};

  for(int ipi=0;ipi<nPi0;ipi++){ 
    for(int iclus=0; iclus<2 ; iclus++){
      Gamma const &clus = (iclus==0)
	? klVec[0].pi0().at(ipi).g1() : klVec[0].pi0().at(ipi).g2();
      cx[ipi*2+iclus] = clus.x();
      cy[ipi*2+iclus] = clus.y();
      for(int i=0; i<clus.clusterIdVec().size(); i++){
	clusterCsiId[nClusterCsi++] = clus.clusterIdVec().at(i);
      }
    }
  }
  

  for(int icsi = 0; icsi < data.CsiNumber; icsi++){
    int id = data.CsiModID[icsi];
    double edep  = data.CsiEne[icsi];
    
    bool find = false;
    for(int i=0; i<nClusterCsi; i++){
      if(clusterCsiId[i]==id){
	find = true;
	break;
      }
    }

    if(find) continue;
    
    double x,y,w;
    CsiMap::getCsiMap()->getXYW(id,x,y,w);

    double dist = sqrt(pow(x-cx[0],2)+pow(y-cy[0],2));
    for(int i = 1; i<2*nPi0; i++){
      double tmpdist = sqrt(pow(x-cx[i],2)+pow(y-cy[i],2));
      if(dist>tmpdist) dist = tmpdist;
    }
    
    if(edep > getCsiThreshold(dist)){
      isVeto = true;
      break;
    }
  }
  if(isVeto) return 1<<DigiReader::CSI;
  return 0;
}

double getCsiThreshold(double distance){
  double constant = 3.73275;
  double slope = 0.0074758;
  double threshold = exp(constant-slope*distance);

  if(threshold>10) return 10;
  if(threshold<1.5) return 1.5;
  return threshold;
}

void   user_cut(E14GNAnaDataContainer &data,std::vector<Klong> const &klVec){
  
  int &iCut = data.CutCondition;
  iCut = 0;
  iCut += shapeCut(klVec);
  iCut += kinematicCut(klVec);
  
  int &vCut = data.VetoCondition;
  vCut = 0;
  vCut += vetoCut(data);
  vCut += csiVetoCut(data,klVec);
}

