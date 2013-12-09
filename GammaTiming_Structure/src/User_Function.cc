#include "User_Function.h"
#include "klong/Klong.h"
#include "klong/RecKlong.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "csimap/CsiMap.h"
#include <iostream>
#include <list>

#include "rec2g/Rec2g.h"
#include "pi0/Pi0.h"
#include "TMath.h"

double         rec_mass2g(std::list<Gamma> const &glist, double recPosition){
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
void           recVtxWithConstM( const Gamma& g1, const Gamma& g2, double Mass,double* recZ, double* recZsig2 )
{
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
std::list<Pi0> recPi0withConstM( std::list<Gamma> glist, double mass )
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
bool           user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList,double &mass , double recPosition){
  //static Rec2g rec2g;        
  // reconstruction 
  
  mass = rec_mass2g( glist,recPosition);  
  //std::cout<< "user_rec : " << mass << std::endl;
  piList = recPi0withConstM(glist,mass);
  if(piList.size()!=1)   return false;
  
  // position correction for angle dependency
  E14GNAnaFunction::getFunction()->correctPosition(piList.front().g1());  
  E14GNAnaFunction::getFunction()->correctPosition(piList.front().g2());  
  correctTime( piList.front().g1() );
  correctTime( piList.front().g2() );

 
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
bool           user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList){
  static Rec2g rec2g;      

  // reconstruction 
  piList = rec2g.recPi0withConstM(glist);
  if(piList.size()!=1)   return false;
  
  // position correction for angle dependency
  E14GNAnaFunction::getFunction()->correctPosition(piList.front());  
  // enrgy correction for angle&energy dependency
  E14GNAnaFunction::getFunction()->correctEnergyWithAngle(piList.front());  
  
  // re-reconstruction with corrected gamma
  std::list<Gamma> glist2;
  glist2.push_back(piList.front().g1());
  glist2.push_back(piList.front().g2());
  piList = rec2g.recPi0withConstM(glist2);
  if(piList.size()!=1)   return false;

  // shape chi2 evaluation 
  E14GNAnaFunction::getFunction()->shapeChi2(piList.front());  

  return true;
}
void           correctTime( Gamma &gam){
  double const PositionCorrection[2] = {6.49003, 0.99254};
  double const CsIRadiationLength    = 18.5;//mm
  if( gam.p3().x()==0 && gam.p3().y()==0 && gam.p3().z()==0){
    return;
  }
  double GammaTime = gam.t();
  double GammaE    = gam.e();
  double L         = CsIRadiationLength*(PositionCorrection[0]+PositionCorrection[1]*log(GammaE/1000.));//mm
  double GammaT0   = GammaTime - L / 300.; // ns
  gam.setTime( GammaT0);
}
//int shapeCut(Pi0 const &pi);
void           user_cut(E14GNAnaDataContainer &data,std::list<Pi0> const &piList){
  Pi0 const &pi = piList.front();

  int &iCut = data.CutCondition;
  iCut = 0;
  //  iCut += shapeCut(pi);
  iCut += standardCut(pi);

  int &vCut = data.VetoCondition;
  vCut = 0;
  vCut += vetoCut(data);
  vCut += csiVetoCut(data,pi);
}
int            standardCut(Pi0 const &pi0){
  
  Gamma const &g1 = pi0.g1();
  Gamma const &g2 = pi0.g2();
  
  int iCut  = 0;

  
  double Egamma_Min = 0.1 * 1000; //GeV to MeV
  double Egamma_Max = 2.0 * 1000; //GeV to MeV
  //  double Egamma_Max = 10.0;
  
  //  double R_Min =175;
  double XY_Min =150;
  double R_Max =850;
  //  double R_Max =1350;  // for step2

  double Zvert_Min =3000;
  double Zvert_Max =5000;

  double Pt_Min =0.13 * 1000; //GeV to MeV
  double Pt_Max =0.25 * 1000; //GeV to MeV

  double Acop_Min =30.; // deg

  double gDist_Min =300.; //mm

  double TotE_Min =0.5 * 1000; // GeV to MeV

  double Etheta_Min =2.5 * 1000  ; //GeV to MeV 
  
  double Eratio_Min = 0.2;

  //Energy Cut
  if( g1.e()  < Egamma_Min || g2.e() < Egamma_Min ||
      g1.e() > Egamma_Max || g2.e() > Egamma_Max)
    iCut += 1<<ENERGY_CUT ;
  
  //Fiducial in CsI
  double R1 = g1.pos().perp();
  double R2 = g2.pos().perp();
  //  if ( R1 < R_Min || R1 > R_Max || R2 < R_Min || R2 > R_Max)  
  if( R1>R_Max ||  (fabs(g1.x())<XY_Min && fabs(g1.y())<XY_Min) ||
      R2>R_Max ||  (fabs(g2.x())<XY_Min && fabs(g2.y())<XY_Min) )
    iCut += 1<<CSI_FIDUCIAL_CUT;
  
  // Vertex Cut
  double RecVertZ = pi0.recZ();
  if ( RecVertZ < Zvert_Min || RecVertZ > Zvert_Max)
    iCut += 1<<VERTEX_CUT ;

  // PT Cut
  double RecPi0Pt = pi0.p3().perp();
  if ( RecPi0Pt < Pt_Min || RecPi0Pt > Pt_Max) 
    iCut += 1<<PT_CUT;

  // Acp_angle  Cut
  double acop;
  double  ct = (g1.x()*g2.x()+g1.y()*g2.y())/R1/R2;
  if(TMath::Abs(ct) >= 1){
    acop = 0.0;
  }
  else{
    acop = TMath::Pi() - TMath::ACos(ct);
    acop = acop * 180. / TMath::Pi();
  }
  
  double Acop_angle = acop;
  if (Acop_angle < Acop_Min) iCut += 1<<COLLINEAR_CUT;
  
  //Two gamma distance
  CLHEP::Hep3Vector gdist = g1.pos()-g2.pos();
  double gamma_dist = gdist.mag();
                              
  if (gamma_dist < gDist_Min) iCut += 1<<DISTANCE_CUT;
    
  //Total energy of two gammas
  double E_tot = g1.e() + g2.e();
  if (E_tot < TotE_Min)  iCut += 1<<E_TOTAL_CUT;
  
  // In_ang cut : Odd-pair cut
  double  In_Ang1 = TMath::ATan(R1/(g1.z()-pi0.recZ()));
  In_Ang1 = In_Ang1 * 180 / TMath::Pi();
  double  Etheta1 = g1.e() * In_Ang1;
  double  In_Ang2 = TMath::ATan(R2/(g2.z()-pi0.recZ()));
  In_Ang2 = In_Ang2 * 180 / TMath::Pi();
  double  Etheta2 = g2.e() * In_Ang2;


  if (Etheta1 < Etheta_Min || Etheta2 < Etheta_Min)
    iCut += 1<<E_THETA_CUT;
  
  // E_ratio_cut
  double Eratio = TMath::Min(g1.e(),g2.e())
    / TMath::Max(g1.e(),g2.e());
  if (Eratio < Eratio_Min)    iCut += 1<<E_RATIO_CUT;

  // Pi0 Kinematic Cut
  if(pi0kine_cut(pi0))  {
    iCut += 1<<PI_KINE_CUT;
  }

  return iCut;  
}
bool           cutline(double x1, double y1, double x2, double y2,double var1, double var2)
{
  double slope = ( y1 - y2 ) / ( x1 - x2 );
  double intercept = y1 - slope * x1 ;
  double y = slope * var1 + intercept;
  if ( var2 > y ) return kTRUE; // above the line
  return kFALSE;
}
bool   pi0kine_cut(Pi0 const &pi0)
{
  double cpi0pt=pi0.p3().perp()*1e-3;//RecPi0Pt*1e-3;
  double cpi0pz=pi0.p3().z()*1e-3;//RecPi0Pz*1e-3;
  double cpi0recz=pi0.recZ()*1e-1;//RecVertZ*1e-1;
  double cpi0e=pi0.e()*1e-3;//RecPi0E*1e-3;

  bool flag = kFALSE;
  double pratio = cpi0pt/cpi0pz;
  if (( cpi0recz < 400 )&&( pratio < 0.1 )) return kTRUE;
  if (( cpi0recz > 400 )&&
      !( cutline(400., 0.1, 500., 0.15, cpi0recz, pratio) )) return kTRUE;
  if ( !(cutline(300., 0.8, 500., 0.4, cpi0recz, cpi0e)) ) return kTRUE;
  if ( (cutline(300., 0.2, 500., 0.34, cpi0recz, pratio)) ) return kTRUE;
  return kFALSE;

}
/*
int shapeCut(Pi0 const &pi){
  double cutVal = 2.5;
  if(pi.g1().chisq() > cutVal || pi.g2().chisq() > cutVal )
    return 1<<SHAPE_CHI2_CUT;
  
  return 0;
}
*/
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
  threshold[DigiReader::NCC]=3;
  threshold[DigiReader::OEV]=3;
  
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
double getCsiThreshold(double distance){
  double constant = 3.73275;
  double slope = 0.0074758;
  double threshold = exp(constant-slope*distance);

  if(threshold>10) return 10;
  if(threshold<1.5) return 1.5;
  return threshold;
}
int    csiVetoCut(E14GNAnaDataContainer const &data,Pi0 const &pi){
  bool isVeto = false;

  for(int icsi = 0; icsi < data.CsiNumber; icsi++){
    int id = data.CsiModID[icsi];
    double edep  = data.CsiEne[icsi];
    
    bool find = false;
    double cx[2]={0},cy[2]={0};
    for(int igam=0; igam<2 ; igam++){
      Gamma const &gam = (igam==0)? pi.g1() : pi.g2();
      cx[igam] = gam.coex();
      cy[igam] = gam.coey();
      for(int i=0; i<gam.clusterIdVec().size()&&!find; i++){
	if(gam.clusterIdVec().at(i)==id)  find = true; 
      }
    }
    
    if(find) continue;    
    double x,y,w;
    CsiMap::getCsiMap()->getXYW(id,x,y,w);
    double dist = std::min(sqrt(pow(x-cx[0],2)+pow(y-cy[0],2)),
			   sqrt(pow(x-cx[1],2)+pow(y-cy[1],2)));

    if(edep > getCsiThreshold(dist)){
      isVeto = true;
      //      break;
    }
  }
  if(isVeto) return 1<<DigiReader::CSI;
  return 0;
}



  
