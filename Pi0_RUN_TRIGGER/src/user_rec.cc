#include "rec2g/Rec2g.h"
#include "pi0/Pi0.h"
#include "gnana/E14GNAnaFunction.h"
#include <iostream>
#include <list>
#include "TMath.h"


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


void 
recVtxWithConstM( const Gamma& g1, const Gamma& g2, double Mass,
		  double* recZ, double* recZsig2 )
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

std::list<Pi0>  
recPi0withConstM( std::list<Gamma> glist, double mass )
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


bool user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList,double &mass , double recPosition){
  //static Rec2g rec2g;        
  // reconstruction 
  
  mass = rec_mass2g( glist,recPosition);  
  //std::cout<< "user_rec : " << mass << std::endl;
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
