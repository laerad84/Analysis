#include "rec2g/Rec2g.h"
#include "pi0/Pi0.h"
#include "gnana/E14GNAnaFunction.h"
#include <iostream>
#include <list>
#include "TMath.h"


bool user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList, double& mass){
  static Rec2g rec2g;      
  CLHEP::Hep3Vector Vtx(0,0,3148);
  std::list<Gamma>::iterator itGamma;  
  CLHEP::Hep3Vector GPos[2];
  GPos[0] = glist.front().pos() - Vtx;
  GPos[1] = glist.back().pos()  - Vtx;
  GPos[0].setMag(glist.front().e());
  GPos[1].setMag(glist.back().e());
  
  double MassSq = (glist.front().e()+ glist.back().e())*(glist.front().e()+ glist.back().e())-(GPos[0]+GPos[1]).mag2();    
  mass   = TMath::Sqrt(MassSq);

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
