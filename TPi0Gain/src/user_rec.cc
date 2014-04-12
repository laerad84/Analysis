#include "rec2g/Rec2g.h"
#include "pi0/Pi0.h"
#include "gnana/E14GNAnaFunction.h"
#include <iostream>
#include <list>


bool user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList){
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
