#include "klong/RecKlong.h"
#include "gnana/E14GNAnaFunction.h"
#include <iostream>
#include <list>


bool user_rec(std::list<Gamma> const &glist, std::vector<Klong> &klVec){
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
