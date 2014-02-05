#ifndef _KlongParamter_H_
#define _KlongParameter_H_
#include "klong/Klong.h"
#include <vector>
#include <list>
#include <iostream>
#include "TMath.h"
#include "TTree.h"

class KlongParameter{
 public:
  KlongParameter();
  ~KlongParameter();
  
  void SetBranchAddress(TTree* trin );
  void Branch( TTree* trout );
  
  void Judge(Klong kl);
  void InitPar();

  double GMaxR;
  double GMinX;
  double GMinY;
  double GMinE;
  double GMaxTDelta;

  double klChisqZ;
  double klpt;
  double klpos[3];
  double klE;
  double pi0ptMax;
  double pi0pt[3];
};




#endif
