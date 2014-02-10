#ifndef USER_FUNCTIONS_H_
#define USER_FUNCTIONS_H_
#include <list>
#include <vector>
#include <iostream>

#include "rec2g/Rec2g.h"
#include "gamma/Gamma.h"
#include "pi0/Pi0.h"
#include "klong/Klong.h"
#include "klong/RecKlong.h"
#include "gnana/E14GNAnaFunction.h"
#include "csimap/CsiMap.h"


#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"

#include "L1TrigCounter.h"
#include "CrateIDHandler.h"

bool User_RecG2(std::list<Gamma> const &glist, std::list<Pi0>& piList);
bool User_RecG4(std::list<Gamma> const &glist, std::vector<Klong>& klVec);
bool User_RecG6(std::list<Gamma> const &glist, std::vector<Klong>& klVec);
double const SpeedOfLight = 299.792458;

void SetGammaTime(Gamma &g);
void SetGammaTime(std::list<Gamma> glist);
double GetWeight(Gamma g);
double GetTiming(Gamma g);
double GetClusterTSigma(Gamma g);
double CalGammaTOF( Klong kl, Gamma g );

class GammaCut{
 public:
  GammaCut();
  ~GammaCut();

  void Branch( TTree* tr );
  void SetBranchAddress( TTree* tr );
  void Decision(Klong kl);
  void Decision(std::vector<Klong> kl);
  void Decision(std::list<Gamma>   g);
  void Decision(Pi0   pi);
  void Reset();
  
  double GMinEnergy;
  double GMinDist; 
  double GMinX;
  double GMinY;
  double GMaxR;
  double GMaxY;
  double GMaxChi;
  double GMaxDeltaT;
  int    GMaxDeltaTID;

};

class CsiCut{
 public:
  CsiCut();
  ~CsiCut();
  void Branch( TTree *tr );
  void SetBranchAddress( TTree* tr );
  void Decision( int  ICsiNumber, int* ICsiID, double* ICsiEne, double* ICsiTime ,double* ICsiSignal, double* ICsiChisq, short* ICsiNDF);
  void Reset();
  void SetCutValue(double wTimeWindow);

  int    CsiNumber;
  int    CsiID[2716];
  double CsiEne[2716];
  double CsiTime[2716];
  double CsiHHTime[2716];
  double CsiSignal[2716];
  double CsiL1TrigCount[20];
  int    CsiL1nTrig;
  double CsiChisq[2716];
  short  CsiNDF[2716];
  short  CsiCrate[2716];
  short  CsiL1[2716];
  short  CsiGB[2716];
  short  CsiPosID[2716];

  TH1D*  hisCsiTime;
  double CsiEventTime;
  double CsiEventTimeSigma;
 private:
  double mWTimeWindow;
  double CsiL1TrigCountThreshold[20];
  CrateIDHandler* CIDHandler;
  CsiMap*         map;
  L1TrigCounter* l1;

};
#endif //USER_FUNCTIONS_H_
