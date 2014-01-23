#include "TTree.h"

class CosmicTriggerTree{
 public:
  CosmicTriggerTree();
  ~CosmicTriggerTree();
  void SetBranchAddress(TTree* trin);
  void Branch(TTree* trOut);
  bool TriggerDecision();
  static int CosmicArr[20];
  static Double_t COSMIC_THRESHOLD[20];
  
  void InitVar();
  int HitUp;
  int HitDn;
  int HitUpCoin;
  int HitDnCoin;
  short UpID[5];
  short DnID[5];
  short UpCoinID[5];
  short DnCoinID[5];

  Int_t    CosmicNumber;
  Double_t CosmicSignal[20];
  Double_t CosmicTime[20];
  Short_t CosmicID[20];
  
};
