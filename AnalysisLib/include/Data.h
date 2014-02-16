#ifndef DATA_H_
#defined DATA_H_
#include "TTree.h"


class Data{
 public:
  
  int RunNo;
  int EventNumber;

  int      CsiNumber;
  double   CsiEne[2716];
  double   CsiID[2716];
  double   CsiTime[2716];
  double   CsiSignal[2716];
  int      CsiL1nTrig;
  double   CsiL1TrigCount[20];

  int      CC03Number;
  Short_t  CC03ID[32];
  Double_t CC03Ene[32];
  Double_t CC03Signal[32];

  int      OEVNumber;
  Short_t  OEVID[44];
  Double_t OEVEne[44];
  Double_t OEVSignal[44];

  int      CVNumber;
  Short_t  CVID[256];
  Double_t CVEne[256];
  Double_t CVSignal[256];

  int      SciNumber;
  Double_t SciEne[5];
  Double_t SciSignal[5];

  int      EtcNumber;
  Double_t EtcEne[20];
  Double_t EtcSignal[20];
  Double_t EtcTime[20];

  int      CosmicNumber;
  Double_t CosmicEne[20];
  Double_t CosmicSignal[20];
  Double_t CosmicTime[20];
  
  Data();
  ~Data();
  void  SetBranchAddress( TTree* tr );
  void  Branch(TTree* tr);
  void  Init();

};


#endif //DATA_H_
