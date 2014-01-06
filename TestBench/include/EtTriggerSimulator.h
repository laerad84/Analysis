#include "TTree.h"

class EtTriggerSimulator {

 public:

  double Et;
  int    EtFlag;
  
  void SetBranchAddress(TTree* tr);
  void Branch(TTree* tr);
  
};


void EtTriggerSimulator::SetBranchAddress(TTree* tr){
  tr->SetBranchAddress("Et",&Et);
  tr->SetBranchAddress("EtFlag",&EtFlag);
}
void EtTriggerSimulator::Branch(TTree* tr){
  tr->Branch("Et",&Et,"Et/D");
  tr->Branch("EtFlag",&EtFlag,"EtFlag/I");
}
void EtTriggerSimulator::EtCalculator( int nArr, double* EArr ){
  Et = 0; 
  for( int i = 0; i< nArr; i++){
    Et += EArr[i];
  }
}
