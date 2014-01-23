#include "CosmicTriggerTree.h"
int CosmicTriggerTree::CosmicArr[20]= {4 ,5 ,2 ,3 ,6 ,7 ,0 ,1 ,12,13,
				       10,11,14,15,8 ,9 ,16,17,18,19};
Double_t CosmicTriggerTree::COSMIC_THRESHOLD[20] = {1000,1000,1000,1000,1000,
						    1000,1000,1000,1000,1000,
						    1000,1000,1000,1000,1000,
						    1000,1000,1000,1000,1000};

CosmicTriggerTree::CosmicTriggerTree(){
  InitVar();
}
CosmicTriggerTree::~CosmicTriggerTree(){
  ;
}
void CosmicTriggerTree::InitVar(){
  HitUp = 0;
  HitDn = 0;
  HitUpCoin = 0;
  HitDnCoin = 0;
}
void CosmicTriggerTree::SetBranchAddress(TTree* trin){
  trin->SetBranchAddress("CosmicNumber",&CosmicNumber);
  trin->SetBranchAddress("CosmicSignal",CosmicSignal);
  trin->SetBranchAddress("CosmicTime",CosmicTime);
  trin->SetBranchAddress("CosmicID",CosmicID);
}
void CosmicTriggerTree::Branch(TTree* trOut ){
  trOut->Branch("CosmicNumber",&CosmicNumber,"CosmicNumber/I");
  trOut->Branch("CosmicID",CosmicID,"CosmicID[CosmicNumber]/S");//CosmicNumber
  trOut->Branch("CosmicSignal",CosmicSignal,"CosmicSignal[CosmicNumber]/D");//CosmicNumber
  trOut->Branch("CosmicTime",CosmicTime,"CosmicTime[CosmicNumber]/D");//CosmicNumber

  trOut->Branch("HitUp",&HitUp,"HitUp/I");
  trOut->Branch("HitDn",&HitDn,"HitDn/I");
  trOut->Branch("HitUpCoin",&HitUpCoin,"HitUpCoin/I");
  trOut->Branch("HitDnCoin",&HitDnCoin,"HitDnCoin/I");
}


bool CosmicTriggerTree::TriggerDecision(){
  InitVar();
  Double_t mSignal[20];
  Double_t mTime[20];
  Short_t mID[20];
  for( int i = 0; i< 20; i++){
    mSignal[i] = 0;
    mTime[i]   = 0; 
    mID[i]     = 0;
  }
  
  for( int i = 0; i< CosmicNumber; i++){
    if( CosmicSignal[i] > 1000 ){
      CosmicID[i] = CosmicArr[CosmicID[i]];    
      mSignal[CosmicID[i]] = CosmicSignal[i];
      mTime[CosmicID[i]]   = CosmicTime[i];
      mID[CosmicID[i]]     = CosmicID[i];
    }
  }
  CosmicNumber  =0;
  for( int i = 0; i< 20; i++){
    if( mSignal[i] > 1000 && mTime[i] > 50 ){
      CosmicSignal[CosmicNumber] = mSignal[i];
      CosmicTime[CosmicNumber]   = mTime[i];
      CosmicID[CosmicNumber]     = mID[i];
      CosmicNumber++;
    }
  }

  for( int iCosmic = 0; iCosmic< 5; iCosmic++){
    if( mSignal[ iCosmic    ] > COSMIC_THRESHOLD[ iCosmic     ] &&
	mSignal[ iCosmic+10 ] > COSMIC_THRESHOLD[ iCosmic +10 ]){
      HitUpCoin |= 1 << iCosmic;
    }
    if( mSignal[ iCosmic    ] > COSMIC_THRESHOLD[ iCosmic     ] ||
	mSignal[ iCosmic+10 ] > COSMIC_THRESHOLD[ iCosmic +10 ]){
      HitUp     |= 1 << iCosmic;
    }
    if( mSignal[ iCosmic+5  ] > COSMIC_THRESHOLD[ iCosmic +5  ] &&
	mSignal[ iCosmic+15 ] > COSMIC_THRESHOLD[ iCosmic +15 ]){
      HitDnCoin |= 1 << iCosmic;
    }
    if( mSignal[ iCosmic+5  ] > COSMIC_THRESHOLD[ iCosmic +5  ] ||
	mSignal[ iCosmic+15 ] > COSMIC_THRESHOLD[ iCosmic +15 ]){
      HitDn     |= 1 << iCosmic;
    }       	
  }
  if( HitUp && HitDn ) { 
    return true;
  }else{
    return false;
  }
}
