#include "Data.h"

Data::Data(){
  Init();
}
Data::~Data(){
  ;
}

void SetBranchAddress( TTree* tr ){
  tr->SetBranchAddress("RunNumber",&RunNumber);
  tr->SetBranchAddress("EventNumber",&EventNumber);

  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiID"    ,CsiID     );//CsiNumber
  tr->SetBranchAddress("CsiSignal",CsiSignal );//CsiNumber
  tr->SetBranchAddress("CsiEne"   ,CsiEne    );//CsiNumber
  tr->SetBranchAddress("CsiTime"  ,CsiTime   );//CsiNumber

  tr->SetBranchAddress("CVNumber",&CVNumber);
  tr->SetBranchAddress("CVID"    ,CVID     );//CVNumber
  tr->SetBranchAddress("CVSignal",CVSignal );//CVNumber
  tr->SetBranchAddress("CVEne"   ,CVEne    );//CVNumber
  tr->SetBranchAddress("CVTime"  ,CVTime   );//CVNumber

  tr->SetBranchAddress("CC03Number",&CC03Number);
  tr->SetBranchAddress("CC03ID"    ,CC03ID     );//CC03Number
  tr->SetBranchAddress("CC03Signal",CC03Signal );//CC03Number
  tr->SetBranchAddress("CC03Ene"   ,CC03Ene    );//CC03Number
  tr->SetBranchAddress("CC03Time"  ,CC03Time   );//CC03Number

  tr->SetBranchAddress("OEVNumber",&OEVNumber);
  tr->SetBranchAddress("OEVID"    ,OEVID     );//OEVNumber
  tr->SetBranchAddress("OEVSignal",OEVSignal );//OEVNumber
  tr->SetBranchAddress("OEVEne"   ,OEVEne    );//OEVNumber
  tr->SetBranchAddress("OEVTime"  ,OEVTime   );//OEVNumber

  tr->SetBranchAddress("SciNumber",&SciNumber);
  tr->SetBranchAddress("SciID"    ,SciID     );//SciNumber
  tr->SetBranchAddress("SciSignal",SciSignal );//SciNumber
  tr->SetBranchAddress("SciEne"   ,SciEne    );//SciNumber
  tr->SetBranchAddress("SciTime"  ,SciTime   );//SciNumber

  tr->SetBranchAddress("EtcNumber",&EtcNumber);
  tr->SetBranchAddress("EtcID"    ,EtcID     );//EtcNumber
  tr->SetBranchAddress("EtcSignal",EtcSignal );//EtcNumber
  tr->SetBranchAddress("EtcEne"   ,EtcEne    );//EtcNumber
  tr->SetBranchAddress("EtcTime"  ,EtcTime   );//EtcNumber

  tr->SetBranchAddress("CosmicNumber",&CosmicNumber);
  tr->SetBranchAddress("CosmicID"    ,CosmicID     );//CosmicNumber
  tr->SetBranchAddress("CosmicSignal",CosmicSignal );//CosmicNumber
  tr->SetBranchAddress("CosmicEne"   ,CosmicEne    );//CosmicNumber
  tr->SetBranchAddress("CosmicTime"  ,CosmicTime   );//CosmicNumber

}

void SetBranchAddress( TTree* tr ){

  tr->Branch("RunNumber",&RunNumber,"RunNumber/I");
  tr->Branch("EventNumber",&EventNumber,"EventNumber/I");

  tr->Branch("CsiNumber",&CsiNumber,"CsiNumber/I");
  tr->Branch("CsiID"    ,CsiID     ,"CsiID[CsiNumber]/I");//CsiNumber
  tr->Branch("CsiSignal",CsiSignal ,"CsiSignal[CsiNumber]/D");//CsiNumber
  tr->Branch("CsiEne"   ,CsiEne    ,"CsiEne[CsiNumber]/D");//CsiNumber
  tr->Branch("CsiTime"  ,CsiTime   ,"CsiTime[CsiNumber]/D");//CsiNumber

  tr->Branch("CVNumber",&CVNumber,"CVNumber/I");
  tr->Branch("CVID"    ,CVID     ,"CVID[CVNumber]/I");//CVNumber
  tr->Branch("CVSignal",CVSignal ,"CVSignal[CVNumber]/D");//CVNumber
  tr->Branch("CVEne"   ,CVEne    ,"CVEne[CVNumber]/D");//CVNumber
  tr->Branch("CVTime"  ,CVTime   ,"CVTime[CVNumber]/D");//CVNumber

  tr->Branch("CC03Number",&CC03Number,"CC03Number/I");
  tr->Branch("CC03ID"    ,CC03ID     ,"CC03ID[CC03Number]/I");//CC03Number
  tr->Branch("CC03Signal",CC03Signal ,"CC03Signal[CC03Number]/D");//CC03Number
  tr->Branch("CC03Ene"   ,CC03Ene    ,"CC03Ene[CC03Number]/D");//CC03Number
  tr->Branch("CC03Time"  ,CC03Time   ,"CC03Time[CC03Number]/D");//CC03Number

  tr->Branch("OEVNumber",&OEVNumber,"OEVNumber/I");
  tr->Branch("OEVID"    ,OEVID     ,"OEVID[OEVNumber]/I");//OEVNumber
  tr->Branch("OEVSignal",OEVSignal ,"OEVSignal[OEVNumber]/D");//OEVNumber
  tr->Branch("OEVEne"   ,OEVEne    ,"OEVEne[OEVNumber]/D");//OEVNumber
  tr->Branch("OEVTime"  ,OEVTime   ,"OEVTime[OEVNumber]/D");//OEVNumber

  tr->Branch("SciNumber",&SciNumber,"SciNumber/I");
  tr->Branch("SciID"    ,SciID     ,"SciID[SciNumber]/I");//SciNumber
  tr->Branch("SciSignal",SciSignal ,"SciSignal[SciNumber]/D");//SciNumber
  tr->Branch("SciEne"   ,SciEne    ,"SciEne[SciNumber]/D");//SciNumber
  tr->Branch("SciTime"  ,SciTime   ,"SciTime[SciNumber]/D");//SciNumber

  tr->Branch("EtcNumber",&EtcNumber,"EtcNumber/I");
  tr->Branch("EtcID"    ,EtcID     ,"EtcID[EtcNumber]/I");//EtcNumber
  tr->Branch("EtcSignal",EtcSignal ,"EtcSignal[EtcNumber]/D");//EtcNumber
  tr->Branch("EtcEne"   ,EtcEne    ,"EtcEne[EtcNumber]/D");//EtcNumber
  tr->Branch("EtcTime"  ,EtcTime   ,"EtcTime[EtcNumber]/D");//EtcNumber

  tr->Branch("CosmicNumber",&CosmicNumber,"CosmicNumber/I");
  tr->Branch("CosmicID"    ,CosmicID     ,"CosmicID[CosmicNumber]/I");//CosmicNumber
  tr->Branch("CosmicSignal",CosmicSignal ,"CosmicSignal[CosmicNumber]/D");//CosmicNumber
  tr->Branch("CosmicEne"   ,CosmicEne    ,"CosmicEne[CosmicNumber]/D");//CosmicNumber
  tr->Branch("CosmicTime"  ,CosmicTime   ,"CosmicTime[CosmicNumber]/D");//CosmicNumber

}

void Init(){
  Int_t nDetector[7]={2716,32,44,256,5,20,20};
  CsiNumber = 0; 
  for( int i = 0; i< 2716; i++){
    CsiID[i] = -1; 
    CsiSignal[i] = 0;
    CsiEne[i] = 0;
    CsiTime[i] = 0;
  }

  CC03Number = 0; 
  for( int i = 0; i< 32; i++){
    CC03ID[i] = -1; 
    CC03Signal[i] = 0;
    CC03Ene[i] = 0;
    CC03Time[i] = 0;
  }

  OEVNumber = 0; 
  for( int i = 0; i< 44; i++){
    OEVID[i] = -1; 
    OEVSignal[i] = 0;
    OEVEne[i] = 0;
    OEVTime[i] = 0;
  }

  CVNumber = 0; 
  for( int i = 0; i< 256; i++){
    CVID[i] = -1; 
    CVSignal[i] = 0;
    CVEne[i] = 0;
    CVTime[i] = 0;
  }

  SciNumber = 0; 
  for( int i = 0; i< 5; i++){
    SciID[i] = -1; 
    SciSignal[i] = 0;
    SciEne[i] = 0;
    SciTime[i] = 0;
  }

  EtcNumber = 0; 
  for( int i = 0; i< 20; i++){
    EtcID[i] = -1; 
    EtcSignal[i] = 0;
    EtcEne[i] = 0;
    EtcTime[i] = 0;
  }

  CosmicNumber = 0; 
  for( int i = 0; i< 20; i++){
    CosmicID[i] = -1; 
    CosmicSignal[i] = 0;
    CosmicEne[i] = 0;
    CosmicTime[i] = 0;
  }

}
