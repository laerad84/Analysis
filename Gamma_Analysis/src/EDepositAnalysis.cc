#include "EDepositAnalysis.h"

void EDepositAnalysis::ConvertPosition(const  Double_t Radius,const  Double_t Theta, const Double_t x, const  Double_t y, Double_t &nx, Double_t &ny ){
  nx = (x+Radius)*TMath::Cos( Theta ) - y*TMath::Sin(Theta);
  ny = (x+Radius)*TMath::Sin( Theta ) - y*TMath::Cos(Theta);
}
const double EDepositAnalysis::Speed_of_Signal = 100.;//[mm/ns]
EDepositAnalysis::EDepositAnalysis(){
  Init();
}
EDepositAnalysis::~EDepositAnalysis(){
  ;
}
void EDepositAnalysis::Init(){
  PrepareDump();
  PrepareIO();
}
int  EDepositAnalysis::PrepareDump(){
  gen    = new PulseGenerator();
  CsIEne = new CsIPoly("CsIEnergy", "CsIEnergy");
  clusterFinder = new ClusterFinder();
  return 0;
}
int  EDepositAnalysis::ResetDump(){
  CsIEne->Reset();
  return 0;
}

void EDepositAnalysis::PrepareIO(){
  chain = new TChain("eventTree00");
  chain->Add(Form("/Volume0/gamma/template_gamma_210MeV_10deg-1E5-0.root"));
  trin  = new EventTree( chain );

  OutputFile = new TFile("Cluster_210MeV_10deg-0.root","recreate");
  OutputTree = new TTree("T", "Clustered data" );

  OutputTree->Branch("nDigi"  ,&nDigi    ,"nDigi/I");
  OutputTree->Branch("CsiID"  ,ID        ,"CsIID[nDigi]/I");//nDigi
  OutputTree->Branch("CsiEne" ,Energy    ,"CsiEne[nDigi]/D");//nDigi
  OutputTree->Branch("CsiTime",SignalTime,"CsiTime[nDigi]/D");//nDigi  
  data = new E14GNAnaDataContainer();
  data->branchOfClusterList( OutputTree );

}
void EDepositAnalysis::DrawEvent(){
  CsIEne->Draw("colz L");
}
  
void EDepositAnalysis::Export(){
  OutputTree->Fill();  
}
void EDepositAnalysis::Close(){
  OutputTree->Write();
  OutputFile->Close();
}
int  EDepositAnalysis::EventProcess( int ievent ){
  trin->GetEntry( ievent );
  ResetDump();
  double TotalEnergy = 0; 
  std::vector<int> CrystalIDVec;
  std::vector<double> CrystalEVec;
  std::vector<int>::iterator itIDVec;
  std::vector<double>::iterator itEVec;
  std::vector<int> hitIDVec;
  std::vector<double> hitZVec;
  std::vector<double> hitEVec;
  std::vector<double> hitTVec; 

  std::list<Cluster> clist;
  Int_t CsIID = -1;
  Double_t CsiPosX, CsiPosY;
  Double_t Radius = gRandom->Rndm()*300+200; 
  Double_t Theta  = gRandom->Rndm()*2*TMath::Pi();
  for( int ihit = 0; ihit < trin->CSI_hits_; ihit++){
    if( trin->CSI_hits_edep[ihit] == 0 ){ continue; }
    TotalEnergy += trin->CSI_hits_edep[ihit];

    ConvertPosition( Radius, Theta, trin->CSI_hits_r_fX[ihit], trin->CSI_hits_r_fY[ihit], CsiPosX, CsiPosY );

    CsIID = CsIEne->Fill( CsiPosX, CsiPosY, trin->CSI_hits_edep[ihit] ) -1;// bin number  -1 = crystalID;
    if( CsIID < 0 ){ continue; }// no Hit on Configured area;
    
    hitIDVec.push_back( CsIID );
    hitZVec.push_back( trin->CSI_hits_r_fZ[ihit] );
    hitEVec.push_back( trin->CSI_hits_edep[ihit] );
    hitTVec.push_back( trin->CSI_hits_time[ihit] + (500 - trin->CSI_hits_r_fZ[ihit])/Speed_of_Signal);
    
    // Decision ID is new or already added // 
    for( itIDVec = CrystalIDVec.begin(); itIDVec != CrystalIDVec.end(); itIDVec++){
      if(( *itIDVec) == CsIID ){ break; }
    }
    if( itIDVec == CrystalIDVec.end() ){
      CrystalIDVec.push_back( CsIID );
    }
  }

  for( itIDVec = CrystalIDVec.begin(); itIDVec != CrystalIDVec.end() ; itIDVec++){
    Double_t Energy = CsIEne->GetBinContent( (*itIDVec) +1 ); 
    CrystalEVec.push_back( Energy );
  }
  
  itIDVec = CrystalIDVec.begin();
  itEVec  = CrystalEVec.begin();    

  // init OutputData // 
  nDigi = 0;
  for( int i = 0; i< 2716; i++){
    ID[i]  = 0;
    Energy[i] = 0.;
    SignalTime[i]  = 0.;
  }

  for(;itIDVec != CrystalIDVec.end() && itEVec != CrystalEVec.end() ; itIDVec++, itEVec++){
    if( *itEVec < 2 ){ continue; }
    std::vector<double> EVec;
    std::vector<double> TVec;
    for( int ihit = 0; ihit < hitZVec.size(); ihit++){
      if( hitIDVec[ihit] == (*itIDVec)){
	EVec.push_back( hitEVec[ihit] );
	TVec.push_back( hitTVec[ihit] );
      }
    }
    TF1* func = gen->GetWaveform( EVec, TVec, 150./13.7 );
    ID[nDigi] = (*itIDVec);
    Energy[nDigi] = func->GetMaximum( 0, 500 )/(150./13.7);
    SignalTime[nDigi] = func->GetMaximumX( 0 ,500 );
    nDigi++;
    delete func;
  }
  clist = clusterFinder->findCluster( nDigi, ID, Energy, SignalTime );
  data->setData( clist );
  // Clustering // 
  return 0;
}
int  EDepositAnalysis::Loop( int Entries ){
  
  int nLoop = trin->fChain->GetEntries();
  int nProcessed = 0; 
  if( Entries > 0  && Entries < nLoop){
    nLoop  = Entries; 
  }
  for( int ievent = 0; ievent < nLoop; ievent++){
    //std::cout << ievent << std::endl;
    ResetDump();
    EventProcess ( ievent );
    Export();
    nProcessed++;
  }
  return nProcessed;
}



