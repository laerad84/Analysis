#include "ConvertTrackData.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "GsimData/GsimGenParticleData.h"
#include "GsimData/GsimTrackData.h"
#include "TMath.h"
#include "TChain.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "rec2g/Rec2g.h"
#include "pi0/Pi0.h"
#include "gamma/Gamma.h"
#include <cstdlib>
const Double_t pdgPi0Mass=134.9766;
int main( int argc, char** argv){
  /*
  TFile* tf = new TFile("CVEtaBgNeutronEta.root");
  TTree* tr = (TTree*)tf->Get("eventTree00");
  */
  TChain* tr = new TChain("RecTree");
  char* name[2]={"Down","Up"};
  int Index = std::atoi(argv[1]);
  tr->Add(Form("Pi0Collection_%s.root",name[Index]));
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );

  TFile* tfout = new TFile(Form("GammaSeed_%s.root",name[Index]),"recreate");
  TTree* trout = new TTree("HEPEvt","HEPEvt");
  int nMax = 2;
  int NHEP;
  int ISTHEP[nMax];
  int IDHEP[nMax];
  int JDAHEP1[nMax];
  int JDAHEP2[nMax];
  double PHEP1[nMax];
  double PHEP2[nMax];
  double PHEP3[nMax];
  double PHEP5[nMax];
  double VHEP1[nMax];
  double VHEP2[nMax];
  double VHEP3[nMax];
  double VHEP4[nMax];
  double PHEP0[nMax];
  double PHEP4[nMax];

  trout->Branch("NHEP",&NHEP,"NHEP/I");
  trout->Branch("IDHEP",IDHEP,"IDHEP[NHEP]/I");//NHEP
  trout->Branch("ISTHEP",ISTHEP,"ISTHEP[NHEP]/I");//NHEP
  trout->Branch("JDAHEP1",JDAHEP1,"JDAHEP1[NHEP]/I");//NHEP
  trout->Branch("JDAHEP2",JDAHEP2,"JDAHEP2[NHEP]/I");//NHEP
  trout->Branch("PHEP1",PHEP1,"PHEP1[NHEP]/D");//NHEP
  trout->Branch("PHEP2",PHEP2,"PHEP2[NHEP]/D");//NHEP
  trout->Branch("PHEP3",PHEP3,"PHEP3[NHEP]/D");//NHEP
  trout->Branch("PHEP5",PHEP5,"PHEP5[NHEP]/D");//NHEP
  trout->Branch("VHEP1",VHEP1,"VHEP1[NHEP]/D");//NHEP
  trout->Branch("VHEP2",VHEP2,"VHEP2[NHEP]/D");//NHEP
  trout->Branch("VHEP3",VHEP3,"VHEP3[NHEP]/D");//NHEP
  trout->Branch("VHEP4",VHEP4,"VHEP4[NHEP]/D");//NHEP
  trout->Branch("PHEP0",PHEP0,"PHEP0[NHEP]/D");//NHEP
  trout->Branch("PHEP4",PHEP4,"PHEP4[NHEP]/D");//NHEP

  for( int i = 0; i< tr->GetEntries(); i++){
    tr->GetEntry(i);    
    NHEP=0;
    if( i% 1000 == 0 ){ std::cout<< i << "/" << tr->GetEntries()<< std::endl;}
    for( int it = 0; it < nMax; it++){
      ISTHEP[it] = 0;
      IDHEP[it] = 0;
      JDAHEP1[it] = 0;
      JDAHEP2[it] = 0;
      PHEP1[it] = 0;
      PHEP2[it] = 0;
      PHEP3[it] = 0;
      PHEP5[it] = 0;
      VHEP1[it] = 0;
      VHEP2[it] = 0;
      VHEP3[it] = 0;
      PHEP0[it] = 0;
      PHEP4[it] = 0;
    }

    std::list<Gamma> glist;
    std::list<Pi0> plist;
    data.getData(glist);
    data.getData(plist);
    if( plist.size() >= 2 ){ continue; }
    std::list<Pi0>::iterator pit=plist.begin();
    Double_t pi0Mass = (*pit).m();

    Double_t Pi0MassPeak[2]={131.678,134.747};
    Double_t Pi0MassPeakSig[2]={6.06779,3.56126};

    if( TMath::Abs(pi0Mass - Pi0MassPeak[Index]) > 2*Pi0MassPeakSig[Index] ){
      continue;
    }
    Double_t Correction = pdgPi0Mass/pi0Mass;
    CLHEP::Hep3Vector gp1 = (*pit).g1().p3();
    CLHEP::Hep3Vector gp2 = (*pit).g1().p3();
    JDAHEP1[NHEP] = 0;
    JDAHEP2[NHEP] = 0;
    ISTHEP[NHEP] = 1;
    IDHEP[NHEP] = 22;
    VHEP1[NHEP] = (*pit).vx();
    VHEP2[NHEP] = (*pit).vy();
    VHEP3[NHEP] = (*pit).vz();
    VHEP4[NHEP] = 0;
    PHEP0[NHEP] = (*pit).g1().p3().mag()/Correction/1000;
    PHEP1[NHEP] = (*pit).g1().p3()[0]/Correction/1000;
    PHEP2[NHEP] = (*pit).g1().p3()[1]/Correction/1000;
    PHEP3[NHEP] = (*pit).g1().p3()[2]/Correction/1000;
    PHEP4[NHEP] = (*pit).g1().p3().mag()/Correction/1000;
    PHEP5[NHEP] = 0;
    NHEP++;
    JDAHEP1[NHEP] = 0;
    JDAHEP2[NHEP] = 0;
    ISTHEP[NHEP] = 1;
    IDHEP[NHEP] = 22;
    VHEP1[NHEP] = (*pit).vx();
    VHEP2[NHEP] = (*pit).vy();
    VHEP3[NHEP] = (*pit).vz();
    VHEP4[NHEP] = 0;
    PHEP0[NHEP] = (*pit).g2().p3().mag()/Correction/1000;
    PHEP1[NHEP] = (*pit).g2().p3()[0]/Correction/1000;
    PHEP2[NHEP] = (*pit).g2().p3()[1]/Correction/1000;
    PHEP3[NHEP] = (*pit).g2().p3()[2]/Correction/1000;
    PHEP4[NHEP] = (*pit).g2().p3().mag()/Correction/1000;
    PHEP5[NHEP] = 0;
    NHEP++;
    if( NHEP== 0){ continue; }
    trout->Fill();
  }

  trout->Write();
  tfout->Close();
  return 0;
}

