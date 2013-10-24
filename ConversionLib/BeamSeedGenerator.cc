#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Vector/ThreeVector.h"
int main( int argc, char** argv ){

  TFile* tf = new TFile("beam01545.root");
  TTree* tr = (TTree*)tf->Get("HEPEvt");

  int nMax = 300;
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


  int ANHEP;
  int AISTHEP[nMax];
  int AIDHEP[nMax];
  int AJDAHEP1[nMax];
  int AJDAHEP2[nMax];
  double APHEP1[nMax];
  double APHEP2[nMax];
  double APHEP3[nMax];
  double APHEP5[nMax];
  double AVHEP1[nMax];
  double AVHEP2[nMax];
  double AVHEP3[nMax];
  double AVHEP4[nMax];
  double APHEP0[nMax];
  double APHEP4[nMax];

  tr->SetBranchAddress("NHEP",&NHEP);
  tr->SetBranchAddress("IDHEP",IDHEP);//NHEP
  tr->SetBranchAddress("ISTHEP",ISTHEP);//NHEP
  tr->SetBranchAddress("JDAHEP1",JDAHEP1);//NHEP
  tr->SetBranchAddress("JDAHEP2",JDAHEP2);//NHEP
  tr->SetBranchAddress("PHEP1",PHEP1);//NHEP
  tr->SetBranchAddress("PHEP2",PHEP2);//NHEP
  tr->SetBranchAddress("PHEP3",PHEP3);//NHEP
  tr->SetBranchAddress("PHEP5",PHEP5);//NHEP
  tr->SetBranchAddress("VHEP1",VHEP1);//NHEP
  tr->SetBranchAddress("VHEP2",VHEP2);//NHEP
  tr->SetBranchAddress("VHEP3",VHEP3);//NHEP
  tr->SetBranchAddress("VHEP4",VHEP4);//NHEP
  tr->SetBranchAddress("PHEP0",PHEP0);//NHEP
  tr->SetBranchAddress("PHEP4",PHEP4);//NHEP

  TFile* tfout = new TFile("beam01545_Populate.root","recreate");
  TTree* trout = new TTree("HEPEvt","HEPEvt");
  trout->Branch("NHEP",&ANHEP,"NHEP/I");
  trout->Branch("ISTHEP",AISTHEP,"ISTHEP[NHEP]/I");//NHEP
  trout->Branch("IDHEP",AIDHEP,"IDHEP[NHEP]/I");//NHEP
  trout->Branch("JDAHEP1",AJDAHEP1,"JDAHEP1[NHEP]/I");//NHEP
  trout->Branch("JDAHEP2",AJDAHEP2,"JDAHEP2[NHEP]/I");//NHEP
  trout->Branch("PHEP1",APHEP1,"PHEP1[NHEP]/D");//NHEP
  trout->Branch("PHEP2",APHEP2,"PHEP2[NHEP]/D");//NHEP
  trout->Branch("PHEP3",APHEP3,"PHEP3[NHEP]/D");//NHEP
  trout->Branch("PHEP5",APHEP5,"PHEP5[NHEP]/D");//NHEP
  trout->Branch("VHEP1",AVHEP1,"VHEP1[NHEP]/D");//NHEP
  trout->Branch("VHEP2",AVHEP2,"VHEP2[NHEP]/D");//NHEP
  trout->Branch("VHEP3",AVHEP3,"VHEP3[NHEP]/D");//NHEP
  trout->Branch("VHEP4",AVHEP4,"VHEP4[NHEP]/D");//NHEP
  trout->Branch("PHEP0",APHEP0,"PHEP0[NHEP]/D");//NHEP
  trout->Branch("PHEP4",APHEP4,"PHEP4[NHEP]/D");//NHEP

  int nIt = 1; 
  for( int i = 0; i< tr->GetEntries(); i++ ){
    //for( int i = 0; i< 100; i++ ){
    tr->GetEntry(i);
    for( int j = 0; j < nIt; j++){
      ANHEP=0;
      for( int it = 0; it < nMax; it++){
	AISTHEP[it] = 0;
	AIDHEP[it] = 0;
	AISTHEP[it] = 0;
	AJDAHEP1[it] = 0;
	AJDAHEP2[it] = 0;
	APHEP1[it] = 0;
	APHEP2[it] = 0;
	APHEP3[it] = 0;
	APHEP5[it] = 0;
	AVHEP1[it] = 0;
	AVHEP2[it] = 0;
	AVHEP3[it] = 0;
	APHEP0[it] = 0;
	APHEP4[it] = 0;
      }
      
      for( int itrack = 0; itrack < NHEP; itrack++){
	AISTHEP[itrack] = ISTHEP[itrack];
	AIDHEP[itrack]  = IDHEP[itrack];
	AJDAHEP1[itrack] = JDAHEP1[itrack];
	AJDAHEP2[itrack] = JDAHEP2[itrack];
	
	if(TMath::Abs(VHEP2[itrack]) < 50 &&TMath::Abs(VHEP1[itrack]) < 48 ){ continue; }
	CLHEP::Hep3Vector p   = CLHEP::Hep3Vector(PHEP1[itrack],PHEP2[itrack],PHEP3[itrack]);
	CLHEP::Hep3Vector vtx = CLHEP::Hep3Vector(VHEP1[itrack],VHEP2[itrack],VHEP3[itrack]);
	double phi = CLHEP::RandFlat::shoot()*2*M_PIl;
	//vtx.rotateZ(phi);
	p.rotateZ(phi);	
	double P0=p.mag();
	double P1=P0 + CLHEP::RandGauss::shoot(0., P0*0.02 );
	
	double Theta=p.theta();
	Theta = Theta + CLHEP::RandGauss::shoot(0., Theta*0.01 );
	
	double Phi=p.phi();
	Phi = Phi + CLHEP::RandGauss::shoot(0., Phi*0.001 );
	
	double px = P1 * sin(Theta) * cos(Phi);
	double py = P1 * sin(Theta) * sin(Phi);
	double pz = P1 * cos(Theta);
	p = CLHEP::Hep3Vector(px, py, pz);

	APHEP1[ANHEP] = px;
	APHEP2[ANHEP] = py; 
	APHEP3[ANHEP] = pz;
	AVHEP1[ANHEP] = vtx[0];
	AVHEP2[ANHEP] = vtx[1];
	AVHEP3[ANHEP] = vtx[2];
	AVHEP4[ANHEP] = 0;
	
	APHEP4[ANHEP] = sqrt(p.mag()*p.mag()+PHEP5[itrack]*PHEP5[itrack])-PHEP5[itrack];
	APHEP5[ANHEP] = PHEP5[itrack];
	APHEP0[ANHEP] = p.mag();
	ANHEP++;
      }
      if( ANHEP > 0 ){
	trout->Fill();
      }
    }
    
  }
  trout->Write();
  tfout->Close();

}
