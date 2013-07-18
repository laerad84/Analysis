#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <list>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "TFile.h"
#include "TTree.h"

#include "TApplication.h"
#include "TSystem.h"
#include "CsIPoly.h"
#include "IDHandler.h"
#include "User_Function.h"

int main( int argc, char** argv ){
  IDHandler* handler = new IDHandler();

  std::string ROOTFILE_SUMUP = std::getenv("ROOTFILE_SUMUP");
  std::string ROOTFILE_WAV   = std::getenv("ROOTFILE_WAV");

  char*  filename[2] = {"SIM","WAV"};
  std::cout<< __PRETTY_FUNCTION__ << std::endl;

  ClusterFinder clusterFinder;
  GammaFinder   gFinder;

  TFile* tfOut = new TFile("Kl_Total_SIMADJ.root","recreate");
  TTree* trOut = new TTree("trKL","Output of Adjusted cluster");

  E14GNAnaDataContainer dataAdj;
  dataAdj.branchOfClusterList( trOut);
  dataAdj.branchOfDigi(trOut);
  dataAdj.branchOfKlong( trOut );


  TFile* tfRoot[2]; 
  TTree* tr[2];
  E14GNAnaDataContainer data[2];
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  
  for( int i = 0; i <2; i++){
    tfRoot[i] = new TFile(Form("Kl_Total_%s.root",filename[i]));
    tr[i] = (TTree*)tfRoot[i]->Get("trKL");
    data[i].setBranchAddress(tr[i]);    
  }

  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  
  CsIPoly* csi = new CsIPoly("csi","csi");
  std::cout<< __PRETTY_FUNCTION__ << std::endl;

  int ifile = 0;
    tr[ifile]->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
    tr[ifile]->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
    trOut->Branch("CsiL1TrigCount",CsiL1TrigCount,"CsIL1TrigCount[20]/D");
    trOut->Branch("CsiL1nTrig",&CsiL1nTrig,"CsiL1nTrig/I");
    for( int ievent = 0; ievent < tr[ifile]->GetEntries(); ievent ++){
      tr[ifile]->GetEntry(ievent);
      dataAdj.reset();
      std::list<Cluster> clist;
      std::list<Gamma>   glist;
      std::vector<Klong> klVec;
      data[ifile].getData(klVec);

      if( klVec.size() == 0 ){ continue; }
      if( CsiL1nTrig< 5 ){ continue; }
      
      Double_t EArr[3000]={0};
      Double_t TArr[3000]={0};      
      for( unsigned int ipi  =0; ipi < klVec[0].pi0().size(); ipi++){
	/////////////////////////////////////////////////////////////////////////////
	/// Search Effect of gap /// 
	/////////////////////////////////////////////////////////////////////////////
	/// Test : Ceneter crystal -2% Crystal around Center crystal 0.5% of Center Crystal ///

	std::vector<int>   vecID = klVec[0].pi0()[ipi].g1().clusterIdVec();
	std::vector<double> vecE  = klVec[0].pi0()[ipi].g1().clusterEVec();
	std::vector<int>::iterator itID = vecID.begin();
	std::vector<double>::iterator itEne = vecE.begin();
	for( ; itID != vecID.end(); itID++,itEne++){
	  double x,y;
	  handler->GetMetricPosition((*itID),x,y);
	  EArr[(*itID)] +=(*itEne)*0.96;      
	  Int_t ID;
	  Double_t delta[4][2] = {{-2.5,0},{2.5,0},{0,-2.5},{0,2.5}};
	  for( int i = 0; i< 4; i++){
	    ID = csi->Fill(x+delta[i][0],y+delta[i][1],1.)-1;
	    EArr[ID]+= (*itEne)*0.01;
	  }
	}
	std::vector<int> vecID1 =  klVec[0].pi0()[ipi].g2().clusterIdVec();
	std::vector<double> vecE1 =  klVec[0].pi0()[ipi].g2().clusterEVec();
	itID  = vecID1.begin();
	itEne = vecE1.begin();
	for( ; itID != vecID1.end(); itID++,itEne++){
	  double x,y;
	  handler->GetMetricPosition((*itID),x,y);
	  EArr[(*itID)] +=(*itEne)*0.98;      
	  Int_t ID;
	  Double_t delta[4][2] = {{-2.5,0},{2.5,0},{0,-2.5},{0,2.5}};
	  for( int i = 0; i< 4; i++){
	    ID = csi->Fill(x+delta[i][0],y+delta[i][1],1.)-1;
	    if( ID <0 || ID >=2716 ){ continue;}
	    EArr[ID]+= (*itEne)*0.005;
	  }
	}
      }
      //Clustering
      Int_t nCrystal=0;
      Double_t EnergyArr[2716]={0};
      Double_t TimeArr[2716]={0};
      Int_t IDArr[2716]={-1};
      
      for( int i = 0; i< 3000; i++){
	if( EArr[i] >3 ){
	  IDArr[nCrystal]    = i;
	  EnergyArr[nCrystal]= EArr[i];
	  nCrystal++;
	}
      }

      std::list<Cluster> newclist;
      std::list<Gamma>   newglist;
      std::vector<Klong> newklVec;
      newclist  = clusterFinder.findCluster( nCrystal, IDArr, EnergyArr, TimeArr);
      gFinder.findGamma( newclist, newglist );
      if( newclist.size() < 6 ){continue; }
      if( newglist.size() != 6){continue; }
      if( user_rec(newglist,newklVec)){
	dataAdj.setData(newclist);
	dataAdj.setData(newglist);
	user_cut(dataAdj,newklVec);
	dataAdj.setData(newklVec);
	trOut->Fill();
      }
    }  
    tfOut->cd();
    trOut->Write();
    tfOut->Close();
}
