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
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
int main( int argc, char** argv){
  TFile* tf[3];
  TTree* tr[3];
  char* filetype[3] ={"3pi0Comp","LaserComp","OldComp"};
  for( int i = 0; i< 3 ;i++){
    tf[i] = new TFile(Form("Kl_Total_3pi0_%s.root",filetype[i]));
    tr[i] = (TTree*)tf[i]->Get("trKL");
  }


  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  E14GNAnaDataContainer data;
  TFile* tfOut = new TFile("TestComp.root","recreate");
  const int nFile = 3; 
  TH1D* hisKLZ[nFile];
  TH1D* hisKLE[nFile];
  for( int i = 0; i< 3; i++){
    hisKLE[i] = new TH1D(Form("hisKLE_%d",i),Form("hisKLE_%s",filetype[i]),100,0,10000);
    hisKLZ[i] = new TH1D(Form("hisKLZ_%d",i),Form("hisKLZ_%s",filetype[i]),70,0,7000);
  }
  
  
  for( int i  =0; i< 3; i++){
    data.setBranchAddress(tr[i]);
    tr[i]->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
    tr[i]->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
    for( int ievent = 0; ievent < tr[i]->GetEntries(); ievent++){
      tr[i]->GetEntry(ievent);
      std::list<Cluster> clist;
      std::list<Gamma>   glist;
      std::vector<Klong> klVec;
      data.getData(clist);
      data.getData(glist);
      data.getData(klVec);
      
      Double_t klptSq  = klVec[0].p3()[0]*klVec[0].p3()[0]+klVec[0].p3()[1]*klVec[0].p3()[2];
      Double_t klMom   = TMath::Sqrt(klVec[0].e()*klVec[0].e()-klVec[0].m()*klVec[0].m());
      //Double_t Ratio   = 1; 
      Double_t Ratio   = 1;//sugarFunc->Eval(klMom)/soltFunc->Eval(klMom);

      if( CsiL1nTrig< 5 ){ continue; }
      if( klVec[0].chisqZ() > 6 ){continue;} 
      bool  bInnerGamma = false;
      bool  bOuterGamma = false;
      // Cut on Gamma // 
      Double_t MinGammaE = 200;
      Int_t nEGamma = 0;      
      if( klVec.size() > 1 ){
	if( klVec[1].chisqZ() - klVec[0].chisqZ() < 6 ){ continue; }
      }
      std::list<Gamma>::iterator git = glist.begin();
      for( int igamma = 0; igamma < 6; igamma++,git++){
	if( (*git).e() > MinGammaE ){ nEGamma++;}
	if( TMath::Abs((*git).x()) < 150 &&
	    TMath::Abs((*git).y()) < 150 ){
	  bInnerGamma = true; 	  
	}
	if( TMath::Sqrt((*git).x()*(*git).x()+(*git).y()*(*git).y()) > 850 ){
	  bOuterGamma = true;
	}
	  if( TMath::Abs((*git).y())>575){ bOuterGamma = true; }
      }

      if( nEGamma < 6 ){ continue; }
      if( bInnerGamma ){ continue; }
      //if( klVec[0].vz() < 5000 && klVec[0].vz()>3000){
	hisKLE[i]->Fill(klVec[0].e());
	//}
      hisKLZ[i]->Fill(klVec[0].vz());
    }    
  }

  tfOut->cd();
  for( int i = 0; i< 3; i++){
    hisKLE[i]->Write();
    hisKLZ[i]->Write();
  }
  tfOut->Close();

}
