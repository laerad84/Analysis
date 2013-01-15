#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include <list>
#include <vector>

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TF1.h"

#include "IDHandler.h"
#include "CsIImage.h"
#include "E14_CLUSTER_BUILDER/E14ReadSumFile.h"
#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include "gamma/GammaFinder.h"

int
main( int argc, char** argv){
  
  gStyle->SetOptStat("neMRIuo");
  gSystem->Load("../E14_ANA_COSMIC/lib/libtest.so");
  IDHandler* handler = new IDHandler(" ../E14_ANA_COSMIC/Data/crystal.txt");
  
  const Int_t nRUN     = 4;
  char* name[4]        ={"Double_et","4Cluster_200MeV","4Cluster_160MeV","3Cluster_370MeV"};
  Int_t RunNumber[4]   ={3969,4089,4094,4097};
  Double_t SECValue[4] ={257727,327229,286614,303859};
  Double_t Spill[4]    ={145,184,162,172}; 


  CsIImage*  image[4];
  for( int i = 0; i< 4; i++){
    image[i] = new CsIImage(handler);
    image[i]->SetTitle(Form("%s",name[i]));
  }
  TApplication* app = new TApplication("app",&argc, argv);
  TH1D* hisKLMass[4];
  TH1D* hisKLMom[4];
  TH1D* hisKLMassCut[4];
  TH1D* hisKLMomCut[4];
  TH1D* TextMassEntries[4];
  TH1D* TextMomEntries[4];
  
  TChain* ch[4];
  for( int i =0 ;i< 4; i++){
    ch[i] = new TChain("Tree");
    ch[i]->Add(Form("klongRootFile/kl%04d.root",RunNumber[i]));
    
    hisKLMass[i] 
      = new TH1D(Form("hisKLMassNoCut_%04d",RunNumber[i]),
		 Form("hisKLMassNoCut_%04d;Mass[MeV];N/5MeV",RunNumber[i]),
		 200,0,1000);

    hisKLMassCut[i]
      = new TH1D(Form("hisKLMassCut_%04d",RunNumber[i]),
		 Form("hisKLMassCut_%04d;Mass[MeV];N/5MeV",RunNumber[i]),
		 200,0,1000);

    hisKLMom[i]
      = new TH1D(Form("hisKLMomNoCut_%04d",RunNumber[i]),
		 Form("hisKLMomNoCut_%04d;Mass[MeV];N/50MeV",RunNumber[i]),
		 160,0,8000);

    hisKLMomCut[i]  
      = new TH1D(Form("hisKLMomCut_%04d",RunNumber[i]),
		 Form("hisKLMomCut_%04d;Mass[MeV];N/50MeV",RunNumber[i]),
		 160,0,8000);

    hisKLMass[i]   ->SetLineColor(i+1);
    hisKLMassCut[i]->SetLineColor(i+1);
    hisKLMom[i]    ->SetLineColor(i+1);
    hisKLMomCut[i] ->SetLineColor(i+1);
  }
  
  E14GNAnaDatacontainer data;
  for( int i = 0; i< 4; i++){
    data.setBranchAddress(ch[i]);
    long nentries = ch[i]->GetEntries();
    for( int ievent = 0; ievent< nentries; ievent++){      
      
      // Analysis event by event
      std::list<Gamma> glist;
      data.getData(glist);
      for( int iCluster = 0; iCluster < data.ClusterNumber ; iCluster++){
	for( int iCrystal = 0; iCrystal < data.ClusterSize[iCLuster]; iCrystal++){
	  image[i]->Fill(data.ClusterCsiId[iCluster][iCrystal]);
	}
      }
    }
  }
    
  TCanvas* can = new TCanvas("can","",1200,600);
  can->Divide(2,1);  
  //  textMassEntries[2]->DrawTextNDC(0.5,0.5,Form("#of Klong(FullCut):%d",(int)hisKLMass[2]->GetEntries()));  
  image->Draw();
  can->SaveAs("Image/KlongTriggerstudy.gif");
  app->Run();
}


