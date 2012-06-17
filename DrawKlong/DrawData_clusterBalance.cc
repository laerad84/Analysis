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
#include "TStyle.h"

#include "IDHandler.h"
#include "CsIImage.h"
#include "E14ReadSumFile.h"
#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include "gamma/GammaFinder.h"

int
main( int argc, char** argv){

  std::cout<<  __LINE__ << std::endl;  
  gStyle->SetOptStat("neMRIuo");
  std::cout<<  __LINE__ << std::endl;  
  IDHandler* handler = new IDHandler();  
  std::cout<<  __LINE__ << std::endl;  
  const Int_t nRUN     = 4;
  std::cout<<  __LINE__ << std::endl;  
  char* name[4]        ={"Double_et","4Cluster_200MeV",
			 "4Cluster_160MeV","3Cluster_370MeV"};
  
  Int_t    RunNumber[4] = {3969,4089,4094,4097};
  Double_t SECValue[4]  = {257727,327229,286614,303859};
  Double_t Spill[4]     = {145,184,162,172};   
  
  std::cout<< __LINE__ << std::endl;
  
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
  
  std::cout<< __LINE__ << std::endl;
  
  TChain* ch[4];
  for( int i =0 ;i< 4; i++){
    ch[i] = new TChain("tro");
    //ch[i]->Add(Form("klongRootFile/kl%04d.root",RunNumber[i]));
    ch[i]->Add(Form("clusterRootFile/cl%04d.root",RunNumber[i]));
    //std::cout<< ch[i]->GetEntries() << std::endl;
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

  Int_t eventID;
  std::cout << __LINE__ << std::endl;
  E14GNAnaDataContainer data[4];


  TH1D* hisTrig[4];
  TH2D* hisrlBalance[4];
  for( int i = 0; i< 4; i++){
    hisTrig[i] = new TH1D(Form("hisTrig_%d",i),"Trig",4,0,4);
    hisrlBalance[i] = new TH2D(Form("hisrlBalance%d",i),"TrigBalance",400,-1,1,120,0,6000);
  }

  for( int i = 0; i< 4; i++){
    data[i].setBranchAddress(ch[i]);    
    //ch[i]->SetBranchAddress("eventID",&eventID);
    long nentries = ch[i]->GetEntries();
    for( int ievent = 0; ievent< nentries; ievent++){      
      //std::cout << eventID << std::endl;
      ch[i]->GetEntry(ievent);      
      /*
      Klong kl;
      data[i].getData(kl);
      //std::cout <<kl << std::endl;     
      if( kl.vz() > 5000 || kl.vz()<3000){
	continue;
      }
      */
      /*
      if(( data[i].CutCondition & 1 )!=0){
	continue;
      }
      */
      
      Int_t leftHit  = 0;
      Int_t rightHit = 0;
      Int_t hit=0;
      Int_t MultiHit = 0;
      Double_t EnergyRight = 0;
      Double_t EnergyLeft  = 0;
      Double_t EnergyTotal = 0;
      

      for( int iCluster = 0;
	   iCluster < data[i].ClusterNumber;
	   iCluster++){	
	if( data[i].ClusterCoePos[iCluster][0] > 0){	    	  
	  EnergyLeft+=data[i].ClusterDepE[iCluster];
	}else{
	  EnergyRight+=data[i].ClusterDepE[iCluster];
	}
	
	if( data[i].ClusterDepE[iCluster] > 200 ){
	  if( data[i].ClusterCoePos[iCluster][0] > 0){	    	  
	    leftHit++;
	  }else{
	    rightHit++;
	  }
	}
	
	/*
	  for( int iCrystal = 0; 
	  iCrystal < data[i].GamClusSize[iCluster];
	  iCrystal++){	  
	  std::cout << data[i].GamClusCsiId[iCluster][iCrystal] << std::endl;
	  }
	*/
      }

      if(leftHit>0){
	hit |= 1;
      }
      if(rightHit>0){
	hit |= 2;
      }
	
      EnergyTotal = EnergyRight + EnergyLeft;
      hisrlBalance[i]->Fill((EnergyLeft -EnergyRight)/EnergyTotal,EnergyTotal); 
      hisTrig[i]->Fill(hit);	
      
    }
  }
  
  TCanvas* can = new TCanvas("can","",1000,1000);
  can->Divide(2,2);  
  // textMassEntries[2]->DrawTextNDC(0.5,0.5,Form("#of Klong(FullCut):%d",(int)hisKLMass[2]->GetEntries()));    
  // hisrlBalance
  // hisTrig
  
  can->cd(1);
  hisTrig[0]->Draw();
  can->cd(2);
  hisTrig[1]->Draw();
  can->cd(3);
  hisTrig[2]->Draw();
  can->cd(4);
  hisTrig[3]->Draw();
  
  can->SaveAs("Image/KlongTriggerstudy.gif");
  app->Run();
}


