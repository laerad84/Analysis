#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TApplication.h"
#include "ClusterTimeReader.h"
#include "TMath.h"
#include "TChain.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"

int main( int argc, char** argv ){

  Int_t Energy = atoi( argv[1] );
  Int_t Theta  = atoi( argv[2] );
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string iFileForm = "%s/Cluster_Time_%dMeV_%ddeg-1E5-%d.root";//ROOTFILE_GAMMACLUS
  std::string oFileForm = "%s/ClusterTimeStructure_SIM_%dMeV_%ddeg.root";//ROOTFILE_GAMMACLUS

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// I/O File Setting ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TChain* ch = new TChain("trCluster");
  for( int iIndex =0; iIndex < 1; iIndex++){
    ch->Add(Form(iFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str(),Energy,Theta,iIndex));
  }
  /*
  TFile* tf = new TFile(Form("%s/Data_All.root",ROOTFILE_GAMMACLUS.c_str()));
  TTree* tr = (TTree*)tf->Get("trCluster");
  */
  ClusterTimeReader* reader = new ClusterTimeReader(ch);
  Int_t nEntries = reader->fChain->GetEntries();   
  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str(),Energy,Theta),"recreate");  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  const int nE = 6;
  const int nTheta = 8; 
  const int nBinsR = 41;
  const int nBinsPhi = 41;
  Double_t width = 12.5;
  Int_t    nBinsRD  = 19;
  Double_t RDMax = width*nBinsRD/2.;
  Double_t RDMin = -1*width*nBinsRD/2.; 

  TProfile2D* profRDT;
  TProfile2D* profRDE;
  TProfile2D* profRDE_NR;
  TProfile2D* profRDE_Low;
  TProfile2D* profRDE_High;
  TProfile2D* profRDE_INV;  

  TProfile2D* profRDT_NR;
  TProfile2D* profRDT_Low;
  TProfile2D* profRDT_High;
  TProfile2D* profRDT_INV;  

  TProfile2D* profRET;
  TProfile2D* profDET;
  TProfile2D* profEET;
  TProfile2D* profEET_D;
  TH2D*     hisRT = new TH2D(Form("hisRT_E_%d_Theta_%d",Energy,Theta),
			     Form("hisRT_E_%d_Theta_%d",Energy,Theta),
			     nBinsRD,RDMin,RDMax,100,-10,10);
  TH2D*     hisDT = new TH2D(Form("hisDT_E_%d_Theta_%d",Energy,Theta),
			     Form("hisDT_E_%d_Theta_%d",Energy,Theta),
			     nBinsRD,RDMin,RDMax,100,-10,10);
  TProfile* profRT;
  TProfile* profDT;
  TProfile* profRE;
  TProfile* profDE;
  TProfile* profRT_NR;
  TProfile* profDT_NR;
  TProfile* profRT_High;
  TProfile* profDT_High;
  TProfile* profDT_Low;
  TProfile* profRT_Low;
  TProfile* profET_D;

  
  profRDT = new TProfile2D(Form("profRDT_E_%d_Theta_%d",Energy,Theta),
			   Form("profRDT_E_%d_Theta_%d",Energy,Theta),
			   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
  profRDT_NR = new TProfile2D(Form("profRDT_NR_E_%d_Theta_%d",Energy,Theta),
			      Form("profRDT_NR_E_%d_Theta_%d",Energy,Theta),
			      nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
  profRDT_Low = new TProfile2D(Form("profRDT_Low_E_%d_Theta_%d",Energy,Theta),
			      Form("profRDT_Low_E_%d_Theta_%d",Energy,Theta),
			      nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
  profRDT_High = new TProfile2D(Form("profRDT_High_E_%d_Theta_%d",Energy,Theta),
			      Form("profRDT_High_E_%d_Theta_%d",Energy,Theta),
			      nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
  profRDT_INV = new TProfile2D(Form("profRDT_INV_E_%d_Theta_%d",Energy,Theta),
			      Form("profRDT_INV_E_%d_Theta_%d",Energy,Theta),
			      nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
  profRDE = new TProfile2D(Form("profRDE_E_%d_Theta_%d",Energy,Theta),
			   Form("profRDE_E_%d_Theta_%d",Energy,Theta),
			   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
  profRDE_NR = new TProfile2D(Form("profRDE_NR_E_%d_Theta_%d",Energy,Theta),
			      Form("profRDE_NR_E_%d_Theta_%d",Energy,Theta),
			      nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
  profRDE_Low = new TProfile2D(Form("profRDE_Low_E_%d_Theta_%d",Energy,Theta),
			      Form("profRDE_Low_E_%d_Theta_%d",Energy,Theta),
			      nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
  profRDE_High = new TProfile2D(Form("profRDE_High_E_%d_Theta_%d",Energy,Theta),
			      Form("profRDE_High_E_%d_Theta_%d",Energy,Theta),
			      nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
  profRDE_INV = new TProfile2D(Form("profRDE_INV_E_%d_Theta_%d",Energy,Theta),
			      Form("profRDE_INV_E_%d_Theta_%d",Energy,Theta),
			      nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);

  
  profRET = new TProfile2D(Form("profRET_E_%d_Theta_%d",Energy,Theta),
			   Form("profRET_E_%d_Theta_%d",Energy,Theta),
			   nBinsRD,RDMin,RDMax,100,0,400);
  profDET = new TProfile2D(Form("profDET_E_%d_Theta_%d",Energy,Theta),
			   Form("profDET_E_%d_Theta_%d",Energy,Theta),
			   nBinsRD,RDMin,RDMax,100,0,400);
  profEET = new TProfile2D(Form("profEET_E_%d_Theta_%d",Energy,Theta),
			   Form("profEET_E_%d_Theta_%d",Energy,Theta),
			   20,0,400,20,0,400);
  profEET_D = new TProfile2D(Form("profEET_D_E_%d_Theta_%d",Energy,Theta),
			     Form("profEET_D_E_%d_Theta_%d",Energy,Theta),
			     20,0,400,20,0,400);
  profET_D = new TProfile(Form("profET_D_E_%d_Theta_%d",Energy,Theta),
			  Form("profET_D_E_%d_Theta_%d",Energy,Theta),
			  20,0,400);
  profRT = new TProfile(Form("profRT_E_%d_Theta_%d",Energy,Theta),
			Form("profRT_E_%d_Theta_%d",Energy,Theta),
			nBinsRD,RDMin,RDMax);  
  profRE = new TProfile(Form("profRE_E_%d_Theta_%d",Energy,Theta),
			Form("profRE_E_%d_Theta_%d",Energy,Theta),
			nBinsRD,RDMin,RDMax);  
  profDT = new TProfile(Form("profDT_E_%d_Theta_%d",Energy,Theta),
			Form("profDT_E_%d_Theta_%d",Energy,Theta),
			nBinsRD,RDMin,RDMax);  
  profDE = new TProfile(Form("profDE_E_%d_Theta_%d",Energy,Theta),
			Form("profDE_E_%d_Theta_%d",Energy,Theta),
			nBinsRD,RDMin,RDMax);  
  profRT_NR = new TProfile(Form("profRT_NR_E_%d_Theta_%d",Energy,Theta),
			   Form("profRT_NR_E_%d_Theta_%d",Energy,Theta),
			   nBinsRD,RDMin,RDMax);  
  profDT_NR = new TProfile(Form("profDT_NR_E_%d_Theta_%d",Energy,Theta),
			   Form("profDT_NR_E_%d_Theta_%d",Energy,Theta),
			   nBinsRD,RDMin,RDMax);  
  profRT_High = new TProfile(Form("profRT_High_E_%d_Theta_%d",Energy,Theta),
			     Form("profRT_High_E_%d_Theta_%d",Energy,Theta),
			     nBinsRD,RDMin,RDMax);  
  profDT_High = new TProfile(Form("profDT_High_E_%d_Theta_%d",Energy,Theta),
			     Form("profDT_High_E_%d_Theta_%d",Energy,Theta),
			     nBinsRD,RDMin,RDMax);  
  profRT_Low = new TProfile(Form("profRT_Low_E_%d_Theta_%d",Energy,Theta),
			    Form("profRT_Low_E_%d_Theta_%d",Energy,Theta),
			    nBinsRD,RDMin,RDMax);  
  profDT_Low = new TProfile(Form("profDT_Low_E_%d_Theta_%d",Energy,Theta),
			     Form("profDT_Low_E_%d_Theta_%d",Energy,Theta),
			     nBinsRD,RDMin,RDMax);  
  ///////////////////////////////////////////////////////////////////////
  /// Cut Condition 
  ///////////////////////////////////////////////////////////////////////
  Double_t OuterRadCut = 500;//mm
  Double_t InnerRadCut = 250;//mm
  Double_t LowEnergyCut= 3;//MeV
  Int_t    ClusterSizeCut = 6;//n
  Double_t RDCutNar = 12.5;//mm
  Double_t RDCutWid = 25.;//mm

  std::cout<< nEntries << std::endl;
  for( int evtIndex = 0; evtIndex < nEntries; evtIndex++){
    reader->GetEntry( evtIndex );
    if( evtIndex %1000 == 0 && evtIndex != 1000){ std::cout<< evtIndex <<"/" << nEntries<< std::endl;}
    //if( reader->nCluster > 1 ){ continue; }
    
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){      
      if( reader->ClusterR[clusterIndex] > OuterRadCut ){continue;}
      if( reader->ClusterR[clusterIndex] < InnerRadCut ){continue;}
      if( reader->nCrystal[clusterIndex] < ClusterSizeCut ){ continue; }

      Double_t RCenterCrystal = 1000;
      Double_t ECenterCrystal = 0;
      Double_t TCenterCrystal = 0;
      Double_t RCenter = 0;
      Double_t DCenter = 0;
      for( int crystalIndex  = 0; crystalIndex < reader->nCrystal[ clusterIndex ]; crystalIndex++){       
	if( reader->CrystalEnergy[clusterIndex][crystalIndex] <= LowEnergyCut){ continue; }	
	
	// Get CenterCrystal Data
	for( int crystalIndex = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	  Double_t PhiinCluster = reader->CrystalPhi[clusterIndex][crystalIndex];
	  Double_t RadinCluster = reader->CrystalR[clusterIndex][crystalIndex];
	  Double_t RinCluster   = RadinCluster*TMath::Cos(PhiinCluster);
	  Double_t DinCluster   = RadinCluster*TMath::Sin(PhiinCluster);
	  Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex];
	  Double_t EinCluster   = reader->CrystalEnergy[clusterIndex][crystalIndex];
	  if( RadinCluster < TMath::Abs(RCenterCrystal)){
	    RCenterCrystal = RadinCluster;
	    ECenterCrystal = EinCluster;
	    TCenterCrystal = TinCluster;
	    RCenter = RinCluster;
	    DCenter = DinCluster;
	  }
	}
      }
      
      /// Fill Data to profile
      if( reader->nCrystal[clusterIndex] > 100 ){ continue; }
      for( int crystalIndex = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){	
	if( reader->CrystalEnergy[clusterIndex][crystalIndex] < LowEnergyCut ){continue;}
	Double_t RadinCluster = reader->CrystalR[clusterIndex][crystalIndex];
	Double_t PhiinCluster = reader->CrystalPhi[clusterIndex][crystalIndex];
	Double_t RinCluster   = RadinCluster*TMath::Cos(PhiinCluster);
	Double_t DinCluster   = RadinCluster*TMath::Sin(PhiinCluster);
	Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TCenterCrystal;
	Double_t EinCluster   = reader->CrystalEnergy[clusterIndex][crystalIndex];
	//if(TinCluster != TinCluster){ continue; }
	if( TMath::IsNaN(TinCluster) != 0){ continue; }
	//std::cout<< RinCluster << " : " << DinCluster << " : " << EinCluster << " : "<< TinCluster << std::endl;
	if( EinCluster == 0 ){continue;}
	if(RadinCluster > RCenterCrystal ){
	  profRDT->Fill(RinCluster,DinCluster,TinCluster);
	  profRDE->Fill(RinCluster,DinCluster,EinCluster);
	  
	  if( EinCluster > ECenterCrystal*0.8 && EinCluster < ECenterCrystal*1.2){
	    profRDT_NR->Fill(RinCluster,DinCluster,TinCluster);
	    profRDE_NR->Fill(RinCluster,DinCluster,EinCluster);
	    if( TMath::Abs(DinCluster) < RDCutWid ){
	      profRT_NR->Fill(RinCluster,TinCluster);
	    }
	    if( TMath::Abs(RinCluster) < RDCutWid ){
	      profDT_NR->Fill(DinCluster,TinCluster);
	    }
	  }

	  if( EinCluster > ECenterCrystal ){
	    profRDT_INV->Fill(RinCluster,DinCluster,TinCluster);
	    profRDE_INV->Fill(RinCluster,DinCluster,TinCluster);
	  }
	  
	  if( ECenterCrystal > reader->ClusterEnergy[crystalIndex]*0.3 ){
	    profRDT_High->Fill(RinCluster,DinCluster,TinCluster);
	    profRDE_High->Fill(RinCluster,DinCluster,EinCluster);
	  }else{
	    profRDT_Low->Fill(RinCluster,DinCluster,TinCluster);
	    profRDE_Low->Fill(RinCluster,DinCluster,EinCluster);
	  }
	  
	  if( TMath::Abs(DinCluster) < RDCutNar ){
	    profRT->Fill(RinCluster,TinCluster);
	    profRE->Fill(RinCluster,EinCluster);
	    hisRT->Fill(RinCluster,TinCluster);
	  }
	  
	  if( TMath::Abs(RinCluster) < RDCutNar ){
	    profDT->Fill(DinCluster,TinCluster);
	    profDE->Fill(DinCluster,EinCluster);
	    hisDT->Fill(DinCluster,TinCluster);
	  }
	  
	  profEET->Fill(ECenterCrystal,EinCluster,TinCluster);
	  if( TMath::Abs(RinCluster) < RDCutNar && TMath::Abs(DinCluster) < RDCutWid){
	    profEET_D->Fill(ECenterCrystal,EinCluster,TinCluster);
	    if(EinCluster < 24){
	      profET_D->Fill(ECenterCrystal,TinCluster);
	    }
	  }
	}

      } 
    }
  }



  tfout->cd();
  hisRT->Write();
  hisDT->Write();
  profET_D->Write();
  profRDT->Write();
  profRDE->Write();
  profRDE_NR->Write();
  profRDE_Low->Write();
  profRDE_High->Write();
  profRDE_INV->Write();  
  
  profRDT_NR->Write();
  profRDT_Low->Write();
  profRDT_High->Write();
  profRDT_INV->Write();  
  
  profRET->Write();
  profDET->Write();
  profEET->Write();
  profEET_D->Write();
  
  
  profRT->Write();
  profDT->Write();
  profRE->Write();
  profDE->Write();
  profRT_NR->Write();
  profDT_NR->Write();
  profRT_High->Write();
  profDT_High->Write();
  profDT_Low->Write();
  profRT_Low->Write();  
  tfout->Close();
  return 0; 
}
