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
#include "TPDF.h"
#include "TText.h"

#include "IDHandler.h"
#include "CsIImage.h"
#include "E14ReadSumFile.h"
#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include "gamma/GammaFinder.h"



int
main( int argc, char** argv){

  std::cout << "KLong Monitor" << std::endl;
  int RunNumber;
  std::string RunName="Klong";

  if( argc >= 2){
    RunNumber  = atoi(argv[1]);
  }else{
    std::cerr << "Arguemnet Errors"<<std::endl;
    return -1;
  }
  if( argc ==3 ){
    RunName = argv[2];
  }
  TFile* tf = new TFile("RunTotal.root","Update");

  std::cout<<  __LINE__ << std::endl;  
  gStyle->SetOptStat("neMRIuo");
  gStyle->SetOptFit(111111111);
  
  IDHandler* handler = new IDHandler();  
  std::cout<< __LINE__ << std::endl;
  
  CsIImage*  image[3];
  image[0] = new CsIImage(handler);
  image[1] = new CsIImage(handler);
  image[2] = new CsIImage(handler);
  image[0]->SetTitle(Form("%s_%d",RunName.c_str(),RunNumber));
  image[1]->SetTitle(Form("%sMassZCut_%d",RunName.c_str(),RunNumber));
  image[2]->SetTitle(Form("%sMassZEneCut_%d",RunName.c_str(),RunNumber));
  
  TApplication* app = new TApplication("app",&argc, argv);
  
  TH1D* hisKLMass[3];
  TH1D* hisKLMom[3];
  TH1D* hisKLPt[3];
  TH1D* hisKLPos[3];
  TH1D* hisKLCoeX[3];
  TH1D* hisKLCoeY[3];
  TH2D* hisKLCoe2D[3];
  TH1D* hisClusPos[3];

  TText* txtEvent[3];
  TText* txtCoeX[3];
  TText* txtCoeY[3];
  
  TH1D* hisTrig;
  TH2D* hisrlBalance;
  
  char* CutName[3] = {"NoCut","MassZ","MassZEne"};
  std::cout<< __LINE__ << std::endl;
  
  TChain* ch;
  ch = new TChain("Tree");
  ch->Add(Form("klongRootFile/kl%04d.root",RunNumber));
  std::cout<< ch->GetEntries() << std::endl;
  for( int i = 0; i< 3; i++){
    txtEvent[i] = new TText(0,0,"");
    txtCoeX[i]  = new TText(0,0,"");
    txtCoeY[i]  = new TText(0,0,"");

    hisKLMass[i]
      = new TH1D(Form("hisKLMass%s_%04d",CutName[i],RunNumber),
		 Form("hisKLMass%s_%04d;Mass[MeV];N/5MeV",CutName[i],RunNumber),
		 300,0,1500);
    hisKLMom[i]
      = new TH1D(Form("hisKLMom%s_%04d",CutName[i],RunNumber),
		 Form("hisKLMom%s_%04d;Mom[MeV];N/50MeV",CutName[i],RunNumber),
		 160,0,8000);  
    hisKLPos[i]
      = new TH1D(Form("hisKLPos%s_%04d",CutName[i],RunNumber),
		 Form("hisKLPos%s_%04d;Pos[mm];N/50mm",CutName[i],RunNumber),
		 160,0,8000);  
    hisKLPt[i]
      = new TH1D(Form("hisKLPt%s_%d",CutName[i],RunNumber),
		 Form("hisKLPt%s_%d",CutName[i],RunNumber),
		 300,0,300);  
    hisKLMass[i]->SetLineColor(i+1);
    hisKLPt[i]  ->SetLineColor(i+1);    
    hisKLMom[i] ->SetLineColor(i+1);
    hisKLPos[i] ->SetLineColor(i+1);
    
    hisKLCoeX[i]
      = new TH1D(Form("hisKLCoeX%s_%d",CutName[i],RunNumber),
		 Form("hisKLCoeX%s_%d",CutName[i],RunNumber),
		 200,-200,200);        
    hisKLCoeY[i]
      = new TH1D(Form("hisKLCoeY%s_%d",CutName[i],RunNumber),
		 Form("hisKLCoeY%s_%d",CutName[i],RunNumber),
		 200,-200,200);
    hisKLCoe2D[i]
      = new TH2D(Form("hisKLCoe2D%s_%d",CutName[i],RunNumber),
		 Form("hisKLCoe2D%s_%d",CutName[i],RunNumber),
		 200,-100,100,200,-100,100);
    hisClusPos[i]
      = new TH1D(Form("hisClusPos%s_%d",CutName[i],RunNumber),
		 Form("hisClusPos%s_%d",CutName[i],RunNumber),
		 80,-1000,1000);
    hisKLCoeX[i]->SetLineColor(i+1);
    hisKLCoeY[i]->SetLineColor(i+1);
  }
  
  hisTrig      = new TH1D(Form("hisTrig_%d",RunNumber),"Trig",4,0,4);
  hisrlBalance = new TH2D(Form("hisrlBalance%d",RunNumber),"TrigBalance",
			  40,-1,1,60,0,6000);
  TF1* func    = new TF1("BalanceFunction","400/(1-TMath::Abs(x))",-1,1);
  Int_t eventID;
  std::cout << __LINE__ << std::endl;
  E14GNAnaDataContainer data;
  data.setBranchAddress(ch);    
  //ch->SetBranchAddress("eventID",&eventID);
  long nentries = ch->GetEntries();
  std::cout << __LINE__ << std::endl;
  for( int ievent = 0; ievent< nentries; ievent++){      
    ch->GetEntry(ievent);      
    Klong kl;
    Double_t klcoex=0;
    Double_t klcoey=0;

    // Data Analysis 
    data.getData(kl);
    for( int pi0Index  = 0; pi0Index < kl.pi0().size(); pi0Index++){
      klcoex += kl.pi0()[pi0Index].g1().x() * kl.pi0()[pi0Index].g1().e();
      klcoex += kl.pi0()[pi0Index].g2().x() * kl.pi0()[pi0Index].g2().e();
      klcoey += kl.pi0()[pi0Index].g1().y() * kl.pi0()[pi0Index].g1().e();
      klcoey += kl.pi0()[pi0Index].g2().y() * kl.pi0()[pi0Index].g2().e();
    }
    klcoex = klcoex/kl.e();
    klcoey = klcoey/kl.e();
    hisKLCoeX[0]->Fill(klcoex);
    hisKLCoeY[0]->Fill(klcoey);
    hisKLCoe2D[0]->Fill(klcoex,klcoey);
    if(kl.e()>1000){
      hisKLMass[0]->Fill(kl.m());
    }
    hisKLMom[0]->Fill(kl.e());
    hisKLPt[0]->Fill(kl.p3().perp());
    hisKLPos[0]->Fill(kl.vz());

    if(( data.CutCondition & 8 ) == 0 && TMath::Abs( kl.m()-497.614) < 10){
      hisKLMass[1]->Fill(kl.m());
      if(kl.e() >1000){
	hisKLMom[1]->Fill(kl.e());
      }
      hisKLCoeX[1]->Fill(klcoex);
      hisKLCoeY[1]->Fill(klcoey);
      hisKLCoe2D[1]->Fill(klcoex,klcoey);
      hisKLPt[1]->Fill(kl.p3().perp());
      hisKLPos[1]->Fill(kl.vz());
      if(( data.CutCondition & 1 ) == 0 ){
	if(kl.e() > 1000){
	  hisKLMass[2]->Fill(kl.m());
	}
	hisKLMom[2]->Fill(kl.e());
	hisKLCoeX[2]->Fill(klcoex);
	hisKLCoeY[2]->Fill(klcoey);
	hisKLCoe2D[2]->Fill(klcoex,klcoey);
	hisKLPt[2]->Fill(kl.p3().perp());
	hisKLPos[2]->Fill(kl.vz());
      }
      
    }
    
    /*
      if(( data.CutCondition & 1 )!=0){
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
	 iCluster < data.GamClusNumber ;
	 iCluster++){	
      image[0]->Fill(data.GamClusCsiId[iCluster][0]);
      if(( data.CutCondition & 8 )== 0 && TMath::Abs(kl.m()-497.614) < 10 ){
	image[1]->Fill(data.GamClusCsiId[iCluster][0]);
	if((data.CutCondition & 1 ) == 0 ){
	  image[2]->Fill(data.GamClusCsiId[iCluster][0]);
	}
      }
      if( data.GamClusCoePos[iCluster][0] > 0){	    	  
	EnergyLeft+=data.GamClusDepE[iCluster];
      }else{
	EnergyRight+=data.GamClusDepE[iCluster];
      }
      
      if( data.GamClusDepE[iCluster] > 200 ){
	if( data.GamClusCoePos[iCluster][0] > 0){	    	  
	  leftHit++;
	}else{
	  rightHit++;
	}
      }
      hisClusPos[0]->Fill(data.GamClusCoePos[iCluster][0]);
      if(( data.CutCondition & 8 ) == 0 && TMath::Abs( kl.m()-497.614) < 10){
	hisClusPos[1]->Fill(data.GamClusCoePos[iCluster][0]);
	if(( data.CutCondition & 1)==0){
	  hisClusPos[2]->Fill(data.GamClusCoePos[iCluster][0]);
	}
      }
    }
    if(leftHit>0){
      hit |= 1;
    }
    if(rightHit>0){
      hit |= 2;
    }
    
    EnergyTotal = EnergyRight + EnergyLeft;
    hisrlBalance->Fill((EnergyLeft -EnergyRight)/EnergyTotal,EnergyTotal); 
    hisTrig->Fill(hit);	
  }
  
  TCanvas* can = new TCanvas("can","",1600,800);
  can->Divide(4,2);  
  TPDF* pdf = new TPDF(Form("Image/KlongDraw%d.pdf",RunNumber),111);
  can->cd(1);
  image[1]->Draw();
  can->cd(2);
  hisrlBalance->Draw("colz");
  func->Draw("same");
  TText* txtBalance = new TText(0,0,"");
  txtBalance->DrawTextNDC(0.5,0.2,
			  Form("RLBalance:%.3f",hisrlBalance->GetMean(1)));
  can->cd(3);
  hisKLPos[0]->Draw();  
  hisKLPos[1]->Draw("same");  
  hisKLPos[2]->Draw("same");  
  can->cd(4);
  //gPad->SetLogy();
  hisKLMass[0]->Draw();
  hisKLMass[1]->Draw("same");
  hisKLMass[2]->Draw("same");
  for( int i = 0; i< 3; i++){
    txtEvent[i]->SetTextSize(0.04);
    txtEvent[i]->SetTextColor(i+1);
    if( i ==2 )txtEvent[i]->SetTextColor(i+2);
    txtEvent[i]->DrawTextNDC(0.5,0.6-0.05*i,Form("#of Klong%s:%d",CutName[i],(int)hisKLMass[i]->GetEntries()));
  }

  can->cd(5);
  gPad->SetLogy();
  hisKLMom[0]->Draw();
  hisKLMom[1]->Draw("same");
  hisKLMom[2]->Draw("same");

  can->cd(6);  
  gPad->SetLogy();
  hisKLPt[0]->Draw();
  hisKLPt[1]->Draw("same");
  hisKLPt[2]->Draw("same");
  
  can->cd(7);
  gPad->SetLogy();
  hisKLCoeX[0]->Draw();
  hisKLCoeX[1]->Draw("same");
  hisKLCoeX[2]->Draw("same");
  for( int i = 0; i< 3; i++){
    txtCoeX[i]->SetTextSize(0.04);
    txtCoeX[i]->SetTextColor(i+1);
    if( i ==2 )txtCoeX[i]->SetTextColor(i+2);
    txtCoeX[i]->DrawTextNDC(0.5,0.6-0.05*i,Form("#of Klong%s:%.2f",CutName[i],hisKLCoeX[i]->GetMean()));
  }
  can->cd(8);
  gPad->SetLogy();
  hisKLCoeY[0]->Draw();
  hisKLCoeY[1]->Draw("same");
  hisKLCoeY[2]->Draw("same");
  for( int i = 0; i< 3; i++){
    txtCoeY[i]->SetTextSize(0.04);
    txtCoeY[i]->SetTextColor(i+1);
    if( i ==2 )txtCoeY[i]->SetTextColor(i+2);
    txtCoeY[i]->DrawTextNDC(0.5,0.6-0.05*i,Form("#of Klong%s:%.2f",CutName[i],hisKLCoeY[i]->GetMean()));
  }
  tf->Cd("");
  for( int i = 0; i<3; i++){
    hisKLPt[i]->Write();
    hisKLCoeX[i]->Write();
    hisKLCoeY[i]->Write();
    hisKLMass[i]->Write();
    hisKLMom[i]->Write();
    hisClusPos[i]->Write();
  }
  can->Update();
  pdf->Close();
  tf->Close();   
  //can->SaveAs("Image/KlongTriggerstudy.gif");
  //app->Run();
}


