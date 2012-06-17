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
#include "TProfile.h"
#include "TMath.h"

#include "IDHandler.h"
#include "CsIImage.h"
#include "E14ReadSumFile.h"
#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include "gamma/GammaFinder.h"



int
main( int argc, char** argv){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  std::cout << "KLong Monitor" << std::endl;
  int RunNumber;
  int RunNumberFinal;
  std::string RunName="Klong";
  std::string RunListFile;
  std::vector<int> RunList;
  RunListFile = argv[1];
  std::ifstream ifs(RunListFile.c_str());
  if( !ifs.is_open()){return -1;}
  while( ifs >> RunNumber){
    RunList.push_back(RunNumber);
  }
  RunNumberFinal = RunList[RunList.size()-1];

  std::cout<< "Number Of Run: " << RunList.size() << std::endl;

  /*
  if( argc ==3 ){
    RunNumber  = atoi(argv[1]);
    RunNumberFinal = atoi(argv[2]);    
  }else{
    std::cerr << "Arguemnet Errors"<<std::endl;
    return -1;
  }
  */


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
  TH2D* hisKLRecPosXZ[3];
  TH2D* hisKLRecPosYZ[3];
  TH2D* hisKLRecPosXY[3];
  TH2D* hisGammaSpectrum2D[3];
  TH1D* hisGammaSpectrum[3];
  TH1D* hisKL = new TH1D("MassDistrib","mass",500,0,1000);

  TText* txtEvent[3];
  TText* txtCoeX[3];
  TText* txtCoeY[3];
  
  TH1D* hisTrig;
  TH2D* hisrlBalance;
  
  char* CutName[3] = {"NoCut","MiddleCut","FullCut"};
  std::cout<< __LINE__ << std::endl;
  
  TChain* ch;
  ch = new TChain("trCalibration");
  /*
  ch = new TChain("Tree");
  for( int i = RunNumber; i<= RunNumberFinal; i++){
    ch->Add(Form("klongRootFile/kl%04d.root",i));
  }
  */

  for( std::vector<int>::iterator it = RunList.begin();
       it != RunList.end();
       ++it){
    ch->Add(Form("klongRootFile/Calibration_%04d_8.root",*it));
  }
  

  //RunNumber = 99999;
  std::cout<< ch->GetEntries() << std::endl;
  for( int i = 0; i< 3; i++){
    txtEvent[i] = new TText(0,0,"");
    txtCoeX[i]  = new TText(0,0,"");
    txtCoeY[i]  = new TText(0,0,"");
    hisGammaSpectrum[i]
      = new TH1D(Form("hisGamma%s_%04d_%04d",
		      CutName[i],RunNumber, RunNumberFinal),
		 Form("hisGamma%s_%04d_%04d;Energy[MeV];N/5MeV",
		      CutName[i],RunNumber, RunNumberFinal),
		 400,0,2000);
    hisGammaSpectrum2D[i]
      = new TH2D(Form("hisGamma2D%s_%04d_%04d",
		      CutName[i],RunNumber, RunNumberFinal),
		 Form("hisGamma2D%s_%04d_%04d;Radious[mm];Energy[Mev]",
		      CutName[i],RunNumber, RunNumberFinal),
		 40,0,1000,400,0,2000);
		 
    hisKLMass[i]
      = new TH1D(Form("hisKLMass%s_%04d_%04d",
		      CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLMass%s_%04d_%04d;Mass[MeV];N/5MeV",
		      CutName[i],RunNumber,RunNumberFinal),
		 400,400,800);
    hisKLMom[i]
      = new TH1D(Form("hisKLMom%s_%04d_%04d",
		      CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLMom%s_%04d_%04d;Mom[MeV];N/50MeV",
		      CutName[i],RunNumber,RunNumberFinal),
		 160,0,8000);  
    hisKLPos[i]
      = new TH1D(Form("hisKLPos%s_%04d_%04d",
		      CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLPos%s_%04d_%04d;Pos[mm];N/50mm",
		      CutName[i],RunNumber,RunNumberFinal),
		 160,0,8000);  
    hisKLRecPosXY[i]
      = new TH2D(Form("hisKLPosXY%s_%04d_%04d",
		      CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLPosXY%s_%04d_%04d;X[mm];Y[mm]",
		      CutName[i],RunNumber,RunNumberFinal),
		 400,-200,200,400,-200,200);  
    hisKLRecPosXZ[i]
      = new TH2D(Form("hisKLPosXZ%s_%04d_%04d",
		      CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLPosXZ%s_%04d_%04d;Z[mm];X[mm]",
		      CutName[i],RunNumber,RunNumberFinal),
		 160,0,8000,400,-200,200);  
    hisKLRecPosYZ[i]
      = new TH2D(Form("hisKLPosYZ%s_%04d_%04d",
		      CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLPosYZ%s_%04d_%04d;Z[mm];Y[mm]",
		      CutName[i],RunNumber,RunNumberFinal),
		 160,0,8000,400,-200,200);  
    hisKLPt[i]
      = new TH1D(Form("hisKLPt%s_%d_%04d",CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLPt%s_%d_%04d",CutName[i],RunNumber,RunNumberFinal),
		 300,0,300);  
    hisKLMass[i]->SetLineColor(i+1);
    hisKLPt[i]  ->SetLineColor(i+1);    
    hisKLMom[i] ->SetLineColor(i+1);
    hisKLPos[i] ->SetLineColor(i+1);
    
    hisKLCoeX[i]
      = new TH1D(Form("hisKLCoeX%s_%d_%d",CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLCoeX%s_%d_%d",CutName[i],RunNumber,RunNumberFinal),
		 200,-200,200);        
    hisKLCoeY[i]
      = new TH1D(Form("hisKLCoeY%s_%d_%d",CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLCoeY%s_%d_%d",CutName[i],RunNumber,RunNumberFinal),
		 200,-200,200);
    hisKLCoe2D[i]
      = new TH2D(Form("hisKLCoe2D%s_%d_%d",CutName[i],RunNumber,RunNumberFinal),
		 Form("hisKLCoe2D%s_%d_%d",CutName[i],RunNumber,RunNumberFinal),
		 200,-100,100,200,-100,100);
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
    hisKLMass[0]->Fill(kl.m());
    hisKLMom[0]->Fill(kl.e());
    hisKLPt[0]->Fill(kl.p3().perp());
    hisKLPos[0]->Fill(kl.vz());
    hisKLRecPosXY[0]->Fill(kl.vx(),kl.vy());
    hisKLRecPosYZ[0]->Fill(kl.vz(),kl.vy());
    hisKLRecPosXZ[0]->Fill(kl.vz(),kl.vx());

    if(( data.CutCondition & (8) ) == 0 && TMath::Abs( kl.m()-497.614) < 10){
    //if( kl.chisqZ()<10){
      hisKLRecPosXY[1]->Fill(kl.vx(),kl.vy());
      hisKLRecPosYZ[1]->Fill(kl.vz(),kl.vy());
      hisKLRecPosXZ[1]->Fill(kl.vz(),kl.vx());
      hisKLMass[1]->Fill(kl.m());
      hisKLMom[1]->Fill(kl.e());
      hisKLCoeX[1]->Fill(klcoex);
      hisKLCoeY[1]->Fill(klcoey);
      hisKLCoe2D[1]->Fill(klcoex,klcoey);
      hisKLPt[1]->Fill(kl.p3().perp());
      hisKLPos[1]->Fill(kl.vz());
      //if(( data.CutCondition & (1|2|4) ) == 0 && (kl.chisqZ() < 10)){
      if(( data.CutCondition & (1|4|8)) == 0){
	hisKL->Fill(kl.m());
      }
      if(( data.CutCondition & (1|2|4|8))==0 && TMath::Abs( kl.m()-497.614) < 10){
	hisKLRecPosXY[2]->Fill(kl.vx(),kl.vy());
	hisKLRecPosYZ[2]->Fill(kl.vz(),kl.vy());
	hisKLRecPosXZ[2]->Fill(kl.vz(),kl.vx());
	hisKLMass[2]->Fill(kl.m());
	hisKLMom[2]->Fill(kl.e());
	hisKLCoeX[2]->Fill(klcoex);
	hisKLCoeY[2]->Fill(klcoey);
	hisKLCoe2D[2]->Fill(klcoex,klcoey);
	hisKLPt[2]->Fill(kl.p3().perp());
	hisKLPos[2]->Fill(kl.vz());	
	Double_t x;
	Double_t y;
	Double_t rad;
	Double_t energy;

	for( int pi0Index  = 0; pi0Index < kl.pi0().size(); pi0Index++){
	  energy = kl.pi0()[pi0Index].g1().e();
	  x      = kl.pi0()[pi0Index].g1().x();
	  y      = kl.pi0()[pi0Index].g1().y();
	  rad    = TMath::Sqrt(x*x + y*y);
	  hisGammaSpectrum[2]->Fill(energy);
	  hisGammaSpectrum2D[2]->Fill(rad,energy);
	  energy = kl.pi0()[pi0Index].g2().e();
	  x      = kl.pi0()[pi0Index].g2().x();
	  y      = kl.pi0()[pi0Index].g2().y();
	  rad    = TMath::Sqrt(x*x + y*y);	  
	  hisGammaSpectrum[2]->Fill(energy);
	  hisGammaSpectrum2D[2]->Fill(rad,energy);
	}

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
  /*
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
  gPad->SetLogy();
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
  
  //can->SaveAs("Image/KlongTriggerstudy.gif");
  pdf->Close();
*/ 
 /*
  TCanvas* canvas = new TCanvas("canvas","",1600,800);
  canvas->Divide(4,2);
  canvas->cd(1);
  hisKLRecPosXY[0]->Draw("colz");
  canvas->cd(2);
  hisKLRecPosXY[2]->Draw("colz");
  canvas->cd(3);
  hisKLRecPosXZ[0]->Draw("colz");
  canvas->cd(4);
  hisKLRecPosXZ[2]->Draw("colz");
  TProfile* proxz = (TProfile*)hisKLRecPosXZ[2]->ProfileX();
  canvas->cd(5);
  hisKLRecPosYZ[0]->Draw("colz");
  canvas->cd(6);
  hisKLRecPosYZ[2]->Draw("colz");
  TProfile* proyz = (TProfile*)hisKLRecPosYZ[2]->ProfileX();
  canvas->cd(7);
  proxz->Draw();
  proxz->Fit("pol1","","",3000,5000);
  canvas->cd(8);
  proyz->Draw();
  proyz->Fit("pol1","","",3000,5000);
*/

  TCanvas* canSpec = new TCanvas("canSpec","",1600,800);
  canSpec->Divide(4,2);
  canSpec->cd(1);
  hisGammaSpectrum[2]->Draw();
  canSpec->cd(2);
  hisGammaSpectrum2D[2]->Draw("colz");
  canSpec->cd(3);
  hisKL->Draw();
  app->Run();
}


