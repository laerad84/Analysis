#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>

void TestSpectrum(){

  TFile* tf[3]; 
  TTree* tr[3];
  char* name[3] = {"SIM","SUM","WAV"};
  
  for( int i = 0; i < 3; i++){
    tf[i] = new TFile(Form("Kl_Total_%s.root",name[i]));
    tr[i] = (TTree*)tf[i]->Get(Form("trKL"));
  }
  
  Double_t KLMass;
  Double_t KLChisq;
  Double_t KLE;
  Double_t KLPos[3];
  
  Double_t GammaE[6];
  Double_t GammaPos[6][3]; 

  for( int i = 0; i< 3; i++){
    tr[i]->SetBranchAddress("KLMass",&KLMass);
    tr[i]->SetBranchAddress("KLChisq",&KLChisq);
    tr[i]->SetBranchAddress("KLE",&KLE);
    tr[i]->SetBranchAddress("KLPos",KLPos);
    tr[i]->SetBranchAddress("GammaE",GammaE);
    tr[i]->SetBranchAddress("GammaPos",GammaPos);
  }

  TH1D* hisKlongPosZ[3];
  TH1D* hisKlongE[3];
  TH1D* hisKlongM[3];
  TH1D* hisKlongEMassTail[3];
  TH1D* hisKlongPosZMassTail[3];
  TH1D* hisGammaE[3];

  for( int i = 0; i< 3; i++){
    hisKlongPosZ[i] = new TH1D(Form("hisKlongPos%d",i),
			       Form("hisKlongPos_%s",name[i]),
			       100,3000,6000);
    hisKlongPosZMassTail[i] = new TH1D(Form("hisKlongPosTail%d",i),
				   Form("hisKlongPosTail_%s",name[i]),
				   100,3000,6000);
    hisKlongE[i] = new TH1D(Form("hisKlongE%d",i),
			       Form("hisKlongE_%s",name[i]),
			       120,0,6000);    
    hisKlongM[i] = new TH1D(Form("hisKlongM%d",i),
			    Form("hisKlongM_%s",name[i]),
			    400,400,800);
    hisKlongEMassTail[i] = new TH1D(Form("hisKlongMassTail%d",i),
				    Form("hisKlongMassTail_%s",name[i]),
				    120,0,6000);
    hisGammaE[i] = new TH1D(Form("hisGammaE%d",i),
				    Form("hishisGammaE_%s",name[i]),
				    120,0,6000);

  }

  TH1D* hisKlongMTotal = new TH1D(Form("hisKlongMTotal"),
				  Form("hisKlongMTotal"),
				  400,400,800);
  TH1D* hisKlongMSlice[125];
  for( int i = 0; i< 125; i++){
    hisKlongMSlice[i] = new TH1D(Form("hisKlongMSlice_%d",i),
				 Form("hisKlongMSlice_%d",i),
				 400,400,800);
  }


  std::cout << "LOOP" << std::endl;
  for( int iFile = 0; iFile < 3; iFile++){
    for( int ievent = 0; ievent < tr[iFile]->GetEntries(); ievent++){

      tr[iFile]->GetEntry(ievent);
      //if( ievent  >= 100000 ){ break ; } 
      
      bool fInnerGamma = false;
      Int_t nGamma = 0;
      for( int igamma = 0; igamma < 6; igamma++){
	if( TMath::Abs(GammaPos[igamma][0]) < 150 &&
	    TMath::Abs(GammaPos[igamma][1]) < 150 ){
	  fInnerGamma = true; 
	}
	if( GammaE[igamma] > 200 ){
	  nGamma++;
	}
      }
      
      if( KLChisq > 5 ){ continue; }
      if( fInnerGamma ){ continue; }
      if( nGamma < 5 ){ continue; }
      if( KLPos[2] > 5000 || KLPos[2] < 3000 ){ continue; } 
      if( KLE >5000 ){ continue; }

      hisKlongPosZ[iFile]->Fill( KLPos[2] );
      int ibinE = hisKlongE[iFile]->Fill( KLE );
      int ibinM = hisKlongM[iFile]->Fill( KLMass );
      /*
      for( int iGamma = 0; iGamma < 6; iGamma++){
	std::cout<< "Gamma" << std::endl;
	hisGammaE[iFile]->Fill(GammaE[iGamma]);
      }
      */
      if( iFile == 0 ){
	hisKlongMSlice[ibinE]->Fill( KLMass );
      }
      if( KLMass > 510 ){ 
	hisKlongEMassTail[iFile]->Fill(KLE);
	hisKlongPosZMassTail[iFile]->Fill(KLPos[2]);
      }
    }
  }

  std::cout << "Draw" << std::endl;
  for( int i = 0; i< 3; i++){
    hisKlongPosZ[i]->Scale( 1/hisKlongPosZ[i]->GetEntries());
    hisKlongE[i]->Scale( 1/hisKlongE[i]->GetEntries());
    hisKlongM[i]->Scale( 1/hisKlongE[i]->GetEntries());
    hisKlongEMassTail[i]->Scale( 1./hisKlongEMassTail[i]->GetEntries());
  }

  for( int i = 0; i < hisKlongE[0]->GetNbinsX(); i++){
    //hisKlongMSlice[i]->Scale( 1/hisKlongE[0]->GetEntries());
    hisKlongMTotal->Add(hisKlongMSlice[i],hisKlongE[0]->GetBinContent(i));
  }
  hisKlongMTotal->Scale(1./hisKlongMTotal->GetEntries());

  
  TCanvas* can = new TCanvas("can","",800,800);
  can->Divide( 2, 2);
  can->cd(1);
  hisKlongPosZ[0]->SetLineColor(1);
  hisKlongPosZ[1]->SetLineColor(2);
  hisKlongPosZ[2]->SetLineColor(3);
  hisKlongPosZ[0]->Draw();
  hisKlongPosZ[1]->Draw("same");
  hisKlongPosZ[2]->Draw("same");
  can->cd(2);
  hisKlongE[0]->SetLineColor(1);
  hisKlongE[1]->SetLineColor(2);
  hisKlongE[2]->SetLineColor(3);
  hisKlongE[0]->Draw();
  hisKlongE[1]->Draw("same");
  hisKlongE[2]->Draw("same");
  can->cd(3);
  hisKlongM[0]->SetLineColor(1);
  hisKlongM[1]->SetLineColor(2);
  hisKlongM[2]->SetLineColor(3);
  hisKlongM[0]->Draw();
  hisKlongM[1]->Draw("same");
  hisKlongM[2]->Draw("same");
  hisKlongMTotal->Draw("same");
  can->cd(4);
  /*
  hisKlongEMassTail[0]->SetLineColor(1);
  hisKlongEMassTail[1]->SetLineColor(2);
  hisKlongEMassTail[2]->SetLineColor(3);
  hisKlongEMassTail[0]->Draw();
  hisKlongEMassTail[1]->Draw("same");
  hisKlongEMassTail[2]->Draw("same");
  */
  /*
  hisGammaE[0]->SetLineColor(1);
  hisGammaE[1]->SetLineColor(2);
  hisGammaE[2]->SetLineColor(3);
  hisGammaE[0]->Draw();
  hisGammaE[1]->Draw("same");
  hisGammaE[2]->Draw("same");
  */
  hisKlongPosZMassTail[0]->SetLineColor(1);
  hisKlongPosZMassTail[1]->SetLineColor(2);
  hisKlongPosZMassTail[2]->SetLineColor(3);
  hisKlongPosZMassTail[0]->Draw();
  hisKlongPosZMassTail[1]->Draw("same");
  hisKlongPosZMassTail[2]->Draw("same");


}
