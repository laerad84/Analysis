void TestGammDistribution(){
  TGraphErrors* gr = new TGraphErrors();

  const Int_t nFiles=23;
  Int_t RunNumber=4503;
  Double_t pi0mass[nFiles]={0};
  TFile* tf[nFiles];
  TTree* tr[nFiles];
  TH1D* hisMass[nFiles];
  TH1D* hisMassAll = new TH1D("MassAll","Mass",600,0,600);
  TH2D* hisGammaDist = new TH2D("GammaDistribution","",
				80,-1000,1000,80,-1000,1000);
  TChain* ch = new TChain("Tree");
  for( int i = 0; i< nFiles; i++){

    tf[i] = new TFile(Form("pi%d.root",RunNumber+i));
    ch->Add(tf[i]->GetName());
    tr[i] = (TTree*)tf[i]->Get("Tree");
    hisMass[i] = new TH1D(Form("Mass%d",RunNumber+i),Form("Mass%d",RunNumber+i),
			  600,0,600);
    
    tr[i]->Project(hisMass[i]->GetName(),
		   "Pi0Mass[0]",
		   "GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5&&(TMath::Abs(GammaPos[0][0])>150||TMath::Abs(GammaPos[0][1])>150)&&(TMath::Abs(GammaPos[1][0])>150||TMath::Abs(GammaPos[1][1])>150)");
    
    tr[i]->Project(hisGammaDist->GetName(),
		   "GammaPos[][1]:GammaPos[][0]",
		   "GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5");

    std::cout<< hisMass[i]->GetEntries() << std::endl;
    hisMass[i]->Fit("gaus","Q","",125,140);
    TF1* gaus = hisMass[i]->GetFunction("gaus");

    std::cout<< i << " : " 
	     << hisMass[i]->GetEntries() << " : "
	     << gaus->GetParameter(1) << " : "
	     << gaus->GetParameter(2) << std::endl;
    pi0mass[i] = gaus->GetParameter(1);
    gr->SetPoint(gr->GetN(),i+RunNumber, gaus->GetParameter(1)/pi0mass[0]);
    gr->SetPointError(gr->GetN()-1,0, gaus->GetParError(1)/pi0mass[0]);
    hisMassAll->Add(hisMass[i]);
  }
  
  for( int i = 0; i< nFiles; i++){
    hisMass[i]->SetLineColor(i+1);
    
  }
  TCanvas* can = new TCanvas("can","",1000,500);
  can->Divide(2,1);
  can->cd(1);
  hisMassAll->Draw();
  /*
  for( int i = 0; i< nFiles; i++){
    hisMass[i]->Draw("same");
  }
  */
  can->cd(2);
  //gr->SetMarkerStyle(4);
  //gr->Draw("AP");
  //hisGammaDist->Draw("colz");
  ch->Draw("GammaPos[][1]:GammaPos[][0]>>hist(80,-1000,1000,80,-1000,1000)",
	   "GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5&&(TMath::Abs(GammaPos[0][0])>150||TMath::Abs(GammaPos[0][1])>150)&&(TMath::Abs(GammaPos[1][0])>150||TMath::Abs(GammaPos[1][1])>150)","colz");
	   
}
  
  
		      
  
