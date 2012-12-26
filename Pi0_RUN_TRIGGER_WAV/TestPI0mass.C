void TestPI0mass(){
  TGraphErrors* gr = new TGraphErrors();

  const Int_t nFiles=10;
  Int_t RunNumber=4503;
  Double_t pi0mass[nFiles]={0};
  TFile* tf[nFiles];
  TTree* tr[nFiles];
  TFile* tfout = new TFile("Pi0noCal.root","recreate");
  TH1D* hisMass[nFiles];
  TH1D* hisMassAll = new TH1D("MassAll","Mass",600,0,600);

  for( int i = 0; i< nFiles; i++){
    tf[i] = new TFile(Form("Pi0_data/Pi0_%d.root",RunNumber+i));
    tr[i] = (TTree*)tf[i]->Get("Tree");
    hisMass[i] = new TH1D(Form("Mass%d",RunNumber+i),Form("Mass%d",RunNumber+i),
			  600,0,600);
    tr[i]->Project(hisMass[i]->GetName(),
		   "Pi0Mass[0]",
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
  hisMassAll->Fit("gaus","","",130,140);

  for( int i = 0; i< nFiles; i++){
    hisMass[i]->Draw("same");
  }


  can->cd(2);
  gr->SetMarkerStyle(4);
  gr->Draw("AP");

  /*
  tfout->cd();
  hisMassAll->Write();
  tfout->Close();
  */
}
  
  
		      
  
