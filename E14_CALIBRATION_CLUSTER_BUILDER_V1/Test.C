void Test(){
  const int nit = 4;
  TFile* tf[nit];
  TTree* tr[nit];
  TH1D* his[nit];
  for( int i = 0; i< nit; i++){
    tf[i] = new TFile(Form("CalibrationData/Calibration_4231_%d.root",i));
    tr[i] = (TTree*)tf[i]->Get("trCalibration");
    his[i] = new TH1D(Form("HisCalibration%d",i),
		      Form("HisCalibration%d;Mass[MeV];N/5MeV",i),80,480,520);
    tr[i]->Project(his[i]->GetName(),
		   "KlongMass[0]",
		   "KlongChisqZ[0]<10");  
  }
  TCanvas* can = new TCanvas("can","",800,800);
  his[nit-1]->Draw();
  for( int i = 0 ; i< nit; i++){
    his[i]->SetLineColor(i+1);
    his[i]->Draw("same");
    std::cout <<his[i]->GetEntries() << " " << his[i]->GetRMS() << std::endl;
  }

  //  his[0]->Fit("gaus","","",480,520);
  //his[1]->Fit("gaus","","",480,520);
}

  
  
