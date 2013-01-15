void CompareZ(){
  
  Int_t RunNumber[2];
  RunNumber[0] = 4267;
  RunNumber[1] = 4551;
  Int_t nFiles[2]={10,10};
  TH1D*   hisPos[2];
  TFile*  tfKaon[2];
  TChain*  chKaon[2];
  
  for( int i = 0; i< 2; i++){
    chKaon[i] = new TChain("Tree");
    for( int ifile= 0; ifile < nFiles[i]; ifile++){ 
      chKaon[i]->Add(Form("klongRootFile/kl%d.root",RunNumber[i]+ifile));
    }
  }

  for( int i = 0; i< 2; i++){
    hisPos[i] = new TH1D(Form("his%d",i),Form("his%d",i),
			 50,0,6000);
    chKaon[i]->Project(hisPos[i]->GetName(),"KlongE[0]","KlongChisqZ[0]<5");
    std::cout<< hisPos[i]->GetMean() << std::endl;
    hisPos[i]->SetLineColor(i+1);
    hisPos[i]->Scale(1/(hisPos[i]->GetEntries()));
  }
  
  hisPos[0]->Draw("");
  hisPos[1]->Draw("same");
  
}
  
