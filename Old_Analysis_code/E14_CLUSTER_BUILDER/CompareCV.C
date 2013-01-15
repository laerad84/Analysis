void CompareCV(){
  
  TH1D* hisKLMass_WO_CV = new TH1D("hisKLmom_WO_CV",
				   "KLmom;KLMom[MeV];N/10MeV",
				   200,400,800);
  TH1D* hisKLMass_W_CV = new TH1D("hisKLmom_W_CV",
				   "KLmom;KLMom[MeV];N/10MeV",
				   200,400,800);
  TChain* ch[2];
  ch[0] = new TChain("Tree");
  ch[1] = new TChain("Tree");
  /*
  ch[0]->Add("reckl3588.root");
  ch[0]->Add("reckl3589.root");
  ch[0]->Add("reckl3590.root");
  ch[1]->Add("reckl3591.root");
  ch[1]->Add("reckl3595.root");
  ch[1]->Add("reckl3596.root");
  */

  ch[0]->Add("klongRootFile/kl4247.root");
  ch[1]->Add("klongRootFile/kl4251.root");


  ch[0]->Project(hisKLMass_WO_CV->GetName(),
		 "KlongMass[0]",
		 "TMath::Abs(KlongMass[0]-498)<10 && (CutCondition & 8) == 0");
  ch[1]->Project(hisKLMass_W_CV->GetName(),
		 "KlongMass[0]",
		 "TMath::Abs(KlongMass[0]-498)<10 && (CutCondition & 8) == 0");


  hisKLMass_WO_CV->Draw();
  hisKLMass_WO_CV->GetYaxis()->SetRangeUser(0,500);
  hisKLMass_W_CV->SetLineColor(2);
  hisKLMass_W_CV->Draw("same");
  std::cout<< hisKLMass_WO_CV->GetEntries() << std::endl;
  std::cout<< hisKLMass_W_CV->GetEntries() << std::endl;

}
    


  
  
