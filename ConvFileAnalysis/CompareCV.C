void CompareCV(){
  
  TH1D* hisKLMass_WO_CV = new TH1D("hisKLmom_WO_CV",
				   "KLmom;KLMom[MeV];N/10MeV",
				   120,0,1200);
  TH1D* hisKLMass_W_CV = new TH1D("hisKLmom_W_CV",
				   "KLmom;KLMom[MeV];N/10MeV",
				   120,0,1200);
  TChain* ch[2];
  ch[0] = new TChain("Tree");
  ch[1] = new TChain("Tree");
  
  ch[0]->Add("reckl3588.root");
  ch[0]->Add("reckl3589.root");
  ch[0]->Add("reckl3590.root");
  ch[1]->Add("reckl3591.root");
  ch[1]->Add("reckl3595.root");
  ch[1]->Add("reckl3596.root");

  ch[0]->Project(hisKLMass_WO_CV->GetName(),"KlongMass[0]","KlongPos[0][2]>3000 && KlongPos[0][2]<50000  && CutCondition & (1+2+4+8)");
  ch[1]->Project(hisKLMass_W_CV->GetName(),"KlongMass[0]","KlongPos[0][2]>3000 && KlongPos[0][2]<50000  && CutCondition & (1+2+4+8)");


  hisKLMass_WO_CV->Draw();
  hisKLMass_W_CV->SetLineColor(2);
  hisKLMass_W_CV->Draw("same");

}
    


  
  
