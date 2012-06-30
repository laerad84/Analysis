void Draw_KL(){

  TChain* ch  = new TChain("Tree");
  for( int i = 3821; i < 3836; i++){
    ch->Add(Form("klongRootFile/kl%04d.root",i));
  }  


  TCanvas* can = new TCanvas("can","",800,800); 
  //ch->Draw("KlongMass[0]","(CutCondition & (1+2+3+4) )==0");  
  ch->Draw("KlongMass[0]","(CutCondition & (1+2+3+4) )==0&& KlongMass[0] >480 && KlongMass[0] < 520");  
 
}
