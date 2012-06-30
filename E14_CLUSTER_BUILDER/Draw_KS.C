void Draw_KS(){

  TChain* ch = new TChain("Tree");
  for( int i = 0; i< 68; i++){
    ch->Add(Form("klongRootFile/ks%d.root",4162+i));
  }
  
  ch->Draw("KlongMass[0]>>his(300,0,600)","CutCondition==0");
  //ch->Draw("KlongPt[0]:KlongMass[0]>>(400,200,600,50,0,20)","(CutCondition&(1|2|4|8))==0","colz");
  gPad->SetLogz();
}
