void testStability(){
  TChain* ch = new TChain("laserEventData");
  
  for( int i = 3644; i< 3714; i++){
    ch->Add(Form("laserdata_%04d.root",i));
  }
  ch->Draw("PINData[3]:T+3600*15","","colz");
}

  
