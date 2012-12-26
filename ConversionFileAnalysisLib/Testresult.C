void Testresult(){
  TFile* tf  =new TFile("TEMPLETE_OUT_COSMIC.root");
  TGraphErrors* gr[2716];
  for( int i = 0; i< 2716; i++){
    gr[i] = (TGraphErrors*)tf->Get(Form("Waveform_Cosmic_%d",i));
  }
  TH1D* hisTail = new TH1D("hisTail","",400,0,1);
  TGraph* testgr = new TGraph(); 
  TCanvas* can = new TCanvas("can","",0,0,1800,900);
  gPad->SetGridx();
  gPad->SetGridy();
  //gr[0]->Draw("AP");
  for( int i = 0; i< 2716; i++){
    if( gr[i]->GetN() < 100){
      continue;
    }

    //gr[i]->SetMarkerStyle(7);
    //gr[i]->Draw("P");
    TSpline3* spl = new TSpline3(Form("spl%d",i),(TGraph*)gr[i]);
    hisTail->Fill(spl->Eval(30));
    testgr->SetPoint( testgr->GetN(), i, spl->Eval(30));
  }
  //hisTail->Draw();
  testgr->Draw("AP");
}

