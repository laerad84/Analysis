void TestResult(){
  
  TFile* tf1 = new TFile("TEMPLETE_OUT_HEIGHT_3pi0RunList_200_400.root");
  TFile* tf2 = new TFile("TEMPLETE_OUT_HEIGHT_3pi0RunList_400_800.root");
  /*
  */

  TPostScript* ps = new TPostScript("Template.ps",111);
  ps->NewPage();
  TCanvas* can = new TCanvas("can","",720,720*TMath::Sqrt(2)*0.9);
  Int_t xdiv =2;
  Int_t ydiv =3; 
  Int_t nDiv = xdiv*ydiv;
  can->Divide( xdiv,ydiv);

  for( int i = 0; i< xdiv*ydiv; i++){
    can->cd(i+1);
    gPad->SetGridx();
    gPad->SetGridy();
  }

  for( int ich  = 0; ich < 2716; ich++){
    can->cd( (ich % nDiv) +1 );
    TGraph* gr1 = (TGraph*)tf1->Get(Form("Waveform_Height_200_400_%d",ich));
    TGraph* gr2 = (TGraph*)tf2->Get(Form("Waveform_Height_400_800_%d",ich));    
    //TH2D*   his = (TH2D*)tf1->Get(Form("TempHist_Height_200_400_%d",ich));
    gr1->SetMarkerStyle(6);
    gr2->SetMarkerStyle(6);    
    gr1->SetMarkerColor(2);		       
    gr2->SetMarkerColor(4);    
    //his->Draw("col");
    gr1->Draw("AP");
    gr2->Draw("P");
    if((ich % nDiv) == nDiv-1){
      can->Update();
      can->Modified();
      ps->NewPage();
    }
  }
  ps->Close();
}
