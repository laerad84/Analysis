void Test(){
  
  TFile* tf = new TFile("SimCalibration_Data/Merge_All.root");

  TH1D* his_All[2716];
  TH1D* his_Fit[2716];
  TH1D* his_Mis[2716];
  for( int i =0 ;i< 2716; ++i){
    his_All[i] = (TH1D*)tf->Get(Form("hisCalib_All_%04d",i)); 
    his_Fit[i] = (TH1D*)tf->Get(Form("hisCalib_Fit_%04d",i));
    his_Mis[i] = (TH1D*)tf->Get(Form("hisCalib_Mis_%04d",i));
  }

  TH1D* hisRatio = new TH1D("hisPurity","purity",50,0,0.25);
  for( int i = 0; i< 2716; ++i){
    if( his_All[i]->GetEntries()>50){
      hisRatio->Fill(his_Mis[i]->GetEntries() / his_All[i]->GetEntries());
    }
  }
  
  hisRatio->Draw();
  gPad->SetGridx();
  gPad->SetGridy();
}
      
  
