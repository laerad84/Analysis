void TestProf(){
  
  TFile* tf = new TFile("SimCalibration_Data/nonFixing/Merge_All.root");
  //TFile* tf = new TFile("SimCalibration_Data/Merge_All.root");
  

  TH2D* his_GammaE_DeltaE = (TH2D*)tf->Get("his_GammaE_DeltaE");
  TGraph* gr = new TGraph();
  
  TProfile* prof = his_GammaE_DeltaE->ProfileX();
  TF1* func0 = new TF1("func0","[0]+[1]*x",0,4000);
  prof->Fit(func0,"","",400,2000);
  
  TF1* func = new TF1("func","([0]+[1]*x)*(1+[2]*TMath::Exp([3]*x))",0,4000);
  func->SetParameter(0,func0->GetParameter(0));
  func->SetParameter(1,func0->GetParameter(1));
  func->SetParLimits(0,func0->GetParameter(0),func0->GetParameter(0));
  func->SetParLimits(1,func0->GetParameter(1),func0->GetParameter(1));
  func->SetParLimits(2,-0.2,0);
  func->SetParLimits(3,-1,0);
  

  prof->Fit(func,"","",100,2000);

  for( int i = 0; i< prof->GetNbinsX() ; i++){
    if( prof->GetBinContent(i+1) == 0) continue;
    Double_t y = prof->GetBinContent(i+1);
    Double_t x = prof->GetBinCenter(i+1); 
    gr->SetPoint( gr->GetN(), x, y); 
  }
  
  //gr->SetMarkerStyle(5);
  func->Draw();
  gr->Draw("p");
  TSpline3* spl = new TSpline3("spl",gr);
  //spl->Draw("same");
  prof->Draw("same");
}
