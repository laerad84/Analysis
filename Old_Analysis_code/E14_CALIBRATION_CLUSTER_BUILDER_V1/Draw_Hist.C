void Draw_Hist(){
  //TFile* tf = new TFile("Calibration_Data/Calibration_With_Tight_Mass_Cut/CalhistListADV_KLRunList_2.txt_11.root");
  TFile* tf = new TFile("Calibration_Data/Calibration_NoChange/CalhistListADV_KLRunList_2.txt_11.root");
  
   //TFile* tf = new TFile("Calibration_Data/CalhistListADV_KLRunList_2.txt_1.root");
  TGraph* gr  = new TGraph();
  
  TH2D* hisEneRMS = new TH2D("his_EneRMS",
			     "Energy vs RMS;Mean Energy of CsI[MeV];RMS of Cal Factor", 
			     20,0,1000,
			     100,0,0.1);
  TH2D* hisRatioRMS = new TH2D("his_RatioRMS",
			       "Ratio vs RMS;Mean Ratio of CsI[MeV];RMS of Cal Factor", 
			       20,0,1,
			       100,0,0.1);
  TH2D* hisEneFactor = new TH2D("his_EneFactor",
				"Factor vs Energy ;Mean Energy of CsI[MeV];Cal Factor", 
				20,0,2000,
				100,0.8,1.2);
  
  TH2D* hisEne[2716];
  TH2D* hisRatio[2716];  
  TH2D* hisEneSum = new TH2D("hisSum","hisSum;[MeV];Cal Factor",
			     20,0,2000,
			     100,0,2);
  for( int i =0 ;i< 2716; ++i ){
    hisEne[i] = (TH2D*)tf->Get(Form("hisCalibrationFactorEne_%d",i));
    hisRatio[i] = (TH2D*)tf->Get(Form("hisCalibrationRatio_%d",i));    
    Double_t Mean_energy = hisEne[i]->GetMean(1);
    Double_t RMS_Factor  = hisEne[i]->GetRMS(2);
    Double_t Mean_Ratio  = hisRatio[i]->GetMean(1);
    Double_t Mean_Factor = hisEne[i]->GetMean(2);

    if( Mean_energy ==0 || RMS_Factor == 0){ continue;}
    gr->SetPoint( gr->GetN() , Mean_energy, RMS_Factor);
    hisEneRMS->Fill(Mean_energy, RMS_Factor);    
    hisRatioRMS->Fill(Mean_Ratio, RMS_Factor);
    hisEneFactor->Fill(Mean_energy, Mean_Factor);

    hisEneSum->Add(hisEne[i]);

  }


  //gr->Draw("AP");
  //TCanvas* can = new TCanvas("can","",800,800);
  /*
  hisEneRMS->Draw("colz");
  TProfile* pro = hisEneRMS->ProfileX();
  gPad->SetGridx();
  gPad->SetGridy();
  pro->Draw("same");
  */

  hisRatioRMS->Draw("colz");
  TProfile* pro = hisRatioRMS->ProfileX();
  gPad->SetGridx();
  gPad->SetGridy();
  pro->Draw("same");

  /*
  hisEneSum->Draw("colz");
  TProfile* pro = hisEneSum->ProfileX();
  pro->GetYaxis()->SetRangeUser(0.8,1.2);
  pro->Draw("same");
  gPad->SetGridx();
  gPad->SetGridy();
  */
}
