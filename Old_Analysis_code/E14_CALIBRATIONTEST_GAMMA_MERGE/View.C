void View(){
  TFile* tf = new TFile("SimCalibration_Data/Merge_All.root");
  
  TH2D* his_KlMass_Corr_All = (TH2D*)tf->Get("his_KlMass_Corr_All");
  TH2D* his_Pi0Mass_Corr_All = (TH2D*)tf->Get("his_Pi0Mass_Corr_All");  
  TH2D* his_ChisqZ_Corr_All  = (TH2D*)tf->Get("his_ChisqZ_Corr_All");
  TH2D* his_GammaE_Corr_All  = (TH2D*)tf->Get("his_GammaE_Corr_All");  
  TH2D* his_CalChisq_Corr_All= (TH2D*)tf->Get("his_CalChisq_Corr_All");
  TH2D* his_deltaPz_recPz    = (TH2D*)tf->Get("his_deltaPz_recPz");

  TCanvas* can = new TCanvas("can","",1200,800);
  can->Divide(3,2);
  can->cd(1);
  his_KlMass_Corr_All->Draw("colz");
  can->SetLogz();
  can->cd(2);
  his_Pi0Mass_Corr_All->Draw("colz");
  can->SetLogz();
  can->cd(3);
  his_ChisqZ_Corr_All->Draw("colz");
  can->SetLogz();
  can->cd(4);
  his_GammaE_Corr_All->Draw("colz");
  can->SetLogz();
  can->cd(5);
  his_deltaPz_recPz->Draw("colz");
  can->SetLogz();
  can->cd(6);
  his_CalChisq_Corr_All->Draw("colz");
  can->SetLogz();
}

  
				      
