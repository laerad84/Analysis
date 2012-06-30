#include <fstream>

void ConfirmCalibrationFactor(){
  gStyle->SetOptStat("neRMI");
  gStyle->SetPalette(1);
  TFile* tf  = new TFile("KLCalibration_Sigma5.root");
  const int nCSI = 2716;  
  TH1D*  hisCalFactor[nCSI];
  int    nCal[nCSI];
  double CalFactor[nCSI];
  double CalRMS[nCSI];
  
  TGraph* grChannel   = new TGraph();
  TH2D* hisCalConfirm  = new TH2D("hisCalFactor","",80,0.8,1.2,80,0.8,1.2);  
  std::ofstream ofs("CalibrationN.txt");
  TH1D* hisAspect     = new TH1D("hisAspect","",400,0,2);
  
  for( int i = 0; i< nCSI; i++){

    hisCalFactor[i] = (TH1D*)tf->Get(Form("his_Calibration_%04d",i));
    nCal[i]         = hisCalFactor[i]->Integral();
    CalFactor[i]    = hisCalFactor[i]->GetMean();
    CalRMS[i]       = hisCalFactor[i]->GetRMS();
    grChannel->SetPoint(grChannel->GetN(),0.95+0.005*(i%10),CalFactor[i]);
    hisCalConfirm->Fill(1/(0.95+0.005*(i%10)),CalFactor[i]);
    hisAspect->Fill(CalFactor[i]*(0.95+0.005*(i%10))/0.975);
    ofs << i << "\t" << CalFactor[i] << "\t"<< nCal[i] << std::endl;

  }
  ofs.close();

  TCanvas* can = new TCanvas("can","",600,600);
  //hisAspect->Draw();
  hisCalConfirm->Draw("colz");
  TProfile* pro = hisCalConfirm->ProfileX();
  pro->SetLineColor(2);
  pro->Draw("same");

}



