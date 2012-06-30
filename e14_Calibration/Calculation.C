#include <fstream>

void Calculation(){
  gStyle->SetOptStat("neRMI");

  TFile* tf  = new TFile("KLCalibration_Sigma5.root");
  const int nCSI = 2716;  
  TH1D*  hisCalFactor[nCSI];
  int    nCal[nCSI];
  double CalFactor[nCSI];
  double CalRMS[nCSI];

  TH1D* hisCalDistrib = new TH1D("hisCalDistrib","",200,0,2);
  TH1D* hisRMSDistrib = new TH1D("hisRMSDistrib","",100,0,0.1);
  TH1D* hisHitDistrib = new TH1D("hisHitDistrib","",100,0,100);
  TGraph* grChannel   = new TGraph();
  std::ofstream ofs("CalibrationN.txt");
  for( int i = 0; i< nCSI; i++){
    hisCalFactor[i] = (TH1D*)tf->Get(Form("his_Calibration_%04d",i));
    nCal[i]         = hisCalFactor[i]->Integral();
    CalFactor[i]    = hisCalFactor[i]->GetMean();
    CalRMS[i]       = hisCalFactor[i]->GetRMS();
    hisCalDistrib->Fill(CalFactor[i]);
    hisHitDistrib->Fill(nCal[i]);
    grChannel->SetPoint(grChannel->GetN(),i,nCal[i]);
    if( nCal[i] != 0){
      hisRMSDistrib->Fill(CalRMS[i]/CalFactor[i]);      
    }
    ofs << i << "\t" << CalFactor[i] << "\t"<< nCal[i] << std::endl;

  }
  ofs.close();

  TCanvas* can = new TCanvas("can","",1200,1200);
  can->Divide(2,2);
  can->cd(1);
  hisCalDistrib->Draw();
  can->cd(2);
  hisRMSDistrib->Draw();
  can->cd(3);
  hisHitDistrib->Draw();
  can->cd(4);
  grChannel->SetMarkerStyle(4);
  grChannel->Draw("AP");

}




