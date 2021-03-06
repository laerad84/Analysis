#include <fstream>

void Calculation_CalibrationFactor(Int_t CalibrationNumber){
  gStyle->SetOptStat("neRMI");
  Int_t RunNumber = 4231;
  TFile* tf  = new TFile(Form("CalibrationData/Calibration_%d_%d.root",RunNumber,CalibrationNumber));
  
  const int nCSI = 2716;  
  TH1D*  hisCalFactor[nCSI];
  int    nCal[nCSI];
  double CalFactor[nCSI];
  double CalRMS[nCSI];
  double preCalFactor[nCSI];
  
  TH1D* hisCalDistrib = new TH1D("hisCalDistrib","",200,0,2);
  TH1D* hisRMSDistrib = new TH1D("hisRMSDistrib","",100,0,0.1);
  TH1D* hisHitDistrib = new TH1D("hisHitDistrib","",100,0,100);
  TGraph* grChannel   = new TGraph();
  
  for( int i = 0; i< nCSI; i++){
    preCalFactor[i] = 1;
  }
  
  if( CalibrationNumber!=0){
    std::ifstream ifs(Form("Data/CalibrationFactor_%d.dat",
			   CalibrationNumber));
    if(!ifs.is_open()){ return;}
    Int_t id;
    Double_t precal;
    while( !ifs.eof() ){
      ifs >> id >> precal;
      preCalFactor[id] = precal;
    }
  }   
  Int_t nextCalNum = CalibrationNumber +1;
  std::ofstream ofs(Form("Data/CalibrationFactor_%d.dat",nextCalNum));

  std::cout<< "Loop" << std::endl;
  for( int i = 0; i< nCSI; i++){
    hisCalFactor[i] = (TH1D*)tf->Get(Form("his_Calibration_%04d",i));
    nCal[i]         = hisCalFactor[i]->Integral();
    
    if( hisCalFactor[i]->Integral() < 20 ){
      continue;
    } 
    
    hisCalFactor[i]->Fit("gaus","Q","");
    TF1* calFunction= hisCalFactor[i]->GetFunction("gaus");
    CalFactor[i]    = hisCalFactor[i]->GetMean();
    CalRMS[i]       = hisCalFactor[i]->GetRMS();
    //CalFactor[i] = calFunction->GetParameter(1);
    //CalRMS[i] = calFunction->GetParameter(2);
    hisCalDistrib->Fill(CalFactor[i]);
    hisHitDistrib->Fill(nCal[i]);
    grChannel->SetPoint(grChannel->GetN(),i,nCal[i]);
    if( nCal[i] != 0){  
      hisRMSDistrib->Fill(CalRMS[i]/CalFactor[i]);
    }
    ofs << i << "\t" << CalFactor[i]*preCalFactor[i] << std::endl;    
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




