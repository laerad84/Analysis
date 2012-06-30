#include <fstream>
void Draw(){
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("NeRMi");
  gSystem->Load("../AnalysisLib/lib/libAnalysisLib.so");
  
  IDHandler* handler = new IDHandler();
  CsIImage* image    = new CsIImage(handler);
  TH1D* hisResidual = new TH1D("his","",100,0.8,1.2);
  
 std::ifstream ifs("Calibration_Data/Calibration_With_CC03_Fiducial_CUT/CalibrationFactorADV_10.dat");
 //std::ifstream ifs("Calibration_Data/Calibration_With_Tight_Mass_Cut/CalibrationFactorADV_10.dat");
  //std::ifstream ifs("Calibration_Data/OldCalibration/CalibrationFactor_10.dat");

  if( !ifs.is_open() ) return;
  Double_t gainList[2716]= {0};
  Double_t gainListSpec[2716]={0};

  
  int ID;
  double gain;
  while(ifs >> ID >> gain){
    gainList[ID] = gain;
  }

  std::ifstream ifsStat("Calibration_Data/Calibration_With_CC03_Fiducial_CUT/CalibrationStaticsADV_10.dat");
  //std::ifstream ifsStat("Calibration_Data/OldCalibration/CalibrationStatics_10.dat");

  if( !ifsStat.is_open()) return ;
  Double_t gainStat[2716] = {0};

  while( ifsStat >> ID >> gain ){
    gainStat[ID] = gain;
  }

  std::ifstream ifs1("Data/calibConstKe3.dat");
  while( !ifs1.eof()){
    int ID;
    double gain;
    ifs1 >> ID >> gain;
    gainListSpec[ID] = gain;
  }

  TGraph*gr = new TGraph();
  for( int i = 0; i<2716; i++){
    
    //image->Fill(i,gainList[i]);    
    gr->SetPoint(gr->GetN(), gainList[i],gainListSpec[i]);
    Double_t x,y;    
    handler->GetMetricPosition(i,x,y);

    if( TMath::Sqrt(x*x+y*y) > 850){continue;}
    if( TMath::Abs(x) < 150 && TMath::Abs(y) < 150 ){ continue;}
    if( TMath::Abs(y) > 550 ){continue;}
    
    if( !gainList[i]==0 && gainStat[i] > 300*0.3 ){
      hisResidual->Fill((gainListSpec[i])/gainList[i]);
      image->Fill(i,gainStat[i]);    
      //hisResidual->Fill((gainListSpec[i]));
    }
    
  }

  TCanvas* can = new TCanvas("can","",1600,800);
  //gr->Draw("AP");
  //image->DrawWithRange("colz",0.8,1.2);
  can->Divide(2,1);
  can->cd(1);
  hisResidual->Fit("gaus","","",0.9,1.1);  
  hisResidual->Draw();
  can->cd(2);
  image->Draw();
}
  
