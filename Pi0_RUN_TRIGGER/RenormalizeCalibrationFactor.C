#include <iostream>
#include <fstream>
void RenormalizeCalibrationFactor(){
  
  std::ifstream ifs( "CalibrationFactor/CalibrationFactorADV_15.dat");
  std::ofstream ofs( "CalibrationFactor/CalibrationFactor_Normalize.dat");
  Int_t ID;
  Double_t calFactor;
  
  Int_t chID[2716];
  Double_t CalibrationFactor[2716] = 1; 
  TH1D* hisCalFactor = new TH1D("hisCalFactor","hisCalFactor",100,0.75,1.25);

  TGraph* gr = new TGraph();
  while( ifs >> ID >> calFactor ){
    hisCalFactor->Fill(calFactor); 
    CalibrationFactor[ID] = calFactor;
    gr->SetPoint( gr->GetN() , ID, calFactor);
  };

  std::cout<< hisCalFactor->GetMean() << std::endl;
  std::cout<< hisCalFactor->GetRMS() << std::endl;
  Double_t Normalizer = hisCalFactor->GetMean();

  for( int i = 0; i< gr->GetN(); i++){
    ofs << (int)(gr->GetX()[i]) << gr->GetY()[i] << std::endl; 
  }
  ifs.close();
  ofs.close();
  hisCalFactor->Draw();
}
