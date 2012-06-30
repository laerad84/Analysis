#include <fstream>
void Test(){
  std::ifstream ifs1("Calibration_Data/CalibrationFactor_10.dat");
  std::ifstream ifs2("Calibration_Data/CalibrationFactorADV_10.dat");
  
  Double_t list1[2716];
  Double_t list2[2716];
  for( int i = 0; i< 2716; ++i ){
    list1[i] = 0;
    list2[i] = 0;
  }
  
  int id;
  double cal;
  while( ifs1 >> id >> cal ){
    list1[id] = cal;
  }
  while( ifs2 >> id >> cal) {
    list2[id] = cal;
  }
  TH1D* his  = new TH1D("his","",100,-0.25,0.25);
  for( int i = 0; i< 2716; ++i){
    if( list1[i] !=0 && list2[i]!= 0){
      his->Fill((list1[i]-list2[i])/list1[i]);
    }
  }

  TCanvas* can = new TCanvas("can","",800,400);
  can->Divide(2,1);
  can->cd(1);
  
  TGraph* gr = new TGraph(2716, list1, list2);
  gr->Draw("AP");
  can->cd(2);
  his->Draw();
  his->Fit("gaus","","",-0.05,0.05);
}

		       
		     
  
