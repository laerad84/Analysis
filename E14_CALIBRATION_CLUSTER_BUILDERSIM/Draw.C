#include <fstream>
void Draw(){

  gSystem->Load("../AnalysisLib/lib/libAnalysisLib.so");

  IDHandler* handler = new IDHandler();
  CsIImage* image    = new CsIImage(handler);
  TH1D* hisResidual = new TH1D("his","",80,-0.2,0.2);
  
  std::ifstream ifs("CalibrationData/OldCalibrationData/CalibrationFactor_8.dat");
  if( !ifs.is_open() ){std::cout << "OLD" << std::endl; return;}
  
  Double_t gainList[2716]= {1};
  Double_t gainListSpec[2716]={1};
  
  while( !ifs.eof() ){
    int ID;
    double gain;
    ifs >> ID >> gain;
    std::cout << ID << " : " <<  gain << std::endl;
    gainList[ID] = gain;
  }

  std::ifstream ifs1("Data/calibConstKe3.dat");
  if( !ifs1.is_open()){ std::cout<<"TEST" <<
  while( !ifs1.eof()){
    int ID;
    double gain;
    ifs1 >> ID >> gain;
    gainListSpec[ID] = gain;
  }

  TGraph*gr = new TGraph();
  for( int i = 0; i<2716; i++){
    image->Fill(i,gainList[i]);    
    gr->SetPoint(gr->GetN(), gainList[i],gainListSpec[i]);
    if( !gainList[i]==0){
      hisResidual->Fill((gainList[i]-gainListSpec[i])/gainList[i]);
    }
  }


  TCanvas* can = new TCanvas("can","",800,800);
  //gr->Draw("AP");
  //image->DrawWithRange("colz",0.8,1.2);
  hisResidual->Fit("gaus","","",-0.03,0.06);
  hisResidual->Draw();
  //image->Draw();
}
  
