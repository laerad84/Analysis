#include <fstream>
#include <iostream>
void testResult(){
  
  std::ifstream ifs0("simGain.txt");
  std::ifstream ifs1("CalibrationN.txt");
  if( !ifs0.is_open() ){return;}
  if( !ifs1.is_open() ){return;}

  TH1D* hisGain = new TH1D("test","",200,0,2);
  TH1D* hisInit = new TH1D("test1","",200,0,2);
  double cal0;
  double cal1;
  double calfactor0[2716]={0};
  double calfactor1[2716]={0};
  
  Int_t  N;
  Int_t ID;
  while( !ifs0.eof() ){
    ifs0 >> ID >> cal0;
    calfactor0[ID] = cal0;
    std::cout<< ID << std::endl;
  }

  while( !ifs1.eof() ){
    ifs1 >> ID >> cal1 >> N;
    std::cout<< ID << std::endl;
    calfactor1[ID] = cal1;
  }
  
  for( int i = 0; i< 2716; i++){
    if(calfactor0[i] ==0 || calfactor1[i] ==0){
      continue;
    }
    hisGain->Fill(calfactor0[i]/calfactor1[i]);
    hisInit->Fill(calfactor0[i]);
  }
  hisGain->Draw();
  hisInit->SetLineColor(2);
  hisInit->Draw("same");

  std::cout<< hisGain->GetRMS()/hisGain->GetMean() << "\t"
	   << hisInit->GetRMS()/hisInit->GetMean() << std::endl;
}
