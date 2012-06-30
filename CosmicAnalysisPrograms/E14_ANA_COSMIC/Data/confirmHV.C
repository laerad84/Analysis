#include <fstream>
#include <iostream>

void confirmHV(){
  gSystem->Load("../lib/libtest.so");
  IDHandler* handler  = new IDHandler("../Data/crystal.txt");
  CsIImage*  image    = new CsIImage(handler);
  
  std::ifstream ifs("CH_HV_LIST.txt");
  
  double Volt[2716]={0};
  
  while( !ifs.eof() ){
    int ch;
    double voltage;

    ifs >> ch >> voltage ;
    std::cout << ch << " \t" << voltage << std::endl;
    Volt[ch] = voltage;    
  }
  std::cout<< "image" << std::endl;
  for( int i = 0; i< 2716; i++){
    image->Fill(i, Volt[i]);
  }
  
  TCanvas* can = new TCanvas("can","",800,800);
  image->Draw();

}
  
    
  
  
