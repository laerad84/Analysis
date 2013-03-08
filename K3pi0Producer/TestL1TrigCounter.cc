#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TApplication.h"
#include "L1TrigCounter.h"

int main( int argc , char** argv ){
  TApplication* app = new TApplication("app",&argc, argv );
  std::cout <<"Define" << std::endl;
  L1TrigCounter l1;
  bool rst = l1.ReadMapFile();
  if( !rst ){std::cout<<"File" << std::endl;}
  std::cout<< "SetThreshold"<< std::endl;

  l1.SetThreshold( 10000 );
  std::cout<<__LINE__<<std::endl;
  l1.Reset();
  std::cout<<" prepare to Fill" << std::endl;
  for( int i = 0; i< 2716; i++){
    l1.Fill( i, i*10);
  }
  std::vector<double> output = l1.GetCount();
  for( int i = 0; i < output.size(); i++){
    std::cout<< i << "\t" << output.at(i) << std::endl;
  }
  std::cout<< l1.GetTriggedCrate().size() << std::endl;
  l1.Reset();
  TCanvas* can = new TCanvas("can","can",1200,600);
  can->Divide(2,1);
  can->cd(1);
  l1.DrawTrig();
  can->cd(2);
  l1.DrawTriggerMap();
  can->Update();
  app->Run();
}
