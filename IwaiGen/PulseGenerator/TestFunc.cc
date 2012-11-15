#include "PulseGenerator.h"

#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TApplication.h"
#include "TCanvas.h"
#include <vector>

int main( int argc, char** argv ){
  
  //TApplication* app = new TApplication("app",&argc,argv);
  PulseGenerator* gen = new PulseGenerator();
  PulseGenerator* gen1 = new PulseGenerator(); 
  TCanvas* can = new TCanvas("can","",800,800);
  std::vector<double> EVector;
  std::vector<double> TVector;
  Double_t meanTime = 0; 
  Double_t Energy   = 0; 

  for( int i = 0; i< 10; i++){
    double tempEnergy = gRandom->Gaus(40, 5);
    if( tempEnergy <= 0 ) continue; 
    EVector.push_back( tempEnergy );
    TVector.push_back( i*1);
    meanTime += tempEnergy*i*1;
    Energy   += tempEnergy;
  }
  
  meanTime = meanTime/Energy;
  can->Divide(2,2);
  for( int i = 0; i< 10000 ; i++){
    //std::cout<< "Generator" << std::endl;
    can->cd(1);
    //gen->Reset();
    TF1* func = gen->GetWaveform( EVector, TVector);
    //gen->Reset();
    TF1* func1 = gen->GetWaveform( Energy,meanTime);
    //gen->Reset();
    TF1* func2 = gen->GetWaveform( Energy, 0);
    /*
    std::cout<< func->GetMaximumX(0, 500 )   << "  " 
	     << func1->GetMaximumX( 0, 500 ) << "  "
	     << func2->GetMaximumX( 0, 500 ) << std::endl;
    */
    delete func;
    delete func1;
    delete func2;
    /*
    func->Draw("C");
    func1->Draw("same C");
    can->cd(2);
    gen->pePDF->Draw();
    can->cd(3);
    gen->pfPDF->Draw();
    can->Update();
    can->Modified();
    getchar();
    */
  }
  
  //app->Run();
}
