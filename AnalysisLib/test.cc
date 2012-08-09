#include <iostream>
#include <fstream>
#include "PeakCompensater.h"
#include "TApplication.h"
#include "TCanvas.h"

int
main( int argc, char** argv ){

  std::cout<< "Test" << std::endl;
  PeakCompensater* comp = new PeakCompensater();
  TApplication* app = new TApplication("app",&argc, argv );
  TCanvas* can  =new TCanvas("can","",800,600);
  comp->m_gr[0]->SetMarkerStyle(7);
  comp->m_gr[1]->SetMarkerStyle(7);
  comp->m_gr[2]->SetMarkerStyle(7);
  comp->m_gr[0]->SetMarkerColor(1);
  comp->m_gr[1]->SetMarkerColor(2);
  comp->m_gr[2]->SetMarkerColor(3);
  
  comp->m_gr[0]->Draw("AP");  
  comp->m_gr[1]->Draw("P");
  comp->m_gr[2]->Draw("P");

  std::cout<< comp->Compensate(0, 10) << std::endl;

  double aspect_old = 1.;
  double aspect_new = 1.;
  std::string ANALIB= std::getenv("ANALYSISLIB");
  std::ofstream ofs(Form("%s/Data/LinearityFunctionRegion.txt",ANALIB.c_str()));
  for( int Output  = 0; Output < 16000; Output++){
    aspect_new = comp->Compensate(0,Output)/((double)Output);
    for( int index = 1; index < 20; index++){
      if( aspect_new > 1+index*0.01 && aspect_old < 1+index*0.01 ){
	std::cout  << Output      << " , " 
	  //<< aspect_old  << " : " 
	  //	 << aspect_new  
		 <<  std::endl;
	ofs << 1 + index * 0.01 << "\t" <<  Output << "\t";
      }
      ofs << "\n";
    }
    aspect_old = aspect_new;
  }
  aspect_old = 1;
  for( int Output  = 0; Output < 16000; Output++){
    aspect_new = comp->Compensate(-1,Output)/((double)Output);
    for( int index = 1; index < 20; index++){
      if( aspect_new > 1+index*0.01 && aspect_old < 1+index*0.01 ){
	std::cout  << Output      << " , " 
	  //<< aspect_old  << " : " 
	  //<< aspect_new  
		 <<  std::endl;
	ofs << 1 + index * 0.01 << "\t" <<  Output << "\t";
      }
      ofs << "\n";
    }
    aspect_old = aspect_new;
  }
  aspect_old = 1;
  for( int Output  = 0; Output < 16000; Output++){
    aspect_new = comp->Compensate(2240,Output)/((double)Output);
    for( int index = 1; index < 20; index++){
      if( aspect_new > 1+index*0.01 && aspect_old < 1+index*0.01 ){
	std::cout  << Output      << " , " 
	  //<< aspect_old  << " : " 
	  //<< aspect_new  
		 <<  std::endl;
	ofs << 1 + index * 0.01 << "\t" <<  Output << "\t";
      }
      ofs << "\n";
    }
    aspect_old = aspect_new;
  }

  ofs.close();

  app->Run();
  return 0;
}
