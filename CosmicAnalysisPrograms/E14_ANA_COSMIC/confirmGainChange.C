#include <iostream>
#include <fstream>

void confirmGainChange(){
  

  gStyle->SetPalette(1);
  gSystem->Load("lib/libtest.so");
  IDHandler* handler = new IDHandler("Data/crystal.txt");
  CsIImage*  image   = new CsIImage(handler);
  
  
  std::ifstream ifs("ANA_Cosmic_FEB_2.txt");
  const int nCH = 2716;
  double gainBefore[nCH]      = {0};
  double gainAfter[nCH]       = {0};
  double gainLaserBefore[nCH] = {0};
  double gainLaserAfter[nCH]  = {0};
  double ratio[nCH]           = {0};
  
  while(!ifs.eof()){
    int    ich, voltage;
    double gain0, gain1;
    ifs >> ich >> voltage >>  gain0 >>  gain1;
    gainBefore[ich] = gain0;
    gainAfter[ich]  = gain1;
    if( gainBefore[ich] != 0 ){
      ratio[ich] = gainAfter[ich]/ gainBefore[ich];
    }
  }


  int beforeChangeRunNumber=3681;
  int afterChangeRunNumber=3700;
  TFile* tfBefore = new TFile(Form("laserROOTFiles/testLaser%04d.root",
				   beforeChangeRunNumber));
  TFile* tfAfter  = new TFile(Form("laserROOTFiles/testLaser%04d.root",
				   afterChangeRunNumber));
  TGraph* gr = new TGraph();
  gr->SetNameTitle("Confirm_Gain_Change","Confirm_Gain_Change;Set Value;Laser Value");
  
  std::cout << "OK" << std::endl;
  for( int i = 0; i< nCH; i++){
    TH1D* hisBefore= NULL;
    TH1D* hisAfter = NULL;
    
    hisBefore = (TH1D*)tfBefore->Get(Form("laser%d",i));
    if( hisBefore == NULL){
      continue;
    }
    hisAfter  = (TH1D*)tfAfter->Get(Form("laser%d",i));
    if( hisAfter == NULL ){
      continue;
    }
    
    // taking Mean of Histogram as LY //
    gainLaserBefore[i] = hisBefore->GetMean();
    gainLaserAfter[i]  = hisAfter->GetMean();
    
    //std::cout <<i << "\t" <<  gainLaserBefore[i] << "\t" << gainLaserAfter[i] << std::endl;
  }
  
  for( int i = 0; i< nCH; i++){
    if( gainBefore[i] != 0 && gainLaserBefore[i] != 0){
      gr->SetPoint(gr->GetN(),gainAfter[i]/gainBefore[i],
		   gainLaserAfter[i]/gainLaserBefore[i]);

      image->Fill (i,1 - (gainAfter[i]/gainBefore[i])/(gainLaserAfter[i]/gainLaserBefore[i]));
      if( TMath::Abs(1 - (gainAfter[i]/gainBefore[i])/(gainLaserAfter[i]/gainLaserBefore[i])) >0.1){
	//std::cout<< i << "\t" <<  (gainAfter[i]/gainBefore[i])/(gainLaserAfter[i]/gainLaserBefore[i])
	//<< std::endl; 
      }
    }else{
      std::cout<< i << "\t" << gainBefore[i]  << "\t" 
	       << gainLaserBefore[i] << std::endl;
    }
  }

  //gr->Draw("AP");
  image->Draw();
}


    

    
				     
    
    

