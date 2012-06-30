#include <fstream>
#include <iostream>

void GetCosmicOut(){
  gSystem->Load("lib/libtest.so");
  IDHandler* handler  = new IDHandler("Data/crystal.txt");
  CsIImage*  image    = new CsIImage(handler);
  image->SetTitle("Distribution On image; x[mm]; y[mm]");
  TH1D* gainHist      = new TH1D("gain","Distribution of gain;Integrated ADC;N/10cnt",500,0,5000);				 
  
  std::ofstream ofs("CH_Cosmic_outList_3599_3600.txt");  
  
  Double_t gain;
  Double_t Sigma;
  Double_t Norm;
  Double_t HistEntries;
  TFile* tf = new TFile("cosmic_3599_3600.root");
  
  for( int i = 0; i< 2716; i++){    
    TH1D* hist  = (TH1D*)tf->Get(Form("his_CH%04d",i));
    TF1*  func  = NULL;
    func        = hist->GetFunction("landau");
    if( func == NULL ){
      continue;
    }
    if( i < 2240){      
      gain = func->GetParameter(1);      
      image->Fill(i, func->GetParameter(1));
      ofs << i << "\t" 
	  << gain << "\t"
	  << func->GetParameter(2) << std::endl; 
    }else{
      gain = func->GetParameter(1)/2;
      image->Fill(i, gain);
      ofs << i << "\t" 
	  << gain << "\t" 
	  << func->GetParameter(2)/2   << std::endl; 
    }
    gainHist->Fill(gain);
  }

  ofs.close();
  TCanvas* can = new TCanvas("can","",1600,800);
  can->Divide(2,1);  
  can->cd(1);
  image->Draw();
  can->cd(2);
  gainHist->Draw();

}
