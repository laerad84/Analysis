#include <iostream>
void  CalFactorGenerator(){
  gStyle->SetOptStat(0);
  const  int nPercent = 16;
  Int_t Width[nPercent]= {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
  Int_t ColorSpace[nPercent] ={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

  TH1D* hisDistribution[nPercent];
  TH1D* hisBack = new TH1D("CalibrationFactor",
			   "Calibration Factor;Calibration Factor",
			   100,0.5,1.5);


  for( int i = 0; i< nPercent; i++){
    hisDistribution[i]=  new TH1D(Form("hisDistribution_%d",i),
				  Form("Sigma:%d\%;CalibrationFactor",Width[i]),
				  100,0.5,1.5);

    std::ofstream ofs(Form("CalFactor_%d.txt",Width[i]));
    for( int ich = 0; ich < 2716; ich++ ){
      Double_t Calfactor = gRandom->Gaus(1,0.01*Width[i]);
      ofs << ich << "\t" << Calfactor<< "\n";
      hisDistribution[i]->Fill(Calfactor);
    }
    ofs.close();
  }
  TCanvas* can = new TCanvas("can","can",800,800);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  hisBack->Draw();
  hisBack->GetYaxis()->SetRangeUser(1,1600);  
  gPad->SetFrameLineWidth(2);
  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->SetFillColor(0);

  for( int i = 0; i< nPercent; i++){
    leg->AddEntry(hisDistribution[i],Form("Sigma:%d\%",Width[i]));
  }

  for( int i = 0; i< nPercent; i++){
    hisDistribution[i]->SetLineColor(ColorSpace[i]);
    hisDistribution[i]->SetFillStyle(3001+i);
    hisDistribution[i]->SetFillColor(ColorSpace[i]);
    hisDistribution[i]->SetLineWidth(2);
    hisDistribution[i]->Draw("same");
  }  
  leg->Draw();
}
