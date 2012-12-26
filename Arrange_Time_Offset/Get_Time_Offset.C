#include <fstream>
#include <iostream>

void Get_Time_Offset(){
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("nerMI");

  TFile* ifile  = new TFile("TimeDistributionPi0.root");
  TH2D*  hisTime= (TH2D*)ifile->Get("hisTime");
  TCanvas* can = new TCanvas("can","", 800,800);
  TH1D*  hisTimeChannel[2716];
  TGraph* grEntries = new TGraph();

  TGraph* grTime = new TGraph();
  TGraph* grSigma= new TGraph(); 
  grEntries->SetMarkerStyle(6);
  grTime->SetMarkerStyle(6);
  grSigma->SetMarkerStyle(6);
  std::ofstream TimeDeltaofs("TimeDelta.txt");

  for( int i = 0; i< 2716; i++){
    hisTimeChannel[i] = new TH1D(Form("his%d",i),"",200,-50,50);
    hisTimeChannel[i] = hisTime->ProjectionY(Form("his%d",i),i+1,i+2);
    grEntries->SetPoint( i, i, hisTimeChannel[i]->GetEntries() ); 
  }
  TPostScript* ps = new TPostScript("TimeDelta.ps",111);
  

  Double_t tempPeakX[2716]     = {-999};
  Double_t tempPeakY[2716]     = {-999};
  Double_t tempRMS[2716]       = {-999};
  Double_t gaussianPeak[2716]  = {-999};
  Double_t gaussianSigma[2716] = {0};
  Double_t Chisquare[2716]     = {-1}; 
  Int_t    Fit[2716]           = {0};
  TF1* func = new TF1("Gaus","[0]*TMath::Gaus(x,[1],[2])",-50,50);
  for( int i = 0; i< 2716; i++){
    tempPeakX[i] = hisTimeChannel[i]->GetBinCenter(hisTimeChannel[i]->GetMaximumBin());
    tempPeakY[i] = hisTimeChannel[i]->GetBinContent(hisTimeChannel[i]->GetMaximumBin());
    tempRMS[i]  = hisTimeChannel[i]->GetRMS()/2;  
    if( hisTimeChannel[i]->GetEntries() != 0 ){
      func->SetParameter( 0, tempPeakY[i] );
      func->SetParameter( 1, tempPeakX[i] );
      //func->SetParameter( 2, tempRMS[i] );
      func->SetParLimits( 0, 0.4*tempPeakY[i], 1.5*tempPeakY[i] );
      func->SetParLimits( 1, tempPeakX[i]-5, tempPeakX[i]+5);
      func->SetParLimits( 2, 0 , tempRMS[i]*2);

      hisTimeChannel[i]->Fit(func->GetName(),"","",
			     tempPeakX[i]-tempRMS[i],
			     tempPeakX[i]+tempRMS[i]);
      can->Update();
      can->Modified();

      gaussianPeak[i] = func->GetParameter(1);
      gaussianSigma[i]= func->GetParameter(2);
      grTime->SetPoint(grTime->GetN(), i, gaussianPeak[i] );
      grSigma->SetPoint( grSigma->GetN(), i , gaussianSigma[i]);
      
      TimeDeltaofs << i << "\t" 
		   << gaussianPeak[i] << "\t" 
		   << gaussianSigma[i]<< "\t" 
		   << func->GetChisquare()/func->GetNDF() << "\t" 
		   << func->GetParameter(0) << "\t" 
		   << hisTimeChannel[i]->GetEntries() << "\n";
    }
  }
  TimeDeltaofs.close();
  ps->NewPage();
  can->Clear();
  can->Divide(2,2);
  can->cd(1);
  grTime->Draw("AP");
  can->cd(2);
  grSigma->Draw("AP");
  can->cd(3);
  grEntries->Draw("AP");
  can->Update();
  can->Modified();
  ps->Close();

}
