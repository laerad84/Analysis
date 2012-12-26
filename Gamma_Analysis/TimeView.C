#include <iostream>
#include "TMath.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"

void TimeView(){
  TFile* tf = new TFile("CrystalTimeEnergy_1.root");
  const int nCrystalE = 18;
  TProfile* profTE[17];
  TProfile* profE[3];
  profE[0] = (TProfile*)tf->Get(Form("profTimeEnergyCrystal_4_0"));
  profE[1] = (TProfile*)tf->Get(Form("profTimeEnergyCrystal_4_1"));
  profE[2] = (TProfile*)tf->Get(Form("profTimeEnergyCrystal_4_2"));
  
  TProfile* profLowE  =(TProfile*)tf->Get(Form("profTimeEnergyLowEnergy"));
  for( int i = 0; i< 17; i++){
    profTE[i] = (TProfile*)tf->Get(Form("profTimeEnergyCrystalMerge_%d",i));
  }
  TGraphErrors* gre = new TGraphErrors();  
  TGraphErrors* grLowE = new TGraphErrors();
  TGraphErrors* grTotal = new TGraphErrors();

  profLowE->Fit("pol0","","",50,75);
  TF1* func0 = (TF1*)profLowE->GetFunction("pol0");
  Double_t LowEOffSet = func0->GetParameter(0);
  TF1* func1 = new TF1("func","[0]+[1]*TMath::Exp([2]*x)",0,75);
        
  Int_t ECrystalArr[nCrystalE] = {0,12.5,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400};
  Double_t OffSet = profTE[3]->GetBinContent(3);
  for( int i = 2; i< 17; i++){
    std::cout<< profTE[i]->GetBinContent(3) << std::endl;
    gre->SetPoint( gre->GetN(), (ECrystalArr[i]+ECrystalArr[i+1])/2,
		   -(profTE[i]->GetBinContent(3)-OffSet));
    gre->SetPointError( gre->GetN()-1, 0, profTE[i]->GetBinError(3));    
  }
  for( int i = 1; i< profLowE->GetNbinsX(); i++){
    if( profLowE->GetBinEntries(i) == 0 ){ continue; }
    grLowE->SetPoint( grLowE->GetN(), profLowE->GetBinCenter(i), profLowE->GetBinContent(i)-LowEOffSet);
    grLowE->SetPointError( grLowE->GetN()-1, 0, profLowE->GetBinError( i));
  }
  func1->SetParameter(0, 0);
  func1->SetParameter(1, grLowE->GetY()[0] );

  grLowE->Fit( func1, "","", 0, 75 );
  TF1* func2 = new TF1("funcAll","[0]+[1]*TMath::Exp([2]*x)+[3]*TMath::Exp([4]*x)",0,400);
  TF1* func3  =new TF1("funcpol2","[0]+[1]*TMath::Exp([2]*x)",0,400);
  TF1* func4 = new TF1("funcTest","[0]+[1]*TMath::Exp([2]*x)+[3]*TMath::Exp([4]*x)",0,400);

  func3->SetParameter(0,0);
  func3->SetParameter(1,0.01);
  func3->SetParameter(2,0.005);
  

  gre->SetMarkerStyle(7);
  gre->Fit(func3,"","",50,400);

  func2->SetParameter(0, (func1->GetParameter(0)+func3->GetParameter(0))/2);
  func2->SetParameter(1, func1->GetParameter(1));
  func2->SetParameter(2, func1->GetParameter(2));
  func2->SetParameter(3, func3->GetParameter(1));
  func2->SetParameter(4, func3->GetParameter(2));
  func2->SetParLimits(1,  func1->GetParameter(1)*0.5, func1->GetParameter(1)*2);
  func2->SetParLimits(2,  func1->GetParameter(2)*2, func1->GetParameter(2)*0.5);
  func2->SetParLimits(3,  func3->GetParameter(1)*0.5, func3->GetParameter(1)*2);
  func2->SetParLimits(4,  func3->GetParameter(2)*0.5, func3->GetParameter(2)*2);


  func4->SetParameter(0, (func1->GetParameter(0)+func3->GetParameter(0))/2);
  func4->SetParameter(1, func1->GetParameter(1));
  func4->SetParameter(2, func1->GetParameter(2));
  func4->SetParameter(3, func3->GetParameter(1));
  func4->SetParameter(4, func3->GetParameter(2)*2);
  func4->SetParLimits(1,  func1->GetParameter(1)*0.5, func1->GetParameter(1)*2);
  func4->SetParLimits(2,  func1->GetParameter(2)*2, func1->GetParameter(2)*0.5);
  func4->SetParLimits(3,  func3->GetParameter(1)*0.5, func3->GetParameter(1)*2);
  func4->SetParLimits(4,  func3->GetParameter(2)*0.5*2, func3->GetParameter(2)*2*2);

  gre->Draw("AP");
  //profLowE->Draw("same ");
  grLowE->Draw("same P");
  for( int i = 0; i< gre->GetN(); i++){
    grTotal->SetPoint( grTotal->GetN(), gre->GetX()[i], gre->GetY()[i]);
    grTotal->SetPointError( grTotal->GetN()-1, 0, gre->GetEY()[i] );
  }
  for( int i = 0; i < grLowE->GetN(); i++){
    if( grLowE->GetX()[i] > 75 ){ continue; }
    grTotal->SetPoint( grTotal->GetN(), grLowE->GetX()[i], grLowE->GetY()[i]);
    grTotal->SetPointError( grTotal->GetN()-1, 0, grLowE->GetEY()[i] );
  }

  TCanvas* can = new TCanvas("can","",800,800);
  can->Divide(1,2);
  can->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  grTotal->Fit(func2, "","",0,400);
  grTotal->Draw("AP");
  //grLowE->Fit(func4,"","",0,400);
  //grLowE->Draw("AP");
  std::cout << func4->Eval(50) << std::endl;
  for( int i = 0; i< 5; i++){
    std::cout<< i << "\t" << func2->GetParameter(i) << "\t" << func2->GetParError(i) << std::endl;
  }
  std::cout << func2->GetMinimumX() << "\t" << func2->GetMinimum() << std::endl;
  TGraphErrors* grLowE = new TGraphErrors();
  TGraphErrors* grTotalDelta = new TGraphErrors();
  for( int i = 0; i< grTotal->GetN();i++){
    grTotalDelta->SetPoint( i, grTotal->GetX()[i], grTotal->GetY()[i]- func2->Eval( grTotal->GetX()[i]));
    grTotalDelta->SetPointError( i, 0, grTotal->GetEY()[i] );
  }
  can->cd(2);
  grTotalDelta->Draw("AP");
  gPad->SetGridx();
  gPad->SetGridy();

}
