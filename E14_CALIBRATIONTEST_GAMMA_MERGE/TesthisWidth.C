void TesthisWidth(){
  
  TFile* tf = new TFile("SimCalibration_Data/Merge_All.root");
  //TFile* tf = new TFile("SimCalibration_Data/nonFixing/Merge_All.root");

  TH2D* his_GammaE_Pi0= (TH2D*)tf->Get("his_GammaE_Pi0_All"); 
  TH2D* his_Pi0Mass_GammaE_All_CUT = (TH2D*)tf->Get("his_Pi0Mass_GammaE_All_CUT");
  TH2D* his_Cal_Pi0Mass = (TH2D*)tf->Get("his_Cal_Pi0Mass");
  TH2D* his_Cal_KLMass  = (TH2D*)tf->Get("his_Cal_KLMass");
  TH2D* his_Cal_GammaE  = (TH2D*)tf->Get("his_Cal_GammaE");
  TH2D* his_Cal_LGammaE  = (TH2D*)tf->Get("his_Cal_LGammaE");
  TH2D* his_GammaE_DeltaE = (TH2D*)tf->Get("his_GammaE_DeltaE");

  TH2D* his_Pi0Mass_GammaE = (TH2D*)tf->Get("his_Pi0Mass_GammaE");
  TH2D* his_Pi0Mass_LGammaE = (TH2D*)tf->Get("his_Pi0Mass_LGammaE");

  TH1D* his_Pi0Mass_GammaE_proy[40];
  TH1D* his_Pi0Mass_LGammaE_proy[40];
  for( int i = 0; i< 40; i++){
    his_Pi0Mass_GammaE_proy[i] = his_Pi0Mass_GammaE->ProjectionY(Form("his_Pi0Mass_GammaE_proy_%d",i), i,i+1);
    his_Pi0Mass_LGammaE_proy[i] = his_Pi0Mass_LGammaE->ProjectionY(Form("his_Pi0Mass_LGammaE_proy_%d",i), i,i+1);
  }
  
  std::cout<< "DEBUG" << std::endl;
  TGraphErrors* gr_Pi0Mass_GammaE_RMS = new TGraphErrors();
  TGraphErrors* gr_Pi0Mass_LGammaE_RMS = new TGraphErrors();
  std::cout<< "DEBUG" << std::endl;
  
  for( int i = 0; i< 40; i++){
    std::cout<< "DEBUG:"<< i << std::endl;
    

    if( his_Pi0Mass_GammaE_proy[i]->GetEntries()>0){
      gr_Pi0Mass_GammaE_RMS->SetPoint(gr_Pi0Mass_GammaE_RMS->GetN(),i*50,
				      his_Pi0Mass_GammaE_proy[i]->GetRMS());
      gr_Pi0Mass_GammaE_RMS->SetPointError(gr_Pi0Mass_GammaE_RMS->GetN()-1,0,
					   his_Pi0Mass_GammaE_proy[i]->GetRMSError());
    }
    std::cout<< "DEBUG:"<< i << std::endl;

    if( his_Pi0Mass_LGammaE_proy[i]->GetEntries()>0){

      gr_Pi0Mass_LGammaE_RMS->SetPoint(gr_Pi0Mass_LGammaE_RMS->GetN(),i*50.,
				       his_Pi0Mass_LGammaE_proy[i]->GetRMS());
      gr_Pi0Mass_LGammaE_RMS->SetPointError(gr_Pi0Mass_LGammaE_RMS->GetN()-1,0,
					    his_Pi0Mass_LGammaE_proy[i]->GetRMSError());
    }
  }
  std::cout<< "DEBUG" << std::endl;
								 

  TH1D* his[100];  
  std::cout<< his_Pi0Mass_GammaE_All_CUT->GetNbinsX() << std::endl;
  std::cout<< his_Pi0Mass_GammaE_All_CUT->GetNbinsY() << std::endl;
  for( int i = 0; i< 100; ++i ){
    //his[i] = his_GammaE_Pi0_All->ProjectionX(Form("his%d",i),i,i+1);
    his[i] = his_Pi0Mass_GammaE_All_CUT->ProjectionY(Form("his%d",i),i, i+1);
  }
  
  //his[10]->Draw();
  TGraph* gr = new TGraph();
  for( int i = 0; i< 50; ++i){
    //his[i]->SetLineColor(i%10);
    //his[i]->Draw("same");
    gr->SetPoint( gr->GetN(), i , his[i]->GetRMS() );
    if( his[i]->GetEntries() != 0)
      his[i]->Scale(1/his[i]->GetBinContent(his[i]->GetMaximumBin()));
  }
  TProfile* prof = his_Cal_GammaE->ProfileX();
  TProfile* prof_GammaE_DeltaE = his_GammaE_DeltaE->ProfileX();
  TGraph*   gr_Aspect = new TGraph();
  for( int i = 1; i < prof->GetNbinsX()+1;i++){    
    Double_t CalDelta = prof->GetBinContent(i);
    Double_t GamDelta = prof_GammaE_DeltaE->GetBinContent(i);
    if( GamDelta !=0)
      gr_Aspect->SetPoint( gr_Aspect->GetN(),i,CalDelta/GamDelta);
  }

  TCanvas * can = new TCanvas("can","",1200,1200);
  can->Divide(3,3);
  can->cd(1);
  
  //gr->SetMarkerStyle(5);
  //gr->Draw("AP");
  his_GammaE_DeltaE->Draw("colz");
  

  can->cd(2);
  prof_GammaE_DeltaE->Draw();
  prof->SetLineColor(2);
  prof->Draw("same");

  

  /*
  his[6]->SetLineColor(2);
  his[6]->Draw();
  for( int i = 0 ; i< 10; i++){
    his[(i+1)]->SetLineColor(i+1);
    his[(i+1)]->Draw("same");
  }

  gPad->SetLogy();
  */
  can->cd(3);
  his_GammaE_Pi0->Draw("colz");
  can->cd(4);
  gr_Aspect->SetMarkerStyle(4);
  gr_Aspect->Draw("AP");
  //his_Pi0Mass_GammaE_All_CUT->Draw("colz");  
  //prof->Draw();

  gPad->SetGridx();
  gPad->SetGridy();

  can->cd(5);
  //his_Cal_Pi0Mass->Draw("colz");
  gr_Pi0Mass_GammaE_RMS->SetMarkerStyle(5);
  gr_Pi0Mass_GammaE_RMS->GetYaxis()->SetRangeUser(0,3.5);
  gr_Pi0Mass_GammaE_RMS->GetXaxis()->SetRangeUser(0,2000);
  gr_Pi0Mass_GammaE_RMS->Draw("AP");
  can->cd(6);
  //his_Cal_KLMass->Draw("colz");
  gr_Pi0Mass_LGammaE_RMS->SetMarkerStyle(5);
  gr_Pi0Mass_LGammaE_RMS->GetYaxis()->SetRangeUser(0,3.5);
  gr_Pi0Mass_LGammaE_RMS->GetXaxis()->SetRangeUser(0,2000);
  gr_Pi0Mass_LGammaE_RMS->Draw("AP");

  can->cd(7);
  his_Cal_GammaE->Draw("colz");
  prof->Draw("same");
  can->cd(8);
  his_Pi0Mass_GammaE->Draw("colz");
  can->cd(9);
  his_Pi0Mass_LGammaE->Draw("colz");

}

