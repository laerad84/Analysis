void TestDrawAll(){
  const int nE        = 6;
  const int nTheta    = 8;

  TH2D* hisRT[nE-1][nTheta-1];
  TH2D* hisDT[nE-1][nTheta-1];
  TProfile* profRT[nE-1][nTheta-1];
  TProfile* profDT[nE-1][nTheta-1];
  Int_t EArr[nE] = {100,200,300,500,800,1300};
  Int_t ThetaArr[nTheta] = {10,15,20,25,30,35,40,45};

  const int nESIM     = 7;
  const int nThetaSIM = 15; 
  TFile* tfSIM[nESIM][nThetaSIM][2];
  TH2D* hisRTSIM[nESIM][nThetaSIM];
  TH2D* hisDTSIM[nESIM][nThetaSIM];
  TProfile* profRTSIM[nESIM][nThetaSIM];
  TProfile* profDTSIM[nESIM][nThetaSIM];
  TH2D* hisRTSIMB[nESIM][nThetaSIM];
  TH2D* hisDTSIMB[nESIM][nThetaSIM];
  TProfile* profRTSIMB[nESIM][nThetaSIM];
  TProfile* profDTSIMB[nESIM][nThetaSIM];
  
  Int_t ESIMArr[nESIM] = { 100, 145, 210, 310, 450, 650, 950};
  Int_t ThetaSIMArr[nThetaSIM] = { 0, 10, 12, 14, 16, 18, 20, 22, 25, 28, 30, 32, 35, 40, 50};
  
  
  TFile* tf = new TFile("ClusterTimeStructure.root");
  for( int iIndex  = 0; iIndex < nE -1; iIndex++){
    for( int jIndex = 0; jIndex < nTheta -1; jIndex++){
      hisRT[iIndex][jIndex] = (TH2D*)tf->Get(Form("hisRT_r_E_%d_%d_Theta_%d_%d",
						  EArr[iIndex],EArr[iIndex+1],
						  ThetaArr[jIndex],ThetaArr[jIndex+1]));
      hisDT[iIndex][jIndex] = (TH2D*)tf->Get(Form("hisDT_r_E_%d_%d_Theta_%d_%d",
						  EArr[iIndex],EArr[iIndex+1],
						  ThetaArr[jIndex],ThetaArr[jIndex+1]));
      profRT[iIndex][jIndex] = hisRT[iIndex][jIndex]->ProfileX();
      profDT[iIndex][jIndex] = hisDT[iIndex][jIndex]->ProfileX();
    }
  }

  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  
  for( int iIndex =  0; iIndex < nESIM; iIndex++){
    for( int jIndex  = 0; jIndex < nThetaSIM; jIndex++){
     tfSIM[iIndex][jIndex][0] = new TFile(Form("%s/ClusterTimeStructure_SIM_%dMeV_%ddeg.root",
					       ROOTFILE_GAMMACLUS.c_str(),ESIMArr[iIndex],ThetaSIMArr[jIndex]));
     hisRTSIM[iIndex][jIndex] = (TH2D*)tfSIM[iIndex][jIndex][0]->Get(Form("hisRT_E_%d_%d_Theta_%d_%d",
									  ESIMArr[iIndex],ESIMArr[iIndex],
									  ThetaSIMArr[jIndex],ThetaSIMArr[jIndex]));
     profRTSIM[iIndex][jIndex] = hisRTSIM[iIndex][jIndex]->ProfileX();
     tfSIM[iIndex][jIndex][1] = new TFile(Form("%s/ClusterTimeStructure_SIM_Back_%dMeV_%ddeg.root",
					       ROOTFILE_GAMMACLUS.c_str(),ESIMArr[iIndex],ThetaSIMArr[jIndex]));
     hisRTSIMB[iIndex][jIndex] = (TH2D*)tfSIM[iIndex][jIndex][1]->Get(Form("hisRT_E_%d_%d_Theta_%d_%d",
									   ESIMArr[iIndex],ESIMArr[iIndex],
									   ThetaSIMArr[jIndex],ThetaSIMArr[jIndex]));
     profRTSIMB[iIndex][jIndex]= hisRTSIMB[iIndex][jIndex]->ProfileX();
    }    
  }
  

  
  TPostScript* ps  = new TPostScript("ClusterTimeStructure.ps",112);
  TCanvas * can = new TCanvas("can","",1400,1000);
  can->Divide(3,3); 
  for( int iIndex = 0; iIndex <nE-1; iIndex++){
    ps->NewPage();
    for( int jIndex = 0; jIndex < nTheta-1; jIndex++ ){
      can->cd( jIndex +1);
      gPad->SetGridx();
      gPad->SetGridy();
      profRT[iIndex][jIndex]->Draw();
    }
    can->Update();
  }
  
  for( int iIndex = 0; iIndex <nE-1; iIndex++){
    ps->NewPage();
    for( int jIndex = 0; jIndex < nTheta-1; jIndex++ ){
      can->cd( jIndex +1);
      gPad->SetGridx();
      gPad->SetGridy();
      profDT[iIndex][jIndex]->Draw();
    }
    can->Update();
  }

  can->Clear();
  can->Divide(4,4);
  for( int iIndex = 0; iIndex <nESIM; iIndex++){
    ps->NewPage();
    for( int jIndex = 0; jIndex < nThetaSIM; jIndex++ ){
      can->cd( jIndex +1);
      gPad->SetGridx();
      gPad->SetGridy();
      profRTSIM[iIndex][jIndex]->SetLineColor(1);
      profRTSIMB[iIndex][jIndex]->SetLineColor(2);
      profRTSIM[iIndex][jIndex]->GetYaxis()->SetRangeUser( -2,2 );
      profRTSIM[iIndex][jIndex]->Draw();
      profRTSIMB[iIndex][jIndex]->Draw("same");      
      
    }
    can->Update();
  }






  ps->Close();
}
