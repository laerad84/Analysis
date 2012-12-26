void TestDraw(){
  const int nE = 6;
  const int nTheta = 8;
  
  TH2D* hisRT[nE-1][nTheta-1];
  TH2D* hisDT[nE-1][nTheta-1];
  TProfile* profRT[nE-1][nTheta-1];
  TProfile* profDT[nE-1][nTheta-1];
  Int_t EArr[nE] = {100,200,300,500,800,1300};
  Int_t ThetaArr[nTheta] = {10,15,20,25,30,35,40,45};
  
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
  
  ps->Close();
}
