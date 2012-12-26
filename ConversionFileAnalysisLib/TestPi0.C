void TestPi0(){
  gStyle->SetOptFit(111111111);
  TChain* tr1 = new TChain("Tree");
  TChain* tr2 = new TChain("Tree");
  for( int i = 4503; i< 4510; i++){
    tr1->Add(Form("/Volume0/ExpData/2012_Feb_Beam/sumup_data/RECONSTRUCT_RESULT_%d.root",i));
    tr2->Add(Form("/home/jwlee/local/Analysis/Pi0_RUN_TRIGGER/test%d.root",i));
  }

  TH1D* his[3];
  for( int i = 0; i< 3; i++){
    his[i] = new TH1D(Form("his%d",i),Form("his%d",i),50,0,200);
  }
  tr1->Project(his[2]->GetName(),"Pi0Mass[0]","GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5&&TimeSigma<2");

  tr1->Project(his[1]->GetName(),"Pi0Mass[0]","GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5");
  tr2->Project(his[0]->GetName(),"Pi0Mass[0]","GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5");


  his[0]->Fit("gaus","","",125,145);
  his[1]->Fit("gaus","","",125,145);
  his[2]->Fit("gaus","","",125,145);
  
  his[0]->Draw();
  his[1]->SetLineColor(2);
  his[1]->Draw("same");
  his[2]->SetLineColor(3);
  his[2]->Draw("same");

}

    
