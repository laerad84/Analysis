void Pi0Mass(){
  gStyle->SetOptStat("");

  TGraphErrors* gr = new TGraphErrors();
  const Int_t nFiles=10;
  Int_t RunNumber=4515;

  TChain* chBeam = new TChain("Tree");
  TChain* chSim  = new TChain("Tree");
  chSim->Add("Pi0SIM.root");
  chBeam->Add("Pi0WAV.root");
  TH1D* hisBeamMass = new TH1D("hisBeamMass","Pi0MassBeamData",100,100,200);
  TH1D* hisSimMass  = new TH1D("hisSimMass" ,"Pi0MassSimData" ,100,100,200);
  TH1D* hisBeamGamma= new TH1D("hisBeamGamma","GammaEne",600,0,2400);
  TH1D* hisSimGamma = new TH1D("hisSimGamma","GammaEne",600,0,2400);
  chBeam->Project(hisBeamMass->GetName(),
		  "Pi0Mass[0]",
		  "GammaE[0]>250&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5");
  chSim->Project(hisSimMass->GetName(),
		 "Pi0Mass[0]",
		 "GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5");
  
  chBeam->Project(hisBeamGamma->GetName(),
		  "GammaE*0.983",
		  "GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5");

  chSim->Project(hisSimGamma->GetName(),
		 "GammaE",
		 "GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5");
  
  TF1* funcBeam = new TF1("funcBeam","gaus",0,600);
  TF1* funcSim  = new TF1("funcSim" ,"gaus",0,600);

  funcBeam->SetLineColor(1);
  funcBeam->SetLineWidth(2);
  funcSim->SetLineColor(2);
  funcSim->SetLineWidth(2);

  hisBeamMass->Scale(1);
  hisSimMass->Scale(2);
  Double_t xMean[2];
  xMean[0] = hisBeamMass->GetBinCenter(hisBeamMass->GetMaximumBin());
  xMean[1] = hisSimMass->GetBinCenter(hisSimMass->GetMaximumBin());
  hisBeamMass->Fit(funcBeam->GetName(),"","",xMean[0]-7,xMean[0]+7);
  hisSimMass->Fit(funcSim->GetName(),"","",xMean[1]-7,xMean[1]+7);
  
  TCanvas* can = new TCanvas("can","",1000,500);
  gPad->SetGridx();
  gPad->SetGridy();

  hisBeamMass->SetLineColor(1);
  hisBeamMass->Draw();
  hisSimMass->SetLineColor(2);
  hisSimMass->Draw("same");

  TText* txt[10];
  for( int i = 0; i< 10; ++i){
    txt[i] = new TText(0.,0.,"");
  }
  TCanvas *can = new TCanvas("canvas","",1000,500);
  canvas->Divide(2,1);
  canvas->cd(1);
  hisBeamGamma->Draw();
  //canvas->cd(2);
  hisSimGamma->Scale(hisBeamGamma->GetEntries()/hisSimGamma->GetEntries());
  hisSimGamma->Draw("same");
}
  
  
		      
  
