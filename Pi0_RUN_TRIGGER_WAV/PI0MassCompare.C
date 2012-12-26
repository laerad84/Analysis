void PI0MassCompare(){

  TGraphErrors* gr = new TGraphErrors();
  const Int_t nFiles=10;
  Int_t RunNumber=4503;

  TChain* chBeam = new TChain("Tree");
  TChain* chSim  = new TChain("Tree");
  std::string SimData= "/group/had/koto/ps/klea/work/jwlee/RootFiles/Pi0Run/Pi0Build/SIM_e14_AL_Target.mac_1000000_%d_Pi.root";
  std::string BeamData = "/home/had/jwlee/local/Analysis/Data/2012_FEB/Pi0_data/Pi0_CAL_%d.root";

  for( int i = 0; i< nFiles; ++i){
    chBeam->Add(Form(BeamData.c_str(),RunNumber+i));
  }
  for( int i = 0; i< 250; ++i ){
    chSim->Add(Form(SimData.c_str(),i));
  }

  TH1D* hisBeamMass = new TH1D("hisBeamMass","Pi0MassBeamData",600,0,600);
  TH1D* hisSimMass  = new TH1D("hisSimMass" ,"Pi0MassSimData" ,600,0,600);
  TH1D* hisBeamGamma= new TH1D("hisBeamGamma","GammaEne",600,0,2400);
  TH1D* hisSimGamma = new TH1D("hisSimGamma","GammaEne",600,0,2400);
  chBeam->Project(hisBeamMass->GetName(),
		  "Pi0Mass[0]",
		   "GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5");
  chSim->Project(hisSimMass->GetName(),
		 "Pi0Mass[0]",
		 "GammaE[0]>250&&GammaE[1]>250&&GammaChi2[0]<5&&GammaChi2[1]<5");
  
  chBeam->Project(hisBeamGamma->GetName(),
		  "GammaE",
		  "GammaChi2[0]<5&&GammaChi2[1]<5");
  chSim->Project(hisSimGamma->GetName(),
		 "GammaE",
		 "GammaChi2[0]<5&&GammaChi2[1]<5");
  



  
  TF1* funcBeam = new TF1("funcBeam","gaus",0,600);
  TF1* funcSim  = new TF1("funcSim" ,"gaus",0,600);

  funcBeam->SetLineColor(1);
  funcBeam->SetLineWidth(1);
  funcSim->SetLineColor(2);
  funcSim->SetLineWidth(1);

  hisBeamMass->Scale(1);
  hisSimMass->Scale(1);
  hisBeamMass->Fit(funcBeam->GetName(),"","",128,140);
  hisSimMass->Fit(funcSim->GetName(),"","",128,140);
  
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
  
  txt[9]->DrawTextNDC(0.5,0.75,Form("Parameter(Simulation Data is Normalized(0.7))"));
  txt[8]->DrawTextNDC(0.5,0.7,Form("Beam(80kevent),Sim(2.5Mevent)"));
  txt[0]->DrawTextNDC(0.5,0.65,Form("Total(Beam):%5.1f",hisBeamMass->GetEntries()));
  txt[1]->DrawTextNDC(0.5,0.6,Form("Peak(Beam):%5.1f",funcBeam->GetParameter(1)));
  txt[2]->DrawTextNDC(0.5,0.55,Form("Sigma(Beam):%5.1f",funcBeam->GetParameter(2)));

  txt[3]->SetTextColor(2);
  txt[4]->SetTextColor(2);
  txt[5]->SetTextColor(2);
  txt[3]->DrawTextNDC(0.5,0.5,Form("Total(Sim):%5.1f",hisSimMass->GetEntries()));
  txt[4]->DrawTextNDC(0.5,0.45,Form("Peak(Sim):%5.1f",funcSim->GetParameter(1)));
  txt[5]->DrawTextNDC(0.5,0.4,Form("Sigma(Sim):%5.1f",funcSim->GetParameter(2)));
  TCanvas *can = new TCanvas("canvas","",1000,500);
  canvas->Divide(2,1);
  canvas->cd(1);
  hisBeamGamma->Draw();
  canvas->cd(2);
  hisSimGamma->Draw();

}
  
  
		      
  
