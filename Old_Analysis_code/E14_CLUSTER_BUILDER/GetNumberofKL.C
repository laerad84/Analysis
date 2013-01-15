void GetNumberofKL(){
  
  int initialRunNumber= 4158;
  int finalRunNumber  = 4624;

  TGraph* grTotal = new TGraph();
  TGraph* grMassZ = new TGraph();
  
  long nTotal=0;
  long nTotalCut=0;
  for( int runNumber  = initialRunNumber;
       runNumber < finalRunNumber;
       runNumber++){
    
    TFile* tf =new TFile(Form("klongRootFile/kl%04d.root",runNumber));
    TTree* tr = (TTree*)tf->Get("Tree");
    if( tr == NULL ){
      continue;
    }
    TH1D* histMass  = new TH1D("histMass","",400,0,2500);
    TH1D* hisMassCut = new TH1D("histMassCut","",400,0,2500);
    tr->Project(histMass->GetName()   ,"KlongMass[0]");
    tr->Project(histMassCut->GetName(),"KlongMass[0]",
		"(CutCondition & 8) == 0 && TMath::Abs(KlongMass[0]-498)<10");
    nTotal    += histMass->GetEntries();
    nTotalCut += histMassCut->GetEntries();
    
    grTotal->SetPoint(grTotal->GetN(),runNumber,nTotal);
    grMassZ->SetPoint(grMassZ->GetN(),runNumber,nTotalCut);
  }
  grTotal->SetMarkerColor(1);
  grTotal->SetMarkerStyle(4);
  grMassZ->SetMarkerColor(2);
  grMassZ->SetMarkerStyle(4);
  
  TText* txt  = new TText(0.2,0.7,Form("TotalKL:%d",nTotal));
  TText* txt1 = new TText(0.2,0.6,Form("KLwithCut:%d",nTotalCut));  
  TCanvas* can = new TCanvas("can","can",800,800);
  gPad->SetGridx();
  gPad->SetGridy();
  grTotal->Draw("AP");
  grMassZ->Draw("P");
  
  txt->DrawTextNDC(0.2,0.7,Form("TotalKL:%d",nTotal));
  txt1->SetTextColor(2);
  txt1->DrawTextNDC(0.2,0.6,Form("KLwithCut:%d",nTotalCut));  
  
}
  
  
