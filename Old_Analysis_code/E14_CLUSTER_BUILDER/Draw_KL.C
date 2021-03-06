

void Draw_KL(int initialRunNumber){
  gStyle->SetOptStat("nerMRIuo");
  TH1D* hisKLMass[3];
  TH1D* hisKLMom[2];
  TText* textMassEntries[3];
  TText* textMomEntries[2];
  
  hisKLMass[0] = new TH1D("hisKLMassWithNOCut",Form("KLMass%06d;Mass[MeV];N/5MeV",initialRunNumber),80,400,800);
  hisKLMass[1] = new TH1D("hisKLMassCutSome",Form("KLMassWithCut%06d;Mass[MeV];N/5MeV",initialRunNumber),80,400,800);
  hisKLMass[2] = new TH1D("hisKLMassCutFull",Form("KLMassWithAllCut%06d;Mass[MeV];N/5MeV",initialRunNumber),80,400,800);
  hisKLMass[0]->SetLineColor(1);
  hisKLMass[1]->SetLineColor(2);
  hisKLMass[2]->SetLineColor(3);
  hisKLMom[0] = new TH1D("hisKLMomWithNOCut","KLMomentum;Momentum[MeV];N/20MeV",300,0,6000);
  hisKLMom[1] = new TH1D("hisKLMomCut0x00F","KLMomentum;Momentum[MeV];N/20MeV",300,0,6000);
  hisKLMom[1]->SetLineColor(2);
  
  
  
  
  TChain* ch  = new TChain("Tree");
  ch->Add(Form("klongRootFile/kl%04d.root",initialRunNumber));

  ch->Project(hisKLMass[0]->GetName(),"KlongMass[0]");
  ch->Project(hisKLMass[1]->GetName(),"KlongMass[0]" ,"(CutCondition & (8))==0 && TMath::Abs(KlongMass[0]-498) < 10");
  ch->Project(hisKLMass[2]->GetName(),"KlongMass[0]" ,"(CutCondition & (1+2+4+8))==0 && TMath::Abs(KlongMass[0]-498)<10");
  ch->Project(hisKLMom[0]->GetName() ,"KlongMom[0][2]");
  ch->Project(hisKLMom[1]->GetName() ,"KlongMom[0][2]","(CutCondition & (1+2+4+8))==0");
  
  TCanvas* can = new TCanvas("can","",1200,600);
  can->Divide(2,1);

  std::cout << hisKLMass[0]->GetEntries() << std::endl;
  std::cout << hisKLMass[1]->GetEntries() << std::endl;
  
  can->cd(1);
  gPad->SetLogy();
  hisKLMass[0]->Draw();
  hisKLMass[1]->Draw("same");
  hisKLMass[2]->Draw("same");
  textMassEntries[0] = new TText(0.5,0.5,"");
  textMassEntries[1] = new TText(0.5,0.5,"");
  textMassEntries[2] = new TText(0.5,0.5,"");
  textMassEntries[0]->SetTextSize(0.04);
  textMassEntries[1]->SetTextSize(0.04);
  textMassEntries[2]->SetTextSize(0.04);
  textMassEntries[0]->SetTextColor(1);
  textMassEntries[1]->SetTextColor(2);
  textMassEntries[2]->SetTextColor(3);
  textMassEntries[0]->DrawTextNDC(0.5,0.6,Form("#of Klong(NoCut):%d",(int)hisKLMass[0]->GetEntries()));
  textMassEntries[1]->DrawTextNDC(0.5,0.55,Form("#of Klong(SumeCut):%d",(int)hisKLMass[1]->GetEntries()));
  textMassEntries[2]->DrawTextNDC(0.5,0.5,Form("#of Klong(FullCut):%d",(int)hisKLMass[2]->GetEntries()));

  can->cd(2);
  gPad->SetLogy();
  hisKLMom[0]->Draw();
  hisKLMom[1]->Draw("same");
  can->SaveAs(Form("Image/KLong_%06d.gif",initialRunNumber));
}


