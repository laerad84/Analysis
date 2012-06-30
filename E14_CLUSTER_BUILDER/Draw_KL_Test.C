void Draw_KL_Test(){

  TChain* ch = new TChain("Tree");
  TChain* ch1 = new TChain("Tree");
  
  TH1D* his  = new TH1D("Klong6g","Klong6g",20,450,550);
  TH1D* his1 = new TH1D("Klong4g","Klong4g",20,450,550);
  TH1D* his2 = new TH1D("Klong4gAll","Klong4gAll",60,250,550);

  for( int i = 0; i< 68; i++){
    ch->Add(Form("klongRootFile/kl%d.root" ,4162+i));
    ch1->Add(Form("klongRootFile/ks%d.root",4162+i));
  }
  ch->Project(his->GetName()  ,"KlongMass[0]","CutCondition==0");
  ch1->Project(his1->GetName(),"KlongMass[0]","CutCondition==0");
  ch1->Project(his2->GetName(),"KlongMass[0]","CutCondition==0");
  
  TF1* func = new TF1("func","gaus(0)+expo(3)",0,550);
  func->SetParameter(1,498);
  func->SetParameter(2,5);
  TF1* func2 = new TF1("func2","gaus(0)",0,550);
  func2->SetParameter(1,498);
  func2->SetParameter(2,5);



  TCanvas* can = new TCanvas("can","",1200,600);
  can->Divide(2,1);
  can->cd(1);
  his2->Fit(func->GetName(),"","",450,550);
  his2->Draw();
  TF1* func1 = new TF1("Test","gaus",450,550);
  func1->SetParameter(0,func->GetParameter(0));
  func1->SetParameter(1,func->GetParameter(1));
  func1->SetParameter(2,func->GetParameter(2));

  can->cd(2);
  his1->SetLineColor(2);  
  his->Draw();
  his->Fit(func2->GetName(),"","",450,550);
  func->Draw("same");
  his1->Draw("same");

  std::cout<< func2->GetParameter(0) << " " 
	   << func->GetParameter(0)  << " " 
	   << func->GetParameter(0)/func2->GetParameter(0)<< std::endl;
  std::cout<< func2->Integral(450,550) << " " 
	   << func1->Integral(450,550)  << " " 
	   << func1->Integral(450,550)/func2->Integral(450,550)
	   << std::endl;
  //ch->Draw("KlongPt[0]:KlongMass[0]>>(400,200,600,50,0,20)","(CutCondition&(1|2|4|8))==0","colz");
  gPad->SetLogz();
  TText* text = new TText(0.5,0.5,"");
  TText* text1 = new TText(0.5,0.5,"");
  text->DrawTextNDC(0.5,0.5,Form("Integral:%2.3lf",func1->Integral(450,550)));
  text1->DrawTextNDC(0.5,0.6,Form("Integral:%2.3lf",func2->Integral(450,550)));
  
}
