void View_with_minute(){

  TChain* ch = new TChain("laserEventData");
  
  Double_t Time;
  Int_t    ADCData[16];
  Double_t Time0;
  Double_t TimeE;
  for( int i = 3644; i< 3750; i++){   
    ch->Add(Form("laserdata_%04d.root",i));
  }
  ch->SetBranchAddress("Time",&Time);
  ch->SetBranchAddress("ADCData",ADCData);
 
  ch->GetEntry(0);
  Time0 = Time;
  ch->GetEntry(ch->GetEntries()-1);
  TimeE = Time;

  const Int_t nMinute = (Int_t)((TimeE - Time0) / 60) + 1;
  TH1D* hisOut[43200];
  TH1D* hisPed[43200];
  for( int i = 0; i< nMinute; i++){
    hisOut[i] = new TH1D(Form("hisOut%04d",i),Form("hisOut%04d",i),500,0,20000);
    hisPed[i] = new TH1D(Form("hisPed%04d",i),Form("hisPed%04d",i),500,0,5000);
  }

  for( int i = 0; i< ch->GetEntries() ; i++){
    ch->GetEntry(i);
    Int_t hisIndex = (Int_t)((Time - Time0) / 60 );
    if( ADCData[1] > 300 ){
      hisOut[hisIndex]->Fill(ADCData[1]);
    }else{
      hisPed[hisIndex]->Fill(ADCData[1]);
    }
  }

  TFile* tf = new TFile("laserOutput.root","recreate");
  TGraphErrors* gre  = new TGraphErrors();
  TGraphErrors* grelog=new TGraphErrors();

  gre->SetNameTitle("LaserPower","LaserPowerDrift;time[minute];PINOut");
  for( int i = 0 ; i< nMinute; i++ ){
    Double_t Output = hisOut[i]->GetMean() - hisPed[i]->GetMean();
    Double_t Outerr = hisOut[i]->GetMeanError();
    gre->SetPoint(i,i,Output);
    gre->SetPointError(i,0,Outerr);
    if( Output > 0 ){
      grelog->SetPoint(grelog->GetN(), i, TMath::Log(Output) );
      //grelog->SetPointError(grelog->GetN()-1, grelog->GetN()-1, TMath::Log(Outerr));
    }
  }
  TCanvas* can = new TCanvas("can","",0,0,1200,600);
  can->Divide(2,1);
  can->cd(1);
  gre->Draw("AP");
  can->cd(2);
  grelog->SetMarkerStyle(4);
  grelog->Draw("AP");

  gre->Write();
  grelog->Write();
  std::cout << TimeE - Time0 << std::endl;
  tf->Close();

}
  
