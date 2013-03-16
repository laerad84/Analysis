void DrawAspect(){
  gStyle->SetOptStat(0);
  TFile* tf = new TFile("DistributionTest.root");
  TH1D* hisKLZAcceptance[2][10];
  TGraphErrors* grAspect[2][10];
  TGraphErrors* grRatio[2][10];
  TH1D* hisRatio[2][10];
  for( int i  =0 ; i< 2; i++){
    for( int j = 0; j< 10; j++ ){
      hisKLZAcceptance[i][j] = (TH1D*)tf->Get(Form("hisKLZAcceptance_%d_%d",i,j));
      hisKLZAcceptance[i][j]->SetLineColor(i+1);
      grAspect[i][j]         = new TGraphErrors();
      grAspect[i][j]->SetNameTitle(Form("grAspect_%d_%d",i,j),Form("grAspect_%d_%d",i,j));
      grAspect[i][j]->SetMarkerStyle(20+i);
      grRatio[i][j]         = new TGraphErrors();
      grRatio[i][j]->SetNameTitle(Form("grRatio_%d_%d",i,j),Form("grRatio_%d_%d",i,j));
      grRatio[i][j]->SetMarkerStyle(20+i);
      hisRatio[i][j] = new TH1D(Form("hisRatio_%d_%d",i,j),Form("hisRatio_%d_%d",i,j),70,0,7000);
      //ihisKlongPosZ[i][j]->Scale(1/hisKlongPosZ[i][j]->GetEntries());
      hisKLZAcceptance[i][j]->SetLineWidth(2);
    }
  }

  TCanvas* can  = new TCanvas("can" ,"can",1000,800);
  can->cd();
  TPad*    pad1= new TPad("pad1","pad1",0,0,1,0.25);
  TPad*    pad2= new TPad("pad2","pad2",0,0.25,1,1);
  pad1->SetMargin(0.1,0.1,0.1,0.0);
  pad2->SetMargin(0.1,0.1,0.005,0.1);  
  pad1->Draw();
  pad2->Draw();
  pad2->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  Double_t ScaleFactor[10]={0};
  
  for( int j = 0; j < 10; j++){
    ScaleFactor[j] = hisKLZAcceptance[1][j]->Integral()/hisKLZAcceptance[0][j]->Integral();
    hisKLZAcceptance[0][j]->Scale(ScaleFactor[j]);
    hisKLZAcceptance[1][j]->SetMarkerStyle(21);
    hisKLZAcceptance[1][j]->SetMarkerColor(2);
    hisRatio[1][j]->SetMarkerStyle(21);
    hisRatio[1][j]->SetMarkerColor(2);
    for( int ibin = 0; ibin< hisKLZAcceptance[1][j]->GetNbinsX();ibin++){
      grAspect[1][j]->SetPoint(ibin,hisKLZAcceptance[1][j]->GetBinCenter(ibin+1),hisKLZAcceptance[1][j]->GetBinContent(ibin+1));
      grAspect[1][j]->SetPointError(ibin,0,TMath::Sqrt(hisKLZAcceptance[1][j]->GetBinContent(ibin+1)));
      if( hisKLZAcceptance[0][j]->GetBinContent(ibin+1)!=0){
	hisRatio[1][j]->SetBinContent( ibin+1,hisKLZAcceptance[1][j]->GetBinContent(ibin+1)/hisKLZAcceptance[0][j]->GetBinContent(ibin+1));
	Double_t Error0 = hisKLZAcceptance[1][j]->GetBinError(ibin+1)/hisKLZAcceptance[0][j]->GetBinContent(ibin+1);
	Double_t Error1 = hisKLZAcceptance[0][j]->GetBinError(ibin+1)/hisKLZAcceptance[0][j]->GetBinContent(ibin+1);

	hisRatio[1][j]->SetBinError( ibin+1, TMath::Sqrt(Error0*Error0+Error1*Error1));
	grRatio[1][j]->SetPoint( grRatio[1][j]->GetN(),hisKLZAcceptance[1][j]->GetBinCenter(ibin+1),hisKLZAcceptance[1][j]->GetBinContent(ibin+1)/hisKLZAcceptance[0][j]->GetBinContent(ibin+1));
	
      }
    }
  }
  for( int i = 0; i< 10 ; i++){
    
    pad2->cd();
    hisKLZAcceptance[1][i]->Draw("ep");
    //grAspect[1][4]->Draw("P");
    hisKLZAcceptance[0][i]->Draw("same");
    pad1->cd();
    gPad->SetGridx();
    gPad->SetGridy();
    hisRatio[1][i]->GetYaxis()->SetRangeUser(0.5,1.5);
    hisRatio[1][i]->GetYaxis()->SetLabelSize(0.1);
    hisRatio[1][i]->GetXaxis()->SetLabelSize(0.1);
    hisRatio[1][i]->Draw("ep");
    
    can->Update();
    can->Modified();
    can->Update();
    can->Print(Form("Aspect_%d.pdf",i));
    getchar();
  }
  //grRatio[1][4]->Draw("AP");
  
  //hisKLZAcceptance[1][4]->Divide(hisKLZAcceptance[0][4]);
  //hisKLZAcceptance[1][4]->SetMarkerStyle(21);
  //hisKLZAcceptance[1][4]->Draw("ep");


}
