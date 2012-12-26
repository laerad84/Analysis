void testDraw(){
  
  Int_t LowLimit = 250;
  Int_t HighLimit = 500;
  //TFile* tf= new TFile("TEMPLATE_HEIGHT_500_1000.root");
  TFile* tf = new TFile(Form( "TEMPLATE_SPLINE_%d_%d.root",LowLimit,HighLimit));
  TH1D*  hisPoints = new TH1D("hisPoints","points",500,0,500);
  TGraph* gr[2716];
  for( int i = 0; i< 2716; i++){

    gr[i] = (TGraph*)tf->Get(Form("Template_graph_%d", i ));
    if( gr[i]->GetN() == 0 ){ continue; }
    hisPoints->Fill( gr[i]->GetN() ); 
    if( gr[i]->GetN() < 450 ) { 
      std::cout << i << " : " 
		<< gr[i]->GetN()  << std::endl;}
  }
  TCanvas* can  = new TCanvas("can","",800,800);
  can->Divide(2,2);
  can->cd(1);
  gr[0]->Draw("AP");
  for( int i = 0; i< 2716; i++){
    if( gr[i]->GetN() == 0 ){ continue;}
    if( i ==  210 || i == 860 ) { continue; }
    gr[i]->Draw("P");
  }
  
  can->cd(2);
  hisPoints->Draw();

}
				
