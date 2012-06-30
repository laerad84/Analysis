void compare(){

  TFile* tf = new TFile("RootData/CosmicOutput.root");
  TGraphErrors* grOld = (TGraphErrors*)tf->Get("grCosmicOutput3381");
  TFile* tf1 = new TFile("test_peak3508.root");
  TGraphErrors* grNew = (TGraphErrors*)tf1->Get("gain");
  
  std::cout<< grOld->GetN() << std::endl;
  std::cout<< grNew->GetN() << std::endl;

  TGraph* gr = new TGraph();
  for( int i = 0; i< 2716; i++){
    if( grOld->GetY()[i] >0 ){
      gr->SetPoint(gr->GetN(), grOld->GetY()[i] ,grNew->GetY()[i]);
    }
  }

  gr->Draw("AP");

}
  
