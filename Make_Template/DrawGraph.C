void DrawGraph(){
  
  TFile* tf = new TFile("TEMPLATE_COMPLETE_GRAPH_3pi0RunList_200_400.root");
  TGraph* gr[2716];
  for( int i = 0; i< 2716; i++){
    gr[i] = (TGraph*)tf->Get(Form("Template_%d",i));
  }
  
  TCanvas* can  = new TCanvas("can","",800,800);
  gr[0]->Draw("AP");
  for( int i = 1; i< 2716; i++){
    gr[i]->Draw("P");
  }
  
}


