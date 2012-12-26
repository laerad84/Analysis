void DrawGraph(){
  
  TFile* tf = new TFile("TEMPLATE_COMPLETE_GRAPH_3pi0RunList_200_400.root");
  TGraph* gr[2716];
  TSpline3 spl[2716];
  for( int i = 0; i< 2716; i++){
    gr[i]  = (TGraph*)tf->Get(Form("Template_%d",i));
    spl[i] = new TSpline3(Form("spl%d",i), gr[i] ); 
  }
  
  TPostScript* ps = new TPostScript("Output.ps",111);
  TCanvas* can  = new TCanvas("can","",800,1024);
  can->Divide(2,3);
  
  Int_t nDiv = 6;
  for( int ich  = 0; ich < 2716; ich++){
    can->cd( (ich % nDiv) +1 );
    gr[i]->SetMarkerStyle(6);    
    gr[i]->SetMarkerColor(2);  
    //his->Draw("col");
    gr[i]->Draw("AP");
    spl[i]->Draw("same");
    if((ich % nDiv) == nDiv-1){
      can->Update();
      can->Modified();
      ps->NewPage();
    }
  }
}


