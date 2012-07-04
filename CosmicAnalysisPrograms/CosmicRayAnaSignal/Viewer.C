void Viewer(){

  TFile* tf = new TFile("TestTimeOutLogic_3897.root");
  TTree* tr = (TTree*)tf->Get("TimeDeltaCosmic");
  Int_t IDFirst;
  Int_t IDSecond;
  Double_t Mean;
  Double_t RMS;
  Double_t Error;
  Int_t Entries;
  tr->SetBranchAddress("IDFirst",&IDFirst);
  tr->SetBranchAddress("IDSecond",&IDSecond);
  tr->SetBranchAddress("Mean",&Mean);
  tr->SetBranchAddress("RMS",&RMS);
  tr->SetBranchAddress("Error",&Error);
  tr->SetBranchAddress("Entries",&Entries);

  TTree* trRead = new TTree("ReadFile","");
  trRead->ReadFile("testNewWORKCompileOffset.txt","CHID/I:Offset/D:OffsetError/D");
  
  int CHID;
  double Offset;
  trRead->SetBranchAddress("CHID",&CHID);
  trRead->SetBranchAddress("Offset",&Offset);

  const int nCH = 2716;
  int CHIDList[nCH]={-1};
  double OffsetList[nCH]={0};
  int CHFlag[nCH] = {0};
  if( trRead->GetEntries() != nCH ){
    std::cout<< "Error" << std::endl; 
    return ;
  }

  for( int i = 0; i< 2716; i++){    
    trRead->GetEntry(i);
    CHIDList[CHID]   = CHID;
    OffsetList[CHID] = Offset;    
  }
  
  TH1D* hisDistribution0 = new TH1D("hisDistribution0","",400,-50,50);
  TH1D* hisDistribution  = new TH1D("hisDistribution","",400,-50,50);



  TH1D* his_before[4];
  TH1D* his_after[4];

  char *Tag[4] = {"All","Small:Small","Small:Large","Large:Large"};

  for( int i = 0; i< 4; i++){
    his_before[i] = new TH1D(Form("his_before_%d",i),Tag[i],400,-50,50);
    his_after[i]  = new TH1D(Form("his_after_%d",i) ,Tag[i],400,-50,50);
    his_before[i]->SetLineColor(i+1);
    his_after[i]->SetLineColor(i+1);
  }

  for( int i = 0; i< tr->GetEntries();i++ ){
    tr->GetEntry(i);
    if( Entries <= 100 || Error >=5 ){ continue; }
    int hisIndex = 0; 
    
    his_before[hisIndex]->Fill(Mean);
    his_after[hisIndex]->Fill( Mean -OffsetList[IDFirst] + OffsetList[IDSecond]);
    
    if( IDFirst < 2240  && IDSecond < 2240 ){
      hisIndex  = 1;
    }else if( IDFirst < 2240 && IDSecond >= 2240 ){
      hisIndex  = 2; 
    }else if( IDFirst >= 2240 && IDSecond >= 2240 ){
      hisIndex  = 3;
    }else{
      continue;
    } 
    
    his_before[hisIndex]->Fill(Mean);
    his_after[hisIndex]->Fill( Mean -OffsetList[IDFirst] + OffsetList[IDSecond]);
    
  }


  TCanvas* can = new TCanvas( "can", "", 800,400 );
  can->Divide(2,1);
  can->cd(1);
  his_before[0]->Draw();
  his_before[1]->Draw("same");
  his_before[2]->Draw("same");
  his_before[3]->Draw("same");
  

  gPad->SetLogy();
  gPad->SetGridx();  
  gPad->SetGridy();  
  can->cd(2);
  his_after[0]->Draw();
  his_after[1]->Draw("same");
  his_after[2]->Draw("same");
  his_after[3]->Draw("same");

  gPad->SetLogy();
  gPad->SetGridx();  
  gPad->SetGridy();  

}
  
  
 
  
