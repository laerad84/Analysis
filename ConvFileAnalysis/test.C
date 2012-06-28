#include <cstdlib>
#include <cstdio>
#include <iostream>

void test(){
  gStyle->SetOptFit(111111111);
  gSystem->Load("/home/jwlee/local/Analysis/DrawLib/lib/libDrawLib.so");

  const int cosmicArr[20] = {0 ,5 ,2 ,1 ,6 ,7 ,4 ,13,12,11,
  			     14,3 ,10,15,8 ,9 ,16,17,18,19};

  //const int cosmicArr[20] = {0 ,3 ,2 ,11,6 ,1 ,4 ,5 ,14,15,
  //12,9 ,8 ,7 ,10,13,16,17,18,19};
  const double OffSet[20]={
    26.6396 ,9.77932 ,7.28189 ,-25.9666,11.6498 ,-10.0956,-15.5716,11.8809 ,13.5989 ,10.4982,
    -2.25069,-22.9987,-19.1044,0       ,-20.5825,-10.2814,-16.9232,11.5399 ,-8.11825, 11.733};
  
  int runNumber = 3949;
  int ch1   = 1599;
  int ch2   = 1600;
  TFile* tf = new TFile(Form("run%d_wav.root", runNumber));
  TTree* tr = (TTree*)tf->Get("WFTree");  
  std::string sumDir= std::getenv("SUMUPFILEDIR");
  
  int CsiNumber;
  int CsiID[4096];
  double CsiHHTiming[4096];
  double CsiSignal[4096];
  int CC03Number;
  int CC03ID[4096];
  double CC03HHTiming[4096];
  double CC03Signal[4096];
  int OEVNumber;
  int OEVID[4096];
  double OEVHHTiming[4096];
  double OEVSignal[4096];

  int CosmicNumber;
  int CosmicID[4096];
  double CosmicHHTiming[4096];
  double CosmicSignal[4096];

  /// Set Address 


  {
  
    tr->SetBranchAddress("CsiNumber",&CsiNumber);
    tr->SetBranchAddress("CsiID",CsiID);//CsiNumber
    tr->SetBranchAddress("CsiTiming",CsiHHTiming);//CsiNumber
    tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
    
    tr->SetBranchAddress("CC03Number",&CC03Number);
    tr->SetBranchAddress("CC03ID",CC03ID);//CC03Number
    tr->SetBranchAddress("CC03Timing",CC03HHTiming);//CC03Number
    tr->SetBranchAddress("CC03Signal",CC03Signal);//CC03Number
    
    tr->SetBranchAddress("OEVNumber",&OEVNumber);
    tr->SetBranchAddress("OEVID",OEVID);//OEVNumber
    tr->SetBranchAddress("OEVTiming",OEVHHTiming);//OEVNumber
    tr->SetBranchAddress("OEVSignal",OEVSignal);//OEVNumber
    
    tr->SetBranchAddress("CosmicNumber",&CosmicNumber);
    tr->SetBranchAddress("CosmicID",CosmicID);//CosmicNumber
    tr->SetBranchAddress("CosmicTiming",CosmicHHTiming);//CosmicNumber
    tr->SetBranchAddress("CosmicSignal",CosmicSignal);//CosmicNumber
    
  }

  CC03_Module* cc03 = new CC03_Module("cc03");
  CsI_Module*   csi = new CsI_Module("csi");  
  OEV_Module*   oev = new OEV_Module("oev");  
  

  TH2D* hisCoincidency[20];
  for( int i = 0; i<20; i++){
    hisCoincidency[i]= new TH2D(Form("coincidency%d",i),"",
				20,0,20,500,-250,250);
  }

  TH1D* hisCoin[20];
  for( int i = 0; i< 20; i++){
    hisCoin[i] = new TH1D(Form("his%d",i),"",20,0,20);
  }

  TGraph* gr[400];
  for( int i  =0; i< 400; i++){
    gr[i] = new TGraph();
  }

  TH1D* hisDelta[20][20];
  for( int i = 0 ; i< 20; i++){
    for( int j = 0; j< 20; j++){
      hisDelta[i][j] = new TH1D(Form("hisDelta_%d_%d",i,j),Form("Delta of Channel %d and %d", i, j), 200,-100,100);
    }
  }

  TH2D* hisHit    = new TH2D("hisHit","",10,0,10,2,0,2);
  TCanvas* canvas = new TCanvas("canvas","",800,400);

  canvas->Divide(2,1);
  double channelTime[2716];
  double channelOut[2716];
  double cosmicTime[20];
  double cosmicOut[20];

  TH1D* hisSum   = new TH1D("hisSum"  ,"",500,0,1000);
  //for( int i = 0; i< tr->GetEntries(); i++){
  for( int i = 0; i< 20000; i++){

    for( int j = 0; j < 2716; j++){
      channelOut[j] = 0;
      channelTime[j] = -1; 
    }

    for( int j = 0; j< 20; j++){
      cosmicOut[j] = 0;
      cosmicTime[j]=-1;
    }

    tr->GetEntry(i); 
    
    if( CsiNumber > 500 ) continue;
    
    for( int j = 0; j< CsiNumber; j++){
      channelOut[CsiID[j]] = CsiSignal[j];
      channelTime[CsiID[j]] = CsiHHTiming[j];
    }
    for( int j = 0; j< CosmicNumber; j++){
      cosmicOut[cosmicArr[CosmicID[j]]] = CosmicSignal[j];
      cosmicTime[cosmicArr[CosmicID[j]]] = CosmicHHTiming[j];
    }
    
    if(cosmicOut[0] > 1000 && cosmicOut[10] > 1000 && cosmicOut[5] > 1000 && cosmicOut[15] > 1000 ){
      hisSum->Fill( cosmicTime[0] + cosmicTime[10] - cosmicTime[5] - cosmicTime[15]);
    }

    for( int iIndex = 0; iIndex < 20; iIndex++ ){
      if( cosmicOut[iIndex] < 1000){continue;}
      for( int jIndex = iIndex+1; jIndex < 20; jIndex++ ){
	if( cosmicOut[jIndex] < 1000 ){ continue; }
	hisDelta[iIndex][jIndex]->Fill( cosmicTime[iIndex] - cosmicTime[jIndex]- (OffSet[iIndex]-OffSet[jIndex]) );
      }
    }    
    //hisHit->Reset();
    for( int j = 0; j< CosmicNumber;j++){
      if( CosmicSignal[j] < 2000 ){ continue; }
      for( int k = 0; k < CosmicNumber; k++){
	if( CosmicSignal[k] < 2000 ){continue;}
	//if( CosmicID[j]==CosmicID[k] ){continue;}
	hisCoincidency[CosmicID[j]]->Fill(CosmicID[k],CosmicHHTiming[j]-CosmicHHTiming[k]);
	hisCoin[CosmicID[j]]->Fill(CosmicID[k]);
	gr[cosmicArr[CosmicID[j]]*20+cosmicArr[CosmicID[k]]]
	  ->SetPoint(gr[cosmicArr[CosmicID[j]]*20+cosmicArr[CosmicID[k]]]->GetN(),
		     CosmicHHTiming[j],CosmicHHTiming[k]-CosmicHHTiming[j]);
	
      }
    }

    
    
    //Draw 
    
    { 
      cc03->Reset();
      csi->Reset();
      oev->Reset();
      
      for( int j = 0; j < CsiNumber; j++){
	csi->Fill( CsiID[j], CsiSignal[j]);
      }
      
      for( int j = 0; j< CC03Number; j++){
	cc03->Fill( CC03ID[j], CC03Signal[j]);
      }
      
      for( int j = 0; j< OEVNumber; j++){
	oev->Fill( OEVID[j],OEVSignal[j]);
      }
      
      for( int j = 0; j< CosmicNumber; j++){
	hisHit->Fill(cosmicArr[CosmicID[j]]%10,cosmicArr[CosmicID[j]]/10);
      }
    }


    /*
    
    canvas->cd(1);
    csi->DrawWithRange(0,16000,"colz");
    oev->DrawWithRange(0,16000,"same colz");
    cc03->DrawWithRange(0,16000,"same colz");
    canvas->cd(2);
    hisHit->Draw("colz");
    canvas->Update();
    canvas->Modified();   
    */
    // getchar();
    
  }  
  
  TPostScript* ps = new TPostScript("test.ps",111);
  TCanvas* can = new TCanvas("can","",1600,2000);
  can->Divide(4,5);
  for( int i = 0; i< 20; i++){
    ps->NewPage();
    for( int j = 0; j< 20; j++){    
      can->cd( j+1 );
      hisDelta[i][j]->Draw();
      std::cout<< i << "\t" 
	       << j << "\t"
	       << hisDelta[i][j]->GetEntries()   << "\t" 
	       << hisDelta[i][j]->GetMean()      <<"\t"
	       << hisDelta[i][j]->GetMeanError() << std::endl;
    }
    can->Update();
    can->Modified();
  }

  ps->NewPage();
  ps->Close();
  TCanvas* canv = new TCanvas("canv","",800,400);
  canv->Divide(2,1);
  canv->cd(1);
  canv->cd(2);
  hisSum->Draw();
}
