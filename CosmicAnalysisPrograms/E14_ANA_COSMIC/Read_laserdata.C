void Read_laserData(int runNumber){

  TFile* tf = new TFile(Form("../sumdata/laserRootFile/laserdata_%04d.root",runNumber));
  
  TTree* tr  =new (TTree*)tf->Get("laserEventfData");
  
  Double_t Time;
  Double_t PINData[4];

  Double_t PINPed[4];
  Double_t PINOut[4];
  Double_t ADC[16];
  Double_t Scaler[16];

  tr->SetBranchAddress("T",&Time);
  tr->SetBranchAddress("PINData",&PINData);
  tr->SetBranchAddress("ADCData",&ADCData);
  tr->SetBranchAddress("Scaler",&Scaler);

  TH1D* hisPINout[4];
  TH1D* hisPINped[4];
  for( int i = 0; i <4; i++){
    hisPINout[i] = new TH1D(Form("hisPIN%d",i),Form("hisPIN%d",i),
			    1000,0,5000);
    hisPINped[i] = new TH1D(Form("hisPINped%d",i),Form("hisPINped%d",i),
			    1000,0,5000);

  }
  
  




  

  
  
  
