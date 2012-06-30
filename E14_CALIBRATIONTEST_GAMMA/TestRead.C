void TestRead(){
  
  static const int nMaxTracks = 400;

  TFile* tf = new TFile("SimRoot/Conv_KL3pi0.mac_1000000_0.root");
  TTree* tr = (TTree*)tf->Get("T");

  Int_t   nTrack;
  Float_t ek[nMaxTracks];
  Double_t v[nMaxTracks][3];//nTrack
 
  tr->SetBranchAddress("nTrack",&nTrack);
  tr->SetBranchAddress("ek",ek);//nTrack
  tr->SetBranchAddress("v",v);//nTrack
  
  for( int i = 0; i< 10; ++i ){
    tr->GetEntry(i);
    std::cout<<"nTrack: "<< nTrack << std::endl;
    for( int j = 0; j < nTrack; ++j){
      std::cout << v[j][0] << "\t";
    }
    std::cout<< std::endl;
  }
}

