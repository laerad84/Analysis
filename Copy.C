void Copy(){

  TFile* tfin = new TFile("/Volume0/ExpData/2012_Feb_Beam/RootFile_cosmic/CosmicResult_20120209.root");

  TFile* tfOut = new TFile("test.root","Recreate");
  TTree* tr   = (TTree*)tfin->Get("GainFitPar");
  TTree* trClone = tr->CopyTree("");
  trClone->Write();
  tfOut->Close();

}

  

  
