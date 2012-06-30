#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
void Statics(){
  
  std::string str = std::getenv("ANALYSISLIB");
  str+= "/lib/libAnalysisLib.so";
  
  gSystem->Load(str.c_str());
  IDHandler* handler = new IDHandler();
  CsIImage* image = new CsIImage(handler);
  
  std::string filename =  "Calibration_Data/KLRunList_2.txt";  
  std::ifstream ifs(filename.c_str());
  std::vector<int> idVec;
  if( ifs.is_open() ){
    Int_t id;
    while ( ifs >> id ){
      idVec.push_back(id);
    }
  }else{ 
    std::cout << "Error" << std::endl;
    return ;
  }
  
  TChain* ch = new TChain("trCalibration");
  std::vector<int>::iterator it;
  for( it = idVec.begin();
       it != idVec.end();
       ++it){
    ch->Add(Form("Calibration_Data/Calibration_NoChange/CalibrationADV_%d_%d.root",*it,0));
  }
  
  TH1D* his 
    = new TH1D("his", "", 2716,0,2716);
  TH2D* his2D 
    = new TH2D("his2D", "" ,30,0,300,200,0,2);
  TH1D* hisKLRaw 
    = new TH1D("hisKlong","Reconstructed Klong Mass;Mass[MeV];N/1Mev",200,400,600);
  TH1D* hisKLUsed
    = new TH1D("hisMass6GammaRec",
	       "Reconstructed Mass with 6 Gamma Event;Mass[MeV];N/1Mev",
	       300,400,700);
  TH1D* hisKLCalibrated 
    = new TH1D("hisKlongCalibrated","Klong;Mass[MeV];N/1Mev",200,400,600);
  

  //ch->Project(his->GetName(),"CorrID","FlagCalibrated==0");
  //ch->Project(his2D->GetName(),"Corr:CorrE","FlagCalibrated==0");
  
  //  ch->Project(hisKLRaw->GetName(),"KlongMass[0]");
  ch->Project(hisKLUsed->GetName(),"KlongMass[0]");//,
  //	      "FlagKL_prefit==0&&(CutCondition&(1+8))==0");
  //ch->Project(hisKLCalibrated->GetName(),"KlongMass[0]",
  //	      "FlagKL_prefit==0&&(CutCondition&(1+8))==0&&nCalibrated>0");

  std::cout<< hisKLRaw->GetEntries() << std::endl;
  std::cout<< hisKLUsed->GetEntries() << std::endl;
  std::cout<< hisKLCalibrated->GetEntries() << std::endl;
  
  hisKLRaw->SetLineColor(1);
  hisKLUsed->SetLineColor(2);
  hisKLCalibrated->SetLineColor(4);
  

  for( int i = 0; i< 2716; ++i ){
    image->Fill(i, his->GetBinContent(i+1));
  }

  /*
  TCanvas* can  = new TCanvas("can","",1600,800);
  can->Divide(2,1);
  can->cd(1);
  image->Draw();
  can->cd(2);
  his2D->Draw("colz");
  TProfile* pro = his2D->ProfileX();
  pro->Draw("same");
  can->SaveAs("Statics.pdf");
  */

  TCanvas* canvas = new TCanvas("canvas","",800,800);
  //hisKLRaw->Draw();
  hisKLUsed->SetLineColor(1);
  hisKLUsed->Draw("");
  gPad->SetGridx();
  gPad->SetGridy();
  //hisKLCalibrated->Draw("same");

}

  
