#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"

#include "TPostScript.h"

//void Merge(int iterationNumber){

int 
main( int argc, char** argv ){
  int iterationNumber  = 10;

  TFile* tf = new TFile("Calibration_Data/MergeOutput.root","recreate");

  TH2D*  his_klpos_delta_All=new TH2D("his_klpos_delta_All","recklz-simklz:recklz",
				      75,0,7500,100,-100,900);
  TH2D*  his_klpos_delta_Fit=new TH2D("his_klpos_delta_Fit","recklz-simklz:recklz",
				      75,0,7500,100,-100,900);
  TH2D*  his_klpos_delta_Mis=new TH2D("his_klpos_delta_Mis","recklz-simklz:recklz",
				      75,0,7500,100,-100,900);
  TH1D*  his_klpos_All = new TH1D("his_klpos_All","recklz",
				  75,0,7500);
  TH1D*  his_klpos_Fit = new TH1D("his_klpos_Fit","recklz",
				  75,0,7500);  
  TH1D*  his_klpos_Mis = new TH1D("his_klpos_Mis","recklz",
				  75,0,7500);

  TH1D*  his_delta_All = new TH1D("his_delta_All","recklz-simklz",100,-100,900);
  TH1D*  his_delta_Fit = new TH1D("his_delta_Fit","recklz-simklz",100,-100,900);
  TH1D*  his_delta_Mis = new TH1D("his_delta_Mis","recklz-simklz",100,-100,900);

  TH2D*  his_klpos_klChisqZ_All = new TH2D("his_klpos_klChisqZ_All","klChisqZ:recklz",
					   75,0,7500,100,0,100);
  TH2D*  his_klpos_klChisqZ_Fit = new TH2D("his_klpos_klChisqZ_Fit","klChisqZ:recklz",
					   75,0,7500,100,0,100);
  TH2D*  his_klpos_klChisqZ_Mis = new TH2D("his_klpos_klChisqZ_Mis","klChisqZ:recklz",
					   75,0,7500,100,0,100);

  TH1D*  his_klChisqZ_All  = new TH1D("his_klChisqZ_All","klChisqZ",100,0,100);
  TH1D*  his_klChisqZ_Fit  = new TH1D("his_klChisqZ_Fit","klChisqZ",100,0,100);
  TH1D*  his_klChisqZ_Mis  = new TH1D("his_klChisqZ_Mis","klChisqZ",100,0,100);
  TH1D*  his_SecklChisqZ_All  = new TH1D("his_SecklChisqZ_All","Second klChisqZ",100,0,100);
  TH1D*  his_SecklChisqZ_Fit  = new TH1D("his_SecklChisqZ_Fit","Second klChisqZ",100,0,100);
  TH1D*  his_SecklChisqZ_Mis  = new TH1D("his_SecklChisqZ_Mis","Second klChisqZ",100,0,100);


  TH1D* his_GammaE_All = new TH1D("his_GammaE_All","GammaE",100,0,2000);
  TH1D* his_GammaE_Fit = new TH1D("his_GammaE_Fit","GammaE",100,0,2000);
  TH1D* his_GammaE_Mis = new TH1D("his_GammaE_Mis","GammaE",100,0,2000);
  
  TH1D* his_klE_All = new TH1D("his_klE_All","klEnergy",80,0,8000);
  TH1D* his_klE_Fit = new TH1D("his_klE_Fit","klEnergy",80,0,8000);
  TH1D* his_klE_Mis = new TH1D("his_klE_Mis","klEnergy",80,0,8000);
  
  TH1D* his_klMass_All = new TH1D("his_klMass_All","klMass",200,400,600);
  TH1D* his_klMass_Fit = new TH1D("his_klMass_Fit","klMass",200,400,600);
  TH1D* his_klMass_Mis = new TH1D("his_klMass_Mis","klMass",200,400,600);

  TH1D* his_klDeltaChisqZ_All = new TH1D("his_klDeltaChisqZ_All","klDeltaChisqZ",100,0,100);
  TH1D* his_klDeltaChisqZ_Fit = new TH1D("his_klDeltaChisqZ_Fit","klDeltaChisqZ",100,0,100);
  TH1D* his_klDeltaChisqZ_Mis = new TH1D("his_klDeltaChisqZ_Mis","klDeltaChisqZ",100,0,100);
  
  /*
  TH1D* his_GammaR_All = new TH1D("his_GammaR_All","R of Gamma Position",100,0,1000);
  TH1D* his_GammaR_Fit = new TH1D("his_GammaR_Fit","R of Gamma Position",100,0,1000);
  TH1D* his_GammaR_Mis = new TH1D("his_GammaR_Mis","R of Gamma Position",100,0,1000);
  */

  std::string Filename = "Calibration_Data/CalibrationADV_%04d_%d.root";
  std::string FileList = "Calibration_Data/KLRunList_2.txt";
  std::vector<int> listVec;
  int runNum;
  std::ifstream ifs(FileList.c_str());
  while( ifs >>runNum){
    listVec.push_back(runNum);
  }

  TChain* ch = new TChain("trCalibration");
  for( int i = 0; i < listVec.size(); ++i){
    ch->Add(Form(Filename.c_str(),listVec[i],iterationNumber));
  }
  std::cout <<ch->GetEntries() << std::endl;
  std::cout<< __LINE__<< std::endl;
  ch->Project(his_klpos_delta_All->GetName(),"KlongPos[0][2]-SimKLVtx[2]:KlongPos[0][2]","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klpos_delta_Fit->GetName(),"KlongPos[0][2]-SimKLVtx[2]:KlongPos[0][2]","KlongFit==1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klpos_delta_Mis->GetName(),"KlongPos[0][2]-SimKLVtx[2]:KlongPos[0][2]","KlongFit!=1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  std::cout<< __LINE__<< std::endl;

  ch->Project(his_delta_All->GetName(),"KlongPos[0][2]-SimKLVtx[2]","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_delta_Fit->GetName(),"KlongPos[0][2]-SimKLVtx[2]","KlongFit==1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_delta_Mis->GetName(),"KlongPos[0][2]-SimKLVtx[2]","KlongFit!=1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  std::cout<< __LINE__<< std::endl;
  
  ch->Project(his_klpos_All->GetName(),"KlongPos[0][2]","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klpos_Fit->GetName(),"KlongPos[0][2]","KlongFit==1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klpos_Mis->GetName(),"KlongPos[0][2]","KlongFit!=1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");

  std::cout<< __LINE__<< std::endl;

  ch->Project(his_klpos_klChisqZ_All->GetName(),"KlongChisqZ[0]:KlongPos[0][2]","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klpos_klChisqZ_Fit->GetName(),"KlongChisqZ[0]:KlongPos[0][2]","KlongFit==1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klpos_klChisqZ_Mis->GetName(),"KlongChisqZ[0]:KlongPos[0][2]","KlongFit!=1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  std::cout<< __LINE__<< std::endl;

  ch->Project(his_klChisqZ_All->GetName(),"KlongChisqZ[0]","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klChisqZ_Fit->GetName(),"KlongChisqZ[0]","KlongFit==1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klChisqZ_Mis->GetName(),"KlongChisqZ[0]","KlongFit!=1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  std::cout<< __LINE__<< std::endl;

  ch->Project(his_SecklChisqZ_All->GetName(),"KlongChisqZ[1]","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_SecklChisqZ_Fit->GetName(),"KlongChisqZ[1]","KlongFit==2&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_SecklChisqZ_Mis->GetName(),"KlongChisqZ[1]","KlongFit!=2&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  std::cout<< __LINE__<< std::endl;

  ch->Project(his_klE_All->GetName(),"KlongE","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klE_Fit->GetName(),"KlongE","KlongFit==1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klE_Mis->GetName(),"KlongE","KlongFit!=1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  std::cout<< __LINE__<< std::endl;
  
  ch->Project(his_GammaE_All->GetName(),"GammaE","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_GammaE_Fit->GetName(),"GammaE","KlongFit==1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_GammaE_Mis->GetName(),"GammaE","KlongFit!=1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  std::cout<< __LINE__<< std::endl;
  ch->Project(his_klMass_All->GetName(),"KlongMass[0]","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klMass_Fit->GetName(),"KlongMass[0]","KlongFit==1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klMass_Mis->GetName(),"KlongMass[0]","KlongFit!=1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  
  ch->Project(his_klDeltaChisqZ_All->GetName(),"KlongChisqZ[1]-KlongChisqZ[0]","GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klDeltaChisqZ_Fit->GetName(),"KlongChisqZ[1]-KlongChisqZ[0]","KlongFit==1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");
  ch->Project(his_klDeltaChisqZ_Mis->GetName(),"KlongChisqZ[1]-KlongChisqZ[0]","KlongFit!=1&&GammaE[0]>200&&GammaE[1]>200&&GammaE[2]>200&&GammaE[3]>200&&GammaE[4]>200&&GammaE[5]>200");

  his_klpos_delta_All->Write();
  his_klpos_delta_Fit->Write();
  his_klpos_delta_Mis->Write();

  his_klpos_All->Write();
  his_klpos_Fit->Write();
  his_klpos_Mis->Write();

  his_delta_All->Write();
  his_delta_Fit->Write();
  his_delta_Mis->Write();
  
  his_klpos_klChisqZ_All->Write();
  his_klpos_klChisqZ_Fit->Write();
  his_klpos_klChisqZ_Mis->Write();
  
  his_klChisqZ_All->Write();
  his_klChisqZ_Fit->Write();
  his_klChisqZ_Mis->Write();

  his_SecklChisqZ_All->Write();
  his_SecklChisqZ_Fit->Write();
  his_SecklChisqZ_Mis->Write();

  his_klE_All->Write();
  his_klE_Fit->Write();
  his_klE_Mis->Write();
  
  his_GammaE_All->Write();
  his_GammaE_Fit->Write();
  his_GammaE_Mis->Write();

  his_klMass_All->Write();
  his_klMass_Fit->Write();
  his_klMass_Mis->Write(); 

  his_klDeltaChisqZ_All->Write();
  his_klDeltaChisqZ_Fit->Write();
  his_klDeltaChisqZ_Mis->Write(); 
  

  TPostScript* ps  = new TPostScript("Merge.ps",111);
  TCanvas*     can = new TCanvas("can","",800,1200);
  can->Divide(2,3);
  ps->NewPage();
  can->cd(1);
  his_klpos_delta_All->Draw("colz");
  can->cd(3);
  his_klpos_delta_Fit->Draw("colz");
  can->cd(5);
  his_klpos_delta_Mis->Draw("colz");
  can->cd(2);
  his_klpos_klChisqZ_All->Draw("colz");
  can->cd(4);
  his_klpos_klChisqZ_Fit->Draw("colz");
  can->cd(6);
  his_klpos_klChisqZ_Mis->Draw("colz");
  
  can->Modified();
  can->Update();
  ps->NewPage();
  
  can->cd(1);
  his_klpos_All->SetLineColor(1);
  his_klpos_Fit->SetLineColor(2);
  his_klpos_Mis->SetLineColor(3);
  his_klpos_All->Draw();
  his_klpos_Fit->Draw("same");
  his_klpos_Mis->Draw("same");
  gPad->SetLogy();

  can->cd(2);
  his_delta_All->SetLineColor(1);
  his_delta_Fit->SetLineColor(2);
  his_delta_Mis->SetLineColor(3);
  his_delta_All->Draw();
  his_delta_Fit->Draw("same");
  his_delta_Mis->Draw("same");
  gPad->SetLogy();

  can->cd(3);
  his_klE_All->SetLineColor(1);
  his_klE_Fit->SetLineColor(2);
  his_klE_Mis->SetLineColor(3);
  his_klE_All->Draw();
  his_klE_Fit->Draw("same");
  his_klE_Mis->Draw("same");
  gPad->SetLogy();

  can->cd(4);
  his_GammaE_All->SetLineColor(1);
  his_GammaE_Fit->SetLineColor(2);
  his_GammaE_Mis->SetLineColor(3);
  his_GammaE_All->Draw();
  his_GammaE_Fit->Draw("same");
  his_GammaE_Mis->Draw("same");
  gPad->SetLogy();
  
  can->cd(5);
  his_klChisqZ_All->SetLineColor(1);
  his_klChisqZ_Fit->SetLineColor(2);
  his_klChisqZ_Mis->SetLineColor(3);
  his_klChisqZ_All->Draw();
  his_klChisqZ_Fit->Draw("same");
  his_klChisqZ_Mis->Draw("same");
  gPad->SetLogy();

  can->cd(6);
  his_SecklChisqZ_All->SetLineColor(1);
  his_SecklChisqZ_Fit->SetLineColor(2);
  his_SecklChisqZ_Mis->SetLineColor(3);
  his_SecklChisqZ_All->Draw();
  his_SecklChisqZ_Fit->Draw("same");
  his_SecklChisqZ_Mis->Draw("same");
  gPad->SetLogy();
 
  can->Modified();
  can->Update();

  ps->NewPage();
  can->cd(1);
  his_klMass_All->SetLineColor(1);
  his_klMass_Fit->SetLineColor(2);
  his_klMass_Mis->SetLineColor(3);
  his_klMass_All->Draw();
  his_klMass_Fit->Draw("same");
  his_klMass_Mis->Draw("same");
  gPad->SetLogy();

  can->cd(2);
  his_klDeltaChisqZ_All->SetLineColor(1);
  his_klDeltaChisqZ_Fit->SetLineColor(2);
  his_klDeltaChisqZ_Mis->SetLineColor(3);
  his_klDeltaChisqZ_All->Draw();
  his_klDeltaChisqZ_Fit->Draw("same");
  his_klDeltaChisqZ_Mis->Draw("same");
  gPad->SetLogy();

  can->Modified();
  can->Update();
  ps->Close();
  tf->Close();
  return 0;
}
    
  
