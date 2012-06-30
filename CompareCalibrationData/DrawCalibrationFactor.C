#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
void DrawCalibrationFactor(){
  gStyle->SetPalette(1);
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("neRMI");

  char* ANALIBDIR = std::getenv("ANALYSISLIB");
  char* libFile   = Form("%s/lib/libAnalysisLib.so",ANALIBDIR);
  gSystem->Load(libFile);
  IDHandler* handler = new IDHandler();
  CsIImage*  image3pi= new CsIImage(handler);
  CsIImage*  imageke3= new CsIImage(handler);
  CsIImage*  imageRat= new CsIImage(handler);
  CsIImage*  imageCUT= new CsIImage(handler);
  CsIImage*  imageSTAT = new CsIImage(handler);

  const int nCSI = 2716;
  double Cal3pi[nCSI]={0};
  double Cal3piStatics[nCSI]={0};
  double Calke3[nCSI]={0};

  std::string Cal3piFile="Data/CalibrationFactor_4.dat";
  std::string Calke3File="Data/ke3CalibConst.dat";
  std::string Cal3piStaticsFile = "Data/CalibrationStatics_4.dat";

  std::ifstream fsCal3pi(Cal3piFile.c_str());
  std::ifstream fsCalke3(Calke3File.c_str());
  std::ifstream fsCal3pistatics(Cal3piStaticsFile.c_str());

  if(!fsCal3pi.is_open()){return;}
  if(!fsCalke3.is_open()){return;}
  if(!fsCal3pistatics.is_open()){return;}
  
  TH1D* histRatio_RAW = new TH1D("Cal_Ratio_RAW","CalRatio;ke3/3pi0",100,0,2);
  TH1D* histRatio_CUT = new TH1D("Cal_Ratio_CUT","CalRatio;ke3/3pi0",100,0,2);
  TH1D* histCal_3pi   = new TH1D("Cal_2pi","Cal_3pi",100,0,2);
  TH1D* histCal_ke3   = new TH1D("Cal_ke3","Cal_ke3",100,0,2);

  int id;
  double cal;
  while( fsCal3pi >> id >> cal ){
    Cal3pi[id] = cal;
    image3pi->Fill(id,cal);
    histCal_3pi->Fill(cal);
  }
  while( fsCalke3 >> id >> cal ){
    Calke3[id] = cal;
    imageke3->Fill(id,cal);
    histCal_ke3->Fill(cal);
  }
  double statics;
  while( fsCal3pistatics >> id >> statics ){
    Cal3piStatics[id] = statics;
    imageSTAT->Fill(id, statics);
  }

  TGraph* gr  = new TGraph();
  TGraph* grStat = new TGraph();
  for( int ich = 0; ich < nCSI; ich++){
    if( Cal3pi[ich] <= 0 || Calke3[ich] <= 0 ){
      continue;
    }
    double x,y; 
    handler->GetMetricPosition( ich , x, y);
    double R = TMath::Sqrt( x*x+y*y);
    double ax = TMath::Abs(x);
    double ay = TMath::Abs(y);
    double calRatio = ((Calke3[ich]-Cal3pi[ich])/Cal3pi[ich]+1);

    gr->SetPoint(gr->GetN(), Cal3pi[ich], Calke3[ich]);
    imageRat->Fill(ich,calRatio);
    histRatio_RAW->Fill(calRatio);
    if( ax < 150 && ay < 150){
      continue;
    }
    if( R > 850 ){
      continue;
    }
    if( ay > 550){
      continue;
    }
    grStat->SetPoint(grStat->GetN(),Cal3piStatics[ich],calRatio);
    histRatio_CUT->Fill(calRatio);
    imageCUT->Fill(ich,calRatio);
  }
  TCanvas* can  = new TCanvas("can","",1200,900);
  can->Divide(4,3);
  can->cd(1);
  gr->SetNameTitle("ke3_VS_3pi","ke3_vs_3pi;3pi0;ke3");
  gr->GetXaxis()->SetRangeUser(0.8,1.4);
  gr->GetYaxis()->SetRangeUser(0.8,1.4);
  gr->Draw("AP");
  can->cd(2);
  grStat->SetNameTitle("(ke3-3pi0)/3pi0","(ke3-3pi0)/3pi0_fiducialRegion;nentries;(ke3-3pi0)/3pi0");
  grStat->Draw("AP");

  can->cd(3);
  //imageRat->DrawWithRange("colz",0.5,1.5);
  imageSTAT->SetTitle("Statics 3pi0");
  imageSTAT->Draw();
  gPad->SetLogz();

  can->cd(5);
  image3pi->SetTitle("3pi0CalibrationFactor");
  image3pi->DrawWithRange("colz",0.8,1.4);
  can->cd(9);
  imageke3->SetTitle("ke3CalibrationFactor");
  imageke3->DrawWithRange("colz",0.8,1.4);

  can->cd(6);
  histCal_3pi->Draw();
  can->cd(10);
  histCal_ke3->Draw();

  can->cd(7);
  imageRat->SetTitle("Ratio");
  imageRat->DrawWithRange("colz",0.8,1.4);
  can->cd(11);
  imageCUT->SetTitle("Ratio_Fiducial_region");
  imageCUT->DrawWithRange("colz",0.8,1.4);


  can->cd(8);
  histRatio_RAW->SetLineColor(1);
  histRatio_RAW->Draw();
  histRatio_RAW->Fit("gaus","","",0.9,1.1); 
  histRatio_RAW->GetFunction("gaus")->SetLineWidth(1);
  can->cd(12);
  histRatio_CUT->SetLineColor(2);
  histRatio_CUT->Fit("gaus","","",0.9,1.1);
  histRatio_CUT->Draw();
  histRatio_CUT->GetFunction("gaus")->SetLineWidth(1);

  
  /*
  histCal_3pi->SetLineColor(4);
  histCal_ke3->SetLineColor(6);
  histCal_3pi->Draw("same");
  histCal_ke3->Draw("same");
  */

}


	 





  
  
  
  
