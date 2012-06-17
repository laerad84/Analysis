#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"

#include "TF2.h"
#include "E14ReadSumFile.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "HoughCsI.h"
#include "Chisq_cosmic.h"

int
main(int argc, char** argv){
  gStyle->SetPalette(1);
  // ARGV[0]  <InputROOTFile> <OutputROOTFile> [<Vision>]
  std::string InputFile;
  std::string OutputFile;
  
  if( argc == 3){
    InputFile = argv[1];
    OutputFile = argv[2];
  }else{
    std::cerr << "ARGUEMENT ERROR" << std::endl;
    return -1;
  }
  
  //
  TApplication* app = new TApplication("app",&argc, argv);
  
 
  E14ReadSumFile* Reader = new E14ReadSumFile();
  Reader->Add(InputFile.c_str());
  
  IDHandler* handler = new IDHandler("Data/crystal.txt");
  HoughCsI*  hough   = new HoughCsI();
  Chisq_cosmic* chi2Cosmic = new Chisq_cosmic();
  
  CsIImage*  image   = new CsIImage(handler);
  CsIImage*  image1  = new CsIImage(handler);
  TCanvas*   can     = new TCanvas("can","",800,0,800,800);
  can->Divide(2,2);
  
  
  long nEntries = Reader->GetEntries();
  std::cout << "TotalEntries::" << nEntries  << std::endl;
  
  TH2D*    hisHough;
  TF1*     funcHough;
  TLine*   line;
  TLine*   lineChi;
  TF2*     func;
  Double_t roh;
  Double_t theta;
  Int_t    HitFlag[N_TOTAL_CSI]={};
  Int_t    HitBit;
  
  for( int iEntry = 0; iEntry < nEntries; iEntry++){
    Reader->GetEntry(iEntry);
    HitBit=0;
    ///AnalysisCode    
    TGraph* gr = new TGraph();
    hough->Reset();
    image->Reset();
    image1->Reset();
    
    for( int iCosmic = 0; iCosmic < Reader->CosmicNumber; iCosmic++){
    }
    std::cout <<"Run Analysis "<< std::endl;
    for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
      if( Reader->CsiEne[idigi] > 3 ){
	double x,y;
	handler->GetMetricPosition(Reader->CsiModID[idigi], x, y);
	image->Fill(Reader->CsiModID[idigi]);
	gr->SetPoint( gr->GetN(), x,y);	
      }
    }    
    if( gr->GetN() > 1000 ){
      for( int idigi = 0; idigi < Reader->CsiNumber; idigi++){
	if( Reader->CsiEne[idigi] < 300 ){
	  std::cout << "ID::" << Reader->CsiModID[idigi] <<"\t" <<  Reader->CsiEne[idigi] << std::endl;
	}
      }
    }
    /*
      Int_t findList[46]=
    {1664,1665,1712,1713,1636,1637,1638,1639,1684,1685,
    1686,1687,1732,1733,1734,1735,1780,1781,1782,1783,
    1828,1829,1830,1831,1876,1877,1878,1879,1924,1925,
    1926,1927,2070,2071,2118,2119,2118,2119,2166,2167,
    2214,2215,2138,2139,2186,2187};
    */
    
    if(hough->CosmicJudgment(gr)){
      roh = hough->GetRoh();
      theta = hough->GetTheta();
      for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
	if( Reader->CsiEne[idigi] > 3 ){
	  double x,y;
	  //std::cout << "ID::" << Reader->CsiModID[idigi] << std::endl;
	  handler->GetMetricPosition(Reader->CsiModID[idigi], x, y);
	  if( hough->CalDistance(x,y) <= 50 ){
	    image1->Fill(Reader->CsiModID[idigi]);
	    gr->SetPoint( gr->GetN(), x,y);
	  }	  
	}
      }    
    }
    
    chi2Cosmic->Reset();
    chi2Cosmic->SetFunction(gr);
    std::cout << hough->GetRoh() << "\t" << hough->GetTheta() << std::endl;
    chi2Cosmic->SetRange(hough->GetRoh(),hough->GetTheta());
    chi2Cosmic->CalChisq();
    
    can->cd(3);
    image->Draw();
    
    hough->Reset();
    if( hough->CosmicJudgment(gr)){
      funcHough = hough->GetFunction();
      funcHough->GetXaxis()->SetRangeUser(-1000,1000);
      funcHough->GetYaxis()->SetRangeUser(-1000,1000);
      funcHough->Draw("same");
      line = hough->GetLine();
      line->SetLineColor(4);
      line->SetLineWidth(2);
      line->Draw("same");
    }
    can->cd(2);
    func = chi2Cosmic->GetFunction();
    func->Draw("colz");
    std::cout<< func->GetExpFormula() << std::endl;
    can->cd(1);
    if( hough->CosmicJudgment(gr)){
      image1->Draw();    
      lineChi = chi2Cosmic->GetLine();
      lineChi->SetLineColor(5);
      lineChi->SetLineWidth(2);
      lineChi->Draw("same");
    }
    can->cd(4);
    hisHough = hough->GetHisHough();
    hisHough->Draw("colz");
    can->Modified();
    can->Update();
    getchar();    
  }
  app->Run();
}
