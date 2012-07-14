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

#include "E14ReadSumFile.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "HoughCsI.h"

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

 
  E14ReadSumFile* Reader = new E14ReadSumFile(0);
  Reader->Add(InputFile.c_str());

  IDHandler* handler = new IDHandler();
  CsIImage*  image   = new CsIImage(handler);
  CsIImage*  image1   = new CsIImage(handler);
  TCanvas*   can     = new TCanvas("can","",800,0,800,800);
  can->Divide(2,2);


  long nEntries = Reader->GetEntries();
  std::cout << "TotalEntries::" << nEntries  << std::endl;
 
  TH2D*   hisHough;
  TF1*    funcHough;
  TLine*  line;

  Double_t roh;
  Double_t theta;
  Int_t    HitFlag[N_TOTAL_CSI]={};
  HoughCsI* hough = new HoughCsI();
  for( int iEntry = 0; iEntry < nEntries; iEntry++){
    Reader->GetEntry(iEntry);
    
    ///AnalysisCode    
    TGraph* gr = new TGraph();
    hough->Reset();
    image->Reset();
    //image1->Reset();
    std::cout <<"Run Analysis "<< std::endl;
    for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
      if( Reader->CsiEne[idigi] > 300 ){
	double x,y;
	//std::cout << "ID::" << Reader->CsiModID[idigi] << " " << Reader->CsiEne[idigi] << std::endl;
	handler->GetMetricPosition(Reader->CsiModID[idigi], x, y);
	//image->Fill(Reader->CsiModID[idigi],Reader->CsiEne[idigi]);
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

    if(hough->CosmicJudgment(gr)){
      roh = hough->GetRoh();
      theta = hough->GetTheta();
      for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
	if( Reader->CsiEne[idigi] > 300 ){
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

    can->cd(1);
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
    gr->GetYaxis()->SetRangeUser(-1000,1000);
    gr->GetXaxis()->SetRangeUser(-1000,1000);
    
    gr->Draw("AP");
    can->cd(3);
    image1->Draw();
    can->cd(4);
    hisHough = hough->GetHisHough();
    hisHough->Draw("colz");
    can->Modified();
    can->Update();
    getchar();
    
  }
  app->Run();
}
