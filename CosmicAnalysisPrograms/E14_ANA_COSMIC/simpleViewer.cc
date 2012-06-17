
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "TH1.h"
#include "E14ReadSumFile.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TAxis.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "TStyle.h"
#include "TCanvas.h"
int main(int argc, char** argv){

  gStyle->SetPalette(1);
  
  std::cout<< __LINE__ << std::endl;
  
  TApplication* app = new TApplication("app",&argc, argv);
  IDHandler* handler = new IDHandler("Data/crystal.txt");
  
  const int nCH = 2716;
  const int nPoint  = 4; 
  E14ReadSumFile* Reader[nPoint];
  CsIImage* image[nPoint];
  
  int RunNumber[nPoint] = {4796,4797,4798,4799};
  std::string filename = "/Volume0/ExpData/2012_Feb_Beam/Sum%d.root";
  std::cout<< __LINE__ << std::endl;
  
  for(int ipoint  =0 ;ipoint < nPoint;++ipoint){
    Reader[ipoint] = new E14ReadSumFile();
    Reader[ipoint]->Add(Form(filename.c_str(),RunNumber[ipoint]));     
    image[ipoint]  = new CsIImage(handler);
  }
  std::cout << "ANALYSIS" << std::endl;

  /*
  E14ReadSumFile* Reader2  = new E14ReadSumFile();
  Reader->Add("/Volume0/ExpData/2012_Feb_Beam/Sum4798.root");
  */

  TH1D* hist[nPoint][nCH];
  for( int ipoint  = 0; ipoint < nPoint ; ++ipoint){
    for( int ich = 0; ich < nCH; ++ich){
      hist[ipoint][ich] = new TH1D(Form("hist_%d_%d",ich,ipoint),
				   Form("hist_%d_%d",ich,ipoint),
				   1600,
				   0,
				   160000);
    }
  }
  
  for( int ipoint  = 0; ipoint < nPoint; ++ipoint){
    for( int ievent = 0; ievent < Reader[ipoint]->GetEntries(); ++ievent){
      Reader[ipoint]->GetEntry(ievent);
      for(int idigi = 0; idigi < Reader[ipoint]->CsiNumber; idigi++){
	  if(Reader[ipoint]->CsiEne[idigi] > 10){
	    hist[ipoint][ Reader[ipoint]->CsiModID[idigi] ]->Fill(Reader[ipoint]->CsiEne[idigi]);
	  }
	}
    }
  }

  TGraph* gr[nPoint];
  for( int i = 0; i< nPoint; ++i){
    gr[i]= new TGraph();
  }
  
  for( int ich = 0; ich < nCH; ich++){
    if( hist[0][ich]->GetMean() > 10){
      image[0]->Fill(ich,hist[0][ich]->GetMean());
    }
  }
 
  for( int ipoint =1; ipoint < nPoint; ipoint++){
    for( int ich = 0; ich < nCH; ich++){
      if(hist[ipoint][ich]->GetMean()> 10 &&
	 hist[0][ich]->GetMean()     > 10 ){
	gr[ipoint-1]->SetPoint(gr[ipoint-1]->GetN(),ich,hist[ipoint][ich]->GetMean()/hist[0][ich]->GetMean());
	image[ipoint]->Fill(ich,hist[ipoint][ich]->GetMean()/hist[0][ich]->GetMean());
      }
    }
  }
  
  for( int ipoint = 0; ipoint < nPoint; ipoint++){
    gr[ipoint]->SetMarkerStyle(5);
    gr[ipoint]->SetMarkerColor(ipoint+1);
  }
  
  TCanvas* can  = new TCanvas("can","",1200,800);
  can->Divide(3,2);
  for( int i = 0; i< nPoint; i++){
    can->cd(i+1);
    if( i==0){
      image[i]->DrawWithRange("colz",0,3000);
    }else{
      image[i]->DrawWithRange("colz",0,3);
    }
  }
  can->cd(6);
  gr[0]->GetYaxis()->SetRangeUser(0,4);
  gr[0]->Draw("AP");
  for( int i = 0; i< nPoint-1; i++){
    gr[i]->Draw("p");
  }  
  app->Run();
}

	  

