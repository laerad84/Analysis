#include <iostream>
#include <fstream>

#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TApplication.h"

//int
//main ( int argc, char** argv){
void confirmCosmicData(){

  //TFile* tf = new TFile("cosmicPeak3599.root");
  TFile* tf = new TFile("cosmic_3599_3600.root");
  TGraphErrors* gain   = (TGraphErrors*)tf->Get("gain");
  TGraphErrors* cosmic = (TGraphErrors*)tf->Get("cosmic");
  
  TH1D* his[2716];
  for( int i = 0; i< 2716; i++){
    his[i] = (TH1D*)tf->Get(Form("his_CH%04d",i));
  }
  TCanvas* can = new TCanvas("can","",800,800);
  can->Draw();
  for( int i = 0; i< 2716; i++){
    if( gain->GetY()[i] < 700 ){
      if(his[i]->Integral()>0){
	his[i]->Draw();
	
	TF1* func  = his[i]->GetFunction("landau");
	std::cout<< i << "\t" <<  func->GetParameter(1) << " : " 
		 << func->GetParameter(2)<< std::endl;
	can->Modified();
	can->Update();
	getchar();
      }
    }
  }
}
      
  
  
