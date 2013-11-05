#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TF1.h"

int main( int argc , char** argv ){

  TF1* func = new TF1("func","gaus",-5,5 );
  TFile* tf = new TFile("GammaTimeShape.root");
  const int nAngle= 9;
  const int nD    = 5;
  TH2D* hisRDT[nAngle][nD];
  TH1D* hisRDTMean[nAngle][nD][3];
  for( int i = 0; i< nAngle; i++){
    for( int j = 0; j< nD; j++){
      hisRDT[i][j] = (TH2D*)tf->Get(Form("hisRDT_%d_%d",i,j));
      hisRDT[i][j]->FitSlicesY(func);
      hisRDTMean[i][j][0] = (TH1D*)gDirectory->FindObject(Form("%s_0",hisRDT[i][j]->GetName()));
      hisRDTMean[i][j][1] = (TH1D*)gDirectory->FindObject(Form("%s_1",hisRDT[i][j]->GetName()));
      hisRDTMean[i][j][2] = (TH1D*)gDirectory->FindObject(Form("%s_2",hisRDT[i][j]->GetName()));
    }
  }


  TFile* tfOut = new TFile("GammaTimeShapeRD.root","recreate");  
  TH2D* hisGammaTimeShape[nAngle];
  for( int i = 0; i< nAngle; i++){
    hisGammaTimeShape[i] = new TH2D(Form("hisGammaTimeShape_%d",i),
				    Form("hisGammaTimeShape_%d;R[mm];D[mm]",i),
				    9,-112.5,112.5,
				    9,-112.5,112.5);
  }

  for( int i = 0; i< nAngle; i++){
    for( int j  =0; j< nD; j++){
      for( int ibin = 1; ibin <= hisRDT[i][j]->GetNbinsX(); ibin++){
	if( hisRDTMean[i][j][0]->GetBinContent(ibin) < 5){ continue; }
	Double_t D = (j)*25+6;
	Double_t R = (ibin-5)*25;
	Double_t TimeValue = hisRDTMean[i][j][1]->GetBinContent(ibin);
	Double_t Error     = hisRDTMean[i][j][2]->GetBinContent(ibin);
	int nBin = hisGammaTimeShape[i]->Fill(R,D,TimeValue);
	hisGammaTimeShape[i]->SetBinError( nBin, Error);
	if( TMath::Abs( D ) > 12.5 ){
	  nBin = hisGammaTimeShape[i]->Fill(R,-D,TimeValue);
	  hisGammaTimeShape[i]->SetBinError( nBin, Error);
	}
      }
    }
  }

  for( int i = 0; i< nAngle; i++){
    hisGammaTimeShape[i]->Write();    
  }
  tfOut->Close();
  return 0;
}
