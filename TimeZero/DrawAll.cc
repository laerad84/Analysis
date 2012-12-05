////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// DrawAll - no arguement //
// Input  : Pi0Out.root << output of TimeZeroTestWithPi0Run >> 
// Output : Pi0Peak.dat , Pi0Out.ps
// GetPeak From histogram and print out peak value and rms to text File. 
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPostScript.h"
#include "TCanvas.h"

//void DrawAll(){

int main( int argc, char** argv ){
  gStyle->SetOptFit(111111111);
  TFile* tf = new TFile("Pi0Out.root");
  TH1D* his[2716];
  for( int i = 0; i< 2716; i++){
    his[i] = (TH1D*)tf->Get(Form("hisTimeDelta%d",i));
  }

  TGraphErrors *grTimeDelta= new TGraphErrors();
  TH1D*        hisTimeDelta= new TH1D("hisTimeDeltaDistrib",
				      "hisTimeDeltaDistrib",
				      100, -37.5 , 12.5);
  TH1D*        hisTimeDeltaSmall = new TH1D("hisTimeDeltaSmall",
					    "hisTimeDeltaSmall",
					    100,-37.5,12.5 );

  Double_t TimeDelta[2716]={0};
  Double_t TimeDeltaSigma[2716]={-1};
  
  std::ofstream ofs("Pi0Peak.dat");
  TPostScript* ps = new TPostScript("Pi0Out.ps",111);
  TCanvas* can = new TCanvas("can","",764,1080);
  can->Divide(2,3);
  for( int i = 0; i< 2716; i++){
    can->cd( i%6+1 );
    if( his[i]->GetEntries() != 0 ){
      if( i == 2237 ){
	his[i]->Fit("gaus","Q","",his[i]->GetBinCenter(his[i]->GetMaximumBin())-5,
		    his[i]->GetBinCenter( his[i]->GetMaximumBin() ) +6 );
	
	his[i]->Draw();
	
      }else{
      his[i]->Fit("gaus","Q","",his[i]->GetBinCenter(his[i]->GetMaximumBin())-5,
		  his[i]->GetBinCenter( his[i]->GetMaximumBin() ) +6 );
      his[i]->Draw();
      }
      TF1* func = his[i]->GetFunction("gaus");
      TimeDelta[ i ]      = func->GetParameter( 1 );
      TimeDeltaSigma[ i ] = func->GetParameter( 2 );
	
      /*
      ofs << i << "\t" 
	  << func->GetParameter( 1 ) << "\t" 
	  << func->GetParameter( 2 ) << "\n";
      */

      grTimeDelta->SetPoint(grTimeDelta->GetN(), i, func->GetParameter( 1 ));
      grTimeDelta->SetPointError( grTimeDelta->GetN() -1 , 0, func->GetParError( 2 ));
      hisTimeDelta->Fill( func->GetParameter(1) );
      if( i < 2240 ){
	hisTimeDeltaSmall->Fill( func->GetParameter(1) );
      }
			    
    }else{
      TimeDelta[i]       = 0;
      TimeDeltaSigma[i]  = -1;
      /*
      ofs << i << "\t" 
	  << 0 << "\t"
	  << 0 << "\n";
      */
    }
    if( (i+1)%6 == 0 ){
      can->Update();
      can->Modified();
      ps->NewPage();
    }
  }

  for( int i = 0; i< 2716; i++){
    ofs << i                                           << "  " 
	<< TimeDelta[i] - hisTimeDeltaSmall->GetMean() << "  " 
	<< TimeDeltaSigma[i]                           << "\n";
  }
  can->Clear();
  can->Divide(1,3);
  can->cd(1);
  grTimeDelta->Draw("AP");
  can->cd(2);
  hisTimeDelta->Draw();
  can->cd(3);
  hisTimeDeltaSmall->Draw();
  can->Update();
  can->Modified();
  ofs.close();
  ps->Close();
}
  
  

  
    
  
			    
