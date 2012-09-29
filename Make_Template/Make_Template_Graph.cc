#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TF1.h"
#include "TMath.h"
#include "TSpline.h"
#include "TApplication.h"
int 
main( int argc , char** argv ){
  TApplication* app = new TApplication( "app", &argc, argv );
  TCanvas* canvas  = new TCanvas("canvas","",800,800);
  TFile* tf1 = new TFile("TEMPLETE_OUT_HEIGHT_3pi0RunList_200_400.root");
  TFile* tfOut = new TFile("TEMPLATE_COMPLETE_GRAPH_3pi0RunList_200_400.root",
			   "RECREATE");

  //TFile* tf2 = new TFile("TEMPLETE_OUT_HEIGHT_3pi0RunList_400_800.root");
  TGraph* grWaveform[2716];
  TGraph* grTemplateGraph[2716];
  TGraph* grTemplateFinal[2716];
  
  TSpline3* splTemplateGraph[2716];

  for( int ich = 0; ich< 2716; ich++){
    grWaveform[ich]      = (TGraph*)tf1->Get(Form("Waveform_Height_200_400_%d",ich));
    //splTemplateGraph[ich]= new TSpline3(Form("Spl%d",ich),grWaveform[ich]);
    grTemplateGraph[ich] = new TGraph();
    grTemplateGraph[ich]->SetName(Form("Template_Graph_%d",ich));
    grTemplateGraph[ich]->SetTitle(Form("Template_Graph_%d",ich));
    grTemplateFinal[ich] = new TGraph();
    grTemplateFinal[ich]->SetName(Form("Template_%d",ich));
    grTemplateFinal[ich]->SetTitle(Form("Template_%d",ich));
  }

  for( int i = 0; i< 2716; i++){
    if( grWaveform[i]->GetN() > 300 && i != 1862){
      Double_t Pedestal= 0.;
      Int_t    nPoint  = 0; 
      Double_t FormerYValue=0;
      Double_t FormerXValue=-150.5;      
      for( int ipoint = 0; ipoint < grWaveform[i]->GetN(); ipoint++){
	if( grWaveform[i]->GetX()[ipoint] > -100 ){ break; }
	Pedestal += grWaveform[i]->GetY()[ipoint];
	nPoint++;
      }
      Pedestal = Pedestal / nPoint; 

      if( grWaveform[i]->GetX()[0] > -149.5){
	grTemplateGraph[i]->SetPoint(0,-149.5,0);
      }      
      for( int ipoint = 0; ipoint < grWaveform[i]->GetN(); ipoint++){
	if( i != 2395 ){
	  if( grWaveform[i]->GetX()[ipoint] < 100){
	    if( TMath::Abs(grWaveform[i]->GetY()[ipoint] - FormerYValue) >
		TMath::Abs(grWaveform[i]->GetX()[ipoint] - FormerXValue)*0.03 ){ continue; }
	  }else{
	    if( TMath::Abs(grWaveform[i]->GetY()[ipoint] - FormerYValue) >
		TMath::Abs(grWaveform[i]->GetX()[ipoint] - FormerXValue)*0.01 ){ continue; }
	  }
	}else{
	  if( grWaveform[i]->GetX()[ipoint] < 100){
	    if( TMath::Abs(grWaveform[i]->GetY()[ipoint] - FormerYValue) >
		TMath::Abs(grWaveform[i]->GetX()[ipoint] - FormerXValue)*0.03 ){ continue; }
	  }else{
	    if( TMath::Abs(grWaveform[i]->GetY()[ipoint] - FormerYValue) >
		TMath::Abs(grWaveform[i]->GetX()[ipoint] - FormerXValue)*0.01 ){ continue; }
	  }

	}
	if( grWaveform[i]->GetX()[ipoint] < -100 ){
	  grTemplateGraph[i]->SetPoint(grTemplateGraph[i]->GetN(), 
				       grWaveform[i]->GetX()[ipoint], 
				       0);
				       //grWaveform[i]->GetY()[ipoint]);
	  FormerXValue = grWaveform[i]->GetX()[ipoint];
	  FormerYValue = grWaveform[i]->GetY()[ipoint];
	  
	}else if( grWaveform[i]->GetX()[ipoint]< 200 ){
	  grTemplateGraph[i]->SetPoint(grTemplateGraph[i]->GetN(), 
				       grWaveform[i]->GetX()[ipoint], 
				       grWaveform[i]->GetY()[ipoint]-Pedestal);
	  FormerXValue = grWaveform[i]->GetX()[ipoint];
	  FormerYValue = grWaveform[i]->GetY()[ipoint];
	}else{
	  break;
	}
      }
            
      grTemplateGraph[i]->Fit("pol1", "Q","",175,200);
      TF1* func = grTemplateGraph[i]->GetFunction("pol1");
      grTemplateGraph[i]->Draw("AP");

      /*
      canvas->Update();
      canvas->Modified();
      getchar();
      */

      if( func->GetParameter(1) > 0 ){
	std::cerr << "Tail is Up !!!! : " << i << std::endl;
      }
      for( int ipoint = 0; ipoint < 100; ipoint++){
	Double_t y = func->Eval( 200.5+ipoint );
	if( y > 0){
	  grTemplateGraph[i]->SetPoint( grTemplateGraph[i]->GetN(),
					200.5+ipoint,y);
	}else{
	  grTemplateGraph[i]->SetPoint( grTemplateGraph[i]->GetN(),
					200.5+ipoint,0);
	}
      }

      splTemplateGraph[i] = new TSpline3(Form("Spl%d",i), grTemplateGraph[i]);
					      
      for( int ipoint =0; ipoint < grTemplateGraph[i]->GetN(); ipoint++){
	grTemplateFinal[i]->SetPoint( ipoint,
				      grTemplateGraph[i]->GetX()[ipoint],
				      (grTemplateGraph[i]->GetY()[ipoint]-Pedestal)/splTemplateGraph[i]->Eval(0));
      }
    }else{
      for( int ipoint = 0; ipoint < 450;  ipoint++){
	grTemplateGraph[i]->SetPoint( ipoint, -149.5 + ipoint, 0);
	grTemplateFinal[i]->SetPoint( ipoint, -149.5 + ipoint , 0);
      }
    }
  }
    

  TPostScript* ps = new TPostScript("TemplateGraph.ps",111);
  ps->NewPage();
  TCanvas* can = new TCanvas("can","",720,720*TMath::Sqrt(2)*0.9);
  Int_t xdiv =2;
  Int_t ydiv =3; 
  Int_t nDiv = xdiv*ydiv;
  can->Divide( xdiv,ydiv);

  for( int i = 0; i< xdiv*ydiv; i++){
    can->cd(i+1);
    gPad->SetGridx();
    gPad->SetGridy();
  }

  for( int ich  = 0; ich < 2716; ich++){
    can->cd( (ich % nDiv) +1 );
    grTemplateGraph[ich]->SetMarkerStyle(6);
    grTemplateGraph[ich]->SetMarkerColor(4);
    grTemplateGraph[ich]->Draw("AP");
    if((ich % nDiv) == nDiv-1){
      can->Update();
      can->Modified();
      ps->NewPage();
    }
  }
  ps->Close();
  for( int ich = 0; ich < 2716; ich++){
    if( grTemplateGraph[ich]->GetListOfFunctions()->GetEntries() != 0 ){
      grTemplateGraph[ich]->GetListOfFunctions()->Delete();
    }
    grTemplateGraph[ich]->Write();
    grTemplateFinal[ich]->Write();
  }

  tfOut->Close();
  app->Run();
}
