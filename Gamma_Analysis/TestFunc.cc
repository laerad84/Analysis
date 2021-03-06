#include "PulseGenerator.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TApplication.h"
#include "TCanvas.h"

#include <vector>

int main( int argc, char** argv ){
  
  TApplication* app = new TApplication("app",&argc,argv);
  PulseGenerator* gen = new PulseGenerator();
  PulseGenerator* gen1 = new PulseGenerator(); 
  TCanvas* can = new TCanvas("can","",900,900);
  
  Double_t EnergyArr[8] ={10,20,40,80,160,320,640,1280};
  TH1D*    hisEnergy[8];
  TH1D*    hisTime[8];
  for( int i = 0; i< 8; i++){
    hisEnergy[i] = new TH1D(Form("his_%d",i),
			    Form("his_%d",i),
			    100,0,2);
    hisTime[i] = new TH1D(Form("hisT_%d",i),Form("hisT_%d",i),100,150,200);
  }
  /*
  gen->pfPDF->Draw();
  can->Update();
  can->Modified();
  getchar();
  gen->pfPDF->SetParameter(0,1);
  can->Update();
  can->Modified();
  getchar();
  */
  for( int iloop = 0; iloop < 8; iloop++){
    std::vector<double> EVector;
    std::vector<double> TVector;
    Double_t meanTime = 0; 
    Double_t Energy   = 0; 
    Double_t EnergySum= 0;
    Energy = EnergyArr[iloop];
    /*
    for( int i = 0; i< 10; i++){
    double tempEnergy = gRandom->Gaus(EnergyArr[iloop]*0.01, 5);
    if( tempEnergy <= 0 ) continue; 
    EVector.push_back( tempEnergy );
    TVector.push_back( i*1);
    meanTime += tempEnergy*i*1;
    Energy   += tempEnergy;
    EnergySum+= tempEnergy;
    }
    meanTime = meanTime/Energy;
    */

    for( int i = 0; i< 250 ; i++){
      //std::cout<< "Generator" << std::endl;
      //gen->Reset();
      //TF1* func = gen->GetWaveform( EVector, TVector);
      double HpE = atof(argv[1]);
      double gnd = 100;
      double nppE=12.7;
      TF1* func = gen->GetWaveform( Energy, meanTime,HpE,nppE);
      //std::cout << func->Eval(func->GetMaximumX(0,500))/(HpE*Energy) << std::endl;
      //std::cout<< gen->pfPDF->GetParameter(0) << std::endl;
      TGraph* grWave = gen->GenSignal(Energy,meanTime,HpE,nppE,gnd);
      grWave->SetMarkerStyle(21);
      grWave->Draw("AP");
      can->Update();
      can->Modified();
      getchar();

      //gen->Reset();
      TF1* func1 = gen->GetWaveform( Energy,meanTime);

      //gen->Reset();
      TF1* func2 = gen->GetWaveform( Energy, 0);
      
      //std::cout<< func->GetMaximumX(0, 500 )   << "  " 
      //<< func1->GetMaximumX( 0, 500 ) << "  "
      //<< func2->GetMaximumX( 0, 500 ) << std::endl;
      
      const Int_t nPoint= 48;
      TGraph* gr = new TGraph(nPoint);
      for( int ipoint = 0 ; ipoint < nPoint; ipoint++){
	double offset = 8*(gRandom->Rndm()*-0.5);
	double tsum   = 0.;
	double d_tmp  = func1->Eval(ipoint*8.-offset)+gRandom->Gaus(0.,2.4);
	gr->SetPoint(ipoint,ipoint*8,d_tmp);
      }      
      /*
      std::cout<< iloop << "\t" 
	       << func->Eval(func->GetMaximumX(0,500))/(150./13.7*1.63968*0.95)/EnergySum 
	       << std::endl;
      */
      /*
      std::cout<< func->Eval(func->GetMaximumX(0,500)) << std::endl;
      getchar();
      */
      //hisEnergy[iloop]->Fill( func->Eval(func->GetMaximumX(0,500))/(150./13.7*1.63968*0.952)/Energy);
      hisEnergy[iloop]->Fill( func->Eval(func->GetMaximumX(0,500))/(HpE*Energy));
      hisTime[iloop]->Fill(func->GetMaximumX(0,500));
      delete func;
      delete func1;
      delete func2;
      delete gr;
    }
  }
  can->Divide(4,4);
  TGraph* grResolution = new TGraph();
  TGraph* grTimeRes  = new TGraph();
  for( int i = 0; i< 8; i++){
    can->cd( i+1 );
    hisEnergy[i]->Fit("gaus");
    TF1* func = (TF1*)hisEnergy[i]->GetFunction("gaus");
    grResolution->SetPoint(grResolution->GetN(),EnergyArr[i],func->GetParameter(2));
    hisEnergy[i]->Draw();

    std::cout << hisEnergy[i]->GetMean() << std::endl;
    /*
    hisEnergy[i]->Fit("gaus","","",
		      hisEnergy[i]->GetMean()-3*hisEnergy[i]->GetRMS(),
		      hisEnergy[i]->GetMean()+3*hisEnergy[i]->GetRMS());
    */
    can->cd(i+9);

    Int_t FitStatus = hisTime[i]->Fit("gaus");
    TF1* tFunc = hisTime[i]->GetFunction("gaus");
    hisTime[i]->Draw();
    std::cout<< FitStatus << std::endl;
    if( FitStatus >= 0 ){ 
      grTimeRes->SetPoint(grTimeRes->GetN(), EnergyArr[i],tFunc->GetParameter(2));    
    }
  }

  can->cd(9);
  can->Update();
  can->Modified();
  TCanvas* canvas  =new TCanvas("canvas","",800,400);
  canvas->Divide(2,1);
  canvas->cd(1);
  grResolution->SetMarkerStyle(21);
  grResolution->Draw("AP");
  canvas->cd(2);
  grTimeRes->SetMarkerStyle(21);
  grTimeRes->Draw("AP");

  app->Run();
}
