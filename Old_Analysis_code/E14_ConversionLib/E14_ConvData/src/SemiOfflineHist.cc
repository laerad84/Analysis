#include "SemiOfflineHist.h"
#include <iostream>
#include <iomanip>

#include <cmath>

#if !defined(__CINT__)
ClassImp(SemiOfflineHist)
#endif

using namespace std;

SemiOfflineHist::SemiOfflineHist()
{
  nFADC = 16;
  nCH = 16;

  TStyle* gStyle = new TStyle();
  gStyle->SetPalette(1);

  cLaserRaw = new TCanvas("cLaserRaw","cLaserRaw",800,600);
  cPedestal = new TCanvas("cPedestal","cPedestal",800,600);

  NDet = 6;
  DetName[0] = "Csi";
  DetName[1] = "CC03";
  DetName[2] = "OEV";
  DetName[3] = "CV";
  DetName[4] = "Cosmic";
  DetName[5] = "Laser";

  char name[256];
  for(int idet=0;idet<NDet;idet++){
    map[idet] = new E14Mapper();
    sprintf(name,"/share/apps/production/2012VME/sumup/map/ch_map_%s.txt",DetName[idet].c_str());
    map[idet]->GetDataFromText( name );
  }

  box = new TBox();
  box -> SetFillColor(3);
  
}

SemiOfflineHist::~SemiOfflineHist()
{
  ;
}

void SemiOfflineHist::InitHist()
{

  hfile = new TFile( Form("/disk/kotodaq/2012Feb/run%d/RootFile/SemiOffline_%d_%d.root",RunNum,RunNum,CrateID),"recreate");

  for(int i=0;i<nFADC;i++){
    for(int j=0;j<nCH;j++){
      hLaserADC[i][j] = new TH1F( Form("hLaserADC_%d_%d",i,j),"",400,0,400000);
      hLaserTime[i][j] = new TH1F( Form("hLaserTime_%d_%d",i,j),"",48,0,48);
    }
  }

  for(int i=0;i<nFADC;i++){
    for(int j=0;j<nCH;j++){
      hPed[i][j] = new TH1F( Form("hPed_%d_%d",i,j),Form("Pedestal for Crate:%d FADC:%d for CH:%d Crystal ID:%d; ADC;",CrateID,i,j,map[0]->GetE14ID(CrateID,i,j)),600,200,800);
    }
  }

  for( int i=0; i<nFADC; i++ ){
    for( int j=0; j<2; j++ ){
      hRaw[i][j] = new TH2F( Form("hRaw_%d_%d",i,j), Form("RawData for Crate:%d FADC:%d CH:%d Crystal ID:%d; nsample; ADC",CrateID,i,j,map[0]->GetE14ID(CrateID,i,j)), 48, 0, 48, 390, 100, 4000);
      hRaw[i][j]->SetMarkerStyle(20);
      hRaw[i][j]->SetMarkerSize(0.5);
    }
  }

  hSpikeHist = new TH1F("hSpikeHist",Form("Spike count for Crate%d; FADCnum*16+ch",CrateID),nFADC*nCH,0,nFADC*nCH);

  sprintf( psFileName1, "/disk/kotodaq/2012Feb/run%d/psfile/LaserRawForCrate_%d_%d.ps", RunNum, RunNum, CrateID);
  cLaserRaw->Print( Form( "%s[", psFileName1 ) );
  cLaserRaw->SetRightMargin(0.2);
  sprintf( psFileName2, "/disk/kotodaq/2012Feb/run%d/psfile/PedestalForCrate_%d_%d.ps", RunNum, RunNum, CrateID);
  cPedestal->Print( Form( "%s[", psFileName2 ) );
  cPedestal->SetRightMargin(0.2);

}

void SemiOfflineHist::Fill()
{
  double SumIntegratedADC = 0;
  for(int i=0;i<nFADC;i++){
    for(int j=0;j<nCH;j++){
      SumIntegratedADC += IntegratedADC[i][j];
    }
  }
  
  for(int i=0;i<nFADC;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<48;k++){
	hRaw[i][j] -> Fill( k, Data[i][j][k] );
      }
    }
  }

  bool LASERTAG = false;
  if( CrateID == 10 ){
    //    if( 5000 < IntegratedADC[6][4] || IntegratedADC[6][4] < -1000 ){
    if( 5000 < IntegratedADC[6][4] ){
      LASERTAG = true;      
    }
  }else{
    //    if( 5000000 < SumIntegratedADC || SumIntegratedADC < -100000 ){
    if( 5000000 < SumIntegratedADC ){
      LASERTAG = true;
    }
  }

  if( LASERTAG ){ // Laser

    for(int i=0;i<nFADC;i++){
      for(int j=0;j<nCH;j++){
	hLaserADC[i][j] -> Fill( IntegratedADC[i][j] );
	//	if( (float)PeakHeight[i][j]-Pedestal[i][j] > 2){
	if( IntegratedADC[i][j] > 200 ){
	  hLaserTime[i][j] -> Fill( PeakTime[i][j] );
	}
      }
    }    

  }else{
    
    if( -10000 < SumIntegratedADC ){

      for(int i=0;i<nFADC;i++){
	for(int j=0;j<nCH;j++){
	  hPed[i][j] -> Fill( Pedestal[i][j] );
	}
      }
      
      SpikeCheck();
    }
  }

}

void SemiOfflineHist::WriteAll()
{
  double EventNum = 0;
  double x[16*16] = {0};
  double xE[16*16] = {0};
  double Mean[16*16] = {0};
  double MeanE[16*16] = {0};

  EventNum = hLaserADC[0][0] -> GetEntries();

  // For Laser ADC
  for(int i=0;i<nFADC;i++){
    for(int j=0;j<nCH;j++){
      x[i*16+j] = i*16+j+0.5;
      Mean[i*16+j]  = hLaserADC[i][j] -> GetMean();
      if( EventNum == 0 ) continue;
      MeanE[i*16+j] = hLaserADC[i][j] -> GetRMS() / sqrt(EventNum);
    }
  }
  gLaserADC = new TGraphErrors(nFADC*16,x,Mean,xE,MeanE);
  gLaserADC -> SetName("gLaserADC");
  gLaserADC -> SetTitle( Form("Integrated ADC value for Laser in run %d, crate %d ; FADCnum*16+ch ; Mean (Error = RMS/#sqrt{%d} )",RunNum,CrateID,(int)EventNum ));

  // For Laser Time
  for(int i=0;i<nFADC;i++){
    for(int j=0;j<nCH;j++){
      x[i*16+j] = i*16+j+0.5;
      Mean[i*16+j]  = hLaserTime[i][j] -> GetMean();
      if( EventNum == 0 ) continue;
      MeanE[i*16+j] = hLaserTime[i][j] -> GetRMS() / sqrt(EventNum);
    }
  }
  for(int i=0;i<nFADC;i++){
    gLaserTime[i] = new TGraphErrors(nCH,x,&Mean[i*16],xE,&MeanE[i*16]);
    gLaserTime[i] -> SetName( Form("gLaserTime%d",i) );
    gLaserTime[i] -> SetTitle( Form("Peak time value for Laser in run %d, crate:%d FADC:%d; ch ; Mean (Error = RMS/#sqrt{%d} )",RunNum,CrateID,i,(int)EventNum ));
    gLaserTime[i] -> SetMaximum(48);
    gLaserTime[i] -> SetMinimum(0);
    gLaserTime[i] -> SetMarkerStyle(20);

  }
    // Writing for psfile
    /*
      cLaserRaw -> cd();
      gLaserADC -> Draw("AP");
      cLaserRaw -> Print( psFileName );
      cLaserRaw -> Clear();
    */
   
  for(int j=0;j<4;j++){
    cLaserRaw -> cd();
    cLaserRaw -> Divide(2,2);
    for(int i=0;i<4;i++){
      //      cLaserRaw -> GetPad(i+1) -> SetLogy();
      //      cLaserRaw -> GetPad(i+1) -> SetLogz();
      cLaserRaw -> cd(i+1);
      gLaserTime[j*4+i] -> Draw("AP");

      for(int k=0;k<16;k++){
	int e14id = 9999;
	for(int l=0;l<3;l++){
	  if( e14id == 9999 ){
	    e14id = map[l]->GetE14ID(CrateID,j*4+i,k);
	  }
	}
	if( e14id == 9999 ){
	  box -> DrawBox(k,0,k+1,48);
	}
      }
      
    }
    cLaserRaw -> Print( psFileName1 );
    cLaserRaw -> Clear();
  }

  for(int j=0;j<2;j++){
    cLaserRaw -> cd();
    cLaserRaw -> Divide(4,4);
    for(int i=0;i<16;i++){
      cLaserRaw -> GetPad(i+1) -> SetLogy();
      cLaserRaw -> GetPad(i+1) -> SetLogz();
      cLaserRaw -> cd(i+1);
      hRaw[i][j] -> Draw("COLTZ");
    }
    cLaserRaw -> Print( psFileName1 );
    cLaserRaw -> Clear();
  }

  cLaserRaw -> cd();
  cLaserRaw -> SetLogy();
  hSpikeHist -> Draw();
  cLaserRaw -> Print( psFileName1 );
  cLaserRaw -> Clear();  

  for(int i=0;i<nFADC;i++){
    cPedestal -> cd();
    cPedestal -> Divide(4,4);
    for(int j=0;j<nCH;j++){
      cPedestal -> GetPad(j+1) -> SetLogy();
      cPedestal -> cd(j+1);
      hPed[i][j] -> Draw();
    }
    cPedestal -> Print( psFileName2 );
    cPedestal -> Clear();
  }

  cLaserRaw->Print( Form( "%s)", psFileName1 ) );
  cPedestal->Print( Form( "%s)", psFileName2 ) );

  // Writing for rootfile
  hfile -> cd();
  gLaserADC -> Write();
  for(int i=0;i<nFADC;i++){
    gLaserTime[i] -> Write();
    for(int j=0;j<nCH;j++){
      hLaserADC[i][j] -> Write();
      hLaserTime[i][j] -> Write();
      hPed[i][j] -> Write();
    }
  }
  hfile -> Close();
  
  hfile = new TFile( Form("/disk/kotodaq/2012Feb/run%d/RootFile/FinishTag_%d_%d.root",RunNum,RunNum,CrateID),"recreate");
  hfile -> Close();
  
  }

void SemiOfflineHist::SpikeCheck()
{
  int old=0, present=0, next=0, average=0;
  bool spikeFlag = false;
  bool spikeEventFlag = false;

  for( int i=0; i<nFADC; i++ ){
    for( int j=0; j<nCH; j++ ){
      for( int isample=1; isample<47; isample++ ){
	old     = Data[i][j][isample-1];
	present = Data[i][j][isample];
	next    = Data[i][j][isample+1];
	average = (old+next)/2;
	if( average!=0 && abs((double)(present-average))/average>0.2 ){
	  spikeFlag = true;
	  spikeEventFlag = true;
	  break;
	}
      }
      if( spikeFlag ){
	if( !SPIKEID[i][j] ){
	  printf("Found spike!!  crate ID=%d : FADC ID=%d : CH ID=%d\n",
		 CrateID,i,j); 
	  SPIKEID[i][j] = true;
	}

	hSpikeHist -> Fill(i*16+j);
	spikeFlag = false;
      }
    } 
    
    if( spikeEventFlag ) SpikeEvent++;
    
  }
}
