#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TMath.h"

TSpline3*  spl;
TSpline3*  splLaser;
TSpline3*  splCsi[2716];
TSpline3*  splLaserTemp;
TSpline3*  splCsiTemp[2716];

double FuncsplCsi( double *x, double *par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double value = spl->Eval(x0-p0)*p1+par[2];
  return value;
}
double FuncsplLaser( double *x, double *par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double value = splLaser->Eval(x0-p0)*p1+p2;
  return value;
}

double FuncCorr( double *x , double *par){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double value = p0+p1*exp(p2*x0);
  return value;
}

double AdjFunc( double* x, double* par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double value = p0 + p1*exp(p2*x0/10.07) + p3*exp(p4*x0/10.07);
  return value;
}

//void DrawLaserWaveform_Fit(){
int main(int argc, char** argv){
  int DivNumber = atoi( argv[1] );
  if( DivNumber >= 8 || DivNumber < 0 ){ return -1; }
  gStyle->SetOptFit(111111111);
  TF1* TimeAdjFunc = new TF1("TimeAdjFunc",AdjFunc, 0, 16000,5);
  Double_t Par[5] = {-0.0905327+0.0112319,1.54915,-0.114423,0.0758477,0.00487457};
  Double_t ParErrors[5] = {0.00245834,0.0263153,0.00188467,0.007594,4.81501e-05};
  TimeAdjFunc->SetParameters(Par);
  TimeAdjFunc->SetParErrors(ParErrors);

  TH2D*     tempLaser;
  TProfile* profLaser;
  TGraph*   grLaser;
  TGraph*   grLaserTemp;

  TH2D*     tempCsi[2716]={NULL};
  TProfile* profCsi[2716]={NULL};
  TGraph*   grCsi[2716]={NULL};
  TGraph*   grCsiTemp[2716]={NULL};
  
  // Set STemplate to hist // 
  TFile* tfLaserTemplate = new TFile("LaserWaveTemplate.root");
  tempLaser = (TH2D*)tfLaserTemplate->Get("hisLaserTemplate");
  profLaser = tempLaser->ProfileX();
  profLaser->SetLineColor(1);
  grLaser     = new TGraph();
  grLaserTemp = new TGraph();
  for( int iBin = 0; iBin < profLaser->GetNbinsX();iBin++){
    grLaserTemp->SetPoint( grLaserTemp->GetN(),profLaser->GetBinCenter(iBin+1),profLaser->GetBinContent(iBin+1));
  }
  splLaserTemp = new TSpline3("splLaserTemp",grLaserTemp);
  for( int iBin = 0; iBin < profLaser->GetNbinsX();iBin++){
    grLaser->SetPoint( grLaser->GetN(),profLaser->GetBinCenter(iBin+1),profLaser->GetBinContent(iBin+1)/splLaserTemp->Eval(0));
  }
  splLaser = new TSpline3("splLaser",grLaser);
  tfLaserTemplate->Close();

  TFile* tfTemplate = new TFile("laserTemplate.root");
  for( int i = 0; i< 2716;i++){
    tempCsi[i]=(TH2D*)tfTemplate->Get(Form("ChannelTemplate_%d",i));
    if( tempCsi[i]->GetEntries() < 4800 ){ continue; }
    profCsi[i]=tempCsi[i]->ProfileX();
    grCsi[i] = new TGraph();
    grCsiTemp[i] = new TGraph();
    for( int iBin = 0; iBin < profCsi[i]->GetNbinsX();iBin++){
      grCsiTemp[i]->SetPoint( grCsiTemp[i]->GetN(),profCsi[i]->GetBinCenter(iBin+1),profCsi[i]->GetBinContent(iBin+1));
    }
    splCsiTemp[i] = new TSpline3(Form("splCsiTemp_%d",i),grCsiTemp[i]);
    for( int iBin = 0; iBin < profCsi[i]->GetNbinsX();iBin++){
      grCsi[i]->SetPoint( grCsi[i]->GetN(),profCsi[i]->GetBinCenter(iBin+1),profCsi[i]->GetBinContent(iBin+1)/splCsiTemp[i]->Eval(0));
    }
    splCsi[i] = new TSpline3(Form("splCsi_%d",i),grCsi[i]);
  }
  tfTemplate->Close();
  
  TFile* tf = new TFile("/Volume0/ExpData/2012_Feb_Beam/RootFile_wav/WaveformExtract_4747.root");
  TTree *tr = (TTree*)tf->Get("WFTree");
  Double_t LaserTrigWF[48];
  Double_t ChannelOutWF[2716][48];
  Double_t LaserSignal[5];
  Double_t LaserTime[5];
  Double_t LaserPedestal[5];

  Int_t CsiNumber;
  Double_t CsiPedestal[2716];
  Double_t CsiTime[2716];
  Double_t CsiSignal[2716];
  Short_t CsiID[2716];
  Int_t EventNo;
  tr->SetBranchAddress("EventNo"      ,&EventNo);
  tr->SetBranchAddress("LaserTrigWF"  ,LaserTrigWF);
  tr->SetBranchAddress("ChannelOutWF" ,ChannelOutWF);
  tr->SetBranchAddress("LaserSignal"  ,LaserSignal);
  tr->SetBranchAddress("LaserTime"    ,LaserTime);
  tr->SetBranchAddress("LaserPedestal",LaserPedestal);
  tr->SetBranchAddress("CsiNumber"    ,&CsiNumber);
  tr->SetBranchAddress("CsiID"        ,CsiID);//CsiNumber
  tr->SetBranchAddress("CsiSignal"    ,CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiTime"      ,CsiTime);//CsiNumber
  tr->SetBranchAddress("CsiPedestal"  ,CsiPedestal);//CsiNumber


  //Int_t CsiChannelID = 9;
  //TCanvas* can = new TCanvas("can","",800,800);
  //can->Divide(2,2);
  
  Int_t nEntries = tr->GetEntries();
  TGraph* grCsiEvent[2716];
  for( int i = 0; i< 2716; i++){
    grCsiEvent[i] = new TGraph();
  }
  Int_t DivStart = nEntries/8*DivNumber;
  Int_t DivEnd   = nEntries/8*(DivNumber+1);
  TGraph* grLaserEvent = new TGraph();
  TF1* funcLaser = new TF1("funcLaser",FuncsplLaser,-100,200,3);
  
  TFile*   tfOut  = new TFile(Form("LaserTime_m100p50_%d.root",DivNumber),"recreate");
  Int_t    EventNumber;
  Double_t LaserT;
  Double_t LaserHeight;
  Int_t    nDigi;
  Int_t    CsiModID[2716];
  Double_t CsiT[2716];
  Double_t CsiHeight[2716];
  TTree*   trout  =new TTree("FitTimeHeight","");
  trout->Branch("EventNumber",&EventNumber,"EventNumber/I");
  trout->Branch("LaserT"     ,&LaserT     ,"LaserT/D");
  trout->Branch("LaserHeight",&LaserHeight,"LaserHeight/D");
  trout->Branch("nDigi"      ,&nDigi      ,"nDigi/I");
  trout->Branch("CsiModID"   ,CsiModID    ,"CsiModID[nDigi]/I");//nDigi
  trout->Branch("CsiT"       ,CsiT        ,"CsiT[nDigi]/D");//nDigi
  trout->Branch("CsiHeight"  ,CsiHeight   ,"CsiHeight[nDigi]/D");//nDigi
  /*
  TH1D* hisDelta[2716];
  TH2D* hisDelta2D[2716];
  TH1D* hisDeltaStd[2716];
  TH2D* hisDeltaMin2D[2716];
  for( int i = 0; i< 2716; i++){
    hisDelta[i]     = new TH1D(Form("hisDelta_%d",i),Form("hisDelta_%d",i),100,0,50);
    hisDelta2D[i]   = new TH2D(Form("hisDelta2D_%d",i),Form("hisDelta2D_%d",i),160,0,16000,500,-10,40);
    hisDeltaStd[i]  = new TH1D(Form("hisDeltaStd_%d",i),Form("hisDeltaStd_%d",i),300,10,40);
    hisDeltaMin2D[i]= new TH2D(Form("hisDeltaMin2D_%d",i),Form("hisDeltaMin2D_%d",i),160,0,1600,500,-10,40);
  }
  */
  //  for( int ievent  = 0; ievent < nEntries; ievent++){
  std::cout<< "Loop" << std::endl;

  for( int ievent  = DivStart; ievent < DivEnd; ievent++){
    if( (ievent%100) == 0 && ievent ){ std::cout<< ievent << "/" << nEntries << std::endl; }

    for( int i = 0; i< 2716; i++){
      CsiModID[i] = -1;
      CsiT[i] = 0; 
      CsiHeight[i] = 0;
    }
    for( int i = 0; i< 2716; i++){
      grCsiEvent[i]->Set(0);
    }
    grLaserEvent->Set(0);
    
    //for( int ievent  = 0; ievent < 20; ievent++){
    tr->GetEntry(ievent);
    EventNumber = EventNo;
    Double_t Minimum = 0xFFFF;
    Double_t Maximum = 0;
    Int_t    MaximumIndex = 0;
    for( int i = 0; i< 48; i++){
      if(LaserTrigWF[i] > Maximum){
	Maximum = LaserTrigWF[i];
	MaximumIndex = i;
      }
    }
    for( int i = 0; i< 48; i++){
      if( i > MaximumIndex ){ break; }
      if( LaserTrigWF[i] < Minimum ){
	Minimum = LaserTrigWF[i];
      }
    }
    if( Maximum - Minimum < 400){ continue; }
    
    Int_t LaserMaxPoint;
    Double_t LaserMax=0;
    
    Int_t    CsiMaxPoint[2716];
    Double_t CsiMax[2716]={0};
    // SetLaser WF & Fit
    for( int iPoint = 0; iPoint < 48; iPoint++){
      grLaserEvent->SetPoint(iPoint,iPoint*8,LaserTrigWF[iPoint]-LaserPedestal[0]);      
      if( LaserMax < LaserTrigWF[iPoint] ){
	LaserMaxPoint = iPoint;
	LaserMax = LaserTrigWF[iPoint];
      }
    }
    funcLaser->SetParameter(1,LaserSignal[0]);
    funcLaser->SetParameter(0,LaserTime[0]);
    funcLaser->SetParLimits(1,LaserSignal[0]*0.8,LaserSignal[0]*1.2);
    funcLaser->SetParLimits(0,(LaserMaxPoint-1)*8,(LaserMaxPoint+1)*8);
    funcLaser->SetParameter(2,0);
    funcLaser->SetParLimits(2,-20,20);
    grLaserEvent->Fit(funcLaser,"Q","",LaserMaxPoint*8-50,LaserMaxPoint*8+20);
    LaserT      = funcLaser->GetParameter(0);
    LaserHeight = funcLaser->GetParameter(1);
    nDigi = 0; 
    for( int i = 0; i< 2716; i++){
      CsiMaxPoint[i] = 0;
      CsiMax[i] = 0.;
    }
    for( int iCh = 0; iCh < CsiNumber; iCh++){
      //std::cout<<ievent << " : " <<  iCh << std::endl;
      //if( tempCsi[CsiID[iCh]] == NULL ){ continue; }
      //if( tempCsi[CsiID[iCh]]->GetEntries() < 4800 ){ continue; }
      if( splCsiTemp[CsiID[iCh]] == NULL ){ continue; }
      if( splCsi[CsiID[iCh]] == NULL ){ continue; }
      spl = splCsi[CsiID[iCh]];
      TF1* funcCsi   = new TF1("funcCsi",FuncsplCsi,-100,200,2);      
      for( int iPoint = 0; iPoint < 48; iPoint++){	
	grCsiEvent[CsiID[iCh]]->SetPoint(iPoint,iPoint*8,ChannelOutWF[CsiID[iCh]][iPoint]-CsiPedestal[iCh]);
	if( CsiMax[CsiID[iCh]] < ChannelOutWF[CsiID[iCh]][iPoint] ){ 
	  CsiMaxPoint[CsiID[iCh]] = iPoint;
	  CsiMax[CsiID[iCh]] = ChannelOutWF[CsiID[iCh]][iPoint];
	}
      }
      funcCsi->SetParameter(1,CsiSignal[iCh]);
      funcCsi->SetParameter(0,CsiTime[iCh]);
      funcCsi->SetParLimits(1,CsiSignal[iCh]*0.8,CsiSignal[iCh]*1.2);
      funcCsi->SetParLimits(0,(CsiMaxPoint[CsiID[iCh]]-1)*8,(CsiMaxPoint[CsiID[iCh]]+1)*8);
      // -70,20 -- initial
      grCsiEvent[CsiID[iCh]]->Fit( funcCsi,"Q","",CsiMaxPoint[CsiID[iCh]]*8-100,CsiMaxPoint[CsiID[iCh]]*8+50);
      CsiModID[nDigi]  = CsiID[iCh];
      CsiT[nDigi]      = funcCsi->GetParameter(0);
      CsiHeight[nDigi] = funcCsi->GetParameter(1);
      nDigi++;
      
      /*
      hisDelta[CsiID[iCh]]->Fill(funcCsi->GetParameter(0)-funcLaser->GetParameter(0));
      if( funcCsi->GetParameter(1) > 200 && funcCsi->GetParameter(1) < 400 ){
	hisDeltaStd[CsiID[iCh]]->Fill(funcCsi->GetParameter(0)-funcLaser->GetParameter(0));
      }
      hisDelta2D[CsiID[iCh]]   ->Fill(funcCsi->GetParameter(1),funcCsi->GetParameter(0)-funcLaser->GetParameter(0));
      hisDeltaMin2D[CsiID[iCh]]->Fill(funcCsi->GetParameter(1),funcCsi->GetParameter(0)-funcLaser->GetParameter(0));
      */
    }    
    grLaserEvent->GetListOfFunctions()->Delete();
    for( int i = 0; i< 2716; i++){
      grCsiEvent[i]->GetListOfFunctions()->Delete();
    }
    trout->Fill();
  }
  std::cout<< "Write to File" << std::endl;
  trout->Write();
  
  /*
  for( int i = 0; i< 2716; i++){
    hisDelta[i]->Write();
    hisDelta2D[i]->Write();
    hisDeltaStd[i]->Write();
    hisDeltaMin2D[i]->Write();
  }
  */
  tfOut->Close();
  std::cout<< "Close File" << std::endl;
  
}
