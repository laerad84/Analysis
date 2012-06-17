#include <ctime>
#include <iostream>
#include <string>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TGraph.h"
void ReturnTime(time_t now, struct tm &tmTime);
const char* ReturnFilename(time_t now, struct tm &tmTime);
const char* ReturnFilename(int year, int month, int day){
  return Form("Test_%04d_%02d_%02d.root",year,month, day);
}


int main(int argc, char** argv ){
  if( argc != 2 ){
    std::cerr << "Please input " <<  argv[0] << " "  << "[RunNumber]" 
	      << std::endl;
  }
    
  std::string StrRunNumber = argv[1];

  TChain* ch = new TChain("KEYTHLEY");
  ///////////////////////////////////////////////////////////////////////////
  // Read Time Info
  ///////////////////////////////////////////////////////////////////////////
  time_t TimeKeythley;
  time_t TimeStartRun;
  time_t TimeStopRun;
  unsigned int iTimeStart;
  unsigned int iTimeStop;
  struct tm  StartTimeinfo;
  struct tm  StopTimeinfo;
  int* DaysPerMonth=NULL;
  int normalYear[12]={31,28,31,30,31,30,31,31,30,31,30,31};
  int leapYear[12]  ={31,29,31,30,31,30,31,31,30,31,30,31};

  TFile* tfRun = new TFile(Form("TrigCrateData/run4329.root",StrRunNumber.c_str()));
  TTree* trRun = (TTree*)tfRun->Get("runTree");

  trRun->SetBranchAddress("StartDate",&iTimeStart);
  trRun->SetBranchAddress("EndDate"  ,&iTimeStop);
  trRun->GetEntry(0);
  TimeStartRun = (time_t)iTimeStart;
  TimeStopRun  = (time_t)iTimeStop;

  // Set TChain 
  ReturnTime(TimeStartRun,StartTimeinfo);
  ReturnTime(TimeStopRun,StopTimeinfo);

  TFile* tfOut = new TFile(Form("RunInfo_Temp_%s.root",StrRunNumber.c_str()),"RECREATE");

  //Assume one run is less than 1day//
  if( (int)StartTimeinfo.tm_year == (int)StopTimeinfo.tm_year &&
      (int)StartTimeinfo.tm_mon  == (int)StopTimeinfo.tm_mon  &&
      (int)StartTimeinfo.tm_mday == (int)StopTimeinfo.tm_mday ){
    ch->Add(Form("EnvironmentData/KOTO_MONITOR_%04d_%02d_%02d.root",
		 (int)StartTimeinfo.tm_year+1900,
		 (int)StartTimeinfo.tm_mon+1,
		 (int)StartTimeinfo.tm_mday));
  }else{
    ch->Add(Form("EnvironmentData/KOTO_MONITOR_%04d_%02d_%02d.root",
		 (int)StartTimeinfo.tm_year+1900, 
		 (int)StartTimeinfo.tm_mon+1,
		 (int)StartTimeinfo.tm_mday));
    ch->Add(Form("EnvironmentData/KOTO_MONITOR_%04d_%02d_%02d.root",
		 (int)StopTimeinfo.tm_year+1900, 
		 (int)StopTimeinfo.tm_mon+1,
		 (int)StopTimeinfo.tm_mday));    
  }	
  
  double KTime;
  double Data[60];
  ch->SetBranchAddress("T",&KTime);
  ch->SetBranchAddress("Data",Data);
  TH1D* hisTemp[60];
  TGraph* grTemp[60];  
  for( int ihist = 0; ihist < 60; ihist++){
    hisTemp[ihist] = new TH1D(Form("hisTemp_CH%02d",ihist),
			      Form("hisTemp_CH%02d;Temp[ºC]",ihist),
			      400,0,100);
    grTemp[ihist] = new TGraph();
    grTemp[ihist]->SetNameTitle(Form("grTemp_CH%02d",ihist),
				Form("grTemp_CH%02d;Time[sec];Temp[ºC]",ihist));
  }
  
  for( int ievent = 0; ievent < ch->GetEntries(); ievent++){
    ch->GetEntry(ievent);
    if(TimeStartRun > KTime ){ continue;}
    if(TimeStopRun  < KTime){ break;}
    for( int ihist = 0; ihist < 60; ihist++){
      if( Data[ihist] < 100 && Data[ihist]>=0 ){
	hisTemp[ihist]->Fill(Data[ihist]);
	grTemp[ihist]->SetPoint(grTemp[ihist]->GetN(),KTime-TimeStartRun,Data[ihist]);
      }
    }
  }  

  for( int ihist = 0; ihist <60; ihist++){
    hisTemp[ihist]->Write();
    grTemp[ihist]->Write();
  }
  tfOut->Close();

  return 0;
}

void ReturnTime(time_t now, struct tm &tmTime){
  struct tm * tmNow;
  tmNow  = localtime(&now);
  tmTime.tm_year = tmNow->tm_year;
  tmTime.tm_mon  = tmNow->tm_mon;
  tmTime.tm_mday = tmNow->tm_mday;
  tmTime.tm_hour = tmNow->tm_hour;
  tmTime.tm_min  = tmNow->tm_min;
  tmTime.tm_sec  = tmNow->tm_sec;
}

const char*  ReturnFilename(time_t now, struct tm &tmTime){
  ReturnTime(now,tmTime);
  std::string filename = Form("Test_%04d_%02d_%02d.root",
			      tmTime.tm_year+1900,
			      tmTime.tm_mon+1,
			      tmTime.tm_mday);
  return filename.c_str();
}
