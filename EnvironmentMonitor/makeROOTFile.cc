
#include <iostream>
#include <fstream>
#include <string>

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TThread.h"
#include "TApplication.h"
#include "TROOT.h"

#include "TFile.h"
#include "TTree.h"


void convertData( char* data, double *temperature ){
  char chtemp[180][30];
  int  ich   = 0; 
  int  icomp = 0;

  for( int cnt  = 1; cnt  < 4096; cnt++){
    if( data[cnt] == ',' ){
      chtemp[ich][icomp] = '\0';
      ich++;
      icomp = 0;
    }else{
      chtemp[ich][icomp] = data[cnt];
      icomp++;
    }
    if(data[cnt] =='\0'){
      break;
    }
  }
  for( int ich  = 0; ich < 60; ich++){
    temperature[ich] = atof(chtemp[ich*3]);
  }
}

void Usage(char* name){  
  std::cout << Form("Usage:: %s year month day", name) << std::endl;
  exit(-1);
}

int main(int argc, char** argv){
//void  makeROOTFile(){
  if( argc != 4 ){
    Usage(argv[0]);
  }

  int year  = atoi(argv[1]);
  int month = atoi(argv[2]); 
  int day   = atoi(argv[3]);
  const int days[12] = {31,29,31,30,31,30,31,31,30,31,30,31};
  if( month >12 || month < 1 || day >days[month-1] || day < 1){
    Usage(argv[0]);
  }
  

  char Kfilename[128];
  char Dfilename[128];
  sprintf(Kfilename,"EnvironmentData/nfs/KEYTHLEY_%04d_%02d_%02d.dat",year, month, day);
  sprintf(Dfilename,"EnvironmentData/nfs/DR230_%04d_%02d_%02d.dat",year, month, day);

  Double_t  KtimeStamp;
  Double_t  DtimeStamp;
  Double_t  Kch[60];
  Double_t  Dch[20];
  int       dummy;
  std::string KData;

  TFile* tf = new TFile(Form("EnvironmentData/KOTO_MONITOR_%04d_%02d_%02d.root",year,month,day),"recreate");
  std::cout << "OPEN DR230" << std::endl;
  TTree* Dtr = new TTree("DR230","DR230");
  Dtr->ReadFile(Dfilename,"T/D:CH00/I:E00/I:Data00/D:CH01/I:E01/I:Data01/D:CH02/I:E02/I:Data02/D");
  Dtr->Write();
  std::cout << "END WRITE DR230" << std::endl;
  std::cout << "OPEN KEYTHLEY" << std::endl;
  TTree* Ktr = new TTree("KEYTHLEY","KEYTHLEY");
  Ktr->Branch("T",&KtimeStamp,"T/D");
  Ktr->Branch("Data",Kch,"Data[60]/D");
  std::ifstream  Kifs(Kfilename);

  if( !Kifs.is_open() ){
    std::cout << "OpenError" << std::endl;
    exit(-1);
  }

  while( !Kifs.eof() ){
    Kifs >> KtimeStamp >> KData >> dummy;
    convertData( const_cast<char*>(KData.c_str()), Kch );
    if( KtimeStamp == 0){break;}
    printf("%09d\n",(int)KtimeStamp );
    //std::cout << KtimeStamp << " : " << Ktr->GetEntries() << std::endl;
    Ktr->Fill();
  }
  Ktr->Write();
  std::cout << "END WRITE KEYTHLEY" << std::endl;

  tf->Close();

}
