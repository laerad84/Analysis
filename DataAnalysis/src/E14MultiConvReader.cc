#ifndef E14MULTICONVREADER__H__
#include "DataAnalysis/E14MultiConvReader.h"
#endif
#include "GeneralMacros.h"

E14MultiConvReader::E14MultiConvReader(int nCrate ){
  NumberOfCrate = nCrate;
  Init();
}

E14MultiConvReader::~E14MultiConvReader(){
  delete IDhandler;
  for( int i = 0; i< NumberOfCrate; ++i){
    delete conv[i];
  }
  delete [] conv;
}

bool E14MultiConvReader::Init(){
  std::cout<< "Init" << std::endl;
  for( int i = 0; i< NumberOfCrate; ++i){
    conv[i] = new E14ConvReader();
    conv[i]->SetBranchAddress();
  }    
  IDhandler = new E14IDHandler();
  return true;
}

int E14MultiConvReader::AddFile(const char* ConvFileDir,int RunNumber){
  std::string ReadFileForm = ConvFileDir;
  ReadFileForm += "/crate%d/run%d_conv.root";  
  int rst = 0;
  for( int i = 0; i< NumberOfCrate; ++i){    
    rst += conv[i]->AddFile(Form(ReadFileForm.c_str(),i,RunNumber));
  }
  return rst;
}

bool E14MultiConvReader::SetAddress(){

  for( int i = 0; i< NumberOfCrate; ++i){
    if(conv[i] ==NULL ){return false;}
    Data[i]          = &(conv[i]->Data[0][0][0]);
    PeakHeight[i]    = &(conv[i]->PeakHeight[0][0]);
    PeakTime[i]      = &(conv[i]->PeakTime[0][0]);
    Pedestal[i]      = &(conv[i]->Pedestal[0][0]);
    IntegratedADC[i] = &(conv[i]->IntegratedADC[0][0]);
    Error[i]         = &(conv[i]->Error[0]);    
  }
  for( int i = 0; i< E14_NCsI; ++i){
    short CrateNum   = IDhandler->GetCrate(i);
    std::cout << CrateNum << std::endl;
    if( CrateNum < 0 ){ CsIData[i] = NULL;}
    else{
      short SlotNum    = IDhandler->GetSlot(i);
      short ChannelNum = IDhandler->GetChannel(i);      
      CsIData[i] = conv[CrateNum]->Data[SlotNum][ChannelNum];
    }
  }
  return true;
}
 
int E14MultiConvReader::GetChainEntry(int ientry){
  int rst = 0; 
  for( int i = 0; i< NumberOfCrate; ++i){
    rst += conv[i]->GetChainEntry(ientry);
  }
  return rst;
}

long E14MultiConvReader::GetChainEntries(){
  int rst[N_MAXIMUM_CRATE_NUMBER]={0};
  for( int i = 0; i< NumberOfCrate; ++i){
    rst[i] = conv[i]->GetChainEntries();
    if( rst[i] != rst[0] ){
      return 0; 
    }    
  }
  return rst[0];
}


