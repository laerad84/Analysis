#include "CrateIDHandler.h"

CrateIDHandler::CrateIDHandler(){
  Set();
}
CrateIDHandler::~CrateIDHandler(){
  ;
}

void CrateIDHandler::Set(){
  std::string ANALYSISLIB = std::getenv("ANALYSISLIB");
  std::cout<< Form("READ MAP FILE : %s/Data/ch_map_CsI_L1.txt",ANALYSISLIB.c_str()) << std::endl;
  std::ifstream ifs(Form("%s/Data/ch_map_CsI_L1.txt",ANALYSISLIB.c_str()));
  if( ifs.is_open() ){
    std::cout<< "File Opened" << std::endl;
  }else{
    std::cout<< "File is not exist" << std::endl;
  }
  for( int i = 0; i< 2716; i++){
    Crate[i]   = -1;
    Slot[i]    = -1;
    Channel[i] = -1;
    L1[i]      = -1;
  }
  Int_t tmpID, tmpC, tmpS, tmpH, tmpL;  
  while( ifs >> tmpID >> tmpC >> tmpS >> tmpH >> tmpL ){
    Crate[tmpID] = tmpC;
    Slot[tmpID]  = tmpS;
    Channel[tmpID] = tmpH;
    L1[tmpID]      = tmpL;
  }
}

Short_t CrateIDHandler::GetCrate(int id){
  if( id < 0 || id >= 2716 ){ 
    return -1;
  }
  return Crate[id];
}
Short_t CrateIDHandler::GetSlot( int id ){
  if( id < 0 || id >= 2716 ){
    return -1;
  }
  return Slot[id];
}
Short_t CrateIDHandler::GetChannel( int id ){
  if( id < 0 || id >= 2716 ){
    return -1;
  }
  return Channel[id];
}
Short_t CrateIDHandler::GetL1( int id ){
  if( id < 0 || id >= 2716 ){
    return -1;
  }
  return L1[id];
}
