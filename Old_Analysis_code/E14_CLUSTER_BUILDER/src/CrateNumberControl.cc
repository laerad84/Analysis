#include "E14_CLUSTER_BUILDER/CrateNumberControl.h"
#include <fstream>
#include <iostream>

CrateNumberControl::CrateNumberControl(){
  std::string dataFilePath=std::getenv("ANALYSISLIB");
  dataFilePath+="/Data/ch_map_CsI_L1.txt";
  
  short CsICH;
  short CrateID;
  short SlotID;
  short FADCCHID;
  short l1ID;
  
  std::ifstream ifs(dataFilePath.c_str());
  while( ifs >> CsICH >> CrateID >> SlotID >> FADCCHID >> l1ID ){    
    CsICHListMap.insert( std::map<short, short>::value_type( CsICH, l1ID ));
    L1ID[CsICH]=l1ID;
    CsICHList[l1ID].push_back(CsICH);
  }
}

CrateNumberControl::~CrateNumberControl(){
}

short CrateNumberControl::GetL1ID( short csiID){
  short l1ID = CsICHListMap[csiID];
  return l1ID;
}

std::vector<short> CrateNumberControl::GetCsIList(short l1ID){
  return CsICHList[l1ID];
}

