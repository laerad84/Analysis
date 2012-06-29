#ifndef E14IDHANDLER__H__
#include "E14IDHandler.h"
#endif //E14IDHANDLER__H__
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <fstream>

#define DEBUG_TEST()\
  std::cout << __FILE__ << " : "<< __LINE__ << " : " << __FUNCTION__ << " : " << std::endl;
 
E14IDHandler::E14IDHandler(){
  DEBUG_TEST();
  Init();
  ReadMap();
}
E14IDHandler::~E14IDHandler(){
}

bool E14IDHandler::Init(){
  DEBUG_TEST();
  for( short icsi = 0 ; icsi < E14_NCsI; ++icsi){
    CsICrateMap[icsi]   = -1;
    CsISlotMap[icsi]    = -1;
    CsIChannelMap[icsi] = -1;
    CsIL1Map[icsi]      = -1; 
  }
  for( short icrate =0; icrate < 15; ++icrate){
    for( short islot  = 0; islot < 20; ++islot){
      for( short ichannel = 0; ichannel < 16; ++ichannel){
	CsIChMap[icrate][islot][ichannel] = -1;
      }
    }
  }
  return true;
}

bool E14IDHandler::ReadMap(){
  DEBUG_TEST();
  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  std::string MapFilename;
  MapFilename += ANALIBDIR.c_str();
  MapFilename += "/Data/ch_map_CsI_L1.txt";
  std::ifstream ifs(MapFilename.c_str());
  if( !ifs.is_open()){
    return false;
  }
  int icsi;
  int icrate;
  int islot;
  int ichannel;
  int il1;
  DEBUG_TEST();
  while( ifs >> icsi >> icrate >> islot >> ichannel >> il1 ){
    if( icrate > 20 || islot >20 || ichannel > 20 || il1 > 20 ){
      continue;
    }
    CsICrateMap[icsi]   = icrate;
    CsISlotMap[icsi]    = islot;
    CsIChannelMap[icsi] = ichannel;
    CsIL1Map[icsi]      = il1;
    CsIChMap[icrate][islot][ichannel] = icsi;
  }

  DEBUG_TEST();
  return true;
}

short E14IDHandler::GetCrate(short CsINo) const {
  return CsICrateMap[CsINo];
}
short E14IDHandler::GetSlot(short CsINo) const {
  return CsISlotMap[CsINo];
}
short E14IDHandler::GetChannel(short CsINo) const {
  return CsIChannelMap[CsINo];
}
short E14IDHandler::GetL1(short CsINo) const {
  return CsIL1Map[CsINo];
}
short E14IDHandler::GetCsICH( short CrateNo, short SlotNo, short ChannelNo)const
{
  return CsIChMap[CrateNo][SlotNo][ChannelNo];
}


			    
