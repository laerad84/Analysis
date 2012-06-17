//============================================================================
// Name        : IDHandler.cpp
// Author      : Takahiko Masuda
// Version     : 1.0
// Copyright   : taka
// Date        : 2010/11/07
//============================================================================
#include "IDHandler.h"
#include <string>
#include <cstring>
#if !defined(__CINT__)
ClassImp(IDHandler)
#endif


IDHandler::IDHandler(){
  std::string dataFilePathStr = std::getenv("ANALYSISLIB");
  dataFilePathStr +="/Data/crystal.txt";
  std::cout << dataFilePathStr << std::endl;
  // Load crystal positions and fills position data to crystalPositions map
  int crystalID=0, xPosition=0, yPosition=0;
  double xMetric=0, yMetric=0;
  char dummy[256] = {0};
  std::map<std::string,int> crystalMap;

  std::ifstream crystalFile( dataFilePathStr.c_str() );
  if( !crystalFile.is_open() ){
    std::cerr << "File don't exist : " << dataFilePathStr << std::endl;    
  }
  crystalFile.getline( dummy, 256); // remove header
  while( crystalFile >> crystalID >> xMetric >> yMetric >> xPosition >> yPosition ){
    crystalMap["crystalID"] = crystalID;
    crystalMap["xMetric"]   = (int)(xMetric*1000);
    crystalMap["yMetric"]   = (int)(yMetric*1000);
    crystalMap["xPosition"] = xPosition;
    crystalMap["yPosition"] = yPosition;

    crystalMaps.push_back( crystalMap );
    //    std::cout << crystalID << std::endl;
  }
  crystalFile.close();

  for( int crystal=0; crystal<numberOfChannels; crystal++ ){
    crystalMaps[crystal].insert( std::pair<std::string,int>("crate",-1) );
    crystalMaps[crystal].insert( std::pair<std::string,int>("module",-1) );
    crystalMaps[crystal].insert( std::pair<std::string,int>("channel",-1) );
  }
}

IDHandler::~IDHandler(){
}

void IDHandler::Show(){
  std::cout << "crystal\txPosition\tyPosition\txMetric\tyMetric" << std::endl;
  for( int index=0; index<numberOfChannels; index++ ){
    std::cout << crystalMaps[index]["crystalID"]     << "\t" 
	      << crystalMaps[index]["xPosition"]     << "\t" << crystalMaps[index]["yPosition"]    << "\t" 
	      << crystalMaps[index]["xMetric"]/1000. << "\t" << crystalMaps[index]["yMetric"]/1000. 
	      << std::endl;
  }
}

bool IDHandler::GetPosition( int crystalID, int& xPosition, int& yPosition ){
  if( !isCrystal(crystalID) ) return false;

  int index = GetIndex( crystalID );

  xPosition = crystalMaps[index]["xPosition"];
  yPosition = crystalMaps[index]["yPosition"];

  return true;
}



bool IDHandler::GetMetricPosition( int crystalID, double& xMetric, double& yMetric ){
  if( !isChannel(crystalID) ) return false;

  int index = GetIndex( crystalID );

  xMetric = crystalMaps[index]["xMetric"]/1000.;
  yMetric = crystalMaps[index]["yMetric"]/1000.;

  return true;
}


const int IDHandler::GetNumberOfCrystals( void ){
  return numberOfCrystals;
}


bool IDHandler::isSmall( int crystalID ){
  if( crystalID < 0 ) return false;
  if( crystalID >= numberOfSmall ) return false;
  return true;
}


bool IDHandler::isLarge( int crystalID ){
  if( crystalID < numberOfSmall ) return false;
  if( crystalID >= numberOfCrystals ) return false;
  return true;
}


bool IDHandler::isCrystal( int crystalID ){
  if( crystalID < 0 ) return false;
  if( crystalID >= numberOfCrystals ) return false;
  return true;
}


bool IDHandler::isCC03( int crystalID ){
   if( crystalID < 3000 ) return false;
   if( crystalID > 3031 ) return false;
   return true;
 }


bool IDHandler::isChannel( int crystalID ){
  if( isCrystal(crystalID) ) return true;
  if( isCC03(crystalID) )    return true;
  return false;
}


int  IDHandler::GetIndex( int crystalID ){
  if( !isChannel( crystalID ) ) return -1;

  // For CsI crystal, crystalID==index
  if( isCrystal( crystalID ) ) return crystalID;
  
  for( int index=numberOfCrystals; index<numberOfChannels; index++ ){
    if( crystalMaps[index]["crystalID"]==crystalID ){
      return index;
    }
  }
  return -1;  // Not found.
}
