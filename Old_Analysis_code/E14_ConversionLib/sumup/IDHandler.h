#ifndef IDHandler_h
#define IDHandler_h

//============================================================================
// Name        : IDHandler.cpp
// Author      : Takahiko Masuda
// Version     : 1.0
// Copyright   : taka
// Date        : 2010/11/07
//============================================================================

 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>



class IDHandler{
  
public:
  //  IDHandler( const char* crystalFileName, const char* FADCFileName );
  IDHandler( const char* crystalFileName);
  ~IDHandler();
  
  bool GetPosition( int crystalID, int& xPosition, int& yPosition);
  bool GetMetricPosition( int crystalID, double& xMetric, double& yMetric );
  int  GetCrystalID( int crate, int module, int channel );
  int  GetIndex( int crystalID );
  bool GetCMC( int crystalID, int& crate, int& module, int& channel );
  int  GetNumberOfCrates( void );
  const int  GetNumberOfCrystals( void );
  void Show( void );
  
  
  bool isSmall( int crystalID );
  bool isLarge( int crystalID );
  bool isCrystal( int crystalID );
  bool isCC03( int crystalID );
  bool isChannel( int crystalID );
  
private:
  int numberOfSmall;
  int numberOfLarge;
  int numberOfCrystals;
  int numberOfCC03;
  int numberOfChannels;
  
  std::vector< std::map<std::string,int> > crystalMaps;
  int CMCMap[6][15][16];

};

#endif
