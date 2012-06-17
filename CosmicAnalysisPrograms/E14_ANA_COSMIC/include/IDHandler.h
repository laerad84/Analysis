#ifndef IDHANDLER__H__
#define IDHANDLER__H__
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include <TROOT.h>
#include <TObject.h>

class IDHandler{
  
public:
  IDHandler( const char* crystalFileName);
  virtual ~IDHandler();
  
  bool GetPosition( int crystalID, int& xPosition, int& yPosition);
  bool GetMetricPosition( int crystalID, double& xMetric, double& yMetric );
  int  GetIndex( int crystalID );
  const int  GetNumberOfCrystals( void );
  void Show( void );
  
  
  bool isSmall( int crystalID );
  bool isLarge( int crystalID );
  bool isCrystal( int crystalID );
  bool isCC03( int crystalID );
  bool isChannel( int crystalID );


private:
  static const int numberOfSmall    = 2240;
  static const int numberOfLarge    = 476;
  static const int numberOfCrystals = 2716;
  static const int numberOfCC03     = 16*2;
  static const int numberOfChannels = 2716 + 16*2;
  
  std::vector< std::map<std::string,int> > crystalMaps;
  int SmallMapXY[48][48];
  int largeMapXY[38][38];
  ClassDef(IDHandler,1)

};

#endif
