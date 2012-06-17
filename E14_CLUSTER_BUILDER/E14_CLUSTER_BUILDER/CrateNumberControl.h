#ifndef CRATENUMBERCONTROL__H__
#define CRATENUMBERCONTROL__H__


#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <cstdlib>
#include <cstdio>

static const short NumberCsI = 2716;
static const short NumberCrate = 20;

class CrateNumberControl{
 private:
  
  short L1ID[NumberCsI];
  std::vector<short> CsICHList[NumberCrate];
  std::map<short,short> CsICHListMap;
  
 public:
  CrateNumberControl();
  ~CrateNumberControl();
  short GetL1ID( short  csiID);
  std::vector<short> GetCsIList(short l1ID);  
};

#endif
