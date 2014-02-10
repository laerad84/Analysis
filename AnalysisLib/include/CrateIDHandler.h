#ifndef CRATEIDHANDLER_H_
#define CRATEIDHANDLER_H_
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include "TROOT.h"

class CrateIDHandler {
 public:
  CrateIDHandler();
  ~CrateIDHandler();
  Short_t Crate[2716];
  Short_t Slot[2716];
  Short_t Channel[2716];
  Short_t L1[2716];

  void Set();
  Short_t GetCrate( int id );
  Short_t GetSlot( int id );
  Short_t GetChannel( int id );
  Short_t GetL1( int id );

};
#endif //CRATEIDHANDLER_H_
