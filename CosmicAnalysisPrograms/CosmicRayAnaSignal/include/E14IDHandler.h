#ifndef E14IDHANDLER__H__
#define E14IDHANDLER__H__
#include "GeneralTypes.h"

class E14IDHandler
{
 private:  
  short CsICrateMap[E14_NCsI];
  short CsISlotMap[E14_NCsI];
  short CsIChannelMap[E14_NCsI];
  short CsIL1Map[E14_NCsI];
  short CsIChMap[15][20][16];
  bool Init();
  bool ReadMap();

 public:

  E14IDHandler();
  ~E14IDHandler();

  short GetCrate(short CsINo)   const ;
  short GetSlot(short CsINo)    const ;
  short GetChannel(short CsINo) const ;
  short GetL1(short CsINo)      const ;
  short GetCsICH(short CrateNo, short SlotNo, short ChannelNo) const ;
};

#endif //E14IDHANDLER__H__
