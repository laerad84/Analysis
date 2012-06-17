#ifndef E14RawData_h
#define E14RawData_h

//includes
#include "TObject.h"
#include <cstring>
/**
 *  @class GsimDigiData
 *  @brief Digitization data.
 *  This class provides raw information.
 */
class E14RawData : public TObject
{
 public:

  E14RawData();
  virtual ~E14RawData();
  
  virtual void   Clear(Option_t* opt="");

  void initializeDataValues();

  void dump(int i, int j);

  std::string getClassName();
  UShort_t nFADC;
  UShort_t nSamples;
 
  Int_t     EventNo;

  /**
   *  @brief Raw data in FADC.
   */
  UShort_t     Data[20][16][48];

  /**
   *  @brief Rough calculation.
   */
  // Pedestal is mean value of first 3 points.
  Float_t     Pedestal[20][16];
  // Hight value of signal peak
  Float_t     PeakHeight[20][16];
  // Time at PeakH
  Float_t     PeakTime[20][16];
  // Integrated ADC value after pedestal subtraction.
  Float_t     ADCSum[20][16];
  // # of photo electron
  Float_t     pe[20][16];

  Int_t EtSum_FADC[20][48];
  //  Long_t BufferLoopNo;

  UShort_t Error[20];
  Int_t    TimeStamp[20];
  UShort_t TrigNo[20];
  Int_t    SpillNo[20];
  UShort_t SlotNo[20];
  Int_t    Compression_flag[20][16];

 public:
  virtual int CalcAll();
  
 private:
  
  virtual int CalcPedestal();
  virtual int CalcPeakHeight();
  virtual int CalcADCsum();
  virtual int CalcEtSumFADC();

 protected:
  std::string name;
  std::string className;
  ClassDef(E14RawData,0)
};

inline std::string E14RawData::getClassName(){
  return className;
}

#endif // E14RawData_h
