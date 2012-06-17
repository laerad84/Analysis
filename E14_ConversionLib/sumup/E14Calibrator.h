#ifndef E14Calibrator_h
#define E14Calibrator_h

//includes
#include <TObject.h>
#include <TSpline.h>

#include <fstream>

class E14Calibrator : public TObject
{
 public:

  E14Calibrator();
  virtual ~E14Calibrator();
  
  void initializeDataValues();
  int dump();

  int SetThisID( int n ){ ThisID = n; return ThisID;};
 
  // To read mapping, choose one of followings.
  void GetDataFromText( char *InFileName );
  void GetDataFromDB( );

  // Data Accessor 1
  int   GetAllNum();
  double CalcEnergy( int ich, double IntegratedADC, short PeakHeight );
  double GetLinearityFactor(double peak);

  // Data Accessor 2 (for interactive)
  void ShowGainInfo( int e14id ); 

 private:
  int ThisID;
  int AllNum;
  bool DataStoreTag;
  float GainMean[4096];
  float GainSigma[4096];

  ClassDef(E14Calibrator,1)  
};

#endif // E14Calibrator_h
