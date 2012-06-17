#ifndef E14Mapper_h
#define E14Mapper_h

//includes
#include <TObject.h>

class E14Mapper : public TObject
{
 public:

  E14Mapper();
  virtual ~E14Mapper();
  
  void initializeDataValues();
  int dump();

  int VerifyData();

  int SetThisID( int n ){ ThisID = n; return ThisID;};
 
  // To read mapping, choose one of followings.
  void GetDataFromText( char *InFileName );
  void GetDataFromDB();

  // Data Accessor 1
  int GetAllNum();
  int GetCrateID( int n ){ return CrateID[n]; }
  int GetFADCID( int n ){ return FADCID[n]; }
  int GetCHID( int n ){ return CHID[n]; }

  // Data Accessor 2 (for interactive)
  void ShowModuleInfo( int e14id ); 
  int ShowE14ID( int crate_id, int fadc_id, int ch_id );

 private:
  int ThisID;
  int AllNum;
  bool DataStoreTag;
  int CrateID[4096];
  int FADCID[4096];
  int CHID[4096];

  ClassDef(E14Mapper,1)  
};

#endif // E14Mapper_h
