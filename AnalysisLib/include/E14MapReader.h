#ifndef E14MAPREADER__H__
#define E14MAPREADER__H__

#include "E14ModuleMap.h"
#include "TTree.h"
#include "TFile.h"
#include <string>
#include <iostream>

#include "Structs.h"
////////////////////////////////////////////////////////////////////////////////////////////////////
// Class Usage :
//  E14MapReader* map = new MapReader("SumupFile-ROOTfile")
//  map->Add("ModuleName: ex> Csi")
//  add need modules......
//  map->Setup(); // After Setup, Add new modules is disable.
//  int CsiMap[4096][3]; and other modules ... // [0]:crate, [1]:FADC, [2]:channel; 
//  int nCsiMod = map->CopyCFCtoMap(CsiMap);
////////////////////////////////////////////////////////////////////////////////////////////////////

class E14MapReader{
 public:
  static const int nMaxModule =  128;
  E14ModuleMap* module[nMaxModule];
 private:
  std::string   name; 
  int  mNmodule;
  bool bSet;
  bool bAdd;

  TFile* mfile;
  TTree* mtree; 
 public:
  E14MapReader(const char*);
  ~E14MapReader();
  void Init();
  bool Add(const char*);
  bool SetMap();
  bool VerifyModuleID( int );
  bool VerifyModuleCHID(int ,int );
  bool GetCFC(const char*, int, int& , int&, int& );
  int  FindModuleNumber( const char* );
  int  GetNmodule() const { return mNmodule; }
  int  CopyCFCtoMap( int, int[4096][3]);
  int  CopyCFCtoMap( const char*,int[4096][3] );
  bool CopyMap( int, struct MapStruct& );
  const char* GetModuleName( int );  
};

#endif
