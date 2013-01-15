#include "E14Mapper.h"
#include <iostream>
#include <iomanip>

E14Mapper::E14Mapper() : TObject()
{

  AllNum = 0;
  DataStoreTag = false;

  ThisID = 9999;

}


E14Mapper::~E14Mapper()
{
  ;
}

void E14Mapper::initializeDataValues(){

  for(int i=0;i<4096;i++){
    CrateID[i] = 9999;
    FADCID[i]  = 9999;
    CHID[i]    = 9999;
  }

}

int E14Mapper::dump(){

  printf("\nDump mapping data\n");
  for(int i=0;i<AllNum;i++){
    printf("e14 ID = %d : ",i);
    printf("Crate ID = %d : ",CrateID[i]);
    printf("FADC ID = %d : ",FADCID[i]);
    printf("CH ID = %d\n",CHID[i]);
  }

  return AllNum;
}

int E14Mapper::GetAllNum(){

  if( DataStoreTag == false ){
    printf("Data is not stored yet... aborting...\n");
    return -1;
  }

  return AllNum;
}

void E14Mapper::GetDataFromText( char *InFileName ){

  int detID = 9999;
  
  initializeDataValues();

  FILE *fp;
  fp = fopen( InFileName ,"rt");
  if( fp == NULL ){
    printf("Mapping file %s is not found.\n", InFileName);
    return;
  }
  printf("Mapping file is open as %s.\n",InFileName);
  
  fscanf(fp,"%d\n",&AllNum);  
  for(int i=0;i<AllNum;i++){
    fscanf(fp,"%d",&detID);
    fscanf(fp,"%d %d %d\n",&CrateID[detID], &FADCID[detID], &CHID[detID] );    
  }
  fclose( fp );
  
  DataStoreTag = true;
}

// Not implement yet
void E14Mapper::GetDataFromDB(){

  if( ThisID == 9999 ){
    printf("Set ID for mapping file. e.g. = 0 for CsI\n");
    return;
  }
  
  DataStoreTag = true;
}

void E14Mapper::ShowModuleInfo( int e14id ){
  
  if( e14id < 0 || AllNum < e14id ){
    printf("e14id should be less than %d\n",e14id);
  }

  printf("e14id = %d : Crate ID = %d : FADC ID = %d : CH ID = %d\n",
	 e14id, CrateID[e14id], FADCID[e14id], CHID[e14id] );
}

int E14Mapper::ShowE14ID( int crate_id, int fadc_id, int ch_id ){
  
  for(int i=0;i<AllNum;i++){
    if( CrateID[i] == crate_id && FADCID[i] == fadc_id && CHID[i] == ch_id){
      printf("e14id = %d\n",i);
      return i;
    }    
  }
  
  printf("Crate ID = %d : FADC ID = %d : CH ID = %d is not found\n",
	 crate_id, fadc_id, ch_id);
  return 0;
}

int E14Mapper::GetE14ID( int crate_id, int fadc_id, int ch_id ){
  
  for(int i=0;i<AllNum;i++){
    if( CrateID[i] == crate_id && FADCID[i] == fadc_id && CHID[i] == ch_id){
      return i;
    }    
  }
  
  return 9999;
}

int E14Mapper::VerifyData(){

  int BadNum = 0;

  if( DataStoreTag == false ){
    printf("Data is not stored yet... aborting...\n");
    return -1;    
  }

  // check range; upto nFADC=20, nCH=16 
  for(int i=0;i<AllNum;i++){
    if( FADCID[i] == 9999 ) continue;    
    if( FADCID[i] < 0 || 20 < FADCID[i] ){
      printf("e14id=%d : FADC ID %d is over range\n",i,FADCID[i]);
      BadNum += 1;
    }
  }
  for(int i=0;i<AllNum;i++){
    if( CHID[i] == 9999 ) continue;    
    if( CHID[i] < 0 || 16 < CHID[i] ){
      printf("e14id=%d : CH ID %d is over range\n",i,CHID[i]);
      BadNum += 1;
    }
  }

  // check duplication
  for(int i=0; i<AllNum; i++){
    if( CrateID[i] == 9999 ) continue;

    int tmpCrateID = CrateID[i];
    int tmpFACID   = FADCID[i];
    int tmpCHID    = CHID[i];

    for(int j=(i+1); j<AllNum; j++){
      if( CrateID[j] == tmpCrateID && FADCID[j] == tmpFACID &&
	  CHID[j] == tmpCHID ){
	printf("e14id=%d is overlapped to ID %d: CrateID=%d FADCID=%d CHID=%d\n",i,j,tmpCrateID,tmpFACID,tmpCHID);
	BadNum += 1;
      }
    }    

  }

  return BadNum;

}
