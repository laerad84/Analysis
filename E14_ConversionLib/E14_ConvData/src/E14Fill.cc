#include "E14Fill.h"
#include <iostream>
#include <iomanip>

#if !defined(__CINT__)
ClassImp(E14Fill)
#endif

using namespace std;

E14Fill::E14Fill()
{  

  //  nCrate = 6;  // for 2010 run
  nCrate = 1;  // for 2012 run
  for(int nc=0;nc<nCrate;nc++){
    INFILEJUDGE[nc] = true;
    e14data[nc] = new E14DataReadVME();
    rawdata[nc] = new E14RawData();
  }
  sol = new  SemiOfflineHist();

  RunNumber = 0;
  ADCThreshold = 500;

  readend = 0;

  //  DebugMode = false;

  /*
  char name[256];
  for(int i=0;i<NDet;i++){
    sprintf(name,"ch_map_%s.txt",DetName[i].c_str());
    map[i].GetDataFromText( name );
    sprintf(name,"calib_%s.txt",DetName[i].c_str());
    calib[i].GetDataFromText( name );
  }
  */
}


E14Fill::~E14Fill()
{
  ;
}

void E14Fill::PrintParameters()
{
  printf("# of crate = %d\n",nCrate);
}

void E14Fill::Fill()
{
  
  FILE *fp;
  for(int cn=0; cn<nCrate; cn++){
    //    sprintf(InFileName,"%s/crate%d/run%d.dat",
    //            id_rootFile.c_str(),cn,RunNumber);
    sprintf(InFileName,"/disk/%s/partition2/raw_data/crate%d/run%04d.dat",
            ic_rootFile.c_str(),CrateID,RunNumber);
    cout << "Input file = " << InFileName << endl;
    sprintf(OutFileName,"/state/partition2/conv_data/crate%d/run%04d_conv.root",
            CrateID,RunNumber);
    cout << "Output file = " << OutFileName << endl;

    fp = fopen( InFileName ,"r");
    if( fp==NULL ){
      INFILEJUDGE[cn] = false;
    }else{
      INFILEJUDGE[cn] = true;
      fclose(fp);
    }
  }

  PrintParameters();

  for(int cn=0; cn<nCrate; cn++){
    if( !INFILEJUDGE[cn] ) continue;

    // input file
    e14data[cn]->FileOpen( InFileName );
    e14data[cn]->HeaderRead();  
    // output file
    rawdata[cn]->InitTree( OutFileName ); 	
  }

  // For Semioffline hist
  sol->SetCrateID( CrateID );
  sol->SetRunNumber( RunNumber );
  sol->InitHist();

  while( readend==0 ){
  //  for(int evtn=0;evtn<10;evtn++){

    for(int cn=0; cn<nCrate; cn++){
      if( !INFILEJUDGE[cn] ) continue;
      int tmpreadend = e14data[cn]->GetBuffer();
      readend += tmpreadend;
    }
    if( readend>0 ) continue;
    
    // Get raw data from .dat file / e14data class to rawdata class
    for(int sn=0; sn<42; sn++){

      for(int cn=0; cn<nCrate; cn++){
	if( !INFILEJUDGE[cn] ) continue;

	rawdata[cn]->initializeDataValues(); 	
	e14data[cn]->DataRead( sn );

	rawdata[cn]->nFADC    = e14data[cn]->GetnFADC();
	rawdata[cn]->nSamples = e14data[cn]->GetnSamples();
	rawdata[cn]->EventNo  = e14data[cn]->GetEventNo();
	if( cn == 0 && (rawdata[cn]->EventNo%5000)==0 ){
	  printf("processing event# %d\n",rawdata[cn]->EventNo);
	}
	
	for(int i=0;i<20;i++){
	  for(int j=0;j<16;j++){
	    for(int k=0;k<48;k++){
	      rawdata[cn]->Data[i][j][k] = e14data[cn]->GetData(i,j,k);
	    }
	    rawdata[cn]->Compression_flag[i][j] = e14data[cn]->GetCompressionFlag(i,j);
	  }
	  rawdata[cn]->Error[i]     = e14data[cn]->GetError(i);
	  rawdata[cn]->TimeStamp[i] = e14data[cn]->GetTimeStamp(i);
	  rawdata[cn]->TrigNo[i]    = e14data[cn]->GetTrigNo(i);
	  rawdata[cn]->SpillNo[i]   = e14data[cn]->GetSpillNo(i);
	  rawdata[cn]->SlotNo[i]    = e14data[cn]->GetSlotNo(i);
	}
	
	rawdata[cn]->CalcAll();
	//	if( sn%2 == 0 ) rawdata[cn]->Fill();  
	rawdata[cn]->Fill();  
	
	// For Semi Offline
	//	sol->nFADC    = e14data[cn]->GetnFADC();
	sol->nSamples = e14data[cn]->GetnSamples();
	sol->EventNo  = e14data[cn]->GetEventNo();
	for(int i=0;i<20;i++){
	  for(int j=0;j<16;j++){
	    sol->IntegratedADC[i][j]  = rawdata[cn]->IntegratedADC[i][j];
	    sol->PeakHeight[i][j] = rawdata[cn]->PeakHeight[i][j];
	    sol->PeakTime[i][j]  = rawdata[cn]->PeakTime[i][j];
	    sol->Pedestal[i][j]  = rawdata[cn]->Pedestal[i][j];
	    for(int k=0;k<48;k++){
	      sol->Data[i][j][k]  = rawdata[cn]->Data[i][j][k];
	    }
	  }
	}
	sol->Fill();

      }

    }

  }
 
  for(int cn=0; cn<nCrate; cn++){
    rawdata[cn]->Write();    
  }
  sol->WriteAll();

  printf("end\n");

  for(int cn=0; cn<nCrate; cn++){
    delete e14data[cn];
    delete rawdata[cn];
  }

  //  return;
  
}

