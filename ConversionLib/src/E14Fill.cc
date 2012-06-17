#include "E14Fill.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>

#if !defined(__CINT__)
ClassImp(E14Fill)
#endif

using namespace std;

E14Fill::E14Fill()
{  
  char* Path = getenv("CONVERSIONLIBDIR");
  nCrate = 6;  // for 2010 run
  for(int nc=0;nc<nCrate;nc++){
    INFILEJUDGE[nc] = true;
  }
  RunNumber = 0;
  ADCThreshold = 500;

  readend = 0;

  //  DebugMode = false;

  NDet = 4;
  DetName[0] = "CSI";
  DetName[1] = "CC03";
  DetName[2] = "OEV";
  DetName[3] = "COSMIC";

  char name[256];
  for(int i=0;i<NDet;i++){
    sprintf(name,"%s/data/ch_map_%s.txt",Path,DetName[i].c_str());
    map[i].GetDataFromText( name );
    sprintf(name,"%s/data/calib_%s.txt",Path,DetName[i].c_str());
    calib[i].GetDataFromText( name );
    std::cout << name << std::endl;
  }

}


E14Fill::~E14Fill()
{
  ;
}

void E14Fill::InitTree()
{
  char name[256];

  sprintf(name,"RootFile/Raw%04d.root",RunNumber);
  hfile = new TFile(name,"RECREATE");
  for(int cn=0;cn<nCrate;cn++){
    if( !INFILEJUDGE[cn] ) continue;
    sprintf(name,"crate%d",cn);
    tree[cn] = new TTree(name,"");
    sprintf(name,"RawData.",cn);
    tree[cn] -> Branch(name,&rawdata[cn]);
  }

  sprintf(name,"RootFile/Out%04d.root",RunNumber);
  gsimFile  = new TFile(name,"RECREATE");
  calibTree = new TTree("calib","Calibration data");
  mapTree   = new TTree("map","Mapping data");
  gsimTree  = new TTree("eventTree00","Calibrated Tree");
  for(int i=0;i<NDet;i++){
    string CurDetName = DetName[i] + ".";
    calibTree -> Branch(CurDetName.c_str(), &calib[i]);    
    mapTree   -> Branch(CurDetName.c_str(), &map[i]);
    gsimTree  -> Branch(CurDetName.c_str(), &gsim[i]);
  }
  for(int i=0;i<NDet;i++){
    string CurDetName = DetName[i] + ".";
  }
  for(int i=0;i<NDet;i++){
    string CurDetName = DetName[i] + ".";
  }
  
}

void E14Fill::WriteTree()
{
  hfile -> cd();
  for(int cn=0;cn<nCrate;cn++){
    if( !INFILEJUDGE[cn] ) continue;
    tree[cn] -> Write();
  }
  hfile     -> Close();
  
  gsimFile  -> cd();
  calibTree -> Fill();
  calibTree -> Write();
  mapTree   -> Fill();
  mapTree   -> Write();
  gsimTree  -> Write();
  gsimFile  -> Close();
}

void E14Fill::PrintParameters()
{
  printf("# of crate = %d\n",nCrate);
}

void E14Fill::Fill()
{
  
  FILE *fp;
  for(int cn=0; cn<nCrate; cn++){
    sprintf(InFileName,"%s/crate%d/run%d.dat",
	    id_rootFile.c_str(),cn,RunNumber);
    cout << InFileName << endl;
    fp = fopen( InFileName ,"r");
    if( fp==NULL ){
      INFILEJUDGE[cn] = false;
    }else{
      INFILEJUDGE[cn] = true;
      fclose(fp);
    }
  }

  InitTree();
  PrintParameters();

  for(int cn=0; cn<nCrate; cn++){
    if( !INFILEJUDGE[cn] ) continue;

    sprintf(InFileName,"%s/crate%d/run%d.dat",
	    id_rootFile.c_str(),cn,RunNumber);
    
    e14data[cn].FileOpen( InFileName );
    e14data[cn].HeaderRead();  
  }

  while( readend==0 ){
  //  for(int evtn=0;evtn<10;evtn++){

    for(int cn=0; cn<nCrate; cn++){
      if( !INFILEJUDGE[cn] ) continue;
      readend += e14data[cn].GetBuffer();
    }
    if( readend != 0 ) continue;    

    // Get raw data from .dat file / e14data class to rawdata class
    for(int sn=0; sn<42; sn++){

      for(int cn=0; cn<nCrate; cn++){
	if( !INFILEJUDGE[cn] ) continue;

	rawdata[cn].initializeDataValues(); 	
	readend = e14data[cn].DataRead( sn );

	rawdata[cn].nFADC    = e14data[cn].GetnFADC();
	rawdata[cn].nSamples = e14data[cn].GetnSamples();
	rawdata[cn].EventNo  = e14data[cn].GetEventNo();
	if( cn == 0 && (rawdata[cn].EventNo%5000)==0 ){
	  printf("processing event# %d\n",rawdata[cn].EventNo);
	}
	
	for(int i=0;i<20;i++){
	  for(int j=0;j<16;j++){
	    for(int k=0;k<48;k++){
	      rawdata[cn].Data[i][j][k] = e14data[cn].GetData(i,j,k);
	    }
	    rawdata[cn].Compression_flag[i][j] = e14data[cn].GetCompressionFlag(i,j);
	  }
	  rawdata[cn].Error[i]     = e14data[cn].GetError(i);
	  rawdata[cn].TimeStamp[i] = e14data[cn].GetTimeStamp(i);
	  rawdata[cn].TrigNo[i]    = e14data[cn].GetTrigNo(i);
	  rawdata[cn].SpillNo[i]   = e14data[cn].GetSpillNo(i);
	  rawdata[cn].SlotNo[i]    = e14data[cn].GetSlotNo(i);
	}
	
	rawdata[cn].CalcAll();
	tree[cn] -> Fill();

      }

      // Calculate Digi      
      //      for(int idet=0; idet<NDet; idet++){
      for(int idet=0; idet<NDet; idet++){
	TClonesArray &DigiArray = *(gsim[idet].digi);
	if(DigiArray.GetEntriesFast()>0) gsim[idet].digi->Clear("C");
	
	GsimDigiData* aDigi=0;
	
	double TotE = 0;
	int iDigi=0;
	for(int ich=0; ich<map[idet].GetAllNum(); ich++){
	  if( map[idet].GetCrateID(ich) == 9999 ) continue;       
	  int CrateID  = map[idet].GetCrateID(ich);
	  int FADCID   = map[idet].GetFADCID(ich);
	  int CHID     = map[idet].GetCHID(ich);
	  float ADCSum = rawdata[CrateID].ADCSum[FADCID][CHID];

	  // ADC threshold : Not impliment in fisrt stage
	  // if( ADCSum < ADCThreshold ) continue;

	  // Calibration : Not impliment in fisrt stage
	  //  float energy = calib[idet].CalcEnergy( ich, ADCSum );
	  //	if( energy < ADCThreshold ) continue;
	  float energy = ADCSum;
	  
	  new (DigiArray[iDigi]) GsimDigiData();
	  aDigi = (GsimDigiData*)DigiArray[iDigi];
	  
	  aDigi->initializeDataValues();
	  aDigi->thisID = iDigi;
	  //aDigi->detID=m_sensitiveDetectorID;
	  aDigi->detID = idet;
	  aDigi->energy = energy;
	  TotE += energy;
	  aDigi->modID  =  ich;
	  //aDigi->mtimeEntry=i;
	  aDigi->mtimeSize=0;
	  iDigi++;
	}
	
	gsim[idet].nDigi = iDigi;
	gsim[idet].totalEnergy = TotE;
      }
      
      gsimTree -> Fill();
	
    }

 }

  WriteTree();
}

///////////////////////////////////////////////////////////////////////

void E14Fill::FillFromConv()
{
  InitTree();

  int Entries[16] = {0};

  for(int cn=0; cn<nCrate; cn++){
    sprintf(InFileName,"%s/crate%d/run%d_conv.root",
	    id_rootFile.c_str(),cn,RunNumber);

    Entries[cn] = e14data[cn].FileOpenConv( InFileName );
  }

  for(int evtn=0; evtn<Entries[0]; evtn++){  // # of event is set to Event[0]


    for(int cn=0; cn<nCrate; cn++){
      rawdata[cn].initializeDataValues(); 	
      e14data[cn].DataReadConv( evtn );

      rawdata[cn].nFADC    = e14data[cn].GetnFADC();
      rawdata[cn].nSamples = e14data[cn].GetnSamples();
      rawdata[cn].EventNo  = e14data[cn].GetEventNo();
      if( cn == 0 && (rawdata[cn].EventNo%5000)==0 ){
	printf("processing event# %d\n",rawdata[cn].EventNo);
      }
      
      for(int i=0;i<20;i++){
	for(int j=0;j<16;j++){
	  for(int k=0;k<48;k++){
	    rawdata[cn].Data[i][j][k] = e14data[cn].GetData(i,j,k);
	  }
	  rawdata[cn].Compression_flag[i][j] = e14data[cn].GetCompressionFlag(i,j);
	}
	rawdata[cn].Error[i]     = e14data[cn].GetError(i);
	rawdata[cn].TimeStamp[i] = e14data[cn].GetTimeStamp(i);
	rawdata[cn].TrigNo[i]    = e14data[cn].GetTrigNo(i);
	rawdata[cn].SpillNo[i]   = e14data[cn].GetSpillNo(i);
	rawdata[cn].SlotNo[i]    = e14data[cn].GetSlotNo(i);
      }
      
      rawdata[cn].CalcAll();
      //      tree[cn] -> Fill();

    }
  }
  
  WriteTree();
}
  
