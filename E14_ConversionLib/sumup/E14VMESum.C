#include "E14VMESum.h"
#include <iostream>
#include <iomanip>

#if !defined(__CINT__)
ClassImp(E14VMESum)
#endif

using namespace std;

E14VMESum::E14VMESum()
{  

  nCrate = 11;  // for Vacuum test
  sprintf( ConvCompName[0], "compute-1-2"); ConvCrateID[0] = 0;
  sprintf( ConvCompName[1], "compute-1-2"); ConvCrateID[1] = 1;
  sprintf( ConvCompName[2], "compute-1-3"); ConvCrateID[2] = 2;
  sprintf( ConvCompName[3], "compute-1-3"); ConvCrateID[3] = 3;
  sprintf( ConvCompName[4], "compute-1-4"); ConvCrateID[4] = 4;
  sprintf( ConvCompName[5], "compute-1-4"); ConvCrateID[5] = 5;
  sprintf( ConvCompName[6], "compute-1-5"); ConvCrateID[6] = 6;
  sprintf( ConvCompName[7], "compute-1-5"); ConvCrateID[7] = 7;
  sprintf( ConvCompName[8], "compute-1-6"); ConvCrateID[8] = 8;
  sprintf( ConvCompName[9], "compute-1-6"); ConvCrateID[9] = 9;
  sprintf( ConvCompName[10], "compute-1-7"); ConvCrateID[10] = 10;

  for(int nc=0;nc<nCrate;nc++){
    eventTree[ConvCrateID[nc]] = new EventTree();
    //    eventTree[CrateID[nc]] -> OpenFile( Form("/disk/%s/partition2/conv_data/crate%d/") );
    INFILEJUDGE[nc] = true;
  }

  RunNumber = 0;
  //  EnergyThreshold = 0;
  EnergyThreshold = 1.0;

  readend = 0;
  //  DebugMode = false;

  NDet = 7;
  DetName[0] = "Csi";
  DetName[1] = "CC03";
  DetName[2] = "OEV";
  DetName[3] = "CV";
  DetName[4] = "Cosmic";
  DetName[5] = "Laser";
  DetName[6] = "Etc";

  char name[256];
  for(int idet=0;idet<NDet;idet++){
    sprintf(name,"/share/apps/production/2012VME/sumup/map/ch_map_%s.txt",DetName[idet].c_str());
    map[idet].GetDataFromText( name );
    sprintf(name,"/share/apps/production/2012VME/sumup/calib/calib_%s.txt",DetName[idet].c_str());
    calib[idet].GetDataFromText( name );
  }

  sol = new SemiOnlinePlot();

}

E14VMESum::~E14VMESum()
{
  for(int nc=0;nc<nCrate;nc++)  delete eventTree[ConvCrateID[nc]];
  //  delete sol;
}

int E14VMESum::SetRunNumber( int num )
{
  RunNumber = num;
  
  FILE *fp;
  char ConvFileName[256] = {0};
  
  for(int nc=0;nc<nCrate;nc++){
    sprintf(ConvFileName,"/disk/%s/partition2/conv_data/crate%d/run%04d_conv.root",ConvCompName[nc],ConvCrateID[nc],RunNumber);
    fp = fopen( ConvFileName ,"r");
    if( fp == NULL){
      printf("%s is not found\n",ConvFileName);
      return 0;
    }else{
      fclose(fp);

      printf("%s is open\n",ConvFileName);      
      eventTree[ConvCrateID[nc]] -> FileOpen(ConvFileName);
    }
  }

  sol -> Init( RunNumber );

  return 1;
}

void E14VMESum::InitTree()
{
  char name[256];

  sprintf(name,"/state/partition2/sumup_data/Sum%04d.root",RunNumber);
  hfile  = new TFile(name,"RECREATE");
  calibTree = new TTree("calib","Calibration data");
  mapTree   = new TTree("map","Mapping data");
  for(int i=0;i<NDet;i++){
    string CurDetName = DetName[i] + ".";
    calibTree -> Branch(CurDetName.c_str(), &calib[i]);    
    mapTree   -> Branch(CurDetName.c_str(), &map[i]);
  }
  tree   = new TTree("T","2012 run tree");

  char bname1[256]={0};
  char bname2[256]={0};
  tree -> Branch("EventNo",&EventNo,"EventNo/I");
  tree -> Branch("SpillNo",&SpillNo,"SpillNo/S");
  tree -> Branch("TimeStamp",&TimeStamp,"TimeStamp/I");
  tree -> Branch("Error",&Error,"Error/S");
  for(int idet=0;idet<NDet;idet++){
    sprintf( bname1, "%sNumber", DetName[idet].c_str() );
    sprintf( bname2, "%sNumber/I", DetName[idet].c_str() );
    tree -> Branch( bname1, &DetNumber[idet], bname2 );

    sprintf( bname1, "%sModID", DetName[idet].c_str() );
    sprintf( bname2, "%sModID[%sNumber]/I", DetName[idet].c_str(), 
	     DetName[idet].c_str() );
    tree -> Branch( bname1, &DetModID[idet][0], bname2 );

    sprintf( bname1, "%sEne", DetName[idet].c_str() );
    sprintf( bname2, "%sEne[%sNumber]/D", DetName[idet].c_str(), 
	     DetName[idet].c_str() );
    tree -> Branch( bname1, &DetEne[idet][0], bname2 );

    sprintf( bname1, "%sIntegratedADC", DetName[idet].c_str() );
    sprintf( bname2, "%sIntegratedADC[%sNumber]/D", DetName[idet].c_str(), 
	     DetName[idet].c_str() );
    tree -> Branch( bname1, &DetIntegratedADC[idet][0], bname2 );

    sprintf( bname1, "%sTime", DetName[idet].c_str() );
    sprintf( bname2, "%sTime[%sNumber]/D", DetName[idet].c_str(), 
	     DetName[idet].c_str() );
    tree -> Branch( bname1, &DetTime[idet][0], bname2 );
  }

}

void E14VMESum::InitStorage()
{
  EventNo = 0;
  TimeStamp = 0;
  Error = 0;
  for(int idet=0; idet<NDet; idet++){
    DetNumber[idet] = 0;
    for(int ich=0; ich<map[idet].GetAllNum(); ich++){
      DetModID[idet][ich] = -9999;
      DetEne[idet][ich]   = -9999;
      DetIntegratedADC[idet][ich]   = -9999;
      DetTime[idet][ich]  = -9999;
    }
  }
}

void E14VMESum::WriteTree()
{
  hfile -> cd();
  calibTree -> Fill();
  calibTree -> Write();
  mapTree -> Fill();
  mapTree -> Write();
  tree -> Write();
  hfile -> Close();
}

void E14VMESum::PrintParameters()
{
  printf("# of crate = %d\n",nCrate);
}

void E14VMESum::Fill()
{
  
  for(int nc=0;nc<nCrate;nc++){
    INFILEJUDGE[nc] = true;
  }  

  InitTree();

  int MinEntry = 99999999;
  //  int MinEntry = 1000;
  for(int nc=0;nc<nCrate;nc++){
    int nentry = eventTree[ConvCrateID[nc]] -> GetEntriesFast();
    printf("Event entry for crate%d = %d\n",ConvCrateID[nc],nentry);
    if( nentry < MinEntry ) MinEntry = nentry;
  }
  printf("%d entries will be summing up.\n",MinEntry);

  //  MinEntry = 500;
  for(int jentry=0;jentry<MinEntry;jentry++){
    if( jentry%5000 == 0 ) printf("Event = %d\n",jentry);
    //    printf("Event = %d\n",jentry);

    InitStorage();

    for(int nc=0;nc<nCrate;nc++){
      eventTree[ConvCrateID[nc]] -> GetEntry(jentry);      
    }

    /*
    bool ERRORFLAG = false; 
    for(int nc=0;nc<nCrate;nc++){
      eventTree[ConvCrateID[nc]] -> GetEntry(jentry);      
      if( eventTree[ConvCrateID[nc]]->Error[0] == 1 ) ERRORFLAG = true;
    }
    //    printf("%d\n",ERRORFLAG);
    if( ERRORFLAG ) Error = 1;
    */

    if( !CheckTimeStampSync() ) Error = 1;
    //    printf("%d\n",Error);

    if( Error != 1){

      sol -> InitValue();
      for(int icrate=0;icrate<10;icrate++){
	for(int ifadc=0;ifadc<16;ifadc++){
	  sol -> SetEtSum(icrate, ifadc, eventTree[icrate]->EtSum_FADC[ifadc]);
	}
      }

      EventNo = eventTree[ConvCrateID[0]] -> EventNo; 
      SpillNo = eventTree[ConvCrateID[0]] -> SpillNo[0];
      TimeStamp = eventTree[ConvCrateID[0]] -> TimeStamp[0];
      
      for(int idet=0; idet<NDet; idet++){
	for(int ich=0; ich<map[idet].GetAllNum(); ich++){
	  if( map[idet].GetCrateID(ich) == 9999 ) continue;       
	  
	  int CrateID  = map[idet].GetCrateID(ich);
	  int FADCID   = map[idet].GetFADCID(ich);
	  int CHID     = map[idet].GetCHID(ich);
	  //	printf("%d : CrateID=%d FADCID=%d CHID=%d\n",ich,CrateID,FADCID,CHID);
	  
	  short peak = (short)eventTree[CrateID]->PeakHeight[FADCID][CHID];
	  double IntADC = (double)eventTree[CrateID]->IntegratedADC[FADCID][CHID];
	  double energy = IntADC;
	  double time = (double)eventTree[CrateID]->PeakTime[FADCID][CHID];
	  sol -> SetIntegratedADC(idet, ich, energy);
	  sol -> SetPeakHeight(idet, ich, peak);
	  
	  //	if( idet == 1 ) printf("%f\n",energy);
	  //    Calibration : Not impliment in fisrt stage; only csi so far...
	  if( idet == 0 ){
	    energy = calib[idet].CalcEnergy( ich, IntADC, peak );
	    if( energy < EnergyThreshold ) continue;
	  }
	  
	  sol -> SetDetEne(idet, ich, energy);
	  
	  DetModID[idet][DetNumber[idet]] = ich;
	  DetEne[idet][DetNumber[idet]] = energy;
	  DetIntegratedADC[idet][DetNumber[idet]] = IntADC;
	  DetTime[idet][DetNumber[idet]] = time;
	  DetNumber[idet]++;
	}
      }  

    }

    tree -> Fill();
    sol -> Fill( jentry );
    sol -> RunDisplay();
    if( !(jentry%1000) )  sol -> EventDisplay( jentry );
    //    if( !(jentry%100) )  sol -> EventDisplay( jentry );
    //    sol -> EventDisplay( jentry );
    
  }
  
  WriteTree();
  sol -> DrawAll();

  return;
}

bool E14VMESum::CheckTimeStampSync()
{
  bool SYNCJUDGE = true;

  // Allow +- 1
  for(int n=1;n<nCrate;n++){
    //    printf("%d %d\n",n,eventTree[ConvCrateID[n-1]]->TimeStamp[0]);
    if( abs(eventTree[ConvCrateID[n-1]]->TimeStamp[0]-eventTree[ConvCrateID[n]]->TimeStamp[0])>2 )  SYNCJUDGE = false;
  }

  if( !SYNCJUDGE ){        
    printf("Time stamp is not syncronized !!\n");
    for(int n=0;n<nCrate;n++){
      printf("TimeStamp in crate#%d = %d\n",
	     n,eventTree[ConvCrateID[n]]->TimeStamp[0]);
    }
    printf("\n");    
  }

  return SYNCJUDGE;

}


