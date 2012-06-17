#include "E14DataRead2010.h"
#include <iostream>
#include <iomanip>

#if !defined(__CINT__)
ClassImp(E14DataRead2010)
#endif

  using namespace std;

E14DataRead2010::E14DataRead2010()
{
  EventNo = -1;
}

E14DataRead2010::~E14DataRead2010()
{
  if( fout ) fclose(fout);
  if( m_file ) m_file -> Close();
}


void E14DataRead2010::initializeDataValues(){

  for(int i=0;i<20;i++){
    for(int j=0;j<16;j++){
      for(int k=0;k<48;k++){
	Data[i][j][k] = 0;
      }
      Compression_flag[i][j] = 0; 
    }
    Error[i]     = 0;
    TimeStamp[i] = 0;
    TrigNo[i]    = 0;
    SpillNo[i]   = 0;
    SlotNo[i]    = 0;
  }
}

///////////////////////////////////////////////////////////
//////   For data get from dat    /////////////////////////
///////////////////////////////////////////////////////////

void E14DataRead2010::FileOpen( char *InFileName )
{

  fout = fopen(InFileName,"rb");
  //  setvbuf(fout,NULL,_IOFBF,512*1024);
  //  setvbuf(fout,NULL,_IOFBF,4096*1024);

  if( !fout ){
    cout << InFileName << " is not found" << endl;
    return;
  }
  cout << "Open " << InFileName << endl;

}

int E14DataRead2010::HeaderRead(){

  //  cout << DebugMode << endl;
  
  if( DebugMode ) printf("Header read is starting \n");
  
  unsigned int now = 0;
  int runID_check = 0;
  int crateID_data = 0;
  //  int nFADC = 0,nSamples = 0,bufSize = 0;

  readend=0;

  nread=fread((void*)&dp,sizeof(unsigned int),1,fout);
  if(nread!=1) return 1;
  
  if( dp==0x12341234 ){
    fread((void *)&runID_check,sizeof(int),1,fout);
    fread((void *)&crateID_data,sizeof(int),1,fout);
    fread((void *)&now,sizeof(unsigned int),1,fout);
    fread((void *)&nFADC,sizeof(int),1,fout);
    fread((void *)&nSamples,sizeof(int),1,fout);
    fread((void *)&bufSize,sizeof(int),1,fout);

    cout << "Run ID = "                << runID_check << endl;
    cout << "Cratea ID = "             << crateID_data<< endl;
    cout << "# of FADC = "             << nFADC       << endl;
    cout << "# of Sample = "           << nSamples    << endl;
    cout << "Buffer size in 1 loop = " << bufSize     << endl;
  }    
  
  if( dp==0x12121212) return 1;
  
  readend=0;
  if( DebugMode ) printf("Header read was succeeded \n");

  return 0;
  
}

int E14DataRead2010::GetBuffer(){
  
  dp = 0x0;
  while( dp != 0xaaaaaaaa ){
    nread=fread((void*)&dp,sizeof(unsigned int),1,fout);
    if( nread==0 ) return 1;
    
    if(dp==0x43214321){
      printf("convert finish \n");
      return 1;
    }
  }

  //  printf("dp = %x\n",dp);
  

  for(int j=0; j<nFADC; j++){
    int ZeroCount = 0;
    
    fread((void*)&bufData[j],bufSize*sizeof(unsigned short),1,fout);
    
    for(int i=0;i<bufSize;i++){
      if( bufData[j][i] == 0 ) ZeroCount++;      
    }

    /*
    for(int i=0;i<10;i++){
      printf("FADC#=%d %x\n",j,bufData[j][i]);
    }
    getchar();
    */

    if( ZeroCount > 260 ){
      for(int i=0;i<6;i++) printf("%x ",bufData[j][i+6]);
      printf("\nZeroCount = %d\n\n",ZeroCount);
      //	for(int i=0;i<6;i++) printf("%x ",bufData[j][i+32000]);
      //	printf("\n");
      getchar();	  
    }      
  }

  return 0;
}


int E14DataRead2010::DataRead( int sn ){

  int adata = 0;
  unsigned int headerdata = 0;
  
  long BufferLoopNo=0;
  
  unsigned int headersample[20][6] = {{0},{0}};
  
  bool b_error_header[20] = {0};  
  bool b_overth=kFALSE;
  
  for(int n=0;n<nFADC;n++) b_error_header[n]=kFALSE;
    
  //header check
  for(int n=0;n<nFADC;n++){	
    if( b_error_header[n] ) continue;	  
    
    //for header
    for(int h=0;h<6;h++){
      //	  if(b_error_header[n])continue;	  
      
      headerdata=bufData[n][6+(16*nSamples+6)*sn+h];
      if(((headerdata&0xc000)>>14)!=3){
	//if(l!=(totalloop-1)) std::cerr<<"not header data in loop:"<<l<<", event"<<m<<"FADC:"<<n<<"data="<<headerdata<<std::endl;
	std::cerr<< "No header data in loop:"<<BufferLoopNo<<","
		 << "event" <<sn
		 << "FADC:" <<n
		 << "data=" <<headerdata<<std::endl;
	b_error_header[n]=kTRUE;
	continue;
      }
      headersample[n][h]=headerdata;
      
    }
  }

  initializeDataValues();
      
  for(int n=0;n<nFADC;n++){	
    if(b_error_header[n])continue;  
    
    //for header
    for(int h=0;h<6;h++){
      if(b_error_header[n])continue;	  
      
      headerdata=bufData[n][6+(16*nSamples+6)*sn+h];
      /*if(((headdata&0xc000)>>14)!3){
	std::cerr<< "No header data in loop:"<< l <<","
	         << "event"                  << m
		 << "FADC:"                  << n
		 << "data="headerdata        << std::endl;
	b_error_header[n]=kTRUE;
	continue;
	}*/
      headersample[n][h]=headerdata&0x3fff;
    }
    
    //for data sample
    for ( int k=0;k<nSamples;k++ ){
      if(b_error_header[n])continue;
      
      for ( int j=0;j<16;j++ ){
	if(b_error_header[n])continue;
	  
	adata=bufData[n][(16*nSamples+6)*sn+12+16*k+j];
	/*
	  if(((adata&0xc000)>>14)!=2){
	  std::cerr<<"not sample data loop:"<<l<<", event"<<m<<"FADC:"<<n<<", ch:"<<j<<",sampletime "<<k<<"data="<<adata<<std::endl;
	  b_error_header=kTRUE;
	  continue;
	  //break;
	  }
	*/
	Data[n][j][k]=adata&0x3fff;
	
	b_overth=true;
	//if((adata&0x3fff)>600A)b_overth=true;
      }
    }
    //      if ( b_overflow ){ continue; }
  }//FADC
  
  
  
  //Sort and convert data      
  for(int n=0;n<nFADC;n++){
      
    if( b_error_header[n] ){
      Error[n]=1;
      continue;
    }
    
    TimeStamp[n] = ((headersample[n][1]&0x8ff)<<16)+(headersample[n][2]<<14)+headersample[n][3];
    TrigNo[n]    = (headersample[n][4]&0x3fc0)>>6;
    SpillNo[n]   = headersample[n][5];
    SlotNo[n]    = headersample[n][4]&0x1f;    
    int compword = (int)((headersample[n][0]&0xff))<<8+(headersample[n][1]&0xff);
    for(int i=0;i<16;i++){ Compression_flag[n][i]=compword; }
    
  }

  TimeStampJudge();
  EventNo++;

  return 0;
}


void E14DataRead2010::TimeStampJudge(){
  
  bool SyncroJudge = true;
  
  // Exact
  /*
    for(int n=1;n<nFADC;n++){
    //      printf("%d : %d\n",n,TimeStamp[n]);
    if( TimeStamp[0] != TimeStamp[n] ) SyncroJudge = false;
    }
  */
  
  // Allow +- 1

  for(int n=1;n<nFADC;n++){
    //      printf("%d : %d\n",n,TimeStamp[n]);
    if( abs(TimeStamp[0]-TimeStamp[n])>2 ) SyncroJudge = false;
  }
  
  if( !SyncroJudge ){
    
    if( DebugMode ){
      printf("Event No.= %d : TimeStamp = %d : Error = %d\n",
	     EventNo, TimeStamp[0], Error[0]);
      for(int n=1;n<nFADC;n++){
	printf("%d : spill No. = %d : TimeStamp = %d\n",
	       n, SpillNo[n], TimeStamp[n]);
	if( TimeStamp[0] != TimeStamp[n] ) 
	  SyncroJudge = false;
      }
    }
    
    for(int n=0;n<20;n++){ Error[n]=1; }
  }
  
}


///////////////////////////////////////////////////////////
//////   For data get from conv    ////////////////////////
///////////////////////////////////////////////////////////

int  E14DataRead2010::FileOpenConv( char *InFileName )
{

  m_file = new TFile( InFileName );

  if( !m_file ){
    cout << InFileName << " is not found" << endl;
    return 0;
  }
  cout << "Open " << InFileName << endl;

  m_tree = (TTree*)m_file->Get("EventTree");
  m_tree -> SetBranchAddress("nFADC",&nFADC);
  m_tree -> SetBranchAddress("nSamples",&nSamples);
  m_tree -> SetBranchAddress("EventNo",&EventNo);
  m_tree -> SetBranchAddress("Data",Data);
  m_tree -> SetBranchAddress("Error",Error);
  m_tree -> SetBranchAddress("TimeStamp",TimeStamp);
  m_tree -> SetBranchAddress("TrigNo",TrigNo);
  m_tree -> SetBranchAddress("SpillNo",SpillNo);
  m_tree -> SetBranchAddress("SlotNo",SlotNo);
  m_tree -> SetBranchAddress("Compression_flag",Compression_flag);

  return (int)( m_tree -> GetEntriesFast() );
}

void E14DataRead2010::DataReadConv( int ievt )
{
  m_tree -> GetEntry( ievt );

  return;
}
