#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>

#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>

//#include "TROOT.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TSystem.h"

#include "E14RawData.h"
#include "E14DataRead2010.h"
#include "E14Mapper.h"
#include "E14Calibrator.h"
#include "E14Fill.h"

int main(int argc, char** argv){
  std::cout <<__LINE__ <<  std::endl;
  /*
    if( argc != 2){
    std::cout << argv[0] << " RunNumber " << std::endl;
    return -1;
  }
  std::cout <<__LINE__ <<  std::endl;

  int RunNum = atoi( argv[1] );
  std::cout << RunNum  <<  std::endl;
  std::cout <<__LINE__ <<  std::endl;
  */
  E14Fill* a = new E14Fill;
  /*
  a.SetDebugMode();
  a.SetnCrate(13);
  std::cout <<__LINE__ <<  std::endl;
  a.SetInputDirectory("/disk/compute-0-1/partition2/raw_data");
  std::cout <<__LINE__ <<  std::endl;
  a.SetRunNumber(RunNum);
  std::cout <<__LINE__ <<  std::endl;
  a.Fill();
  */
}
