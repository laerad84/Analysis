#include "Environment.h"

void PrintEnvironment(){

  std::cout << "SUMUPFILEDIR  :" << sumFileDir     << std::endl;
  std::cout << "CONVFILEDIR   :" << convFileDir    << std::endl;
  std::cout << "RAWFILEDIR    :" << rawFileDir     << std::endl;
  std::cout << "EXPCALFILEDIR :" << expcalFileDir  << std::endl;
  std::cout << "SIMCALFILEDIR :" << simcalFileDir  << std::endl;
  std::cout << "WAVEANAFILEDIR:" << waveAnaFileDir << std::endl; 
  return;
}

unsigned int GetEnvironment(){
  unsigned int Flag = 0;  
  /*
  sumFileDir  = std::getenv("SUMUPFILEDIR");
  convFileDir = std::getenv("CONVFILEDIR");
  rawFileDir  = std::getenv("RAWFILEDIR");
  expcalFileDir  = std::getenv("EXPCALFILEDIR");
  simcalFileDir  = std::getenv("SIMCALFILEDIR");
  */
  strcpy(sumFileDir    ,std::getenv("SUMUPFILEDIR"));
  strcpy(convFileDir   ,std::getenv("CONVFILEDIR"));
  strcpy(rawFileDir    ,std::getenv("RAWFILEDIR"));
  strcpy(expcalFileDir ,std::getenv("EXPCALFILEDIR"));
  strcpy(simcalFileDir ,std::getenv("SIMCALFILEDIR"));
  strcpy(waveAnaFileDir,std::getenv("WAVEANAFILEDIR"));
  /*
  if( sumFileDir.size() == 0){
    std::cerr << "There is no Definition : SUMUPFILEDIR" << std::endl;
    Flag |= 1;
  }
  if( convFileDir.size() == 0 ){
    std::cerr << "There is no Definition : CONVFILEDIR" << std::endl;
    Flag |= 2;
  }
  if( rawFileDir.size() == 0 ){
    std::cerr << "There is no Definition : RAWFILEDIR" << std::endl;
    Flag |=4;
  }
  if( expcalFileDir.size() == 0 ){
    std::cerr << "There is no Definition : EXPCALFILEDIR" << std::endl; 
    Flag |= 8; 
  }
  if( simcalFileDir.size() == 0 ){
    std::cerr << "There is no Definition : SIMCALFILEDIR" << std::endl;
    Flag |= 16;
  }
  */
  return Flag;
}

