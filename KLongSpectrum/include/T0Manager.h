#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include <fstream>
#include <iostream>
#include "TROOT.h"

class T0Manager{
 public:
  double m_T0[2716];

 public:
  T0Manager();
  ~T0Manager();
  bool     ReadFile(char* filename);
  double   GetT0Offset(int id);
  void     Reset();
  void     PrintOffset(){
    for( int i = 0; i < 2716; i++){
      std::cout<< m_T0[i] << std::endl;
    }
  }
 private:
  void Normalizer();
  
};
