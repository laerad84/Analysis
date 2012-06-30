#include <iostream>
#include <fstream>

void test(){
  
  gSystem->Load("lib/libtest.so");
  IDHandler* handler = new IDHandler("Data/crystal.txt");
  std::ifstream ifs0("HV_MAP_2012_02_04_AfterAdjusted.txt");
  std::ifstream ifs1("HV_MAP_2012_02_04.txt");

  
  int ID;
  int Volt;
  
  int VoltList0[2716] = {0};
  int VoltList1[2716] = {0};
  int IDList[2716]    = {0};
  while( !ifs0.eof()){
    ifs0 >> ID >> Volt;
    VoltList0[ID] = Volt;
  }
  
  while( !ifs1.eof()){
    ifs1 >> ID >> Volt;
    VoltList1[ID] = Volt;
  }
  
  std::ofstream ofs("HV_MAP_2012_02_04_HALF.txt");
  for( int i = 0; i< 2716; i++){
    if( i >=2240     && i< 2240+120 ){continue;}
    if( i >=2716-120 && i< 2716     ){continue;}
    
    double x,y;
    handler->GetMetricPosition(i, x,y);
    if( x < 0 ){      
      ofs << i << "\t" << VoltList0[i] << std::endl;
    }else{
      ofs << i << "\t" << VoltList1[i] << std::endl;
    }
  } 
  ofs.close();
}  
    
