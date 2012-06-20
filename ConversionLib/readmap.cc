//#include "E14ReadMapper.h"
#include "E14ModuleMap.h"
#include "E14MapReader.h"
#include "TTree.h"
#include "TFile.h"
#include <string>

#include <iostream>

struct MapStruct {
  int modID;
  int nMod;
  int Map[4096][3];
  std::string name;
};

int
main( int argc, char** argv ){

  E14MapReader* map = new E14MapReader("Sum4557.root");
  map->Add("Csi");
  map->Add("CC03");
  map->Add("OEV");
  map->Add("Cosmic");
  map->Add("Laser");
  map->Add("CV");
  map->Add("Etc");
  map->SetMap();

  std::cout << " Setup Map is done. " << std::endl;
  std::cout<< map->GetNmodule() << std::endl; 

  struct MapStruct Csimap;
  Csimap.nMod = map->CopyCFCtoMap(0,Csimap.Map);
  Csimap.name = map->GetModuleName(0);
  std::cout << Csimap.name  << std::endl;
  struct MapStruct CC03map;
  CC03map.nMod = map->CopyCFCtoMap("CC03",Csimap.Map);
  CC03map.name ="CC03";
  struct MapStruct OEVmap;
  OEVmap.nMod = map->CopyCFCtoMap("OEV",OEVmap.Map );
  OEVmap.name = "OEV";
  struct MapStruct Cosmicmap;
  Cosmicmap.nMod = map->CopyCFCtoMap("Comsic",Cosmicmap.Map);
  Cosmicmap.name = "Cosmic";
  
  struct MapStruct Lasermap;
  Lasermap.nMod = map->CopyCFCtoMap("Laser", Lasermap.Map);
  Lasermap.name = "Laser";
  /*
  for( int i = 0; i <  map->module[0]->GetAllNum(); i++){
    int a= -1;
    int b= -1;
    int c= -1;
    map->module[0]->GetCFC(i, a,b,c);
    std::cout<< i << " : " << a << " : " << b << " : " << c << std::endl;
  }
  */
}
  
